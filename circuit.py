from bfv.crt import CRTModuli, CRTPolynomial
from bfv.bfv import BFV, RLWE
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div
from bfv.utils import adjust_negative_coefficients
from random import randint
import copy

def main(): 

    # SETUP
    # q is the ciphertext modulus, which is the product of all qis
    # t is the plaintext modulus
    # n is the degree of the cyclotomic polynomial which is the denominator of the polynomial ring
    # sigma is the standard deviation of the discrete Gaussian distribution
    qis = [
        1152921504606584833,
        1152921504598720513,
        1152921504597016577,
        1152921504595968001,
        1152921504595640321,
        1152921504593412097,
        1152921504592822273,
        1152921504592429057,
        1152921504589938689,
        1152921504586530817,
        1152921504585547777,
        1152921504583647233,
        1152921504581877761,
        1152921504581419009,
        1152921504580894721,
    ]
    crt_moduli = CRTModuli(qis)
    t = 65537
    n = 1024
    sigma = 3.2
    discrete_gaussian = DiscreteGaussian(sigma)
    rlwe_q = RLWE(n, crt_moduli.q, t, discrete_gaussian)
    rlwe_qis = []
    for qi in crt_moduli.qis:
        rlwe_qis.append(RLWE(n, qi, t, discrete_gaussian))

    bfv_rq = BFV(rlwe_q)
    bfv_rqis = [BFV(rlwe_qi) for rlwe_qi in rlwe_qis]

    # s is the secret key polynomial, sampled from the ternary distribution {-1, 0, 1}
    # e is the error polynomial, sampled from the discrete Gaussian distribution
    # m is the message polynomial, sampled from the plaintext space
    s = bfv_rq.rlwe.SampleFromTernaryDistribution()
    e = bfv_rq.rlwe.SampleFromErrorDistribution()
    m = bfv_rq.rlwe.Rt.sample_polynomial()

    inputs = []

    # perform secret key encryption on m in RNS basis and extract the input for the circuit
    for i in range(len(crt_moduli.qis)):
        # a is the polynomial sampled from the uniform distribution in the RNS basis
        a = bfv_rqis[i].rlwe.Rq.sample_polynomial()
        (ciphertext, k0, k1) = bfv_rqis[i].SecretKeyEncrypt(
            s, a, m, e, crt_moduli.q
        )
        input = {
            "ciphertext": ciphertext,
            "k0": k0,
            "k1": k1,
        }
        inputs.append(input)

    # Sanity check. Get the representation of the ciphertext in the Q basis and perform decryption and assert that the decryption is equal to the message polynomial
    rqi_polynomials_ct0 = []
    rqi_polynomials_ct1 = []
    for input in inputs:
        rqi_polynomials_ct0.append(input["ciphertext"][0])
        rqi_polynomials_ct1.append(input["ciphertext"][1])

    rq_polynomial_ct0 = CRTPolynomial.from_rqi_polynomials_to_rq_polynomial(
            rqi_polynomials_ct0, n, crt_moduli)
    
    rq_polynomial_ct1 = CRTPolynomial.from_rqi_polynomials_to_rq_polynomial(
            rqi_polynomials_ct1, n, crt_moduli)
    
    # Perform decryption on the ciphertext in the Q basis
    dec = bfv_rq.SecretKeyDecrypt(s, (rq_polynomial_ct0, rq_polynomial_ct1), e)

    # Assert that the decryption is equal to the message polynomial
    assert dec == m

    # Each loop simulates a circuit in a different RNS basis
    for i, input in enumerate(inputs):

        # PRECOMPUTATION PHASE - performed outside the circuit

        # The polynomials have now coefficients in the range (-qi/2, qi/2]
        # the circuit works with field elements in Zp, where p is the prime field of the circuit
        # In order to express the polynomials in Zp, we represent the coefficients in the range [0, qi).
        # Every coefficient in the range [0, qi/2] is normally represented in the range [0, qi/2], 
        # Every coefficient in the range (-qi/2, 0) is represented in the range (qi/2, qi) using the `adjust_negative_coefficients` function
        minus_ai = input["ciphertext"][1]
        ai = adjust_negative_coefficients(Polynomial([-1]) * minus_ai, qis[i])
        si = adjust_negative_coefficients(s, qis[i])
        ei = adjust_negative_coefficients(e, qis[i])
        k0i = Polynomial([input["k0"]])
        k1i = adjust_negative_coefficients(input["k1"], qis[i])
        cyclo = [1] + [0] * (n - 1) + [1]
        cyclo = Polynomial(cyclo)

        # ai * si + ei + ki = vi (this is ct0 before reduction)
        vi_part1 = ai * si
        vi_part2 = ei
        vi_part3 = k0i * k1i
        vi = vi_part1 + vi_part2 + vi_part3

        # ai * si + ei + ki = ti mod Rqi
        ti = input["ciphertext"][0]
        ti = adjust_negative_coefficients(ti, qis[i])

        # assert that vi = ti mod Rqi
        vi_clone = copy.deepcopy(vi)
        # mod Rqi means that we need to:
        # - reduce the coefficients of vi_clone by the cyclotomic polynomial
        # - reduce the coefficients of vi_clone by the modulus
        vi_clone.reduce_coefficients_by_cyclo(cyclo.coefficients)
        vi_clone_reduced = [] 
        for coeff in vi_clone.coefficients:
            vi_clone_reduced.append(coeff % qis[i])
        assert Polynomial(vi_clone_reduced) == ti

        # Calculate R2
        # divide ti - vi by the cyclotomic polynomial over Zqi to get R2
        num = ti + (Polynomial([-1]) * vi)
        # reduce the coefficients of num by the modulus qi
        num.reduce_coefficients_by_modulus(qis[i])
        (quotient, rem) = poly_div(num.coefficients, cyclo.coefficients)
        # assert that the remainder is zero
        assert rem == []
        r2 = Polynomial(quotient)

        # Assert that ti - vi = R2 * cyclo mod Zqi
        lhs = ti + (Polynomial([-1]) * vi)
        rhs = r2 * cyclo
        # reduce the coefficients of lhs by the modulus qi
        lhs.reduce_coefficients_by_modulus(qis[i])
        assert lhs == rhs
        
        # Calculate R1
        # divide ti - vi - R2 * cyclo by the modulus qi to get R1
        num = ti + (Polynomial([-1]) * vi) + Polynomial([-1]) * (r2 * cyclo)
        (quotient, rem) = poly_div(num.coefficients, [qis[i]])
        # assert that the remainder is zero
        assert rem == []
        r1 = Polynomial(quotient)

        # Assert that ti = vi + r1 * qi + r2 * cyclo mod Zp
        lhs = ti
        rhs = vi + (r1 * Polynomial([qis[i]])) + (r2 * cyclo)

        # remove the leading zeroes from the rhs 
        for j in range(len(rhs.coefficients)):
            if rhs.coefficients[j] != 0:
                rhs.coefficients = rhs.coefficients[j:]
                break

        assert lhs == rhs

        # PHASE 1
        # Assign si as input to the circuit (private)
        # Assign ei as input to the circuit (private)
        # Assign k1i as input to the circuit (private)
        # Assign P1 as input to the circuit (private)
        # Assign P2 as input to the circuit (private)

        # Commit to phase 1 witness and fetch alpha
        # For experiment, just generate a random alpha
        alpha = randint(0, 100)

        # Precomputation phase 
        # Evaluate ai(alpha), cyclo(alpha), ti(alpha)
        ai_alpha = ai.evaluate(alpha)
        cyclo_alpha = cyclo.evaluate(alpha)
        ti_alpha = ti.evaluate(alpha)

        # PHASE 2
        # Assign ai(alpha) as input to the circuit (public)
        # Assign K0 as input to the circuit (public)
        # Assign cyclo(alpha) as input to the circuit (public)
        # Assign ti(alpha) as input to the circuit (public)
        # Assign qi as input to the circuit (public)

        # Prove that ti - vi - R2 * cyclo = R1 * qi mod Zp by evaluating LHS and RHS at alpha

        # Constrain the evaluation of s(alpha)
        # Constrain the evaluation of e(alpha)
        # Constrain the evaluation of k1i(alpha)
        # Constrain the evaluation of R1(alpha)
        # Constrain the evaluation of R2(alpha)
        s_alpha = si.evaluate(alpha)
        e_alpha = ei.evaluate(alpha)
        k1i_alpha = k1i.evaluate(alpha)
        r1_alpha = r1.evaluate(alpha)
        r2_alpha = r2.evaluate(alpha)

        # Compute vi(alpha) inside the circuit
        # vi(alpha) = ai(alpha) * s(alpha) + e(alpha) + (k1i(alpha) * k0i)
        vi_alpha = ai_alpha * s_alpha + e_alpha + (k1i_alpha * k0i.coefficients[0])

        # sanity check
        assert vi_alpha == vi.evaluate(alpha)

        # Assert that ti = vi + r1 * qi + r2 * cyclo mod Zp
        lhs = ti_alpha
        rhs = vi_alpha + (r1_alpha * qis[i]) + (r2_alpha * cyclo_alpha)
        assert lhs == rhs
    
        # TODO: add range check for private inputs
main()