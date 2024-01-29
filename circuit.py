from bfv.crt import CRTModuli, CRTPolynomial
from bfv.bfv import BFV, RLWE
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div
from random import randint
import copy

def main(): 
    
    '''
    ENCRYPTION PHASE - performed outside the circuit

    - `qis` are the small moduli used in the CRT basis, which product is the ciphertext modulus q
    - `t` is the plaintext modulus
    - `n` is the degree of the cyclotomic polynomial which is the denominator of the polynomial ring
    - `sigma` is the standard deviation of the discrete Gaussian distribution
    - `s` is the secret key polynomial, sampled from the ternary distribution {-1, 0, 1}
    - `e` is the error polynomial, sampled from the discrete Gaussian distribution
    - `m` is the message polynomial, sampled from the plaintext space
    '''

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

    # Create a BFV instance for the ciphertext modulus q and a list of BFV instances for each small moduli qis
    bfv_rq = BFV(rlwe_q)
    bfv_rqis = [BFV(rlwe_qi) for rlwe_qi in rlwe_qis]

    s = bfv_rq.rlwe.SampleFromTernaryDistribution()
    e = bfv_rq.rlwe.SampleFromErrorDistribution()
    m = bfv_rq.rlwe.Rt.sample_polynomial()

    inputs = []

    # Perform encryption of m in each CRT basis
    for i in range(len(crt_moduli.qis)):

        # a is a polynomial sampled from the uniform distribution Rqi
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

    # Sanity check. Get the representation of the ciphertext in the Q basis and perform decryption and check if the decryption is equal to `m`
    rqi_ct0_polynomials = []
    rqi_ct1_polynomial = []
    for input in inputs:
        rqi_ct0_polynomials.append(input["ciphertext"][0])
        rqi_ct1_polynomial.append(input["ciphertext"][1])

    rq_ct0_polynomial = CRTPolynomial.from_rqi_polynomials_to_rq_polynomial(
            rqi_ct0_polynomials, n, crt_moduli)
    
    rq_ct1_polynomial = CRTPolynomial.from_rqi_polynomials_to_rq_polynomial(
            rqi_ct1_polynomial, n, crt_moduli)
    
    rq_ciphertext_recovered = (rq_ct0_polynomial, rq_ct1_polynomial)
    
    # Perform decryption on the ciphertext in the Q basis
    dec = bfv_rq.SecretKeyDecrypt(s, rq_ciphertext_recovered, e)

    assert dec == m

    # Each round of the loop simulates a proof generation phase for a different CRT basis
    for i, input in enumerate(inputs):

        '''
        SETUP PHASE - performed outside the circuit

        For each CRT basis, we need to compute the polynomials R1 and R2 (check this doc for more details: https://hackmd.io/@gaussian/HJ8DYyjPp)
        '''

        ai = Polynomial([-1]) * input["ciphertext"][1]
        si = s
        ei = e
        k0i = Polynomial([input["k0"]])
        k1i = input["k1"]
        cyclo = [1] + [0] * (n - 1) + [1]
        cyclo = Polynomial(cyclo)

        # ai * si + ei + ki = vi (this is ct0 before reduction in the Rqi ring)
        vi = ai * si + ei + k0i * k1i

        # ai * si + ei + ki = ti mod Rqi = ct0
        ti = input["ciphertext"][0]

        # assert that vi = ti mod Rqi
        vi_clone = copy.deepcopy(vi)
        # mod Rqi means that we need to:
        # - reduce the coefficients of vi_clone by the cyclotomic polynomial
        # - reduce the coefficients of vi_clone by the modulus
        vi_clone.reduce_coefficients_by_cyclo(cyclo.coefficients)
        vi_clone.reduce_coefficients_by_modulus(qis[i])
        assert vi_clone == ti

        # Calculate R2
        # divide ti - vi by the cyclotomic polynomial over Zqi to get R2
        num = ti + (Polynomial([-1]) * vi)
        # reduce the coefficients of num by the modulus qi 
        num.reduce_coefficients_by_modulus(qis[i])
        (quotient, rem) = poly_div(num.coefficients, cyclo.coefficients)
        # assert that the remainder is zero
        assert rem == []
        r2 = Polynomial(quotient)
        # assert that the degree of R2 is equal to n - 2
        assert len(r2.coefficients) - 1 == n - 2

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

        '''
        CIRCUIT - PHASE 1 - ASSIGNMENT PHASE

        In this phase, the private inputs are assigned to the circuit. These are the polynomials si, ei, k1i, R1, R2. 
        We also assign the public inputs qi, t, k0i and b. N is a constant of the circuit.
        '''

        # ... Perform assignment here and expose public inputs...

        '''
        CIRCUIT - PHASE 1 - RANGE CHECK OF PRIVATE POLYNOMIALS 

        In this phase, the coefficients of the private polynomials are checked to be in the correct range.
        '''

        # ... Perform range check here ...

        # sanity check. The coefficients of ti should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = ((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in ti.coefficients)

        # sanity check. The coefficients of ai should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = ((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in ai.coefficients)

        # constraint. The coefficients of si should be in the range [-1, 0, 1]
        assert all(coeff in [-1, 0, 1] for coeff in si.coefficients)

        # sanity check. The coefficients of ai * si should be in the range $[-N \cdot \frac{q_i - 1}{2}, N \cdot \frac{q_i - 1}{2}]$
        bound = ((qis[i] - 1) / 2) * n
        res = ai * si
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of ei should be in the range [-B, B] where B is the upper bound of the discrete Gaussian distribution
        b = int(discrete_gaussian.z_upper)
        assert all(coeff >= -b and coeff <= b for coeff in ei.coefficients)

        # sanity check. The coefficients of ai * si + ei should be in the range $- (N \cdot \frac{q_i - 1}{2} + B), N \cdot \frac{q_i - 1}{2} + B]$
        bound = ((qis[i] - 1) / 2) * n + b
        res = ai * si + ei
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of R2 should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = ((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in r2.coefficients)

        # sanity check. The coefficients of R2 * cyclo should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = ((qis[i] - 1) / 2)
        res = r2 * cyclo
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of k1i should be in the range [-(t-1)/2, (t-1)/2]
        bound = ((t - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in k1i.coefficients)

        # sanity check. The coefficients of k1i * k0i should be in the range $[-\frac{t - 1}{2} \cdot |K_i^{0}|, \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = ((t - 1) / 2) * abs(k0i.coefficients[0])
        res = k1i * k0i
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # sanity check. The coefficients of vi (ai * si + ei + k1i * k0i) should be in the range $[- (N \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|), N \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = ((qis[i] - 1) / 2) * n + b + ((t - 1) / 2) * abs(k0i.coefficients[0])
        assert all(coeff >= -bound and coeff <= bound for coeff in vi.coefficients)

        # sanity check. The coefficients of ti - vi should be in the range $ [- ((N+1) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|), (N+1) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = ((qis[i] - 1) / 2) * (n + 1) + b + ((t - 1) / 2) * abs(k0i.coefficients[0])
        sub = ti + (Polynomial([-1]) * vi)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # sanity check. The coefficients of ti - vi - R2 * cyclo should be in the range $[- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|), (N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = ((qis[i] - 1) / 2) * (n + 2) + b + ((t - 1) / 2) * abs(k0i.coefficients[0])
        sub = ti + (Polynomial([-1]) * vi) + Polynomial([-1]) * (r2 * cyclo)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # constraint. The coefficients of (ti - vi - R2 * cyclo) / qi = R1 should be in the range $[\frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|)}{q_i}, \frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|}{q_i}]$
        bound = (((qis[i] - 1) / 2) * (n + 2) + b + ((t - 1) / 2) * abs(k0i.coefficients[0])) / qis[i]
        assert all(coeff >= -bound and coeff <= bound for coeff in r1.coefficients)

        '''
        CIRCUIT - END OF PHASE 1 - WITNESS COMMITMENT 

        At the end of phase 1, the witness is committed and a challenge is generated using the Fiat-Shamir heuristic.
        '''

        # For the sake of simplicity, we generate a random challenge here
        alpha = randint(0, 1000)

        '''
        CIRCUIT - PHASE 2 - ASSIGNMENT PHASE
        The public inputs are assigned to the circuit. These are the polynomials ai, cyclo, ti evaluated at alpha.
        '''

        # The evaluation of ai_alpha, cyclo_alpha, ti_alpha is performed outside the circuit
        ai_alpha = ai.evaluate(alpha)
        cyclo_alpha = cyclo.evaluate(alpha)
        ti_alpha = ti.evaluate(alpha)

        # ... Perform assignment here and expose expose public inputs ...

        '''
        CIRCUIT - PHASE 2 - CORRECT ENCRIPTION CONSTRAINT

        We need to prove that ti = vi + r1 * qi + r2 * cyclo mod Zp.
        We do that by proving that LHS(alpha) = RHS(alpha) for a random alpha according to Scwhartz-Zippel lemma.
        '''

        s_alpha = si.evaluate(alpha)
        e_alpha = ei.evaluate(alpha)
        k1i_alpha = k1i.evaluate(alpha)
        r1_alpha = r1.evaluate(alpha)
        r2_alpha = r2.evaluate(alpha)

        lhs = ti_alpha
        rhs = ai_alpha * s_alpha + e_alpha + (k1i_alpha * k0i.coefficients[0]) + (r1_alpha * qis[i]) + (r2_alpha * cyclo_alpha)
        assert lhs == rhs

main()