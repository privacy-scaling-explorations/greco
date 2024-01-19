from bfv.crt import CRTModuli
from bfv.bfv import BFV, RLWE
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div
from utils import SecretKeyEncrypt, adjust_negative_coefficients
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

    rlwe = RLWE(n, crt_moduli.q, t, discrete_gaussian)
    bfv = BFV(rlwe)

    # s is the secret key polynomial, sampled from the ternary distribution {-1, 0, 1}
    # e is the error polynomial, sampled from the discrete Gaussian distribution
    # m is the message polynomial, sampled from the plaintext space
    s = bfv.rlwe.SampleFromTernaryDistribution()
    e = bfv.rlwe.SampleFromErrorDistribution()
    m = bfv.rlwe.Rt.sample_polynomial()

    inputs = []

    # perform secret key encryption on m in RNS basis and extract the input for the circuit
    for qi in qis:
        bfv_rqi = BFV(RLWE(n, qi, t, discrete_gaussian))
        # a is the polynomial sampled from the uniform distribution in the RNS basis
        # TODO: add test for this
        a = bfv_rqi.rlwe.Rq.sample_polynomial()
        (ciphertext, k0, k1) = SecretKeyEncrypt(
            a, s, e, m, t, crt_moduli.q, qi
        )
        input = {
            "ciphertext": ciphertext,
            "k0": k0,
            "k1": k1,
        }
        inputs.append(input)

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

        # ti = ct0
        # ti = vi mod Rqi
        # ti = vi - P1 * cyclo mod Zqi
        # ti = vi - P1 * cyclo - P2*qi mod Zp
        ti = input["ciphertext"][0]
        ti = adjust_negative_coefficients(ti, qis[i])

        # assert that ti = vi mod Rqi
        vi_clone = copy.deepcopy(vi)
        # mod Rqi means that we need to:
        # - reduce the coefficients of vi_clone by the cyclotomic polynomial
        # - reduce the coefficients of vi_clone by the modulus
        vi_clone.reduce_coefficients_by_cyclo(cyclo.coefficients)
        vi_clone_reduced = [] 
        for coeff in vi_clone.coefficients:
            vi_clone_reduced.append(coeff % qis[i])
        assert Polynomial(vi_clone_reduced) == ti

        # Calculate P1 
        # ti = vi - P1 * cyclo mod Zqi
        # P1 * cyclo + ti = vi mod Zqi
        # vi mod Zqi / cyclo = (quotient, remainder) where quotient = P1 and the remainder = ti  

        # reduce the coefficients of vi by the modulus qi
        vi_mod_zqi = []
        for coeff in vi.coefficients:
            vi_mod_zqi.append(coeff % qis[i])     
        (quotient, rem_1) = poly_div(vi_mod_zqi, cyclo.coefficients)
        p1 = Polynomial(quotient)
        # adjust the coefficient of the remainder to be inside the field Zqi
        rem_adj = adjust_negative_coefficients(Polynomial(rem_1), qis[i])
        # assert that the rem_adj is equal to ti
        assert rem_adj == ti

        # assert that vi_mod_zqi = P1 * cyclo + rem_1
        lhs = vi_mod_zqi
        rhs = p1 * cyclo + Polynomial(rem_1)
        assert Polynomial(lhs) == rhs

        # assert that vi_mod_zqi = P1 * cyclo + ti mod Zqi
        # We need to further reduce the coefficients of RHS by the modulus qi because, after the operation `adjust_negative_coefficients` operated on ti, the coefficients of rhs are no longer in the range [0, qi)
        lhs = vi_mod_zqi    
        rhs = p1 * cyclo + ti
        rhs_reduced = []
        for coeff in rhs.coefficients:
            rhs_reduced.append(coeff % qis[i])
        assert lhs == rhs_reduced        

        # Calculate P2
        # ti = vi - P1 * cyclo - P2*qi mod Zp
        # P2*qi + P1 * cyclo + ti = vi mod Zp
        # The operation reducing each coefficient of poly `a` by the modulus `qi` can be represented as: a = b * qi + r
        # where b is the quotient and r is the remainder
        # in this case, we reduce the coefficients of vi mod Zp by the modulus `qi`, the quotient is P2 and the remainder is ti + P1 * cyclo

        # vi mod Zp is equal to vi according to the assumptions of the circuits
        (quotient, rem_2) = poly_div(vi.coefficients, [qis[i]])
        p2 = Polynomial(quotient)

        # assert that rem_2 = rem_1 + P1 * cyclo
        lhs = Polynomial(rem_2)
        rhs = Polynomial(rem_1) + p1 * cyclo
        # remove leading zeros from lhs 
        while lhs.coefficients[0] == 0:
            lhs.coefficients.pop(0)
        # remove leading zeros from rhs
        while rhs.coefficients[0] == 0:
            rhs.coefficients.pop(0)
        assert lhs == rhs

        # assert that rem_2 = ti + P1 * cyclo mod Zqi
        # We need to further reduce the coefficients of RHS by the modulus qi because, after the operation of `adjust_negative_coefficients` operated on ti, the coefficients of rhs are no longer in the range [0, qi)
        lhs = Polynomial(rem_2)
        rhs = ti + p1 * cyclo
        rhs_reduced = []
        for coeff in rhs.coefficients:
            rhs_reduced.append(coeff % qis[i])
        # remove leading zeros from lhs
        while lhs.coefficients[0] == 0:
            lhs.coefficients.pop(0)
        # remove leading zeros from rhs
        while rhs_reduced[0] == 0:
            rhs_reduced.pop(0)
        assert lhs == Polynomial(rhs_reduced)

        # assert that vi mod Zp = P2 * qi + rem_2
        lhs = vi
        rhs = p2 * Polynomial([qis[i]]) + Polynomial(rem_2)
        assert lhs == rhs

        # assert that vi mod Zp = P2 * qi + [ti + P1 * cyclo mod Zqi]
        # We need to further reduce the coefficients of [ti + P1 * cyclo mod Zqi] by the modulus qi because, after the operation of `adjust_negative_coefficients` operated on ti, the coefficients of[ti + P1 * cyclo] are no longer the range [0, qi)
        lhs = vi
        rhs = ti + p1 * cyclo
        rhs_reduced = []
        for coeff in rhs.coefficients:
            rhs_reduced.append(coeff % qis[i])
        rhs_reduced = Polynomial(rhs_reduced) + p2 * Polynomial([qis[i]])
        assert lhs == rhs_reduced

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

        # Prove that vi = P2 * qi + P1 * cyclo + ti mod Zp by evaluating LHS and RHS at alpha
        # Prove that ai * si + ei + k1i * k0 = P2 * qi + P1 * cyclo + ti mod Zp 

        # Constrain the evaluation of s(alpha)
        # Constrain the evaluation of e(alpha)
        # Constrain the evaluation of k1i(alpha)
        # Constrain the evaluation of P1(alpha)
        # Constrain the evaluation of P2(alpha)
        s_alpha = si.evaluate(alpha)
        e_alpha = ei.evaluate(alpha)
        k1i_alpha = k1i.evaluate(alpha)
        p1_alpha = p1.evaluate(alpha)
        p2_alpha = p2.evaluate(alpha)

        # Compute vi(alpha) inside the circuit
        # vi(alpha) = ai(alpha) * s(alpha) + e(alpha) + (k1i(alpha) * k0i)
        vi_alpha = ai_alpha * s_alpha + e_alpha + (k1i_alpha * k0i.coefficients[0])
        
        # Compute P1(alpha) * cyclo(alpha) + ti_alpha
        p1_alpha_times_cyclo_alpha_plus_ti_alpha = p1_alpha * cyclo_alpha + ti_alpha

        # Compute P2(alpha) * qi inside the circuit
        p2_alpha_times_qi = p2_alpha * qis[i]

        # Assert that vi(alpha) = P2(alpha) * qi + ti(alpha) + P1(alpha) * cyclo(alpha)
        # TODO: explain why we need to reduce the coefficients of the RHS by the modulus qi
        lhs = vi_alpha % qis[i]
        rhs = (p2_alpha_times_qi + p1_alpha_times_cyclo_alpha_plus_ti_alpha) % qis[i]
        
        assert lhs == rhs

        # TODO: add range check for private inputs
main()