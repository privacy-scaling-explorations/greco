from bfv.crt import CRTModuli
from bfv.bfv import BFV, RLWE
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div
from utils import SecretKeyEncrypt, adjust_negative_coefficients
from random import randint
import copy


def main():
    # Setup
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
    secret_key = bfv.rlwe.SampleFromTernaryDistribution()
    e = bfv.rlwe.SampleFromErrorDistribution()
    message = bfv.rlwe.Rt.sample_polynomial()

    inputs = []

    # perform secret key encryption on m in RNS basis and extract the input for the circuit
    for qi in qis:
        bfv_rqi = BFV(RLWE(n, qi, t, discrete_gaussian))
        a = bfv_rqi.rlwe.Rq.sample_polynomial()
        (ciphertext, k0, k1) = SecretKeyEncrypt(
            a, secret_key, e, message, t, crt_moduli.q, qi
        )
        input = {
            "ciphertext": ciphertext,
            "k0": k0,
            "k1": k1,
        }
        inputs.append(input)

    # Each loop simulates a circuit in a different RNS basis
    for i, input in enumerate(inputs):

        # Precomputation phase - performed outside the circuit
        # TODO : add explanation about negative coefficients
        ai = Polynomial(input["ciphertext"][1].coefficients)
        ai = Polynomial([-1]) * ai
        ai = adjust_negative_coefficients(ai, qis[i])
        si = adjust_negative_coefficients(secret_key, qis[i])
        ei = adjust_negative_coefficients(e, qis[i])
        k0i = adjust_negative_coefficients(Polynomial([input["k0"]]), qis[i])
        k1i = adjust_negative_coefficients(input["k1"], qis[i])
        cyclo = [1] + [0] * (n - 1) + [1]
        cyclo = Polynomial(cyclo)

        # ai * si + ei + ki = vi (this is ct0 before reduction)
        vi_part1 = ai * si
        vi_part2 = ei
        vi_part3 = k0i * k1i
        vi = vi_part1 + vi_part2 + vi_part3

        # assert that vi_part1 does not have negative coefficients
        for coeff in vi_part1.coefficients:
            assert coeff >= 0

        # assert that vi_part2 does not have negative coefficients
        for coeff in vi_part2.coefficients:
            assert coeff >= 0

        # assert that vi_part3 does not have negative coefficients
        for coeff in vi_part3.coefficients:
            assert coeff >= 0

        # assert that vi does not have negative coefficients
        for coeff in vi.coefficients:
            assert coeff >= 0

        # ti = ct0
        # ti = ai * si + ei + ki mod Rqi
        # ti = ai * si + ei + ki - P1 * cyclo mod Zqi
        # ti = ai * si + ei + ki - P1 * cyclo - P2*qi mod Zp
        ti = Polynomial(input["ciphertext"][0].coefficients)
        ti = adjust_negative_coefficients(ti, qis[i])

        # assert that ti = vi mod Rqi
        vi_clone = copy.deepcopy(vi)
        # mod Rqi means that we need to:
        # - reduce the coefficients of vi by the cyclotomic polynomial
        # - reduce the coefficients of vi by the modulus
        vi_clone.reduce_coefficients_by_cyclo(cyclo.coefficients)
        vi_clone_reduced = [] 
        for coeff in vi_clone.coefficients:
            vi_clone_reduced.append(coeff % qis[i])
        assert Polynomial(vi_clone_reduced) == ti

        # Calculate P1
        # (vi mod Zqi) / cyclo = (quotient, remainder) where quotient = P1 and the remainder is equal to ti        
        # reduce the coefficients of vi by the modulus
        vi_mod_zqi = []
        for coeff in vi.coefficients:
            vi_mod_zqi.append(coeff % qis[i])
        vi_mod_zqi = Polynomial(vi_mod_zqi)

        # assert that vi_mod_zqi does not have negative coefficients
        for coeff in vi_mod_zqi.coefficients:
            assert coeff >= 0

        (quotient, rem) = poly_div(vi_mod_zqi.coefficients, cyclo.coefficients)
        rem_adj = adjust_negative_coefficients(Polynomial(rem), qis[i])
        p1 = Polynomial(quotient)
        # assert that the remainder is equal to ti
        assert rem_adj == ti

        # assert that vi_mod_zqi = P1 * cyclo + rem
        lhs = vi_mod_zqi
        rhs = p1 * cyclo + Polynomial(rem)
        assert lhs == rhs

        # assert that P1 does not have negative coefficients
        for coeff in p1.coefficients:
            assert coeff >= 0

        # Calculate P2
        # reducing each coefficient of poly `a` by the modulus `qi` can be represented as: a = b * qi + r
        # where b is the quotient and r is the remainder
        # vi - P1 * cyclo / qi = (quotient, remainder) where quotient = P2 and the remainder is equal to ti
        p1_times_cyclo = p1 * cyclo
        minus_p1_times_cyclo = Polynomial([-1]) * p1_times_cyclo
        minus_p1_times_cyclo = adjust_negative_coefficients(minus_p1_times_cyclo, qis[i])
        num = vi + minus_p1_times_cyclo
        (quotient, rem) = poly_div(num.coefficients, [qis[i]])
        p2 = Polynomial(quotient)
        # assert that P2 does not have negative coefficients
        for coeff in p2.coefficients:
            assert coeff >= 0
        rem_adj = adjust_negative_coefficients(Polynomial(rem), qis[i])
        assert rem_adj == ti

        # assert that vi - P1 * cyclo = P2 * qi + rem_adj
        lhs = vi + minus_p1_times_cyclo
        rhs = p2 * Polynomial([qis[i]]) + rem_adj
        assert lhs == rhs

        # Commit to phase 1 witness and fetch alpha
        # For experiment, just generate a random alpha
        alpha = randint(0, 100)

        # PHASE 1
        # Assign si as input to the circuit (private)
        # Assign ei as input to the circuit (private)
        # Assign k1i as input to the circuit (private)
        # Assign P1 as input to the circuit (private)
        # Assign P2 as input to the circuit (private)

        # Precomputation phase 
        # Evaluate ai(alpha), cyclo(alpha), ti(alpha)
        ai_alpha = ai.evaluate(alpha)
        cyclo_alpha = cyclo.evaluate(alpha)
        ti_alpha = ti.evaluate(alpha)

        # PHASE 2
        # Assign ai(alpha) as input to the circuit (public)
        # Assign K0i as input to the circuit (public)
        # Assign cyclo(alpha) as input to the circuit (public)
        # Assign ti(alpha) as input to the circuit (public)

        # Prove that vi - P1 * cyclo = P2 * qi + ti mod Zp by evaluating LHS and RHS at alpha
        # ti(alpha), cyclo(alpha), qi are public inputs. Also ai(alpha), which will be necessary to compute ti(alpha) is a public input

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
        # vi(alpha) = ai(alpha) * s(alpha) + e(alpha) + k1i(alpha)
        vi_alpha = ai_alpha * s_alpha + e_alpha + k1i_alpha

        # sanity check
        assert vi_alpha == vi.evaluate(alpha)
        
        # Compute P1(alpha) * cyclo(alpha) inside the circuit
        # TODO
        minus_p1_times_cyclo_alpha = minus_p1_times_cyclo.evaluate(alpha)

        # Compute P2(alpha) * qi inside the circuit
        p2_alpha_times_qi = p2_alpha * qis[i]

        # Assert that vi(alpha) - P1(alpha) * cyclo(alpha) = P2(alpha) * qi + ti(alpha)
        lhs = vi_alpha + minus_p1_times_cyclo_alpha
        rhs = p2_alpha_times_qi + ti_alpha
        assert lhs == rhs

main()
