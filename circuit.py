from bfv.crt import CRTModuli
from bfv.bfv import BFV, RLWE
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div, get_centered_remainder
from utils import SecretKeyEncrypt
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
            "a": a,
        }
        inputs.append(input)

    # Each loop simulates a circuit in RNS basis
    for i, input in enumerate(inputs):
        # Precomputation phase
        ai = input["a"].coefficients
        ai = [int(x) for x in ai]
        si = secret_key.coefficients
        ei = e.coefficients
        k0i = int(input["k0"])
        k1i = input["k1"].coefficients
        k1i = [int(x) for x in k1i]
        cyclo = [1] + [0] * (len(ai) - 1) + [1]

        # ai * si + ei + ki = ti (this is ct0 before reduction)
        ti_part1 = Polynomial(ai) * Polynomial(si)
        ti_part2 = Polynomial(ei)
        ti_part3 = Polynomial([k0i]) * Polynomial(k1i)
        ti = ti_part1 + ti_part2 + ti_part3

        # vi = ct0
        # vi = ai * si + ei + ki mod Rqi
        # vi = ai * si + ei + ki - P1 * cyclo mod Zqi
        # vi = ai * si + ei + ki - P1 * cyclo - P2*qi mod Zp

        # assert that vi = ti mod Rqi
        ti_clone = Polynomial(copy.deepcopy(ti.coefficients))
        # mod Zqi means that we need to:
        # - reduce the coefficients of ti by the cyclotomic polynomial
        # - reduce the coefficients of ti by the modulus
        ti_clone.reduce_coefficients_by_cyclo(cyclo)
        ti_clone.reduce_coefficients_by_modulus(qis[i])
        assert ti_clone.coefficients == input["ciphertext"][0].coefficients

        # Calculate P1
        # (ti mod Zqi) / cyclo = (quotient, remainder) where quotient = P1 and the remainder is equal to vi
        ti_mod_zqi = Polynomial(copy.deepcopy(ti.coefficients))
        ti_mod_zqi.reduce_coefficients_by_modulus(qis[i])
        (quotient, rem) = poly_div(ti_mod_zqi.coefficients, cyclo)
        # Get the centered remainder representation from each coefficient of rem
        rem = [get_centered_remainder(x, qis[i]) for x in rem]
        # Assert that the remainder is equal to vi
        assert rem == input["ciphertext"][0].coefficients
        # Get the centered remainder representation from each coefficient of quotient
        quotient = [get_centered_remainder(x, qis[i]) for x in quotient]
        p1 = Polynomial(quotient)

        # Calculate P2
        # reducing each coefficient of poly `a` by the modulus `qi` can be represented as: a = qi * b + r
        # where b is the quotient and r is the remainder
        # ti - P1 * cyclo / qi = (quotient, remainder) where quotient = P2 and the remainder is equal to vi
        p1_times_cyclo = Polynomial(p1.coefficients) * Polynomial(cyclo)
        minus_p1_times_cyclo = Polynomial([-1]) * p1_times_cyclo
        num = Polynomial(ti.coefficients) + minus_p1_times_cyclo
        (quotient, rem) = poly_div(num.coefficients, [qis[i]])
        # Get the centered remainder representation from each coefficient of rem
        rem = [get_centered_remainder(x, qis[i]) for x in rem]
        # Assert that the remainder is equal to vi
        assert rem == input["ciphertext"][0].coefficients
        # Get the centered remainder representation from each coefficient of quotient
        quotient = [get_centered_remainder(x, qis[i]) for x in quotient]
        p2 = Polynomial(quotient)
        # Commit to phase 1 witness and fetch alpha
        # For experiment, just generate a random alpha
        alpha = randint(0, 100)

        # vi + P1 * cyclo + P2 * qi = ti (note: this is failing atm)
        lhs = Polynomial(input["ciphertext"][0].coefficients) + p1_times_cyclo + Polynomial(p2.coefficients) * Polynomial([qis[i]])

        # PHASE 1
        # Assign si as input to the circuit (private)
        # Assign ei as input to the circuit (private)
        # Assign k1i as input to the circuit (private)
        # Assign P1 as input to the circuit (private)
        # Assign P2 as input to the circuit (private)

        # Precomputation phase 
        # Evaluate ai(alpha), cyclo(alpha), vi(alpha)
        ai_alpha = Polynomial(ai).evaluate(alpha)
        cyclo_alpha = Polynomial(cyclo).evaluate(alpha)
        vi_alpha = Polynomial(input["ciphertext"][0].coefficients).evaluate(alpha)

        # PHASE 2
        # Assign ai(alpha) as input to the circuit (public)
        # Assign K0i as input to the circuit (public)
        # Assign cyclo(alpha) as input to the circuit (public)
        # Assign vi(alpha) as input to the circuit (public)

        # Prove that ti - P1 * cyclo - P2 * qi = vi mod Zp by evaluating LHS and RHS at alpha
        # vi(alpha), cyclo(alpha), qi are public inputs. Also ai(alpha), which will be necessary to compute ti(alpha) is a public input

        # Constrain the evaluation of s(alpha)
        # Constrain the evaluation of e(alpha)
        # Constrain the evaluation of k1i(alpha)
        # Constrain the evaluation of P1(alpha)
        # Constrain the evaluation of P2(alpha)
        s_alpha = Polynomial(si).evaluate(alpha)
        e_alpha = Polynomial(e.coefficients).evaluate(alpha)
        k1i_alpha = Polynomial(k1i).evaluate(alpha)
        p1_alpha = Polynomial(p1.coefficients).evaluate(alpha)
        p2_alpha = Polynomial(p2.coefficients).evaluate(alpha)

        # Compute ti(alpha) inside the circuit
        # ti(alpha) = ai(alpha) * s(alpha) + e(alpha) + k1i(alpha)
        ti_alpha = ai_alpha * s_alpha + e_alpha + k1i_alpha

        # sanity check
        assert ti_alpha == ti.evaluate(alpha)
        
        # Compute P1(alpha) * cyclo(alpha) inside the circuit
        pi_alpha_times_cyclo_alpha = p1_alpha * cyclo_alpha

        # Compute P2(alpha) * qi inside the circuit
        pi_alpha_times_qi = p2_alpha * qis[i]

        # Assert that vi(alpha) + P1(alpha) * cyclo(alpha) + P2(alpha) * qi(alpha) = ti(alpha)
        # Assert that coefficients of the matrix are in the expected range

main()
