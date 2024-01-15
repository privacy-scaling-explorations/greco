from bfv.crt import CRTModuli
from bfv.bfv import BFV, RLWE
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div, get_centered_remainder
from utils import SecretKeyEncrypt


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
        si = secret_key.coefficients
        ei = e.coefficients
        k0i = input["k0"]
        k1i = input["k1"].coefficients
        cyclo = [1] + [0] * (len(ai) - 1) + [1]

        # ai * si + ei + ki = ti (this is ct0 before reduction)
        ti_part1 = Polynomial(ai) * Polynomial(si)
        ti_part2 = Polynomial(ei)
        ti_part3 = Polynomial([k0i]) * Polynomial(k1i)
        ti = ti_part1 + ti_part2 + ti_part3

        # vi = ct0
        # vi = ai * si + ei + ki mod Rqi
        # vi = ai * si + ei + ki - P1(cyclo) mod Zqi
        # vi = ai * si + ei + ki - P1(cyclo) - P2*qi mod Zp

        # assert that vi = t1 mod Rqi
        ti_clone = Polynomial(ti.coefficients)
        # mod Zqi means that we need to:
        # - reduce the coefficients of ti by the cyclotomic polynomial
        # - reduce the coefficients of ti by the modulus
        ti_clone.reduce_coefficients_by_cyclo(cyclo)
        ti_clone.reduce_coefficients_by_modulus(qis[i])
        assert ti_clone.coefficients == input["ciphertext"][0].coefficients

        # Calculate P1
        # (ti mod Zqi) / cyclo = (quotient, remainder) where quotient = P1 and the remainder is equal to vi
        ti_mod_zqi = Polynomial(ti.coefficients)
        ti_mod_zqi.reduce_coefficients_by_modulus(qis[i])
        (quotient, rem) = poly_div(ti_mod_zqi.coefficients, cyclo)
        # Get the centered remainder representation from each coefficient of rem
        rem = [get_centered_remainder(x, qis[i]) for x in rem]
        assert rem == input["ciphertext"][0].coefficients

        # PHASE 1
        # Assign Ai as input to the circuit (public)
        # Assign Si as input to the circuit (private)
        # Assign Ei as input to the circuit (private)
        # Assign K0i as input to the circuit (public)
        # Assign K1i as input to the circuit (public)
        # Assign cyclo as input to the circuit (public)

        # Pass P1 as input to the circuit (private)
        # Pass P2 as input to the circuit (private)
        # Assert that coefficients of the matrix are in the expected range
        # Commit to phase 1 and fetch alpha


# PHASE 2
# Fetch alpha
# Pass A(alpha) as input to the circuit (public)
# Pass Ki^0(alpha) as input to the circuit (public)
# Pass cyclo(alpha) as input to the circuit (public)
# Pass Ti(alpha) as input to the circuit (public)
# Assert that U(IV) = Ti(alpha)


main()
