from bfv.crt import CRTModuli
from bfv.bfv import BFV, RLWE
from bfv.discrete_gauss import DiscreteGaussian


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

    # ciphertexts is a list of ciphertexts, each ciphertext is a tuple (ct0, ct1) generated via secret key encryption in RNS basis
    ciphertexts = []

    # perform secret key encryption on m in RNS basis
    for qi in qis:
        bfv_rqi = BFV(RLWE(n, qi, t, discrete_gaussian))
        ciphertext = bfv_rqi.SecretKeyEncrypt(secret_key, message, e, crt_moduli.q)
        ciphertexts.append(ciphertext)

    # PHASE 1
    # Pass Si as input to the circuit (private)
    # Pass Ei as input to the circuit (private)
    # Pass Ki^1 as input to the circuit (private)
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
