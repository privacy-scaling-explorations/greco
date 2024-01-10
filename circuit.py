from bfv.crt import Q
from bfv.bfv import BFV
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import poly_mul, poly_add, Polynomial
from utils import mod_inverse


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
    q = Q(qis)
    t = 65537
    n = 1024
    sigma = 3.2
    discrete_gaussian = DiscreteGaussian(sigma)

    bfv = BFV(n, q.q, t, discrete_gaussian)

    # s is the secret key polynomial, sampled from the ternary distribution {-1, 0, 1}
    # e is the error polynomial, sampled from the discrete Gaussian distribution
    # m is the message polynomial, sampled from the plaintext space
    s = bfv.SampleFromTernaryDistribution()
    e = bfv.SampleFromErrorDistribution()
    m = bfv.Rt.sample_polynomial()

    # k1 is equal to m * q mod t
    k1 = poly_mul(m.coefficients, [q.q])
    k1 = Polynomial(k1, bfv.Rt)

    # ciphertexts is a list of ciphertexts, each ciphertext is a tuple (ct0, ct1) generated via secret key encryption in RNS basis
    ciphertexts = []

    # perform secret key encryption on m in RNS basis
    for qi in qis:
        bfv_i = BFV(n, qi, t, discrete_gaussian)

        # a_i is a polynomial sampled from Rqi (Ring of polynomials with coefficients in Zqi modulo the cyclotomic polynomial)
        a_i = bfv_i.Rq.sample_polynomial()

        # k0 is equal to the negative of the multiplicative inverse of t modulo qi
        k0_i = mod_inverse(t, qi) * (-1)

        # Perform encryption : ct0 = a_i * s + e + k0 * k1 mod qi | ct1 = a_i mod qi
        a_i_times_s = poly_mul(a_i.coefficients, s.coefficients)
        a_i_times_s_plus_e = poly_add(a_i_times_s, e.coefficients)
        k0_times_k1_i = poly_mul([k0_i], k1.coefficients)
        ct0_i = poly_add(a_i_times_s_plus_e, k0_times_k1_i)
        ct0_i = Polynomial(ct0_i, bfv_i.Rq)
        ciphertext_i = (ct0_i, a_i)

        ciphertexts.append(ciphertext_i)

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


# Questions
# - Error accumulation when working with `crt`

# RLWE!
# add message encoding!
