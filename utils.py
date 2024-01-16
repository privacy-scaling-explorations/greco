from bfv.polynomial import Polynomial, PolynomialRing


def SecretKeyEncrypt(
    a: Polynomial,
    s: Polynomial,
    e: Polynomial,
    m: Polynomial,
    t: int,
    q: int,
    qi: int,
) -> tuple[tuple[Polynomial, Polynomial], int, Polynomial]:
    """
    Encrypt a given message m with a given secret key s.

    Parameters:
    - a: polynomial sampled from the distribution χ Uniform.
    - s: secret key polynomial sampled from the distribution χ Ternary.
    - e: error polynomial sampled from the distribution χ Error.
    - m: message polynomial. This must be a polynomial in Rt.
    - t: plaintext modulus
    - q: ciphertext modulus
    - qi: RNS modulus

    Returns:
    - ciphertext: Generated ciphertext.
    - k0 : equals to -t^{-1}
    - k1 : equals to [QM]t
    """

    # Compute the ciphertext.
    # k^{0} = -t^{-1} namely the multiplicative inverse of t modulo qi
    k0 = mod_inverse(t, qi) * (-1)

    # k^{1} = [QM]t namely the scaled message polynomial
    k1 = Polynomial([q]) * m

    # reduce k^{1} in Rt
    rt = PolynomialRing(len(a.coefficients), qi)
    k1.reduce_in_ring(rt)

    # a * s
    mul = a * s

    # b = a*s + e.
    b = mul + e

    # ct_0 = b + k^{0}k^{1}
    ct_0 = b + (Polynomial([k0]) * k1)

    # ct_0 will be in Rqi
    rqi = PolynomialRing(len(a.coefficients), qi)
    ct_0.reduce_in_ring(rqi)

    # ct_1 = -a
    ct_1 = a * Polynomial([-1])

    ciphertext = (ct_0, ct_1)

    return (ciphertext, k0, k1)


def extended_gcd(a, b):
    """
    Computes the greatest common divisor of a and b.
    Returns a tuple (g, x, y) such that a*x + b*y = g = gcd(a, b).
    """
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = extended_gcd(b % a, a)
        return (g, x - (b // a) * y, y)


def mod_inverse(t, q):
    """
    Computes the multiplicative inverse of t modulo q.
    Returns the inverse, or raises an exception if it doesn't exist.
    """
    g, x, _ = extended_gcd(t, q)
    if g != 1:
        raise ValueError("The multiplicative inverse does not exist")
    else:
        return x % q

def adjust_negative_coefficients(poly : Polynomial, modulus: int) -> Polynomial:
    """
    Adjust the coefficients of the polynomial to be positive.
    """
    return Polynomial([(modulus + coeff if coeff < 0 else coeff) for coeff in poly.coefficients])
