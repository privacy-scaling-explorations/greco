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
