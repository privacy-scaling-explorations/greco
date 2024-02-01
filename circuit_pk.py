from bfv.crt import CRTModuli
from bfv.bfv import BFVCrt
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div
from bfv.utils import mod_inverse
from random import randint
import copy

def main(): 
    
    '''
    ENCRYPTION PHASE - performed outside the circuit
    In this phase, we encryption operation is performed. Later, the circuit will be used to prove that the encryption is correct.
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
    bfv_crt = BFVCrt(crt_moduli, n, t, discrete_gaussian)

    # Perform encryption of m in each CRT basis
    s = bfv_crt.SecretKeyGen()
    e = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    ais = []
    for i in range(len(qis)):
        ai = bfv_crt.bfv_q.rlwe.SampleFromTernaryDistribution()
        ais.append(ai)
    pub_keys = bfv_crt.PublicKeyGen(s, e, ais)
    m = bfv_crt.bfv_q.rlwe.Rt.sample_polynomial()
    e0 = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    e1 = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    u = bfv_crt.bfv_q.rlwe.SampleFromTernaryDistribution()

    ciphertexts = bfv_crt.PubKeyEncrypt(pub_keys, m, e0, e1, u)

    # Sanity check for valid decryption    
    message_prime = bfv_crt.Decrypt(s, ciphertexts)

    assert m == message_prime

    # k1 = [QM]t namely the scaled message polynomial
    k1 = Polynomial([crt_moduli.q]) * m
    k1.reduce_coefficients_by_modulus(t)

    # Each round of the loop simulates a proof generation phase for a different CRT basis
    for i, ciphertext in enumerate(ciphertexts):

        '''
        SETUP PHASE - performed outside the circuit

        For each CRT basis, we need to compute the polynomials R1, R2, P1 and P2 (check this doc for more details: https://hackmd.io/@gaussian/HJ8DYyjPp)
        '''

        pk0i = pub_keys[i][0]
        pk1i = pub_keys[i][1]
        cyclo = [1] + [0] * (n - 1) + [1]
        cyclo = Polynomial(cyclo)

        # k0i = -t^{-1} namely the multiplicative inverse of t modulo qi
        k0i = mod_inverse(t, crt_moduli.qis[i]) * (-1)
        k0i = Polynomial([k0i])

        # pk0i * u + e0 + k0i * k1 = vi (this is ct0 before reduction in the Rqi ring)
        vi = pk0i * u + e0 + k0i * k1
        assert(len(vi.coefficients) - 1 == 2 * n - 2)

        # pk0i * u + e + k0i * k1 = ti mod Rqi = ct0
        ti = ciphertext[0]

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

        # pk1i * u + e1 = yi (this is ct1 before reduction in the Rqi ring)
        yi = pk1i * u + e1
        assert(len(yi.coefficients) - 1 == 2 * n - 2)

        # pk1i * u + e1 = ji mod Rqi = ct1
        ji = ciphertext[1]

        # assert that yi = ji mod Rqi
        yi_clone = copy.deepcopy(yi)
        # mod Rqi means that we need to:
        # - reduce the coefficients of yi_clone by the cyclotomic polynomial
        # - reduce the coefficients of yi_clone by the modulus
        yi_clone.reduce_coefficients_by_cyclo(cyclo.coefficients)
        yi_clone.reduce_coefficients_by_modulus(qis[i])
        assert yi_clone == ji

        # Calculate P2
        # divide ji - yi by the cyclotomic polynomial over Zqi to get P2
        num = ji + (Polynomial([-1]) * yi)
        # reduce the coefficients of num by the modulus qi 
        num.reduce_coefficients_by_modulus(qis[i])
        (quotient, rem) = poly_div(num.coefficients, cyclo.coefficients)
        # assert that the remainder is zero
        assert rem == []
        p2 = Polynomial(quotient)
        # assert that the degree of P2 is equal to n - 2
        assert len(p2.coefficients) - 1 == n - 2

        # Assert that ji - yi = P2 * cyclo mod Zqi
        lhs = ji + (Polynomial([-1]) * yi)
        rhs = p2 * cyclo
        # reduce the coefficients of lhs by the modulus qi
        lhs.reduce_coefficients_by_modulus(qis[i])
        assert lhs == rhs 
        
        # Calculate P1
        # divide ji - yi - P2 * cyclo by the modulus qi to get P1
        num = ji + (Polynomial([-1]) * yi) + Polynomial([-1]) * (p2 * cyclo)
        (quotient, rem) = poly_div(num.coefficients, [qis[i]])
        # assert that the remainder is zero
        assert rem == []
        p1 = Polynomial(quotient)

        # Assert that ti = vi + r1 * qi + r2 * cyclo mod Zp
        lhs = ji
        rhs = yi + (p1 * Polynomial([qis[i]])) + (p2 * cyclo)

        # remove the leading zeroes from the rhs 
        for j in range(len(rhs.coefficients)):
            if rhs.coefficients[j] != 0:
                rhs.coefficients = rhs.coefficients[j:]
                break

        assert lhs == rhs

        '''
        CIRCUIT - PHASE 1 - ASSIGNMENT PHASE

        In this phase, the private inputs are assigned to the circuit. These are the polynomials u, e0, e1, k1, R1, R2, P1 and P2. 
        We also assign the public inputs qi, t, k0i and b. N is a constant of the circuit.
        '''

        # ... Perform assignment here and expose public inputs...

        p = 21888242871839275222246405745257275088548364400416034343698204186575808495617

        # Every assigned value must be an element of the field Zp. Negative values `-z` are assigned as `p - z`
        u_assigned = assign_to_circuit(u, p)
        e0_assigned = assign_to_circuit(e0, p)
        e1_assigned = assign_to_circuit(e1, p)
        k1_assigned = assign_to_circuit(k1, p)
        r1_assigned = assign_to_circuit(r1, p)
        r2_assigned = assign_to_circuit(r2, p)
        p1_assigned = assign_to_circuit(p1, p)
        p2_assigned = assign_to_circuit(p2, p)

        '''
        CIRCUIT - PHASE 1 - RANGE CHECK OF PRIVATE POLYNOMIALS 

        In this phase, the coefficients of the private polynomials are checked to be in the correct range.
        '''

        # ... Perform range check here ...

        # sanity check. The coefficients of ti should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in ti.coefficients)

        # sanity check. The coefficients of pk0i should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in pk0i.coefficients)

        # constraint. The coefficients of u should be in the range [-1, 0, 1]
        assert all(coeff in [-1, 0, 1] for coeff in u.coefficients)
        # After the circuit assignement, the coefficients of u_assigned must be in [0, 1, p - 1]
        assert all(coeff in [0, 1, p - 1] for coeff in u_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of u_assigned to be in [0, 1, 2] (the normalization is constrained inside the circuit)
        u_normalized = Polynomial([(coeff + 1) % p for coeff in u_assigned.coefficients])
        assert all(coeff in [0, 1, 2] for coeff in u_normalized.coefficients)
        
        # sanity check. The coefficients of pk0i * u should be in the range $[-N \cdot \frac{q_i - 1}{2}, N \cdot \frac{q_i - 1}{2}]$
        bound = int((qis[i] - 1) / 2) * n
        res = pk0i * u
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of e0 should be in the range [-B, B] where B is the upper bound of the discrete Gaussian distribution
        b = int(discrete_gaussian.z_upper)
        assert all(coeff >= -b and coeff <= b for coeff in e0.coefficients)
        # After the circuit assignement, the coefficients of e0_assigned must be in [0, B] or [p - B, p - 1] (the normalization is constrained inside the circuit)
        assert all(coeff in range(0, b+1) or coeff in range(p - b, p) for coeff in e0_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of e0_assigned to be in [0, 2B]
        e0_normalized = Polynomial([(coeff + b) % p for coeff in e0_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*b for coeff in e0_normalized.coefficients)

        # sanity check. The coefficients of pk0i * u + e0 should be in the range $- (N \cdot \frac{q_i - 1}{2} + B), N \cdot \frac{q_i - 1}{2} + B]$
        bound = int((qis[i] - 1) / 2) * n + b
        res = pk0i * u + e0
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of R2 should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in r2.coefficients)
        # After the circuit assignement, the coefficients of r2_assigned must be in [0, bound] or [p - bound, p - 1] (the normalization is constrained inside the circuit)
        assert all(coeff in range(0, int(bound) + 1) or coeff in range(p - int(bound), p) for coeff in r2_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of r2_assigned to be in [0, 2*bound]
        r2_normalized = Polynomial([(coeff + int(bound)) % p for coeff in r2_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*bound for coeff in r2_normalized.coefficients)

        # sanity check. The coefficients of R2 * cyclo should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        res = r2 * cyclo
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of k1 should be in the range [-(t-1)/2, (t-1)/2]
        bound = int((t - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in k1.coefficients)
        # After the circuit assignement, the coefficients of k1_assigned must be in [0, bound] or [p - bound, p - 1] (the normalization is constrained inside the circuit)
        assert all(coeff in range(0, int(bound) + 1) or coeff in range(p - int(bound), p) for coeff in k1_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of k1_assigned to be in [0, 2*bound]
        k1_normalized = Polynomial([(coeff + int(bound)) % p for coeff in k1_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*bound for coeff in k1_normalized.coefficients)

        # sanity check. The coefficients of k1 * k0i should be in the range $[-\frac{t - 1}{2} \cdot |K_i^{0}|, \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = int((t - 1) / 2) * abs(k0i.coefficients[0])
        res = k1 * k0i
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # sanity check. The coefficients of vi (pk0i * u + e0 + k0i * k1) should be in the range $[- (N \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|), N \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = int((qis[i] - 1) / 2) * n + b + int((t - 1) / 2) * abs(k0i.coefficients[0])
        assert all(coeff >= -bound and coeff <= bound for coeff in vi.coefficients)

        # sanity check. The coefficients of ti - vi should be in the range $ [- ((N+1) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|), (N+1) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = int((qis[i] - 1) / 2) * (n + 1) + b + int((t - 1) / 2) * abs(k0i.coefficients[0])
        sub = ti + (Polynomial([-1]) * vi)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # sanity check. The coefficients of ti - vi - R2 * cyclo should be in the range $[- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|), (N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = ((qis[i] - 1) / 2) * (n + 2) + b + ((t - 1) / 2) * abs(k0i.coefficients[0])
        sub = ti + (Polynomial([-1]) * vi) + Polynomial([-1]) * (r2 * cyclo)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # constraint. The coefficients of (ti - vi - R2 * cyclo) / qi = R1 should be in the range $[\frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|)}{q_i}, \frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|}{q_i}]$
        bound = (int((qis[i] - 1) / 2) * (n + 2) + b + int((t - 1) / 2) * abs(k0i.coefficients[0])) / qis[i]
        # round bound to the nearest integer
        bound = int(bound)
        assert all(coeff >= -bound and coeff <= bound for coeff in r1.coefficients)
        # After the circuit assignement, the coefficients of r1_assigned must be in [0, bound] or [p - bound, p - 1]
        assert all(coeff in range(0, int(bound) + 1) or coeff in range(p - int(bound), p) for coeff in r1_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of r1_assigned to be in [0, 2*bound]
        r1_normalized = Polynomial([(coeff + int(bound)) % p for coeff in r1_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*bound for coeff in r1_normalized.coefficients)
        
        # sanity check. The coefficients of ji should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in ji.coefficients)

        # sanity check. The coefficients of pk1i should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in pk1i.coefficients)
        
        # sanity check. The coefficients of pk1i * u should be in the range $[-N \cdot \frac{q_i - 1}{2}, N \cdot \frac{q_i - 1}{2}]$
        bound = int((qis[i] - 1) / 2) * n
        res = pk1i * u
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of e1 should be in the range [-B, B] where B is the upper bound of the discrete Gaussian distribution
        b = int(discrete_gaussian.z_upper)
        assert all(coeff >= -b and coeff <= b for coeff in e1.coefficients)
        # After the circuit assignement, the coefficients of e1_assigned must be in [0, B] or [p - B, p - 1] (the normalization is constrained inside the circuit)
        assert all(coeff in range(0, b+1) or coeff in range(p - b, p) for coeff in e1_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of e1_assigned to be in [0, 2B]
        e1_normalized = Polynomial([(coeff + b) % p for coeff in e1_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*b for coeff in e1_normalized.coefficients)

        # sanity check. The coefficients of pk1i * u + e1 (=yi) should be in the range $- (N \cdot \frac{q_i - 1}{2} + B), N \cdot \frac{q_i - 1}{2} + B]$
        bound = int((qis[i] - 1) / 2) * n + b
        res = pk1i * u + e1
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)
        assert all(coeff >= -bound and coeff <= bound for coeff in yi.coefficients)
    
        # constraint. The coefficients of p2 should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in p2.coefficients)
        # After the circuit assignement, the coefficients of p2_assigned must be in [0, bound] or [p - bound, p - 1] (the normalization is constrained inside the circuit)
        assert all(coeff in range(0, int(bound) + 1) or coeff in range(p - int(bound), p) for coeff in p2_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of p2_assigned to be in [0, 2*bound]
        p2_normalized = Polynomial([(coeff + int(bound)) % p for coeff in p2_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*bound for coeff in p2_normalized.coefficients)

        # sanity check. The coefficients of P2 * cyclo should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        res = p2 * cyclo
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # sanity check. The coefficients of ji - yi should be in the range $ [- ((N+1) \cdot \frac{q_i - 1}{2} + B), (N+1) \cdot \frac{q_i - 1}{2} + B]$
        bound = int((qis[i] - 1) / 2) * (n + 1) + b
        sub = ji + (Polynomial([-1]) * yi)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # sanity check. The coefficients of ji - yi - P2 * cyclo should be in the range $[- ((N+2) \cdot \frac{q_i - 1}{2} + B), (N+2) \cdot \frac{q_i - 1}{2} + B]$
        bound = ((qis[i] - 1) / 2) * (n + 2) + b
        sub = ji + (Polynomial([-1]) * yi) + Polynomial([-1]) * (p2 * cyclo)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # constraint. The coefficients of (ji - yi - P2 * cyclo) / qi = P1 should be in the range $[\frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B)}{q_i}, \frac{(N+2) \cdot \frac{q_i - 1}{2} + B}{q_i}]$
        bound = (int((qis[i] - 1) / 2) * (n + 2) + b) / qis[i]
        # round bound to the nearest integer
        bound = int(bound)
        assert all(coeff >= -bound and coeff <= bound for coeff in p1.coefficients)
        # After the circuit assignement, the coefficients of p1_assigned must be in [0, bound] or [p - bound, p - 1]
        assert all(coeff in range(0, int(bound) + 1) or coeff in range(p - int(bound), p) for coeff in p1_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of p1_assigned to be in [0, 2*bound]
        p1_normalized = Polynomial([(coeff + int(bound)) % p for coeff in p1_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*bound for coeff in p1_normalized.coefficients)        
        
        '''
        CIRCUIT - END OF PHASE 1 - WITNESS COMMITMENT 

        At the end of phase 1, the witness is committed and a challenge is generated using the Fiat-Shamir heuristic.
        '''

        # For the sake of simplicity, we generate a random challenge here
        alpha = randint(0, 1000)

        '''
        CIRCUIT - PHASE 2 - ASSIGNMENT PHASE
        The public inputs are assigned to the circuit. These are the polynomials pk0i, pk1i, cyclo, ti, ct1 evaluated at alpha.
        '''

        # The evaluation of pk0i_alpha, cyclo_alpha, ti_alpha, ji_alpha is performed outside the circuit
        pk0i_alpha = pk0i.evaluate(alpha) 
        pk1i_alpha = pk1i.evaluate(alpha)
        cyclo_alpha = cyclo.evaluate(alpha)
        ti_alpha = ti.evaluate(alpha)
        ji_alpha = ciphertext[1].evaluate(alpha)

        # ... Perform assignment here and expose expose public inputs ...

        '''
        CIRCUIT - PHASE 2 - CORRECT ENCRYPTION CONSTRAINT

        We need to prove that ti = vi + r1 * qi + r2 * cyclo mod Zp.
        We need to prove that ji = yi + p1 * qi + p2 * cyclo mod Zp.
        We do that by proving that LHS(alpha) = RHS(alpha) for a random alpha according to Scwhartz-Zippel lemma.
        '''

        u_alpha = u_assigned.evaluate(alpha)
        e0_alpha = e0_assigned.evaluate(alpha)
        k1_alpha = k1_assigned.evaluate(alpha)
        r1_alpha = r1_assigned.evaluate(alpha)
        r2_alpha = r2_assigned.evaluate(alpha)

        lhs = ti_alpha 
        rhs = (pk0i_alpha * u_alpha + e0_alpha + (k1_alpha * k0i.coefficients[0]) + (r1_alpha * qis[i]) + (r2_alpha * cyclo_alpha))

        assert lhs % p == rhs % p

        p1_alpha = p1_assigned.evaluate(alpha)
        p2_alpha = p2_assigned.evaluate(alpha)
        e1_alpha = e1_assigned.evaluate(alpha)

        lhs = ji_alpha
        rhs = (pk1i_alpha * u_alpha + e1_alpha) + (p1_alpha * qis[i]) + (p2_alpha * cyclo_alpha)

        assert lhs % p == rhs % p

def assign_to_circuit(poly: Polynomial, p: int) -> Polynomial:
    '''
    This function takes a polynomial and returns its coefficients in the field Zp
    `poly` is the polynomial to be assigned to the circuit
    `p` is the field modulus
    '''
    assigned_coefficients = []
    for coeff in poly.coefficients:
        if coeff < 0:
            coeff = coeff % p
        if coeff > p:
            coeff = coeff % p
        assigned_coefficients.append(coeff)

    return Polynomial(assigned_coefficients)

main()