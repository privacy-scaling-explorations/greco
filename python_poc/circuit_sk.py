from bfv.crt import CRTModuli
from bfv.bfv import BFVCrt
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div
from bfv.utils import mod_inverse
from random import randint
import copy
from utils import assign_to_circuit

def main(): 
    
    '''
    ENCRYPTION PHASE - performed outside the circuit

    In this phase, the encryption operation is performed. Later, the circuit will be used to prove that the encryption was performed correctly.
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
    for i in range(len(crt_moduli.qis)):
        ais.append(bfv_crt.bfv_qis[i].rlwe.Rq.sample_polynomial())

    m = bfv_crt.bfv_q.rlwe.Rt.sample_polynomial()

    ciphertexts = bfv_crt.SecretKeyEncrypt(s, ais, e, m)

    # Sanity check for valid decryption    
    message_prime = bfv_crt.Decrypt(s, ciphertexts)

    assert m == message_prime

    # k1 = [QM]t namely the scaled message polynomial
    k1 = Polynomial([crt_moduli.q]) * m
    k1.reduce_coefficients_by_modulus(t)

    r2s = []
    r1s = []

    # Each round of the loop simulates a proof generation phase for a different CRT basis
    for i, ciphertext in enumerate(ciphertexts):

        '''
        SETUP PHASE - performed outside the circuit

        For each CRT basis, we need to compute the polynomials R1 and R2 (check this doc for more details: https://hackmd.io/@gaussian/HJ8DYyjPp)
        '''

        ct0i = ciphertext[0]
        ai = Polynomial([-1]) * ciphertext[1]
        cyclo = [1] + [0] * (n - 1) + [1]
        cyclo = Polynomial(cyclo)

        # k0i = -t^{-1} namely the multiplicative inverse of t modulo qi
        k0i = mod_inverse(t, crt_moduli.qis[i]) * (-1)
        k0i = Polynomial([k0i])

        # ai * s + e + k0i * k1 = hat(ct0i) (this is ct0i before reduction in the Rqi ring)
        ct0i_hat = ai * s + e + k0i * k1
        assert(len(ct0i_hat.coefficients) - 1 == 2 * n - 2)

        # ai * s + e + k0i * k1 = ct0i mod Rqi
        # assert that ct0i_hat = ct0i mod Rqi
        ct0i_hat_clone = copy.deepcopy(ct0i_hat)
        # mod Rqi means that we need to:
        # - reduce the coefficients of ct0i_hat_clone by the cyclotomic polynomial
        # - reduce the coefficients of ct0i_hat_clone by the modulus
        ct0i_hat_clone.reduce_coefficients_by_cyclo(cyclo.coefficients)
        ct0i_hat_clone.reduce_coefficients_by_modulus(qis[i])
        assert ct0i_hat_clone == ct0i

        # Calculate R2
        # divide ct0i - ct0i_hat by the cyclotomic polynomial over Zqi to get R2
        num = ct0i + (Polynomial([-1]) * ct0i_hat)
        # reduce the coefficients of num by the modulus qi 
        num.reduce_coefficients_by_modulus(qis[i])
        (quotient, rem) = poly_div(num.coefficients, cyclo.coefficients)
        # assert that the remainder is zero
        assert rem == []
        r2 = Polynomial(quotient)
        # assert that the degree of R2 is equal to n - 2
        assert len(r2.coefficients) - 1 == n - 2

        # Assert that ct0i - ct0i_hat = R2 * cyclo mod Zqi
        lhs = ct0i + (Polynomial([-1]) * ct0i_hat)
        rhs = r2 * cyclo
        # reduce the coefficients of lhs by the modulus qi
        lhs.reduce_coefficients_by_modulus(qis[i])
        assert lhs == rhs 
        
        # Calculate R1
        # divide ct0i - ct0i_hat - R2 * cyclo by the modulus qi to get R1
        num = ct0i + (Polynomial([-1]) * ct0i_hat) + Polynomial([-1]) * (r2 * cyclo)
        (quotient, rem) = poly_div(num.coefficients, [qis[i]])
        # assert that the remainder is zero
        assert rem == []
        r1 = Polynomial(quotient)
        # assert that the degree of R1 is 2n - 2
        assert len(r1.coefficients) - 1 == 2 * n - 2

        # Assert that ct0i = ct0i_hat + r1 * qi + r2 * cyclo mod Zp
        lhs = ct0i
        rhs = ct0i_hat + (r1 * Polynomial([qis[i]])) + (r2 * cyclo)

        # remove the leading zeroes from rhs until the length of rhs.coefficients is equal to n
        while len(rhs.coefficients) > n and rhs.coefficients[0] == 0:
            rhs.coefficients.pop(0)

        assert lhs == rhs

        '''
        CIRCUIT - PHASE 1 - ASSIGNMENT PHASE

        In this phase, the polynomials belonging to the matrix S are assigned to the circuit.
        We also assign the public inputs qi, t, k0i and b, while N is a constant of the circuit.
        '''

        # ... Perform assignment here and expose public inputs...

        p = 21888242871839275222246405745257275088548364400416034343698204186575808495617

        # Every assigned value must be an element of the field Zp. Negative values `-z` are assigned as `p - z`
        s_assigned = assign_to_circuit(s, p)
        e_assigned = assign_to_circuit(e, p)
        k1_assigned = assign_to_circuit(k1, p)
        r1_assigned = assign_to_circuit(r1, p)
        r2_assigned = assign_to_circuit(r2, p)

        '''
        CIRCUIT - PHASE 1 - RANGE CHECK OF PRIVATE POLYNOMIALS 

        In this phase, the coefficients of the private polynomials are checked to be in the correct range.
        '''

        # ... Perform range check here ...

        # sanity check. The coefficients of ct0i should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in ct0i.coefficients)

        # sanity check. The coefficients of ai should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in ai.coefficients)

        # constraint. The coefficients of s should be in the range [-1, 0, 1]
        assert all(coeff in [-1, 0, 1] for coeff in s.coefficients)
        # After the circuit assignement, the coefficients of s_assigned must be in [0, 1, p - 1]
        assert all(coeff in [0, 1, p - 1] for coeff in s_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of s_assigned to be in [0, 1, 2] (the normalization is constrained inside the circuit)
        s_normalized = Polynomial([(coeff + 1) % p for coeff in s_assigned.coefficients])
        assert all(coeff in [0, 1, 2] for coeff in s_normalized.coefficients)
        
        # sanity check. The coefficients of ai * s should be in the range $[-N \cdot \frac{q_i - 1}{2}, N \cdot \frac{q_i - 1}{2}]$
        bound = int((qis[i] - 1) / 2) * n
        res = ai * s
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of e should be in the range [-B, B] where B is the upper bound of the discrete Gaussian distribution
        b = int(discrete_gaussian.z_upper)
        assert all(coeff >= -b and coeff <= b for coeff in e.coefficients)
        # After the circuit assignement, the coefficients of e_assigned must be in [0, B] or [p - B, p - 1] (the normalization is constrained inside the circuit)
        assert all(coeff in range(0, b+1) or coeff in range(p - b, p) for coeff in e_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of e_assigned to be in [0, 2B]
        e_normalized = Polynomial([(coeff + b) % p for coeff in e_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*b for coeff in e_normalized.coefficients)

        # sanity check. The coefficients of ai * s + e should be in the range $- (N \cdot \frac{q_i - 1}{2} + B), N \cdot \frac{q_i - 1}{2} + B]$
        bound = int((qis[i] - 1) / 2) * n + b
        res = ai * s + e
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

        # sanity check. The coefficients of ct0i_hat (ai * s + e + k1 * k0i) should be in the range $[- (N \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|), N \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = int((qis[i] - 1) / 2) * n + b + int((t - 1) / 2) * abs(k0i.coefficients[0])
        assert all(coeff >= -bound and coeff <= bound for coeff in ct0i_hat.coefficients)

        # sanity check. The coefficients of ct0i - ct0i_hat should be in the range $ [- ((N+1) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|), (N+1) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = int((qis[i] - 1) / 2) * (n + 1) + b + int((t - 1) / 2) * abs(k0i.coefficients[0])
        sub = ct0i + (Polynomial([-1]) * ct0i_hat)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # sanity check. The coefficients of ct0i - ct0i_hat - R2 * cyclo should be in the range $[- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|), (N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|]$
        bound = ((qis[i] - 1) / 2) * (n + 2) + b + ((t - 1) / 2) * abs(k0i.coefficients[0])
        sub = ct0i + (Polynomial([-1]) * ct0i_hat) + Polynomial([-1]) * (r2 * cyclo)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # constraint. The coefficients of (ct0i - ct0i_hat - R2 * cyclo) / qi = R1 should be in the range $[\frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_i^{0}|)}{q_i}, \frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_i^{0}|}{q_i}]$
        bound = (int((qis[i] - 1) / 2) * (n + 2) + b + int((t - 1) / 2) * abs(k0i.coefficients[0])) / qis[i]
        # round bound to the nearest integer
        bound = int(bound)
        assert all(coeff >= -bound and coeff <= bound for coeff in r1.coefficients)
        # After the circuit assignement, the coefficients of r1_assigned must be in [0, bound] or [p - bound, p - 1]
        assert all(coeff in range(0, int(bound) + 1) or coeff in range(p - int(bound), p) for coeff in r1_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of r1_assigned to be in [0, 2*bound]
        r1_normalized = Polynomial([(coeff + int(bound)) % p for coeff in r1_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*bound for coeff in r1_normalized.coefficients)

        '''
        CIRCUIT - END OF PHASE 1 - WITNESS COMMITMENT 

        At the end of phase 1, the witness is committed and a challenge is generated using the Fiat-Shamir heuristic.
        '''

        # For the sake of simplicity, we generate a random challenge here
        alpha = randint(0, 1000)

        '''
        CIRCUIT - PHASE 2 - ASSIGNMENT PHASE
        The public inputs are assigned to the circuit. These are the polynomials ai, cyclo, ct0i evaluated at alpha.
        '''

        # The evaluation of ai_alpha, cyclo_alpha, ct0i_alpha is performed outside the circuit
        ai_alpha = ai.evaluate(alpha) 
        cyclo_alpha = cyclo.evaluate(alpha)
        ct0i_alpha = ct0i.evaluate(alpha)

        # ... Perform assignment here and expose expose public inputs ...

        '''
        CIRCUIT - PHASE 2 - CORRECT ENCRYPTION CONSTRAINT

        We need to prove that ct0i = ct0i_hat + r1 * qi + r2 * cyclo mod Zp.
        We do that by proving that LHS(alpha) = RHS(alpha) for a random alpha according to Scwhartz-Zippel lemma.
        '''

        s_alpha = s_assigned.evaluate(alpha)
        e_alpha = e_assigned.evaluate(alpha)
        k1_alpha = k1_assigned.evaluate(alpha)
        r1_alpha = r1_assigned.evaluate(alpha)
        r2_alpha = r2_assigned.evaluate(alpha)

        lhs = ct0i_alpha 
        rhs = (ai_alpha * s_alpha + e_alpha + (k1_alpha * k0i.coefficients[0]) + (r1_alpha * qis[i]) + (r2_alpha * cyclo_alpha))

        assert lhs % p == rhs % p

        r1s.append(r1_assigned)
        r2s.append(r2_assigned)

    # Parse the inputs into a JSON format such this can be used as input for the (real) circuit
    json_input = {
        "s": [str(coef) for coef in s_assigned.coefficients],
        "e": [str(coef) for coef in e_assigned.coefficients],
        "k1": [str(coef) for coef in k1_assigned.coefficients],
        "r2s": [[str(coef) for coef in r2.coefficients] for r2 in r2s],
        "r1s": [[str(coef) for coef in r1.coefficients] for r1 in r1s],
    }

    # write the inputs to a a json file
    import json
    with open('./src/data/sk_enc_input.json', 'w') as f:
        json.dump(json_input, f)


main()