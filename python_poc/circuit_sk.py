from bfv.crt import CRTModuli
from bfv.bfv import BFVCrt
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div
from bfv.utils import mod_inverse
from random import randint
import copy
from utils import assign_to_circuit
import argparse
import json

def main(args): 
    
    '''
    ENCRYPTION PHASE - performed outside the circuit

    In this phase, the encryption operation is performed. Later, the circuit will be used to prove that the encryption was performed correctly.
    '''

    n = args.n
    qis = args.qis
    qis = json.loads(qis)
    t = args.t

    crt_moduli = CRTModuli(qis)
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

    # `k0is` are the negative multiplicative inverses of t modulo each qi. These are normalized to be in the range [0, p), where p is the modulus of the prime field of the circuit
    k0is = []

    # `r2is` are the polynomials r2i for each CRT basis. The coefficients are normalized to be in the range [0, p), where p is the modulus of the prime field of the circuit
    r2is = []

    # `r1is` are the polynomials r1i for each CRT basis. The coefficients normalized to be in the range [0, p), where p is the modulus of the prime field of the circuit
    r1is = []

    # `r1_bounds` are the bounds for the coefficients of r1i for each CRT basis
    r1_bounds = []

    # `r2_bounds` are the bounds for the coefficients of r2i for each CRT basis
    r2_bounds = []

    # `ct0is` are the polynomials ct0i for each CRT basis. The coefficients normalized to be in the range [0, p), where p is the modulus of the prime field of the circuit
    ct0is = []

    # `ais` are the polynomials ai for each CRT basis. The coefficients normalized to be in the range [0, p), where p is the modulus of the prime field of the circuit
    ais= []

    p = 21888242871839275222246405745257275088548364400416034343698204186575808495617

    # Each round of the loop simulates a proof generation phase for a different CRT basis
    for i, ciphertext in enumerate(ciphertexts):

        '''
        SETUP PHASE - performed outside the circuit

        For each CRT basis, we need to compute the polynomials r1i and r2i (check this doc for more details: https://hackmd.io/@gaussian/HJ8DYyjPp)
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

        # Calculate r2i
        # divide ct0i - ct0i_hat by the cyclotomic polynomial over Zqi to get r2i
        num = ct0i + (Polynomial([-1]) * ct0i_hat)
        # reduce the coefficients of num by the modulus qi 
        num.reduce_coefficients_by_modulus(qis[i])
        (quotient, rem) = poly_div(num.coefficients, cyclo.coefficients)
        # assert that the remainder is zero
        assert rem == []
        r2i = Polynomial(quotient)
        # assert that the degree of r2i is equal to n - 2
        assert len(r2i.coefficients) - 1 == n - 2

        # Assert that ct0i - ct0i_hat = r2i * cyclo mod Zqi
        lhs = ct0i + (Polynomial([-1]) * ct0i_hat)
        rhs = r2i * cyclo
        # reduce the coefficients of lhs by the modulus qi
        lhs.reduce_coefficients_by_modulus(qis[i])
        assert lhs == rhs 
        
        # Calculate r1i
        # divide ct0i - ct0i_hat - r2i * cyclo by the modulus qi to get r1i
        num = ct0i + (Polynomial([-1]) * ct0i_hat) + Polynomial([-1]) * (r2i * cyclo)
        (quotient, rem) = poly_div(num.coefficients, [qis[i]])
        # assert that the remainder is zero
        assert rem == []
        r1i = Polynomial(quotient)
        # assert that the degree of r1i is 2n - 2
        assert len(r1i.coefficients) - 1 == 2 * n - 2

        # Assert that ct0i = ct0i_hat + r1i * qi + r2i * cyclo mod Zp
        lhs = ct0i
        rhs = ct0i_hat + (r1i * Polynomial([qis[i]])) + (r2i * cyclo)

        # remove the leading zeroes from rhs until the length of rhs.coefficients is equal to n
        while len(rhs.coefficients) > n and rhs.coefficients[0] == 0:
            rhs.coefficients.pop(0)

        assert lhs == rhs

        '''
        CIRCUIT - PHASE 0 - ASSIGNMENT PHASE

        In this phase, all the polynomials are assigned to the circuit. Namely:
        * polynomial ai and the scalars qi and k0i from the matrix Ui
        * polynomials s, e, k, r1i, r2i from the matrix Si
        * polynomials ct0i

        The cyclotomic polynomial is not assigned to the circuit, as this is not an input but a constant parameter.

        ai, qi, k0i and ct0i are exposed as public inputs, while the other polynomials are kept private.
        '''

        # ... Perform assignment here and expose public inputs...

        # Every assigned value must be an element of the field Zp. Negative values `-z` are assigned as `p - z`
        ai_assigned = assign_to_circuit(ai, p)
        k0i_assigned = assign_to_circuit(k0i, p).coefficients[0]
        qi_assigned = assign_to_circuit(Polynomial([qis[i]]), p).coefficients[0]
        ct0i_assigned = assign_to_circuit(ct0i, p)
        s_assigned = assign_to_circuit(s, p)
        e_assigned = assign_to_circuit(e, p)
        k1_assigned = assign_to_circuit(k1, p)
        r1i_assigned = assign_to_circuit(r1i, p)
        r2i_assigned = assign_to_circuit(r2i, p)

        '''
        CIRCUIT - PHASE 0 - RANGE CHECK OF PRIVATE POLYNOMIALS 

        In this phase, the coefficients of the private polynomials from matrix Si are checked to be in the correct range.
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

        # constraint. The coefficients of r2i should be in the range [-(qi-1)/2, (qi-1)/2]
        r2i_bound = int((qis[i] - 1) / 2)
        r2_bounds.append(r2i_bound)
        assert all(coeff >= -r2i_bound and coeff <= r2i_bound for coeff in r2i.coefficients)
        # After the circuit assignement, the coefficients of r2i_assigned must be in [0, r2i_bound] or [p - r2i_bound, p - 1] (the normalization is constrained inside the circuit)
        assert all(coeff in range(0, int(r2i_bound) + 1) or coeff in range(p - int(r2i_bound), p) for coeff in r2i_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of r2i_assigned to be in [0, 2*r2i_bound]
        r2i_normalized = Polynomial([(coeff + int(r2i_bound)) % p for coeff in r2i_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*r2i_bound for coeff in r2i_normalized.coefficients)

        # sanity check. The coefficients of r2i * cyclo should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        res = r2i * cyclo
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of k1 should be in the range [-(t-1)/2, (t-1)/2]
        k1_bound = int((t - 1) / 2)
        assert all(coeff >= -k1_bound and coeff <= k1_bound for coeff in k1.coefficients)
        # After the circuit assignement, the coefficients of k1_assigned must be in [0, k1_bound] or [p - k1_bound, p - 1] (the normalization is constrained inside the circuit)
        assert all(coeff in range(0, int(k1_bound) + 1) or coeff in range(p - int(k1_bound), p) for coeff in k1_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of k1_assigned to be in [0, 2*k1_bound]
        k1_normalized = Polynomial([(coeff + int(k1_bound)) % p for coeff in k1_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*k1_bound for coeff in k1_normalized.coefficients)

        # sanity check. The coefficients of k1 * k0i should be in the range $[-\frac{t - 1}{2} \cdot |K_{0,i}|, \frac{t - 1}{2} \cdot |K_{0,i}|]$
        bound = int((t - 1) / 2) * abs(k0i.coefficients[0])
        res = k1 * k0i
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # sanity check. The coefficients of ct0i_hat (ai * s + e + k1 * k0i) should be in the range $[- (N \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|), N \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|]$
        bound = int((qis[i] - 1) / 2) * n + b + int((t - 1) / 2) * abs(k0i.coefficients[0])
        assert all(coeff >= -bound and coeff <= bound for coeff in ct0i_hat.coefficients)

        # sanity check. The coefficients of ct0i - ct0i_hat should be in the range $ [- ((N+1) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|), (N+1) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|]$
        bound = int((qis[i] - 1) / 2) * (n + 1) + b + int((t - 1) / 2) * abs(k0i.coefficients[0])
        sub = ct0i + (Polynomial([-1]) * ct0i_hat)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # sanity check. The coefficients of ct0i - ct0i_hat - r2i * cyclo should be in the range $[- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|), (N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|]$
        bound = ((qis[i] - 1) / 2) * (n + 2) + b + ((t - 1) / 2) * abs(k0i.coefficients[0])
        sub = ct0i + (Polynomial([-1]) * ct0i_hat) + Polynomial([-1]) * (r2i * cyclo)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # constraint. The coefficients of (ct0i - ct0i_hat - r2i * cyclo) / qi = r1i should be in the range $[\frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|)}{q_i}, \frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|}{q_i}]$
        r1i_bound = (int((qis[i] - 1) / 2) * (n + 2) + b + int((t - 1) / 2) * abs(k0i.coefficients[0])) / qis[i]
        # round bound to the nearest integer
        r1i_bound = int(r1i_bound)
        r1_bounds.append(r1i_bound)
        assert all(coeff >= -r1i_bound and coeff <= r1i_bound for coeff in r1i.coefficients)
        # After the circuit assignement, the coefficients of r1i_assigned must be in [0, r1i_bound] or [p - r1i_bound, p - 1]
        assert all(coeff in range(0, int(r1i_bound) + 1) or coeff in range(p - int(r1i_bound), p) for coeff in r1i_assigned.coefficients)
        # To perform a range check with a smaller lookup table, we normalize the coefficients of r1i_assigned to be in [0, 2*r1i_bound]
        r1i_normalized = Polynomial([(coeff + int(r1i_bound)) % p for coeff in r1i_assigned.coefficients])
        assert all(coeff >= 0 and coeff <= 2*r1i_bound for coeff in r1i_normalized.coefficients)

        '''
        CIRCUIT - END OF PHASE 0 - WITNESS COMMITMENT 

        At the end of phase 0, the witness is committed and a challenge is generated using the Fiat-Shamir heuristic.
        '''

        # For the sake of simplicity, we generate a random challenge here
        alpha = randint(0, 1000)

        '''
        CIRCUIT - PHASE 1 - ASSIGNMENT PHASE

        Challenge alpha is fetched from phase 0 and used to evaluate the public polynomials at alpha.
        These polynomials are evaluate at alpha outside the circuit and the evaluation is assigned to the circuit and exposed as public inputs.
        Since the polynomials are public from phase 0, the evaluation at alpha doesn't need to be constrained inside the circuit,
        but can safely be performed (and verified) outside the circuit

        The public inputs are assigned to the circuit. These are:
        * ai_alpha
        * cyclo_alpha
        * ct0i_alpha
        '''

        # The evaluation of ai_alpha, cyclo_alpha, ct0i_alpha is performed outside the circuit

        ai_alpha = ai_assigned.evaluate(alpha) 
        cyclo_alpha = cyclo.evaluate(alpha)
        ct0i_alpha = ct0i_assigned.evaluate(alpha)

        # ... Perform assignment here and expose expose public inputs ...

        '''
        CIRCUIT - PHASE 1 - CORRECT ENCRYPTION CONSTRAINT

        We need to prove that ct0i = ct0i_hat + r1i * qi + r2i * cyclo mod Zp.
        We do that by proving that LHS(alpha) = RHS(alpha) for a random alpha according to Scwhartz-Zippel lemma.
        '''

        s_alpha = s_assigned.evaluate(alpha)
        e_alpha = e_assigned.evaluate(alpha)
        k1_alpha = k1_assigned.evaluate(alpha)
        r1i_alpha = r1i_assigned.evaluate(alpha)
        r2i_alpha = r2i_assigned.evaluate(alpha)

        lhs = ct0i_alpha 
        rhs = (ai_alpha * s_alpha + e_alpha + (k1_alpha * k0i_assigned) + (r1i_alpha * qi_assigned) + (r2i_alpha * cyclo_alpha))

        assert lhs % p == rhs % p

        # The verifier shoud fetch alpha and check that ai_alpha, cyclo_alpha, ct0i_alpha were generated correctly
        ai_alpha == ai_assigned.evaluate(alpha)
        cyclo_alpha == cyclo.evaluate(alpha)
        ct0i_alpha == ct0i_assigned.evaluate(alpha)

        r1is.append(r1i_assigned)
        r2is.append(r2i_assigned)
        ais.append(ai_assigned)
        ct0is.append(ct0i_assigned)
        k0is.append(k0i_assigned)

    # Parse the inputs into a JSON format such this can be used as input for the (real) circuit
    json_input = {
        "s": [str(coef) for coef in s_assigned.coefficients],
        "e": [str(coef) for coef in e_assigned.coefficients],
        "qis": [str(qi) for qi in qis],
        "k1": [str(coef) for coef in k1_assigned.coefficients],
        "k0is": [str(k0i) for k0i in k0is],
        "r2is": [[str(coef) for coef in r2i.coefficients] for r2i in r2is],
        "r1is": [[str(coef) for coef in r1i.coefficients] for r1i in r1is],
        "ais": [[str(coef) for coef in ai.coefficients] for ai in ais],
        "ct0is": [[str(coef) for coef in ct0i.coefficients] for ct0i in ct0is],
    }

    # write the inputs to a json file
    with open(args.output, 'w') as f:
        json.dump(json_input, f)

    print(f"pub const N: usize = {n};")
    print(f"pub const E_BOUND: u64 = {b};")
    print(f"pub const R1_BOUNDS: [u64; {len(r1_bounds)}] = [{', '.join(map(str, r1_bounds))}];")
    print(f"pub const R2_BOUNDS: [u64; {len(r2_bounds)}] = [{', '.join(map(str, r2_bounds))}];")
    print(f"pub const K1_BOUND: u64 = {k1_bound};")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate rust constants and json inputs for BFV zk proof of secret key encryption circuit"
    )
    parser.add_argument(
        "-n", type=int, required=True, help="Degree of f(x), must be a power of 2."
    )
    parser.add_argument(
        "-qis", type=str, required=True, help="List of qis such that qis[i] is the modulus of the i-th CRT basis of the modulus q of the ciphertext space."
    )
    parser.add_argument(
        "-t", type=int, required=True, help="Modulus t of the plaintext space."
    )
    parser.add_argument(
        "-output", type=str, required=True, default="./src/data/sk_enc_input.json", help="Output file name."
    )

    args = parser.parse_args()
    main(args)
