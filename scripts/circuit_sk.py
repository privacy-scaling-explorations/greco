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
    ENCRYPTION PHASE - performed outside the circuit.
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

    ctis = bfv_crt.SecretKeyEncrypt(s, ais, e, m)

    # Sanity check for valid decryption    
    message_prime = bfv_crt.Decrypt(s, ctis)

    assert m == message_prime

    # k1 = [QM]t namely the scaled message polynomial
    k1 = Polynomial([crt_moduli.q]) * m
    k1.reduce_coefficients_by_modulus(t)

    # `p` is the modulus of the prime field of the circuit
    p = 21888242871839275222246405745257275088548364400416034343698204186575808495617

    # `r2is` are the polynomials r2i for each i-th CRT basis.
    r2is = []

    # `r1is` are the polynomials r1i for each i-th CRT basis.
    r1is = []

    # `k0is` are the negative multiplicative inverses of t modulo each qi.
    k0is = []

    # `ct0is` are the polynomials ct0i for each CRT basis. 
    ct0is = []

    # `ct0is_hat` are the polynomials ct0i_hat for each CRT basis.
    ct0is_hat = []

    '''
    SETUP PHASE - performed outside the circuit
    For each CRT basis, we need to compute the polynomials r1i and r2i (check this doc for more details: https://hackmd.io/@gaussian/r1W98Kqqa)
    '''

    cyclo = [1] + [0] * (n - 1) + [1]
    cyclo = Polynomial(cyclo)

    for i, cti in enumerate(ctis):

        ct0i = cti[0]
        ai = Polynomial([-1]) * cti[1]

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

        r2is.append(r2i)
        r1is.append(r1i)
        k0is.append(k0i)
        ct0is.append(ct0i)
        ct0is_hat.append(ct0i_hat)

    # `r1_bounds` are the bounds for the coefficients of r1i for each CRT basis
    r1_bounds = []

    # `r2_bounds` are the bounds for the coefficients of r2i for each CRT basis
    r2_bounds = []

    '''
    CIRCUIT - PHASE 0
    '''

    # Every assigned value must be an element of the field Zp. Negative coefficients `-z` are assigned as `p - z`
    s_assigned = assign_to_circuit(s, p)
    e_assigned = assign_to_circuit(e, p)
    k1_assigned = assign_to_circuit(k1, p)

    r1is_assigned = []
    r2is_assigned = []

    for i in range(len(ctis)):
        r1i_assigned = assign_to_circuit(r1is[i], p)
        r2i_assigned = assign_to_circuit(r2is[i], p)
        r1is_assigned.append(r1i_assigned)
        r2is_assigned.append(r2i_assigned)


    # For the sake of simplicity, we generate a random challenge here
    gamma = randint(0, 1000)

    '''
    CIRCUIT - PHASE 1
    '''

    # Every assigned value must be an element of the field Zp. Negative coefficients `-z` are assigned as `p - z`
    ais_at_gamma_assigned = []
    ct0is_at_gamma_assigned = []
    qis_assigned = []
    k0is_assigned = []

    for i in range(len(ctis)):
        ai_at_gamma = ais[i].evaluate(gamma)
        ai_at_gamma_assigned = assign_to_circuit(Polynomial([ai_at_gamma]), p).coefficients[0]
        ais_at_gamma_assigned.append(ai_at_gamma_assigned)

        ct0i_at_gamma = ctis[i][0].evaluate(gamma)
        ct0i_at_gamma_assigned = assign_to_circuit(Polynomial([ct0i_at_gamma]), p).coefficients[0]
        ct0is_at_gamma_assigned.append(ct0i_at_gamma_assigned)

        qis_assigned.append(qis[i])

        k0i_assigned = assign_to_circuit(k0is[i], p).coefficients[0]
        k0is_assigned.append(k0i_assigned)

    cyclo_at_gamma = cyclo.evaluate(gamma)
    cyclo_at_gamma_assigned = assign_to_circuit(Polynomial([cyclo_at_gamma]), p).coefficients[0]

    # constraint. The coefficients of s should be in the range [-1, 0, 1]
    assert all(coeff in [-1, 0, 1] for coeff in s.coefficients)
    # After the circuit assignement, the coefficients of s_assigned must be in [0, 1, p - 1]
    assert all(coeff in [0, 1, p - 1] for coeff in s_assigned.coefficients)
    # To perform a range check with a smaller lookup table, we shift the coefficients of s_assigned to be in [0, 1, 2] (the shift operation is constrained inside the circuit)
    s_shifted = Polynomial([(coeff + 1) % p for coeff in s_assigned.coefficients])
    assert all(coeff in [0, 1, 2] for coeff in s_shifted.coefficients)

    # constraint. The coefficients of e should be in the range [-B, B] where B is the upper bound of the discrete Gaussian distribution
    b = int(discrete_gaussian.z_upper)
    assert all(coeff >= -b and coeff <= b for coeff in e.coefficients)
    # After the circuit assignement, the coefficients of e_assigned must be in [0, B] or [p - B, p - 1]
    assert all(coeff in range(0, b+1) or coeff in range(p - b, p) for coeff in e_assigned.coefficients)
    # To perform a range check with a smaller lookup table, we shift the coefficients of e_assigned to be in [0, 2B] (the shift operation is constrained inside the circuit)
    e_shifted = Polynomial([(coeff + b) % p for coeff in e_assigned.coefficients])
    assert all(coeff >= 0 and coeff <= 2*b for coeff in e_shifted.coefficients)

    # constraint. The coefficients of k1 should be in the range [-(t-1)/2, (t-1)/2]
    k1_bound = int((t - 1) / 2)
    assert all(coeff >= -k1_bound and coeff <= k1_bound for coeff in k1.coefficients)
    # After the circuit assignement, the coefficients of k1_assigned must be in [0, k1_bound] or [p - k1_bound, p - 1] 
    assert all(coeff in range(0, int(k1_bound) + 1) or coeff in range(p - int(k1_bound), p) for coeff in k1_assigned.coefficients)
    # To perform a range check with a smaller lookup table, we shift the coefficients of k1_assigned to be in [0, 2*k1_bound] (the shift operation is constrained inside the circuit)
    k1_shifted = Polynomial([(coeff + int(k1_bound)) % p for coeff in k1_assigned.coefficients])
    assert all(coeff >= 0 and coeff <= 2*k1_bound for coeff in k1_shifted.coefficients)

    for i in range(len(ctis)):
        # sanity check. The coefficients of ct0i should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in ct0is[i].coefficients)

        # sanity check. The coefficients of ai should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        assert all(coeff >= -bound and coeff <= bound for coeff in ais[i].coefficients)

        # sanity check. The coefficients of ai * s should be in the range $[-N \cdot \frac{q_i - 1}{2}, N \cdot \frac{q_i - 1}{2}]$
        bound = int((qis[i] - 1) / 2) * n
        res = ais[i] * s
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # sanity check. The coefficients of ai * s + e should be in the range $- (N \cdot \frac{q_i - 1}{2} + B), N \cdot \frac{q_i - 1}{2} + B]$
        bound = int((qis[i] - 1) / 2) * n + b
        res = ais[i] * s + e
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # constraint. The coefficients of r2i should be in the range [-(qi-1)/2, (qi-1)/2]
        r2i_bound = int((qis[i] - 1) / 2)
        r2_bounds.append(r2i_bound)
        assert all(coeff >= -r2i_bound and coeff <= r2i_bound for coeff in r2is[i].coefficients)
        # After the circuit assignement, the coefficients of r2i_assigned must be in [0, r2i_bound] or [p - r2i_bound, p - 1] 
        assert all(coeff in range(0, int(r2i_bound) + 1) or coeff in range(p - int(r2i_bound), p) for coeff in r2is_assigned[i].coefficients)
        # To perform a range check with a smaller lookup table, we shift the coefficients of r2i_assigned to be in [0, 2*r2i_bound] (the shift operation is constrained inside the circuit)
        r2i_shifted = Polynomial([(coeff + int(r2i_bound)) % p for coeff in r2is_assigned[i].coefficients])
        assert all(coeff >= 0 and coeff <= 2*r2i_bound for coeff in r2i_shifted.coefficients)

        # sanity check. The coefficients of r2i * cyclo should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        res = r2is[i] * cyclo
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # sanity check. The coefficients of k1 * k0i should be in the range $[-\frac{t - 1}{2} \cdot |K_{0,i}|, \frac{t - 1}{2} \cdot |K_{0,i}|]$
        bound = int((t - 1) / 2) * abs(k0is[i].coefficients[0])
        res = k1 * k0is[i]
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # sanity check. The coefficients of ct0i_hat (ai * s + e + k1 * k0i) should be in the range $[- (N \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|), N \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|]$
        bound = int((qis[i] - 1) / 2) * n + b + int((t - 1) / 2) * abs(k0is[i].coefficients[0])
        assert all(coeff >= -bound and coeff <= bound for coeff in ct0is_hat[i].coefficients)

        # sanity check. The coefficients of ct0i - ct0i_hat should be in the range $ [- ((N+1) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|), (N+1) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|]$
        bound = int((qis[i] - 1) / 2) * (n + 1) + b + int((t - 1) / 2) * abs(k0is[i].coefficients[0])
        sub = ct0is[i] + (Polynomial([-1]) * ct0is_hat[i])
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # sanity check. The coefficients of ct0i - ct0i_hat - r2i * cyclo should be in the range $[- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|), (N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|]$
        bound = ((qis[i] - 1) / 2) * (n + 2) + b + ((t - 1) / 2) * abs(k0is[i].coefficients[0])
        sub = ct0is[i]  + (Polynomial([-1]) * ct0is_hat[i]) + Polynomial([-1]) * (r2is[i] * cyclo)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # constraint. The coefficients of (ct0i - ct0i_hat - r2i * cyclo) / qi = r1i should be in the range $[\frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|)}{q_i}, \frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|}{q_i}]$
        r1i_bound = (int((qis[i] - 1) / 2) * (n + 2) + b + int((t - 1) / 2) * abs(k0is[i].coefficients[0])) / qis[i]
        # round bound to the nearest integer
        r1i_bound = int(r1i_bound)
        r1_bounds.append(r1i_bound)
        assert all(coeff >= -r1i_bound and coeff <= r1i_bound for coeff in r1is[i].coefficients)
        # After the circuit assignement, the coefficients of r1i_assigned must be in [0, r1i_bound] or [p - r1i_bound, p - 1]
        assert all(coeff in range(0, int(r1i_bound) + 1) or coeff in range(p - int(r1i_bound), p) for coeff in r1is_assigned[0].coefficients)
        # To perform a range check with a smaller lookup table, we shift the coefficients of r1i_assigned to be in [0, 2*r1i_bound] (the shift operation is constrained inside the circuit)
        r1i_shifted = Polynomial([(coeff + int(r1i_bound)) % p for coeff in r1is_assigned[0].coefficients])
        assert all(coeff >= 0 and coeff <= 2*r1i_bound for coeff in r1i_shifted.coefficients)

        s_at_gamma_assigned = s_assigned.evaluate(gamma)
        e_at_gamma_assigned = e_assigned.evaluate(gamma)
        k1_at_gamma_assigned = k1_assigned.evaluate(gamma)

        for i in range(len(ctis)):
            r1i_gamma_assigned = r1is_assigned[i].evaluate(gamma)
            r2i_gamma_assigned = r2is_assigned[i].evaluate(gamma)
            lhs = ct0is_at_gamma_assigned[i]
            rhs = (ais_at_gamma_assigned[i] * s_at_gamma_assigned + e_at_gamma_assigned + (k1_at_gamma_assigned * k0is_assigned[i]) + (r1i_gamma_assigned * qis_assigned[i]) + (r2i_gamma_assigned * cyclo_at_gamma_assigned))
            
            assert lhs % p == rhs % p

        '''
        VERIFICATION PHASE
        '''

        cyclo_at_gamma = cyclo.evaluate(gamma)
        cyclo_at_gamma_assigned_expected = assign_to_circuit(Polynomial([cyclo_at_gamma]), p).coefficients[0]
        assert cyclo_at_gamma_assigned == cyclo_at_gamma_assigned_expected

        for i in range(len(ctis)):
            ai_gamma = ais[i].evaluate(gamma)
            ai_gamma_assigned_expected = assign_to_circuit(Polynomial([ai_gamma]), p).coefficients[0]
            assert ais_at_gamma_assigned[i] == ai_gamma_assigned_expected

            ct0i_gamma = ctis[i][0].evaluate(gamma)
            ct0i_gamma_assigned_expected = assign_to_circuit(Polynomial([ct0i_gamma]), p).coefficients[0]
            assert ct0is_at_gamma_assigned[i] == ct0i_gamma_assigned_expected

            assert qis[i] == qis_assigned[i]

            k0i_assigned_expected = assign_to_circuit(k0is[i], p).coefficients[0]
            assert k0is_assigned[i] == k0i_assigned_expected

    # ais and ct0is need to be parsed such that their coefficients are in the range [0, p - 1]
    # we don't call them assigned because they are never assigned to the circuit
    ais_in_p = [assign_to_circuit(ai, p) for ai in ais]
    ct0is_in_p = [assign_to_circuit(ct0i, p) for ct0i in ct0is]
    
    # Parse the inputs into a JSON format such this can be used as input for the (real) circuit
    json_input = {
        "s": [str(coef) for coef in s_assigned.coefficients],
        "e": [str(coef) for coef in e_assigned.coefficients],
        "qis": [str(qi) for qi in qis_assigned],
        "k1": [str(coef) for coef in k1_assigned.coefficients],
        "k0is": [str(k0i) for k0i in k0is_assigned],
        "r2is": [[str(coef) for coef in r2i.coefficients] for r2i in r2is_assigned],
        "r1is": [[str(coef) for coef in r1i.coefficients] for r1i in r1is_assigned],
        "ais": [[str(coef) for coef in ai_in_p.coefficients] for ai_in_p in ais_in_p],
        "ct0is": [[str(coef) for coef in ct0i_in_p.coefficients] for ct0i_in_p in ct0is_in_p],
    }

    # write the inputs to a json file
    with open(args.output, 'w') as f:
        json.dump(json_input, f)

    print(f"/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.")
    print(f"pub const N: usize = {n};")
    print(f"/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2")
    print(f"pub const E_BOUND: u64 = {b};")
    print(f"/// The coefficients of the plynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.")
    print(f"pub const S_BOUND: u64 = {1};")
    print(f"/// The coefficients of the polynomials `r1is` should exist in the interval `[-R1_BOUND[i], R1_BOUND[i]]` where `R1_BOUND[i]` is equal to `(qi-1)/2`")
    print(f"pub const R1_BOUNDS: [u64; {len(r1_bounds)}] = [{', '.join(map(str, r1_bounds))}];")
    print(f"/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to $\\frac{{(N+2) \\cdot \\frac{{q_i - 1}}{{2}} + B + \\frac{{t - 1}}{{2}} \\cdot |K_{{0,i}}|}}{{q_i}}$")
    print(f"pub const R2_BOUNDS: [u64; {len(r2_bounds)}] = [{', '.join(map(str, r2_bounds))}];")
    print(f"/// The coefficients of `k1` should exist in the interval `[-K1_BOUND, K1_BOUND]` where `K1_BOUND` is equal to `(t-1)/2`")
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
