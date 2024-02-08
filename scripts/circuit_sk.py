from bfv.crt import CRTModuli
from bfv.bfv import BFVCrt
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div
from bfv.utils import mod_inverse
from random import randint
import copy
from utils import assign_to_circuit, count_advice_cells_needed_for_poly_range_check, print_advice_cells_info
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

    # initiate counters for the number of advice cells needed for each constraint phase
    phase_0_assignment_advice_cell_count = 0
    phase_1_assignment_advice_cell_count = 0
    phase_1_range_check_advice_cell_count = 0
    phase_1_eval_at_gamma_constraint_advice_cell_count = 0
    phase_1_encryption_constraint_advice_cell_count = 0

    '''
    CIRCUIT - PHASE 0 - ASSIGNMENT
    '''

    # Every assigned value must be an element of the field Zp. Negative coefficients `-z` are assigned as `p - z`
    s_assigned = assign_to_circuit(s, p)    
    e_assigned = assign_to_circuit(e, p)
    k1_assigned = assign_to_circuit(k1, p)

    phase_0_assignment_advice_cell_count += len(s_assigned.coefficients)
    phase_0_assignment_advice_cell_count += len(e_assigned.coefficients)
    phase_0_assignment_advice_cell_count += len(k1_assigned.coefficients)

    r1is_assigned = []
    r2is_assigned = []

    for i in range(len(ctis)):
        r1i_assigned = assign_to_circuit(r1is[i], p)
        r2i_assigned = assign_to_circuit(r2is[i], p)
        r1is_assigned.append(r1i_assigned)
        r2is_assigned.append(r2i_assigned)

        phase_0_assignment_advice_cell_count += len(r1i_assigned.coefficients)
        phase_0_assignment_advice_cell_count += len(r2i_assigned.coefficients)

    # For the sake of simplicity, we generate a random challenge here
    gamma = randint(0, 1000)

    '''
    CIRCUIT - PHASE 1 - ASSIGNMENT
    '''

    # Every assigned value must be an element of the field Zp. Negative coefficients `-z` are assigned as `p - z`
    ais_at_gamma_assigned = []
    ct0is_at_gamma_assigned = []
    qi_constants = []
    k0i_constants = []

    for i in range(len(ctis)):
        ai_at_gamma = ais[i].evaluate(gamma)
        ai_at_gamma_assigned = assign_to_circuit(Polynomial([ai_at_gamma]), p).coefficients[0]
        ais_at_gamma_assigned.append(ai_at_gamma_assigned)

        phase_1_assignment_advice_cell_count += 1

        ct0i_at_gamma = ctis[i][0].evaluate(gamma)
        ct0i_at_gamma_assigned = assign_to_circuit(Polynomial([ct0i_at_gamma]), p).coefficients[0]
        ct0is_at_gamma_assigned.append(ct0i_at_gamma_assigned)

        phase_1_assignment_advice_cell_count += 1

        qi_constants.append(qis[i])

        k0i_constant = assign_to_circuit(k0is[i], p).coefficients[0]
        k0i_constants.append(k0i_constant)

    cyclo_at_gamma = cyclo.evaluate(gamma)
    cyclo_at_gamma_assigned = assign_to_circuit(Polynomial([cyclo_at_gamma]), p).coefficients[0]

    phase_1_assignment_advice_cell_count += 1

    '''
    CIRCUIT - PHASE 1 - RANGE CHECK
    '''

    lookup_bits = 8

    s_bound = 1
    # constraint. The coefficients of s should be in the range [-1, 0, 1]
    assert all(coeff >= -s_bound and coeff <= s_bound for coeff in s.coefficients)
    # After the circuit assignement, the coefficients of s_assigned must be in [0, 1, p - 1]
    assert all(coeff in range(0, s_bound+1) or coeff in range(p - s_bound, p) for coeff in s_assigned.coefficients)
    # To perform a range check with a smaller lookup table, we shift the coefficients of s_assigned to be in [0, 1, 2] (the shift operation is constrained inside the circuit)
    s_shifted = Polynomial([(coeff + 1) % p for coeff in s_assigned.coefficients])
    assert all(coeff >= 0 and coeff <= 2*s_bound for coeff in s_shifted.coefficients)

    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(s_assigned, 2*s_bound + 1, lookup_bits)
    
    # constraint. The coefficients of e should be in the range [-B, B] where B is the upper bound of the discrete Gaussian distribution
    b = int(discrete_gaussian.z_upper)
    assert all(coeff >= -b and coeff <= b for coeff in e.coefficients)
    # After the circuit assignement, the coefficients of e_assigned must be in [0, B] or [p - B, p - 1]
    assert all(coeff in range(0, b+1) or coeff in range(p - b, p) for coeff in e_assigned.coefficients)
    # To perform a range check with a smaller lookup table, we shift the coefficients of e_assigned to be in [0, 2B] (the shift operation is constrained inside the circuit)
    e_shifted = Polynomial([(coeff + b) % p for coeff in e_assigned.coefficients])
    assert all(coeff >= 0 and coeff <= 2*b for coeff in e_shifted.coefficients)

    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(e_assigned, 2*b + 1, lookup_bits)
    
    # constraint. The coefficients of k1 should be in the range [-(t-1)/2, (t-1)/2]
    k1_bound = int((t - 1) / 2)
    assert all(coeff >= -k1_bound and coeff <= k1_bound for coeff in k1.coefficients)
    # After the circuit assignement, the coefficients of k1_assigned must be in [0, k1_bound] or [p - k1_bound, p - 1] 
    assert all(coeff in range(0, int(k1_bound) + 1) or coeff in range(p - int(k1_bound), p) for coeff in k1_assigned.coefficients)
    # To perform a range check with a smaller lookup table, we shift the coefficients of k1_assigned to be in [0, 2*k1_bound] (the shift operation is constrained inside the circuit)
    k1_shifted = Polynomial([(coeff + int(k1_bound)) % p for coeff in k1_assigned.coefficients])
    assert all(coeff >= 0 and coeff <= 2*k1_bound for coeff in k1_shifted.coefficients)

    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(k1_assigned, 2*k1_bound + 1, lookup_bits)

    s_at_gamma_assigned = s_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(s_assigned.coefficients) * 2 - 1

    e_at_gamma_assigned = e_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(e_assigned.coefficients) * 2 - 1

    k1_at_gamma_assigned = k1_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(k1_assigned.coefficients) * 2 - 1

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

        phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(r2is_assigned[i], 2*r2i_bound + 1, lookup_bits)

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
        assert all(coeff in range(0, int(r1i_bound) + 1) or coeff in range(p - int(r1i_bound), p) for coeff in r1is_assigned[i].coefficients)
        # To perform a range check with a smaller lookup table, we shift the coefficients of r1i_assigned to be in [0, 2*r1i_bound] (the shift operation is constrained inside the circuit)
        r1i_shifted = Polynomial([(coeff + int(r1i_bound)) % p for coeff in r1is_assigned[i].coefficients])
        assert all(coeff >= 0 and coeff <= 2*r1i_bound for coeff in r1i_shifted.coefficients)

        phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(r1is_assigned[i], 2*r1i_bound + 1, lookup_bits)

        '''
        CIRCUIT - PHASE 1 - EVALUATION AT GAMMA CONSTRAINT
        '''

        r1i_gamma_assigned = r1is_assigned[i].evaluate(gamma)
        r2i_gamma_assigned = r2is_assigned[i].evaluate(gamma)

        phase_1_eval_at_gamma_constraint_advice_cell_count += len(r1is_assigned[i].coefficients) * 2 - 1
        phase_1_eval_at_gamma_constraint_advice_cell_count += len(r2is_assigned[i].coefficients) * 2 - 1

        '''
        CIRCUIT - PHASE 1 - CORRECT ENCRYPTION CONSTRAINT
        '''

        lhs = ct0is_at_gamma_assigned[i]
        rhs = (ais_at_gamma_assigned[i] * s_at_gamma_assigned + e_at_gamma_assigned + (k1_at_gamma_assigned * k0i_constants[i]) + (r1i_gamma_assigned * qi_constants[i]) + (r2i_gamma_assigned * cyclo_at_gamma_assigned))
        phase_1_encryption_constraint_advice_cell_count += 16

        assert lhs % p == rhs % p

        '''
        VERIFICATION PHASE
        '''

        cyclo_at_gamma = cyclo.evaluate(gamma)
        cyclo_at_gamma_assigned_expected = assign_to_circuit(Polynomial([cyclo_at_gamma]), p).coefficients[0]
        assert cyclo_at_gamma_assigned == cyclo_at_gamma_assigned_expected

        ai_gamma = ais[i].evaluate(gamma)
        ai_gamma_assigned_expected = assign_to_circuit(Polynomial([ai_gamma]), p).coefficients[0]
        assert ais_at_gamma_assigned[i] == ai_gamma_assigned_expected

        ct0i_gamma = ctis[i][0].evaluate(gamma)
        ct0i_gamma_assigned_expected = assign_to_circuit(Polynomial([ct0i_gamma]), p).coefficients[0]
        assert ct0is_at_gamma_assigned[i] == ct0i_gamma_assigned_expected

        assert qis[i] == qi_constants[i]

        k0i_assigned_expected = assign_to_circuit(k0is[i], p).coefficients[0]
        assert k0i_constants[i] == k0i_assigned_expected

    total_advice_cell_count = phase_0_assignment_advice_cell_count + phase_1_assignment_advice_cell_count + phase_1_range_check_advice_cell_count + phase_1_eval_at_gamma_constraint_advice_cell_count + phase_1_encryption_constraint_advice_cell_count

    print_advice_cells_info(total_advice_cell_count, phase_0_assignment_advice_cell_count, phase_1_assignment_advice_cell_count, phase_1_range_check_advice_cell_count, phase_1_eval_at_gamma_constraint_advice_cell_count, phase_1_encryption_constraint_advice_cell_count)
    # ais and ct0is need to be parsed such that their coefficients are in the range [0, p - 1]
    # we don't call them assigned because they are never assigned to the circuit
    ais_in_p = [assign_to_circuit(ai, p) for ai in ais]
    ct0is_in_p = [assign_to_circuit(ct0i, p) for ct0i in ct0is]
    
    # Parse the inputs into a JSON format such this can be used as input for the (real) circuit
    json_input = {
        "s": [str(coef) for coef in s_assigned.coefficients],
        "e": [str(coef) for coef in e_assigned.coefficients],
        "k1": [str(coef) for coef in k1_assigned.coefficients],
        "r2is": [[str(coef) for coef in r2i.coefficients] for r2i in r2is_assigned],
        "r1is": [[str(coef) for coef in r1i.coefficients] for r1i in r1is_assigned],
        "ais": [[str(coef) for coef in ai_in_p.coefficients] for ai_in_p in ais_in_p],
        "ct0is": [[str(coef) for coef in ct0i_in_p.coefficients] for ct0i_in_p in ct0is_in_p],
    }

    with open(args.output_constants, 'w') as f:
        f.write(f"/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.\n")
        f.write(f"pub const N: usize = {n};\n")
        f.write(f"/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2\n")
        f.write(f"pub const E_BOUND: u64 = {b};\n")
        f.write(f"/// The coefficients of the polynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.\n")
        f.write(f"pub const S_BOUND: u64 = {1};\n")
        f.write(f"/// The coefficients of the polynomials `r1is` should exist in the interval `[-R1_BOUND[i], R1_BOUND[i]]` where `R1_BOUND[i]` is equal to `(qi-1)/2`\n")
        f.write(f"pub const R1_BOUNDS: [u64; {len(r1_bounds)}] = [{', '.join(map(str, r1_bounds))}];\n")
        f.write(f"/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to $\\frac{{(N+2) \\cdot \\frac{{q_i - 1}}{{2}} + B + \\frac{{t - 1}}{{2}} \\cdot |K_{{0,i}}|}}{{q_i}}$\n")
        f.write(f"pub const R2_BOUNDS: [u64; {len(r2_bounds)}] = [{', '.join(map(str, r2_bounds))}];\n")
        f.write(f"/// The coefficients of `k1` should exist in the interval `[-K1_BOUND, K1_BOUND]` where `K1_BOUND` is equal to `(t-1)/2`\n")
        f.write(f"pub const K1_BOUND: u64 = {k1_bound};\n")
        qis_str = ', '.join(f'"{q}"' for q in qi_constants)
        f.write(f"/// List of scalars `qis` such that `qis[i]` is the modulus of the i-th CRT basis of `q` (ciphertext space modulus)\n")
        f.write(f"pub const QIS: [&str; {len(qi_constants)}] = [{qis_str}];\n")
        k0is_str = ', '.join(f'"{k0i}"' for k0i in k0i_constants)
        f.write(f"/// List of scalars `k0is` such that `k0i[i]` is equal to the negative of the multiplicative inverses of t mod qi.\n")
        f.write(f"pub const K0IS: [&str; {len(k0i_constants)}] = [{k0is_str}];\n")

    with open(args.output_input, 'w') as f:
        json.dump(json_input, f)

    # Initialize a structure to hold polynomials with zero coefficients. This will be used at key generation.
    json_input_zeroes = {
        "s": ["0" for _ in s_assigned.coefficients],
        "e": ["0" for _ in e_assigned.coefficients],
        "k1": ["0" for _ in k1_assigned.coefficients],
        "r2is": [["0" for _ in r2i.coefficients] for r2i in r2is_assigned],
        "r1is": [["0" for _ in r1i.coefficients] for r1i in r1is_assigned],
        "ais": [["0" for _ in ai_in_p.coefficients] for ai_in_p in ais_in_p],
        "ct0is": [["0" for _ in ct0i_in_p.coefficients] for ct0i_in_p in ct0is_in_p],
    }

    original_output_path = args.output_input
    path_parts = original_output_path.rsplit('.', 1)
    zeroed_output_path = f"{path_parts[0]}_zeroes.{path_parts[1]}"

    with open(zeroed_output_path, 'w') as f:
        json.dump(json_input_zeroes, f)


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
        "-output_input", type=str, required=True, help="Path for the output json file containing the inputs for the circuit."
    )
    parser.add_argument(
        "-output_constants", type=str, required=True, help="Path for the output rust file containing the constants for the circuit."
    )

    args = parser.parse_args()
    main(args)
