from bfv.polynomial import Polynomial
import math

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

def count_advice_cells_needed_for_poly_range_check(poly: Polynomial, bound: int, lookup_bits: int) -> int:
    '''
    This function takes a polynomial and a bound and returns the number of advice cells needed for a range check
    `poly` is the polynomial to be checked
    `bound` is the upper bound for the range check
    `lookup_bits` is the number of bits used for the lookup table
    '''

    count = 0

    # 4 advice cells for each coefficient needed for the shift addition operation
    count += 4 * len(poly.coefficients)

    # further advice cells for range check inside `check_less_than_safe`
    bound_bits = bound.bit_length()
    range_bits = math.ceil(bound_bits / lookup_bits) * lookup_bits
    num_limbs = math.ceil(range_bits / lookup_bits)

    if num_limbs > 1:
        # 1 + (3 * (num_limbs - 1)) advice cells
        count += (1 + (3 * (num_limbs - 1))) * len(poly.coefficients)
    else:
        # count is not updated if num_limbs is 1
        pass

    # 7 advice cells required for the `check_less_than` constraint inside `check_less_than_safe`
    count += 7 * len(poly.coefficients)

    # the check_less_than_advice_cells constraint also performs a range check on the check_cell in range_bits
    if num_limbs > 1:
        # 1 + (3 * (num_limbs - 1)) advice cells for the check_less_than_advice_cells constraint
        count += (1 + (3 * (num_limbs - 1))) * len(poly.coefficients)
    else:
        # count is not updated if num_limbs is 1
        pass

    return count

def print_advice_cells_info(total_advice_cell_count, phase_0_count, phase_1_assignment_count, phase_1_range_check_count, phase_1_eval_at_gamma_count, phase_1_encryption_constraint_count):
    print("Halo2 Circuit Profile:")
    print(f"Total Advice Cells Needed: {total_advice_cell_count}")
    
    print("\nPhase 0 - Assignment:")
    print(f" - Count: {phase_0_count}, Percentage: {(phase_0_count / total_advice_cell_count) * 100:.2f}%")

    print("\nPhase 1 - Assignment:")
    print(f" - Count: {phase_1_assignment_count}, Percentage: {(phase_1_assignment_count / total_advice_cell_count) * 100:.2f}%")

    print("\nPhase 1 - Range Check:")
    print(f" - Count: {phase_1_range_check_count}, Percentage: {(phase_1_range_check_count / total_advice_cell_count) * 100:.2f}%")

    print("\nPhase 1 - Evaluation at Gamma Constraint:")
    print(f" - Count: {phase_1_eval_at_gamma_count}, Percentage: {(phase_1_eval_at_gamma_count / total_advice_cell_count) * 100:.2f}%")

    print("\nPhase 1 - Correct Encryption Constraint:")
    print(f" - Count: {phase_1_encryption_constraint_count}, Percentage: {(phase_1_encryption_constraint_count / total_advice_cell_count) * 100:.2f}%")
