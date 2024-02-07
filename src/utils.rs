use halo2_base::utils::ScalarField;

/// Evaluate a polynomial at a point in the field. The polynomial is given in coefficients form starting from the highest degree term.
pub fn evaluate_poly<F: ScalarField>(coeffs: &[F], x: F) -> F {
    let mut acc = coeffs[0];
    for coeff in coeffs.iter().skip(1) {
        acc = acc * x + coeff;
    }
    acc
}
