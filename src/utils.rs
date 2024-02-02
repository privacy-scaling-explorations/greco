use axiom_eth::halo2curves::bn256::Fr;
use halo2_base::utils::ScalarField;

/// Evaluate a polynomial at a point in the field. The polynomial is given in coefficients form starting from the highest degree term.
pub fn evaluate_poly<F: ScalarField>(coeffs: &[F], x: F) -> F {
    let mut acc = coeffs[0];
    for coeff in coeffs.iter().skip(1) {
        acc = acc * x + coeff;
    }
    acc
}

/// Perform polynomial multiplication between `a` and `b`
/// The polynomials are given in coefficients form starting from the highest degree term.
pub fn mul(a: Vec<Fr>, b: Vec<Fr>) -> Vec<Fr> {
    let deg_a = a.len() - 1;
    let deg_b = b.len() - 1;
    assert_eq!(deg_a, deg_b);

    let deg_c = deg_a + deg_b;

    // initialize the output polynomial with zeroes
    let mut c = vec![Fr::from(0_u64); deg_c + 1];

    // perform polynomial multiplication
    for i in 0..=deg_a {
        for j in 0..=deg_b {
            c[i + j] += a[i] * b[j]
        }
    }

    c
}
