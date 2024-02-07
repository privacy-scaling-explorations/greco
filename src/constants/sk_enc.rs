/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.
pub const N: usize = 1024;
/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2
pub const E_BOUND: u64 = 19;
/// The coefficients of the plynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.
pub const S_BOUND: u64 = 1;
/// The coefficients of the polynomials `r1is` should exist in the interval `[-R1_BOUND[i], R1_BOUND[i]]` where `R1_BOUND[i]` is equal to `(qi-1)/2`
pub const R1_BOUNDS: [u64; 15] = [
    1321, 12139, 1692, 1530, 19009, 17587, 3417, 15539, 24450, 19013, 24041, 5934, 31437, 16662,
    15909,
];
/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to $\frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|}{q_i}$
pub const R2_BOUNDS: [u64; 15] = [
    576460752303292416,
    576460752299360256,
    576460752298508288,
    576460752297984000,
    576460752297820160,
    576460752296706048,
    576460752296411136,
    576460752296214528,
    576460752294969344,
    576460752293265408,
    576460752292773888,
    576460752291823616,
    576460752290938880,
    576460752290709504,
    576460752290447360,
];
/// The coefficients of `k1` should exist in the interval `[-K1_BOUND, K1_BOUND]` where `K1_BOUND` is equal to `(t-1)/2`
pub const K1_BOUND: u64 = 32768;
