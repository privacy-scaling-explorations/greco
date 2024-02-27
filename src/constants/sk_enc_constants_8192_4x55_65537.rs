/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.
pub const N: usize = 8192;
/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2
pub const E_BOUND: u64 = 19;
/// The coefficients of the polynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.
pub const S_BOUND: u64 = 1;
/// The coefficients of the polynomials `r1is` should exist in the interval `[-R1_BOUND[i], R1_BOUND[i]]` where `R1_BOUND[i]` is equal to `(qi-1)/2`
pub const R1_BOUNDS: [u64; 4] = [30180, 5459, 4259, 13207];
/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to $\frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|}{q_i}$
pub const R2_BOUNDS: [u64; 4] = [16838097769646402, 16838097769646404, 16838097769646404, 16838097769646404];
/// The coefficients of `k1` should exist in the interval `[-K1_BOUND, K1_BOUND]` where `K1_BOUND` is equal to `(t-1)/2`
pub const K1_BOUND: u64 = 32768;
/// List of scalars `qis` such that `qis[i]` is the modulus of the i-th CRT basis of `q` (ciphertext space modulus)
pub const QIS: [&str; 4] = ["33676195539292805", "33676195539292807", "33676195539292809", "33676195539292811"];
/// List of scalars `k0is` such that `k0i[i]` is equal to the negative of the multiplicative inverses of t mod qi.
pub const K0IS: [&str; 4] = ["26806019389021282", "1400755436472713", "167001290115052", "9363376949158393"];
