use axiom_eth::rlc::{
    chip::RlcChip,
    circuit::{builder::RlcCircuitBuilder, instructions::RlcCircuitInstructions, RlcCircuitParams},
};
use halo2_base::{
    gates::{circuit::BaseCircuitParams, GateInstructions, RangeChip, RangeInstructions},
    utils::ScalarField,
    AssignedValue,
    QuantumCell::Constant,
};
use serde::Deserialize;

use crate::utils::evaluate_poly;

const K: usize = 18;

/// Helper function to define the parameters of the RlcCircuit
pub fn test_params() -> RlcCircuitParams {
    RlcCircuitParams {
        base: BaseCircuitParams {
            k: K,
            num_advice_per_phase: vec![1, 1],
            num_fixed: 1,
            num_lookup_advice_per_phase: vec![1, 0],
            lookup_bits: Some(8),
            num_instance_columns: 0,
        },
        num_rlc_columns: 1,
    }
}

/// `BfvSkEncryptionCircuit` is a circuit that checks the correct formation of a ciphertext resulting from BFV secret key encryption
///
/// TODO: update the documentation
#[derive(Deserialize)]
pub struct BfvSkEncryptionCircuit {
    s: Vec<String>,
    e: Vec<String>,
    k1: Vec<String>,
    r2s: Vec<Vec<String>>,
    r1s: Vec<Vec<String>>,
}

/// `N` is the degree of the cyclotomic polynomial defining the rings Rq and Rt.
/// We define it as a constant because it should be known at compile time
const N: usize = 1024;

/// `B` is the upper bound of the gaussian distribution used to sample the error polynomial `e`
/// We define it as a constant because it should be known at compile time
const B: u64 = 19;

/// `qis` is a vector representing the moduli used for the CRT decomposition of the modulus `q`, which is the modulus of the ring Rq
/// We define it as a constant because it should be known at compile time
const QIS: [u64; 15] = [
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
];

/// `t` is the plaintext modulus used for the BFV encryption scheme
/// We define it as a constant because it should be known at compile time
const T: u64 = 65537;

const R1_0_BOUND: u64 = 1321;

/// Payload returned by the first phase of the circuit to be reused in the second phase
pub struct Payload<F: ScalarField> {
    s: Vec<String>,
    e: Vec<String>,
    k1: Vec<String>,
    r2_0: Vec<String>,
    r1_0: Vec<String>,
    s_assigned: Vec<AssignedValue<F>>,
    e_assigned: Vec<AssignedValue<F>>,
    k1_assigned: Vec<AssignedValue<F>>,
    r2_0_assigned: Vec<AssignedValue<F>>,
    r1_0_assigned: Vec<AssignedValue<F>>,
}

impl<F: ScalarField> RlcCircuitInstructions<F> for BfvSkEncryptionCircuit {
    type FirstPhasePayload = Payload<F>;

    // Phase 0
    // * Assign the secret polynomials from the matrix S to the circuit
    // * Assign the public inputs to the circuit
    // * Perform range checks on the coefficients of the secret polynomials
    fn virtual_assign_phase0(
        &self,
        builder: &mut RlcCircuitBuilder<F>,
        range: &RangeChip<F>,
    ) -> Self::FirstPhasePayload {
        let ctx = builder.base.main(0);

        let mut s_assigned = vec![];
        let mut e_assigned = vec![];
        let mut k1_assigned = vec![];
        let mut r2_0_assigned = vec![];
        let mut r1_0_assigned = vec![];

        // assign polynomial s to the witness

        for i in 0..N {
            let val = F::from_str_vartime(&self.s[i]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            s_assigned.push(coeff_assigned);
        }

        // assign polynomial e to the witness

        for i in 0..N {
            let val = F::from_str_vartime(&self.e[i]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            e_assigned.push(coeff_assigned);
        }

        // assign polynomial k1 to the witness

        for i in 0..N {
            let val = F::from_str_vartime(&self.k1[i]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            k1_assigned.push(coeff_assigned);
        }

        // assign polynomial r2[0] to the witness. Note that the degree of R2 is N - 2

        for i in 0..(N - 1) {
            let val = F::from_str_vartime(&self.r2s[0][i]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            r2_0_assigned.push(coeff_assigned);
        }

        // assign polynomial r1[0] to the witness. Note that the degree of R1 is 2N - 2
        for i in 0..(2 * N - 1) {
            let val = F::from_str_vartime(&self.r1s[0][i]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            r1_0_assigned.push(coeff_assigned);
        }

        // perform range check on the coefficients of `s`

        let constant_one = Constant(F::from(1));
        let bound_s = 3;

        for i in 0..N {
            let normalized_coeff = range.gate().add(ctx, s_assigned[i], constant_one);
            // TODO: we can make this a bit more efficient by not having to reassign b as a constant every time
            range.check_less_than_safe(ctx, normalized_coeff, bound_s);
        }

        // perform range check on the coefficients of `e`

        let e_bound = B;
        let e_bound_constant = Constant(F::from(B));

        for i in 0..N {
            let normalized_coeff = range.gate().add(ctx, e_assigned[i], e_bound_constant);
            // TODO: we can make this a bit more efficient by not having to reassign b as a constant every time
            range.check_less_than_safe(ctx, normalized_coeff, 2 * e_bound);
        }

        // perform range check on the coefficients of `k1`

        let k1_bound = (T - 1) / 2;
        let k1_bound_constant = Constant(F::from(k1_bound));

        for i in 0..N {
            let normalized_coeff = range.gate().add(ctx, k1_assigned[i], k1_bound_constant);
            range.check_less_than_safe(ctx, normalized_coeff, 2 * k1_bound);
        }

        // perform range check on the coefficients of `r2[0]`. Note that the degree of R2 is N - 2

        let r2_0_bound = QIS[0] - 1 / 2;
        let r2_0_bound_constant = Constant(F::from(r2_0_bound));

        for i in 0..(N - 1) {
            let normalized_coeff = range.gate().add(ctx, r2_0_assigned[i], r2_0_bound_constant);
            range.check_less_than_safe(ctx, normalized_coeff, 2 * r2_0_bound);
        }

        // perform range check on the coefficients of `r1[0]`. Note that the degree of R1 is 2N - 2

        let r1_0_bound = R1_0_BOUND;
        let r1_0_bound_constant = Constant(F::from(r1_0_bound));

        for i in 0..(2 * N - 1) {
            let normalized_coeff = range.gate().add(ctx, r1_0_assigned[i], r1_0_bound_constant);
            range.check_less_than_safe(ctx, normalized_coeff, 2 * r1_0_bound);
        }

        Payload {
            s: self.s.clone(),
            e: self.e.clone(),
            k1: self.k1.clone(),
            r2_0: self.r2s[0].clone(),
            r1_0: self.r1s[0].clone(),
            s_assigned,
            e_assigned,
            k1_assigned,
            r2_0_assigned,
            r1_0_assigned,
        }
    }

    // Phase 1
    // * Evaluate polynomials of matrix S at the point `gamma`
    // * Assert that `a(gamma) * b(gamma) = c(gamma)`
    fn virtual_assign_phase1(
        builder: &mut RlcCircuitBuilder<F>,
        _: &RangeChip<F>,
        rlc: &RlcChip<F>,
        payload: Self::FirstPhasePayload,
    ) {
        let Payload {
            s,
            e,
            k1,
            r2_0,
            r1_0,
            s_assigned,
            e_assigned,
            k1_assigned,
            r2_0_assigned,
            r1_0_assigned,
        } = payload;

        let (ctx_gate, ctx_rlc) = builder.rlc_ctx_pair();
        let gamma = *rlc.gamma();

        let rlc_trace_s = rlc.compute_rlc_fixed_len(ctx_rlc, s_assigned);
        let rlc_trace_e = rlc.compute_rlc_fixed_len(ctx_rlc, e_assigned);
        let rlc_trace_k1 = rlc.compute_rlc_fixed_len(ctx_rlc, k1_assigned);
        let rlc_trace_r2_0 = rlc.compute_rlc_fixed_len(ctx_rlc, r2_0_assigned);
        let rlc_trace_r1_0 = rlc.compute_rlc_fixed_len(ctx_rlc, r1_0_assigned);

        // assert the correctness of the RLC
        let s_poly = s
            .iter()
            .map(|s| F::from_str_vartime(s).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_s_eval_at_gamma = evaluate_poly(&s_poly, gamma);

        assert_eq!(rlc_trace_s.rlc_val.value(), &expected_poly_s_eval_at_gamma);

        let e_poly = e
            .iter()
            .map(|e| F::from_str_vartime(e).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_e_eval_at_gamma = evaluate_poly(&e_poly, gamma);

        assert_eq!(rlc_trace_e.rlc_val.value(), &expected_poly_e_eval_at_gamma);

        let k1_poly = k1
            .iter()
            .map(|k1| F::from_str_vartime(k1).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_k1_eval_at_gamma = evaluate_poly(&k1_poly, gamma);

        assert_eq!(
            rlc_trace_k1.rlc_val.value(),
            &expected_poly_k1_eval_at_gamma
        );

        let r2_0_poly = r2_0
            .iter()
            .map(|r2_0| F::from_str_vartime(r2_0).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_r2_0_eval_at_gamma = evaluate_poly(&r2_0_poly, gamma);

        assert_eq!(
            rlc_trace_r2_0.rlc_val.value(),
            &expected_poly_r2_0_eval_at_gamma
        );

        let r1_0_poly = r1_0
            .iter()
            .map(|r1_0| F::from_str_vartime(r1_0).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_r1_0_eval_at_gamma = evaluate_poly(&r1_0_poly, gamma);

        assert_eq!(
            rlc_trace_r1_0.rlc_val.value(),
            &expected_poly_r1_0_eval_at_gamma
        );
    }
}

#[cfg(test)]
mod test {

    use std::{fs::File, io::Read};

    use super::{test_params, K};
    use crate::sk_encryption_circuit::BfvSkEncryptionCircuit;
    use axiom_eth::rlc::{circuit::builder::RlcCircuitBuilder, utils::executor::RlcExecutor};
    use halo2_base::{
        gates::circuit::CircuitBuilderStage,
        halo2_proofs::{dev::MockProver, halo2curves::bn256::Fr},
    };

    #[test]
    pub fn test_sk_enc_valid() {
        // 1. Define the inputs of the circuit
        let file_path = "src/data/sk_enc_input.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let circuit = serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

        // 2. Build the circuit for MockProver
        let params = test_params();
        let mut mock_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0).use_params(params);
        mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let circuit = RlcExecutor::new(mock_builder, circuit);

        // 3. Run the mock prover. The circuit should be satisfied
        MockProver::run(K as u32, &circuit, vec![])
            .unwrap()
            .assert_satisfied();
    }
}
