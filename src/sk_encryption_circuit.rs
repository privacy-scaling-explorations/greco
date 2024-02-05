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

use crate::constants::sk_enc::{E_BOUND, K1_BOUND, N, R1_BOUNDS, R2_BOUNDS};
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
    k0is: Vec<String>,
    qis: Vec<String>,
    r2s: Vec<Vec<String>>,
    r1s: Vec<Vec<String>>,
    ais: Vec<Vec<String>>,
    ct0is: Vec<Vec<String>>,
}

/// Payload returned by the first phase of the circuit to be reused in the second phase
pub struct Payload<F: ScalarField> {
    s: Vec<String>,
    e: Vec<String>,
    k1: Vec<String>,
    k0_0: String,
    q0: String,
    r2_0: Vec<String>,
    r1_0: Vec<String>,
    ai_0: Vec<String>,
    ct0_0: Vec<String>,
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
    // TODO: explain
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
        let e_bound_constant = Constant(F::from(E_BOUND));

        for i in 0..N {
            let normalized_coeff = range.gate().add(ctx, e_assigned[i], e_bound_constant);
            // TODO: we can make this a bit more efficient by not having to reassign b as a constant every time
            range.check_less_than_safe(ctx, normalized_coeff, 2 * E_BOUND);
        }

        // perform range check on the coefficients of `k1`

        let k1_bound_constant = Constant(F::from(K1_BOUND));

        for i in 0..N {
            let normalized_coeff = range.gate().add(ctx, k1_assigned[i], k1_bound_constant);
            range.check_less_than_safe(ctx, normalized_coeff, 2 * K1_BOUND);
        }

        // perform range check on the coefficients of `r2_0`. Note that the degree of r2 is N - 2

        let r2_0_bound_constant = Constant(F::from(R2_BOUNDS[0]));

        for i in 0..(N - 1) {
            let normalized_coeff = range.gate().add(ctx, r2_0_assigned[i], r2_0_bound_constant);
            range.check_less_than_safe(ctx, normalized_coeff, 2 * R2_BOUNDS[0]);
        }

        // perform range check on the coefficients of `r1_0`. Note that the degree of r1 is 2N - 2

        let r1_0_bound_constant = Constant(F::from(R1_BOUNDS[0]));

        for i in 0..(2 * N - 1) {
            let normalized_coeff = range.gate().add(ctx, r1_0_assigned[i], r1_0_bound_constant);
            range.check_less_than_safe(ctx, normalized_coeff, 2 * R1_BOUNDS[0]);
        }

        Payload {
            s: self.s.clone(),
            e: self.e.clone(),
            k1: self.k1.clone(),
            k0_0: self.k0is[0].clone(),
            q0: self.qis[0].clone(),
            r2_0: self.r2s[0].clone(),
            r1_0: self.r1s[0].clone(),
            ai_0: self.ais[0].clone(),
            ct0_0: self.ct0is[0].clone(),
            s_assigned,
            e_assigned,
            k1_assigned,
            r2_0_assigned,
            r1_0_assigned,
        }
    }

    // Phase 1
    // TODO: explain
    fn virtual_assign_phase1(
        builder: &mut RlcCircuitBuilder<F>,
        range: &RangeChip<F>,
        rlc: &RlcChip<F>,
        payload: Self::FirstPhasePayload,
    ) {
        let Payload {
            s,
            e,
            k1,
            k0_0,
            q0,
            r2_0,
            r1_0,
            ai_0,
            ct0_0,
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

        let ai_0_poly = ai_0
            .iter()
            .map(|ai_0| F::from_str_vartime(ai_0).unwrap())
            .collect::<Vec<_>>();

        // Assign ai_0(gamma) to the circuit
        let ai_0_eval_at_gamma = evaluate_poly(&ai_0_poly, gamma);
        let assigned_poly_ai_0_eval_at_gamma = ctx_gate.load_witness(ai_0_eval_at_gamma);

        // cyclo poly is equal to x^N + 1
        let mut cyclo_poly = vec![F::from(0); N + 1];
        cyclo_poly[0] = F::from(1);
        cyclo_poly[N] = F::from(1);

        // Assign cyclo(gamma) to the circuit
        let cyclo_eval_at_gamma = evaluate_poly(&cyclo_poly, gamma);
        let assigned_poly_cyclo_eval_at_gamma = ctx_gate.load_witness(cyclo_eval_at_gamma);

        let ct0_0_poly = ct0_0
            .iter()
            .map(|ct0_0| F::from_str_vartime(ct0_0).unwrap())
            .collect::<Vec<_>>();

        // Assign ct0_0(gamma) to the circuit
        let ct0_0_eval_at_gamma = evaluate_poly(&ct0_0_poly, gamma);
        let assigned_poly_ct0_0_eval_at_gamma = ctx_gate.load_witness(ct0_0_eval_at_gamma);

        // Prove that LHS(gamma) = RHS(gamma)
        // LHS = ct0_0(gamma)
        // RHS = ai_0(gamma) * s(gamma) + e(gamma) + k1(gamma) * k0i_0 + r1_0(gamma) * qis_0 + r2_0(gamma) * cyclo(gamma)

        let gate = range.gate();

        // rhs_0 = ai_0(gamma) * s(gamma)
        let rhs_0 = gate.mul(
            ctx_gate,
            assigned_poly_ai_0_eval_at_gamma,
            rlc_trace_s.rlc_val,
        );

        // rhs_1 = e(gamma)
        let rhs_1 = rlc_trace_e.rlc_val;

        let k0_0_assigned = ctx_gate.load_witness(F::from_str_vartime(&k0_0).unwrap());

        // rhs_2 = k1(gamma) * k0i_0
        let rhs_2 = gate.mul(ctx_gate, rlc_trace_k1.rlc_val, k0_0_assigned);

        let q0_assigned = ctx_gate.load_witness(F::from_str_vartime(&q0).unwrap());

        // rhs_3 = r1_0(gamma) * qis_0
        let rhs_3 = gate.mul(ctx_gate, rlc_trace_r1_0.rlc_val, q0_assigned);

        // rhs_4 = r2_0(gamma) * cyclo(gamma)
        let rhs_4 = gate.mul(
            ctx_gate,
            rlc_trace_r2_0.rlc_val,
            assigned_poly_cyclo_eval_at_gamma,
        );

        let rhs_01 = gate.add(ctx_gate, rhs_0, rhs_1);
        let rhs_23 = gate.add(ctx_gate, rhs_2, rhs_3);
        let rhs_0123 = gate.add(ctx_gate, rhs_01, rhs_23);
        let rhs = gate.add(ctx_gate, rhs_0123, rhs_4);
        let lhs = assigned_poly_ct0_0_eval_at_gamma;

        let res = gate.sub(ctx_gate, lhs, rhs);
        println!("res: {:?}", res.value());
        gate.assert_is_const(ctx_gate, &res, &F::from(0));
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
