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
/// # Parameters:
/// * `s`: secret polynomial, sampled from ternary distribution. The coefficients are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
/// * `e`: error polynomial, sampled from discrete Gaussian distribution. The coefficients are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
/// * `qis`: list of qis such that qis[i] is the modulus of the i-th CRT basis of the modulus q of the ciphertext space.
/// * `k1`: scaled message polynomial. The coefficients are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
/// * `k0is`: list of the negative of the multiplicative inverses of t mod qis[i]. The values are normalized to be in range [0, p) where p is the modulus of the prime field of the circuit
/// * `r2is`: list of r2i polynomials for each i-th CRT basis . The coefficients are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
/// * `r1is`: list of r1i polynomials for each CRT i-th CRT basis. The coefficients are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
/// * `ais`: list of ai polynomials for each CRT i-th CRT basis. The coefficients are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
/// * `ct0is`: list of ct0i (first component of the ciphertext cti) polynomials for each CRT i-th CRT basis. The coefficients are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
#[derive(Deserialize)]
pub struct BfvSkEncryptionCircuit {
    s: Vec<String>,
    e: Vec<String>,
    qis: Vec<String>,
    k1: Vec<String>,
    k0is: Vec<String>,
    r2is: Vec<Vec<String>>,
    r1is: Vec<Vec<String>>,
    ais: Vec<Vec<String>>,
    ct0is: Vec<Vec<String>>,
}

/// Payload returned by the first phase of the circuit to be reused in the second phase
pub struct Payload<F: ScalarField> {
    s: Vec<String>,
    e: Vec<String>,
    k1: Vec<String>,
    r2_0: Vec<String>,
    r1_0: Vec<String>,
    a_0: Vec<String>,
    ct0_0: Vec<String>,
    q_0_assigned: AssignedValue<F>,
    k0_0_assigned: AssignedValue<F>,
    s_assigned: Vec<AssignedValue<F>>,
    e_assigned: Vec<AssignedValue<F>>,
    k1_assigned: Vec<AssignedValue<F>>,
    r2_0_assigned: Vec<AssignedValue<F>>,
    r1_0_assigned: Vec<AssignedValue<F>>,
}

impl<F: ScalarField> RlcCircuitInstructions<F> for BfvSkEncryptionCircuit {
    type FirstPhasePayload = Payload<F>;

    /// ## Phase 0
    /// In this phase, all the polynomials are assigned to the circuit. Namely:
    /// * polynomial `a_0` and the scalars `q_0` and `k0_0` from the matrix U_0
    /// * polynomials `s`, `e`, `k`, `r1_0`, `r2_0` from the matrix S_0
    /// * polynomials `ct0_0`
    /// The cyclotomic polynomial is not assigned to the circuit, as this is not an input but a constant parameter.
    /// `a_0`, `q_0`, `k0_0` and `ct0_0` are exposed as public inputs, while the other polynomials are kept private.
    ///
    /// Furhtermore:
    /// * Perform range checks on the coefficients of the secret polynomials
    fn virtual_assign_phase0(
        &self,
        builder: &mut RlcCircuitBuilder<F>,
        range: &RangeChip<F>,
    ) -> Self::FirstPhasePayload {
        let ctx = builder.base.main(0);

        let mut a_0_assigned = vec![];
        let mut s_assigned = vec![];
        let mut e_assigned = vec![];
        let mut k1_assigned = vec![];
        let mut r2_0_assigned = vec![];
        let mut r1_0_assigned = vec![];
        let mut ct0_0_assigned = vec![];

        // assign polynomial a_0 to the witness
        for j in 0..N {
            let val = F::from_str_vartime(&self.ais[0][j]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            a_0_assigned.push(coeff_assigned);
        }

        // assign q_0 to the witness
        let q_0 = F::from_str_vartime(&self.qis[0]).unwrap();
        let q_0_assigned = ctx.load_witness(q_0);

        // assign k0_0 to the witness
        let k0_0 = F::from_str_vartime(&self.k0is[0]).unwrap();
        let k0_0_assigned = ctx.load_witness(k0_0);

        // assign polynomial s to the witness

        for j in 0..N {
            let val = F::from_str_vartime(&self.s[j]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            s_assigned.push(coeff_assigned);
        }

        // assign polynomial e to the witness

        for j in 0..N {
            let val = F::from_str_vartime(&self.e[j]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            e_assigned.push(coeff_assigned);
        }

        // assign polynomial k1 to the witness

        for j in 0..N {
            let val = F::from_str_vartime(&self.k1[j]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            k1_assigned.push(coeff_assigned);
        }

        // assign polynomial r2is[0] to the witness. Note that the degree of r2is[0] is N - 2

        for j in 0..(N - 1) {
            let val = F::from_str_vartime(&self.r2is[0][j]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            r2_0_assigned.push(coeff_assigned);
        }

        // assign polynomial r1is[0] to the witness. Note that the degree of r1is[0] is 2N - 2
        for j in 0..(2 * N - 1) {
            let val = F::from_str_vartime(&self.r1is[0][j]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            r1_0_assigned.push(coeff_assigned);
        }

        // assign polynomial ct0is[0] to the witness
        for j in 0..N {
            let val = F::from_str_vartime(&self.ct0is[0][j]).unwrap();
            let coeff_assigned = ctx.load_witness(val);
            ct0_0_assigned.push(coeff_assigned);
        }

        // TODO: expose a_0_assigned, q_0_assigned, k0_0_assigned and ct0_0_assigned as public inputs

        // perform range check on the coefficients of `s`

        let constant_one = Constant(F::from(1));
        let bound_s = 3;

        for j in 0..N {
            let normalized_coeff = range.gate().add(ctx, s_assigned[j], constant_one);
            // TODO: we can make this a bit more efficient by not having to reassign b as a constant every time
            range.check_less_than_safe(ctx, normalized_coeff, bound_s);
        }

        // perform range check on the coefficients of `e`
        let e_bound_constant = Constant(F::from(E_BOUND));

        for j in 0..N {
            let normalized_coeff = range.gate().add(ctx, e_assigned[j], e_bound_constant);
            // TODO: we can make this a bit more efficient by not having to reassign b as a constant every time
            range.check_less_than_safe(ctx, normalized_coeff, 2 * E_BOUND);
        }

        // perform range check on the coefficients of `k1`

        let k1_bound_constant = Constant(F::from(K1_BOUND));

        for j in 0..N {
            let normalized_coeff = range.gate().add(ctx, k1_assigned[j], k1_bound_constant);
            // TODO: we can make this a bit more efficient by not having to reassign b as a constant every time
            range.check_less_than_safe(ctx, normalized_coeff, 2 * K1_BOUND);
        }

        // perform range check on the coefficients of `r2_0`

        let r2_0_bound_constant = Constant(F::from(R2_BOUNDS[0]));

        for j in 0..(N - 1) {
            let normalized_coeff = range.gate().add(ctx, r2_0_assigned[j], r2_0_bound_constant);
            // TODO: we can make this a bit more efficient by not having to reassign b as a constant every time
            range.check_less_than_safe(ctx, normalized_coeff, 2 * R2_BOUNDS[0]);
        }

        // perform range check on the coefficients of `r1_0`

        let r1_0_bound_constant = Constant(F::from(R1_BOUNDS[0]));

        for j in 0..(2 * N - 1) {
            let normalized_coeff = range.gate().add(ctx, r1_0_assigned[j], r1_0_bound_constant);
            // TODO: we can make this a bit more efficient by not having to reassign b as a constant every time
            range.check_less_than_safe(ctx, normalized_coeff, 2 * R1_BOUNDS[0]);
        }

        Payload {
            s: self.s.clone(),
            e: self.e.clone(),
            k1: self.k1.clone(),
            r2_0: self.r2is[0].clone(),
            r1_0: self.r1is[0].clone(),
            a_0: self.ais[0].clone(),
            ct0_0: self.ct0is[0].clone(),
            q_0_assigned,
            k0_0_assigned,
            s_assigned,
            e_assigned,
            k1_assigned,
            r2_0_assigned,
            r1_0_assigned,
        }
    }

    /// ## Phase 1
    /// * Fetch challenge `gamma` from phase 0
    /// * Assign public inputs to the circuit: `a_0(gamma)`, `cyclo(gamma)`, `ct0_0(gamma)`. Since the polynomials are public from phase 0, the evaluation at gamma doesn't need to be constrained inside the circuit,
    /// but can safely be performed (and verified) outside the circuit
    /// * Enforce the evaluation of the RLC for the polynomials `s`, `e`, `k1`, `r2_0`, `r1_0` at `gamma`
    /// * Enforce that `ct0_0(gamma) = a_0(gamma) * s(gamma) + e(gamma) + k1(gamma) * k0_0 + r1_0(gamma) * q_0 + r2_0(gamma) * cyclo(gamma)`
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
            r2_0,
            r1_0,
            a_0,
            ct0_0,
            q_0_assigned,
            k0_0_assigned,
            s_assigned,
            e_assigned,
            k1_assigned,
            r2_0_assigned,
            r1_0_assigned,
        } = payload;

        let (ctx_gate, ctx_rlc) = builder.rlc_ctx_pair();
        let gamma = *rlc.gamma();

        let a_0_poly = a_0
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect::<Vec<_>>();

        // Assign a_0(gamma) to the circuit
        let a_0_eval_at_gamma = evaluate_poly(&a_0_poly, gamma);
        let a_0_eval_at_gamma_assigned = ctx_gate.load_witness(a_0_eval_at_gamma);

        // TODO: expose a_0_eval_at_gamma_assigned as a public input

        // cyclo poly is equal to x^N + 1
        let mut cyclo_poly = vec![F::from(0); N + 1];
        cyclo_poly[0] = F::from(1);
        cyclo_poly[N] = F::from(1);

        // Assign cyclo(gamma) to the circuit
        let cyclo_eval_at_gamma = evaluate_poly(&cyclo_poly, gamma);
        let cyclo_eval_at_gamma_assigned = ctx_gate.load_witness(cyclo_eval_at_gamma);

        // TODO: expose cyclo_eval_at_gamma_assigned as a public input

        let ct0_0_poly = ct0_0
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect::<Vec<_>>();

        // Assign ct0_0(gamma) to the circuit
        let ct0_0_eval_at_gamma = evaluate_poly(&ct0_0_poly, gamma);
        let ct0_0_eval_at_gamma_assigned = ctx_gate.load_witness(ct0_0_eval_at_gamma);

        // TODO: expose ct0_0_eval_at_gamma_assigned as a public input

        let rlc_trace_s = rlc.compute_rlc_fixed_len(ctx_rlc, s_assigned);
        let rlc_trace_e = rlc.compute_rlc_fixed_len(ctx_rlc, e_assigned);
        let rlc_trace_k1 = rlc.compute_rlc_fixed_len(ctx_rlc, k1_assigned);
        let rlc_trace_r2_0 = rlc.compute_rlc_fixed_len(ctx_rlc, r2_0_assigned);
        let rlc_trace_r1_0 = rlc.compute_rlc_fixed_len(ctx_rlc, r1_0_assigned);

        // SANITY CHECK
        // assert the correctness of the RLC for s, e, k1, r2_0, r1_0
        let s_poly = s
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_s_eval_at_gamma = evaluate_poly(&s_poly, gamma);

        assert_eq!(rlc_trace_s.rlc_val.value(), &expected_poly_s_eval_at_gamma);

        let e_poly = e
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_e_eval_at_gamma = evaluate_poly(&e_poly, gamma);

        assert_eq!(rlc_trace_e.rlc_val.value(), &expected_poly_e_eval_at_gamma);

        let k1_poly = k1
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_k1_eval_at_gamma = evaluate_poly(&k1_poly, gamma);

        assert_eq!(
            rlc_trace_k1.rlc_val.value(),
            &expected_poly_k1_eval_at_gamma
        );

        let r2_0_poly = r2_0
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_r2_0_eval_at_gamma = evaluate_poly(&r2_0_poly, gamma);

        assert_eq!(
            rlc_trace_r2_0.rlc_val.value(),
            &expected_poly_r2_0_eval_at_gamma
        );

        let r1_0_poly = r1_0
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect::<Vec<_>>();

        let expected_poly_r1_0_eval_at_gamma = evaluate_poly(&r1_0_poly, gamma);

        assert_eq!(
            rlc_trace_r1_0.rlc_val.value(),
            &expected_poly_r1_0_eval_at_gamma
        );

        // SANITY CHECK OVER

        // Prove that LHS(gamma) = RHS(gamma)
        // LHS = ct0_0(gamma)
        // RHS = a_0(gamma) * s(gamma) + e(gamma) + k1(gamma) * k0_0 + r1_0(gamma) * q_0 + r2_0(gamma) * cyclo(gamma)

        let gate = range.gate();

        // rhs_0 = a_0(gamma) * s(gamma)
        let rhs_0 = gate.mul(ctx_gate, a_0_eval_at_gamma_assigned, rlc_trace_s.rlc_val);

        // rhs_1 = e(gamma)
        let rhs_1 = rlc_trace_e.rlc_val;

        // rhs_2 = k1(gamma) * k0_0
        let rhs_2 = gate.mul(ctx_gate, rlc_trace_k1.rlc_val, k0_0_assigned);

        // rhs_3 = r1_0(gamma) * q_0
        let rhs_3 = gate.mul(ctx_gate, rlc_trace_r1_0.rlc_val, q_0_assigned);

        // rhs_4 = r2_0(gamma) * cyclo(gamma)
        let rhs_4 = gate.mul(
            ctx_gate,
            rlc_trace_r2_0.rlc_val,
            cyclo_eval_at_gamma_assigned,
        );

        let rhs_01 = gate.add(ctx_gate, rhs_0, rhs_1);
        let rhs_23 = gate.add(ctx_gate, rhs_2, rhs_3);
        let rhs_0123 = gate.add(ctx_gate, rhs_01, rhs_23);
        let rhs = gate.add(ctx_gate, rhs_0123, rhs_4);
        let lhs = ct0_0_eval_at_gamma_assigned;

        let res = gate.sub(ctx_gate, lhs, rhs);
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
        halo2_proofs::{
            dev::{FailureLocation, MockProver, VerifyFailure},
            halo2curves::bn256::Fr,
            plonk::{Any, SecondPhase},
        },
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

        let rlc_circuit = RlcExecutor::new(mock_builder, circuit);

        // 3. Run the mock prover. The circuit should be satisfied
        MockProver::run(K as u32, &rlc_circuit, vec![])
            .unwrap()
            .assert_satisfied();
    }

    #[test]
    pub fn test_sk_enc_invalid_range() {
        // 1. Define the inputs of the circuit
        let file_path = "src/data/sk_enc_input.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let mut circuit = serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

        // 2. Invalidate the circuit by setting the value of a coefficient of the polynomial `s` to be out of range
        let out_of_range_s = "2";

        circuit.s[0] = out_of_range_s.to_string();

        // 3. Build the circuit for MockProver
        let params = test_params();
        let mut mock_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0).use_params(params);
        mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let rlc_circuit = RlcExecutor::new(mock_builder, circuit);

        // 4. Run the mock prover
        let invalid_prover = MockProver::run(K as u32, &rlc_circuit, vec![]).unwrap();

        // 5. Assert that the circuit is not satisfied
        // In particular, it should fail the range check enforced in the first phase and the equality check in the second phase
        assert_eq!(
            invalid_prover.verify(),
            Err(vec![
                VerifyFailure::Lookup {
                    name: "lookup".to_string(),
                    lookup_index: 0,
                    location: FailureLocation::OutsideRegion { row: 8196 }
                },
                VerifyFailure::Permutation {
                    column: (Any::Fixed, 1).into(),
                    location: FailureLocation::InRegion {
                        region: (2, "base+rlc phase 1").into(),
                        offset: 0
                    }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 35 }
                }
            ])
        );
    }

    #[test]
    pub fn test_sk_enc_invalid_polys() {
        // 1. Define the inputs of the circuit
        let file_path = "src/data/sk_enc_input.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let mut circuit = serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

        // 2. Invalidate the circuit by setting a different `s` polynomial
        let invalid_s = vec!["1".to_string(); 1024];

        circuit.s = invalid_s;

        // 3. Build the circuit for MockProver
        let params = test_params();
        let mut mock_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0).use_params(params);
        mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let rlc_circuit = RlcExecutor::new(mock_builder, circuit);

        // 4. Run the mock prover
        let invalid_prover = MockProver::run(K as u32, &rlc_circuit, vec![]).unwrap();

        // 5. Assert that the circuit is not satisfied
        // In particular, it should only fail the equality check in the second phase
        assert_eq!(
            invalid_prover.verify(),
            Err(vec![
                VerifyFailure::Permutation {
                    column: (Any::Fixed, 1).into(),
                    location: FailureLocation::InRegion {
                        region: (2, "base+rlc phase 1").into(),
                        offset: 0
                    }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 35 }
                }
            ])
        );
    }
}
