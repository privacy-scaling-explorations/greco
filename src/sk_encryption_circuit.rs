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
const TEST_K: usize = 22;

/// Helper function to define the parameters of the RlcCircuit
pub fn test_params() -> RlcCircuitParams {
    RlcCircuitParams {
        base: BaseCircuitParams {
            k: TEST_K,
            num_advice_per_phase: vec![1, 1],
            num_fixed: 1,
            num_lookup_advice_per_phase: vec![0, 1],
            lookup_bits: Some(8),
            num_instance_columns: 0,
        },
        num_rlc_columns: 1,
    }
}

/// `BfvSkEncryptionCircuit` is a circuit that checks the correct formation of a ciphertext resulting from BFV secret key encryption
/// All the polynomials coefficients and scalars are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
///
/// # Parameters:
/// * `s`: secret polynomial, sampled from ternary distribution.
/// * `e`: error polynomial, sampled from discrete Gaussian distribution.
/// * `qis`: list of scalars qis such that qis[i] is the modulus of the i-th CRT basis of the modulus q of the ciphertext space.
/// * `k1`: scaled message polynomial.
/// * `k0is`: list of the scalars equal negative of the multiplicative inverses of t mod qis[i].
/// * `r2is`: list of r2i polynomials for each i-th CRT basis .
/// * `r1is`: list of r1i polynomials for each CRT i-th CRT basis.
/// * `ais`: list of ai polynomials for each CRT i-th CRT basis.
/// * `ct0is`: list of ct0i (first component of the ciphertext cti) polynomials for each CRT i-th CRT basis.
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
    s_assigned: Vec<AssignedValue<F>>,
    e_assigned: Vec<AssignedValue<F>>,
    qis: Vec<String>,
    k1_assigned: Vec<AssignedValue<F>>,
    k0is: Vec<String>,
    r2is_assigned: Vec<Vec<AssignedValue<F>>>,
    r1is_assigned: Vec<Vec<AssignedValue<F>>>,
    ais: Vec<Vec<String>>,
    ct0is: Vec<Vec<String>>,
}

impl<F: ScalarField> RlcCircuitInstructions<F> for BfvSkEncryptionCircuit {
    type FirstPhasePayload = Payload<F>;

    /// ## Phase 0
    /// In this phase, the polynomials for each matrix `Si` are assigned to the circuit. Namely:
    /// * polynomials `s`, `e`, `k1` are assigned only once as common to each `Si` matrix
    /// * polynomials `r1i`, `r2i` are assigned for each `Si` matrix
    ///
    /// At the end of phase 0, the witness is committed and hashed and a challenge is extracted (Fiat-Shamir heuristic).
    fn virtual_assign_phase0(
        &self,
        builder: &mut RlcCircuitBuilder<F>,
        _: &RangeChip<F>,
    ) -> Self::FirstPhasePayload {
        let ctx = builder.base.main(0);

        let mut s_assigned = vec![];
        let mut e_assigned = vec![];
        let mut k1_assigned = vec![];
        let mut r2is_assigned = vec![];
        let mut r1is_assigned = vec![];

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

        // assign polynomials r2is to the witness. The degree of r2i is N - 2
        for z in 0..self.r2is.len() {
            let mut vec = vec![];
            for j in 0..(N - 1) {
                let val = F::from_str_vartime(&self.r2is[z][j]).unwrap();
                let coeff_assigned = ctx.load_witness(val);
                vec.push(coeff_assigned);
            }
            r2is_assigned.push(vec);
        }

        // assign polynomials r1is to the witness. The degree of r1i is 2N - 2
        for z in 0..self.r1is.len() {
            let mut vec = vec![];
            for j in 0..(2 * N - 1) {
                let val = F::from_str_vartime(&self.r1is[z][j]).unwrap();
                let coeff_assigned = ctx.load_witness(val);
                vec.push(coeff_assigned);
            }
            r1is_assigned.push(vec);
        }

        Payload {
            s_assigned,
            e_assigned,
            qis: self.qis.clone(),
            k1_assigned,
            k0is: self.k0is.clone(),
            r2is_assigned,
            r1is_assigned,
            ais: self.ais.clone(),
            ct0is: self.ct0is.clone(),
        }
    }

    /// ## Phase 1
    ///
    /// ### ASSIGNMENT
    /// * Fetch challenge `gamma` from phase 0
    /// * Assign evaluations to the circuit: `ai(gamma)`, `ct0_i(gamma)` for each i-th CRT basis.
    /// * Assign `cyclo(gamma)` to the circuit
    /// * Assign scalars `q_i` and `k0_i` to the circuit for each i-th CRT basis
    /// * Expose `ai(gamma)`, `ct0_i(gamma)`, `q_i`, `k0_i` as public inputs for each i-th CRT basis
    /// * Expose `cyclo(gamma)` as public input
    ///
    /// Since the polynomials are public from phase 0, the evaluation at gamma doesn't need to be constrained inside the circuit,
    /// but can safely be performed (and verified) outside the circuit
    ///
    /// ### RANGE CHECK OF PRIVATE POLYNOMIALS
    /// The coefficients of the private polynomials from each i-th matrix `Si` are checked to be in the correct range.
    /// * polynomials `s`, `e`, `k1` are range checked only once as common to each `Si` matrix
    /// * polynomials `r1i`, `r2i` are range checked for each `Si` matrix
    /// 
    /// Negative coefficients `-z` are assigned as `p - z` to the circuit. For example `-1` is assigned as `p - 1`.
    /// Performing the range check on such large coefficients is not efficient are requires large lookup tables.
    /// To avoid this, we shift the coefficients (both negative and positive) by a constant to make them positive and then perform the range check.
    ///
    /// ### CORRECT ENCRYPTION CONSTRAINT
    /// We need to prove that `ct0i = ct0i_hat + r1i * qi + r2i * cyclo` mod p for each i-th CRT basis.
    /// Where `ct0i_hat = ai * s + e + k1 * k0i`
    /// We do that by proving that `LHS(gamma) = RHS(gamma)` according to Scwhartz-Zippel lemma.
    ///
    /// * Constrain the evaluation of the polynomials `s`, `e`, `k1` at gamma. Need to be performed only once as common to each `Si` matrix
    /// * Constrain the evaluation of the polynomials `r1i`, `r2i` at gamma for each `Si` matrix
    /// * Constrain that LHS(gamma) = RHS(gamma) for each i-th CRT basis
    fn virtual_assign_phase1(
        builder: &mut RlcCircuitBuilder<F>,
        range: &RangeChip<F>,
        rlc: &RlcChip<F>,
        payload: Self::FirstPhasePayload,
    ) {
        let Payload {
            s_assigned,
            e_assigned,
            qis,
            k1_assigned,
            k0is,
            r2is_assigned,
            r1is_assigned,
            ais,
            ct0is,
        } = payload;

        // ASSIGNMENT PHASE

        let (ctx_gate, ctx_rlc) = builder.rlc_ctx_pair();
        let gamma = *rlc.gamma();

        let mut ais_at_gamma_assigned = vec![];
        let mut ct0is_at_gamma_assigned = vec![];

        // TODO: Do everything in a single loop
        for ai in &ais {
            let ai_poly = ai
                .iter()
                .map(|coeff| F::from_str_vartime(coeff).unwrap())
                .collect::<Vec<_>>();

            let ai_at_gamma = evaluate_poly(&ai_poly, gamma);
            let ai_at_gamma_assigned = ctx_gate.load_witness(ai_at_gamma);
            ais_at_gamma_assigned.push(ai_at_gamma_assigned);
        }

        for ct0i in &ct0is {
            let ct0i_poly = ct0i
                .iter()
                .map(|coeff| F::from_str_vartime(coeff).unwrap())
                .collect::<Vec<_>>();

            let ct0i_at_gamma = evaluate_poly(&ct0i_poly, gamma);
            let ct0i_at_gamma_assigned = ctx_gate.load_witness(ct0i_at_gamma);
            ct0is_at_gamma_assigned.push(ct0i_at_gamma_assigned);
        }

        // cyclo poly is equal to x^N + 1
        let mut cyclo_poly = vec![F::from(0); N + 1];
        cyclo_poly[0] = F::from(1);
        cyclo_poly[N] = F::from(1);

        // Assign cyclo(gamma) to the circuit
        let cyclo_at_gamma = evaluate_poly(&cyclo_poly, gamma);
        let cyclo_at_gamma_assigned = ctx_gate.load_witness(cyclo_at_gamma);

        let mut qis_assigned = vec![];
        let mut k0is_assigned = vec![];

        for qi in &qis {
            let qi_val = F::from_str_vartime(qi).unwrap();
            let qi_assigned = ctx_gate.load_witness(qi_val);
            qis_assigned.push(qi_assigned);
        }

        for k0i in &k0is {
            let k0i_val = F::from_str_vartime(k0i).unwrap();
            let k0i_assigned = ctx_gate.load_witness(k0i_val);
            k0is_assigned.push(k0i_assigned);
        }

        // TODO: expose ais_at_gamma_assigned, ct0is_at_gamma_assigned, cyclo_at_gamma_assigned, qis_assigned, k0is_assigned as public inputs

        // perform range check on the coefficients of `s`

        let s_bound_constant = Constant(F::from(1));
        let s_bound = 3;

        for j in 0..N {
            let shifted_coeff = range.gate().add(ctx_gate, s_assigned[j], s_bound_constant);
            // TODO: we can make this a bit more efficient by not having to reassign the constant every time
            range.check_less_than_safe(ctx_gate, shifted_coeff, s_bound);
        }

        // perform range check on the coefficients of `e`
        let e_bound_constant = Constant(F::from(E_BOUND));

        for j in 0..N {
            let shifted_coeff = range.gate().add(ctx_gate, e_assigned[j], e_bound_constant);
            // TODO: we can make this a bit more efficient by not having to reassign the constant every time
            range.check_less_than_safe(ctx_gate, shifted_coeff, 2 * E_BOUND);
        }

        // perform range check on the coefficients of `k1`

        let k1_bound_constant = Constant(F::from(K1_BOUND));

        for j in 0..N {
            let shifted_coeff = range
                .gate()
                .add(ctx_gate, k1_assigned[j], k1_bound_constant);
            // TODO: we can make this a bit more efficient by not having to reassign the constant every time
            range.check_less_than_safe(ctx_gate, shifted_coeff, 2 * K1_BOUND);
        }

        // perform range check on the coefficients of `r2is`
        for z in 0..r2is_assigned.len() {
            let r2i_bound_constant = Constant(F::from(R2_BOUNDS[z]));
            for j in 0..(N - 1) {
                let shifted_coeff =
                    range
                        .gate()
                        .add(ctx_gate, r2is_assigned[z][j], r2i_bound_constant);
                // TODO: we can make this a bit more efficient by not having to reassign the constant every time
                range.check_less_than_safe(ctx_gate, shifted_coeff, 2 * R2_BOUNDS[z]);
            }
        }

        // perform range check on the coefficients of `r1is`
        for z in 0..r1is_assigned.len() {
            let r1i_bound_constant = Constant(F::from(R1_BOUNDS[z]));
            for j in 0..(2 * N - 1) {
                let shifted_coeff =
                    range
                        .gate()
                        .add(ctx_gate, r1is_assigned[z][j], r1i_bound_constant);
                // TODO: we can make this a bit more efficient by not having to reassign the constant every time
                range.check_less_than_safe(ctx_gate, shifted_coeff, 2 * R1_BOUNDS[z]);
            }
        }

        // CORRECT ENCRYPTION CONSTRAINT

        // Constrain the evaluation of the polynomials `s`, `e`, `k1` at gamma
        let rlc_trace_s = rlc.compute_rlc_fixed_len(ctx_rlc, s_assigned);
        let rlc_trace_e = rlc.compute_rlc_fixed_len(ctx_rlc, e_assigned);
        let rlc_trace_k1 = rlc.compute_rlc_fixed_len(ctx_rlc, k1_assigned);

        let s_at_gamma = rlc_trace_s.rlc_val;
        let e_at_gamma = rlc_trace_e.rlc_val;
        let k1_at_gamma = rlc_trace_k1.rlc_val;

        let gate = range.gate();

        // For each `i` Prove that LHS(gamma) = RHS(gamma)
        // LHS = ct0i(gamma)
        // RHS = ai(gamma) * s(gamma) + e(gamma) + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)
        for z in 0..ct0is.len() {
            let rlc_trace_r1_i = rlc.compute_rlc_fixed_len(ctx_rlc, r1is_assigned[z].clone());
            let rlc_trace_r2_i = rlc.compute_rlc_fixed_len(ctx_rlc, r2is_assigned[z].clone());

            let r1i_at_gamma = rlc_trace_r1_i.rlc_val;
            let r2i_at_gamma = rlc_trace_r2_i.rlc_val;

            // rhs_0 = ai(gamma) * s(gamma)
            let rhs_0 = gate.mul(ctx_gate, ais_at_gamma_assigned[z], s_at_gamma);

            // rhs_1 = e(gamma)
            let rhs_1 = e_at_gamma;

            // rhs_2 = k1(gamma) * k0i
            let rhs_2 = gate.mul(ctx_gate, k1_at_gamma, k0is_assigned[z]);

            // rhs_3 = r1i(gamma) * qi
            let rhs_3 = gate.mul(ctx_gate, r1i_at_gamma, qis_assigned[z]);

            // rhs_4 = r2i(gamma) * cyclo(gamma)
            let rhs_4 = gate.mul(ctx_gate, r2i_at_gamma, cyclo_at_gamma_assigned);

            let rhs_01 = gate.add(ctx_gate, rhs_0, rhs_1);
            let rhs_23 = gate.add(ctx_gate, rhs_2, rhs_3);
            let rhs_0123 = gate.add(ctx_gate, rhs_01, rhs_23);
            let rhs = gate.add(ctx_gate, rhs_0123, rhs_4);
            let lhs = ct0is_at_gamma_assigned[z];

            let res = gate.sub(ctx_gate, lhs, rhs);
            // TODO: we can make this a bit more efficient by not having to reassign the constant every time
            gate.assert_is_const(ctx_gate, &res, &F::from(0));
        }
    }
}

#[cfg(test)]
mod test {

    use std::{fs::File, io::Read};

    use super::test_params;
    use crate::{constants::sk_enc::R1_BOUNDS, sk_encryption_circuit::BfvSkEncryptionCircuit};
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
        let circuit_params = test_params();
        let mut mock_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0)
                .use_params(circuit_params.clone());
        mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let rlc_circuit = RlcExecutor::new(mock_builder, circuit);

        // 3. Run the mock prover. The circuit should be satisfied
        MockProver::run(
            circuit_params.base.k.try_into().unwrap(),
            &rlc_circuit,
            vec![],
        )
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

        // 2. Invalidate the circuit by setting the value of a coefficient of the polynomial `r1is[0]` to be out of range
        let out_of_range_coeff = R1_BOUNDS[0] + 1;
        circuit.r1is[0][0] = out_of_range_coeff.to_string();

        // 3. Build the circuit for MockProver
        let circuit_params = test_params();
        let mut mock_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0)
                .use_params(circuit_params.clone());
        mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let rlc_circuit = RlcExecutor::new(mock_builder, circuit);

        // 4. Run the mock prover
        let invalid_prover = MockProver::run(
            circuit_params.base.k.try_into().unwrap(),
            &rlc_circuit,
            vec![],
        )
        .unwrap();

        // 5. Assert that the circuit is not satisfied
        // In particular, it should fail the range check enforced in the second phase for the first coefficient of r1is[0] and the equality check in the second phase for the 0-th basis
        assert_eq!(
            invalid_prover.verify(),
            Err(vec![
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 892172 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 892182 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475591 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475595 }
                },
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
        let circuit_params = test_params();
        let mut mock_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0)
                .use_params(circuit_params.clone());
        mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let rlc_circuit = RlcExecutor::new(mock_builder, circuit);

        // 4. Run the mock prover
        let invalid_prover = MockProver::run(
            circuit_params.base.k.try_into().unwrap(),
            &rlc_circuit,
            vec![],
        )
        .unwrap();

        // 5. Assert that the circuit is not satisfied
        // In particular, it should fail the equality check (LHS=RHS) in the second phase for each i-th CRT basis
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
                    location: FailureLocation::OutsideRegion { row: 1475591 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475595 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475627 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475631 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475663 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475667 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475699 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475703 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475735 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475739 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475771 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475775 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475807 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475811 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475843 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475847 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475879 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475883 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475915 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475919 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475951 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475955 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475987 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475991 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1476023 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1476027 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1476059 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1476063 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1476095 }
                },
            ])
        );
    }
}
