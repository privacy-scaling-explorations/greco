use axiom_eth::rlc::{
    chip::RlcChip,
    circuit::{builder::RlcCircuitBuilder, instructions::RlcCircuitInstructions, RlcCircuitParams},
};
use halo2_base::{
    gates::{circuit::BaseCircuitParams, GateInstructions, RangeChip, RangeInstructions},
    utils::ScalarField,
};
use serde::Deserialize;

use crate::constants::sk_enc::{E_BOUND, K1_BOUND, N, R1_BOUNDS, R2_BOUNDS, S_BOUND};
use crate::poly::{Poly, PolyAssigned};
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
    s_assigned: PolyAssigned<F>,
    e_assigned: PolyAssigned<F>,
    qis: Vec<String>,
    k1_assigned: PolyAssigned<F>,
    k0is: Vec<String>,
    r2is_assigned: Vec<PolyAssigned<F>>,
    r1is_assigned: Vec<PolyAssigned<F>>,
    ais: Vec<Vec<String>>,
    ct0is: Vec<Vec<String>>,
}

impl<F: ScalarField> RlcCircuitInstructions<F> for BfvSkEncryptionCircuit {
    type FirstPhasePayload = Payload<F>;

    /// #### Phase 0
    /// In this phase, the polynomials for each matrix $S_i$ are assigned to the circuit. Namely:
    /// * polynomials `s`, `e`, `k1` are assigned to the witness table. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * polynomials `r1i`, `r2i` are assigned to the witness table for each $S_i$ matrix.

    /// Witness values are element of the finite field $\mod{p}$. Negative coefficients $-z$ are assigned as field elements $p - z$.

    /// At the end of phase 0, the witness generated so far is interpolated into a polynomial and committed by the prover. The hash of this commitment is used as challenge and will be used as a source of randomness $\gamma$ in Phase 1. This feature is made available by Halo2 [Challenge API](https://hackmd.io/@axiom/SJw3p-qX3).
    fn virtual_assign_phase0(
        &self,
        builder: &mut RlcCircuitBuilder<F>,
        _: &RangeChip<F>,
    ) -> Self::FirstPhasePayload {
        let ctx = builder.base.main(0);

        let s = Poly::<F>::new(self.s.clone());
        let s_assigned = PolyAssigned::new(ctx, s);

        let e = Poly::<F>::new(self.e.clone());
        let e_assigned = PolyAssigned::new(ctx, e);

        let k1 = Poly::<F>::new(self.k1.clone());
        let k1_assigned = PolyAssigned::new(ctx, k1);

        let mut r2is_assigned = vec![];
        let mut r1is_assigned = vec![];

        for z in 0..self.ct0is.len() {
            let r2i = Poly::<F>::new(self.r2is[z].clone());
            let r2i_assigned = PolyAssigned::new(ctx, r2i);
            r2is_assigned.push(r2i_assigned);

            let r1i = Poly::<F>::new(self.r1is[z].clone());
            let r1i_assigned = PolyAssigned::new(ctx, r1i);
            r1is_assigned.push(r1i_assigned);
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

    /// #### Phase 1

    /// In this phase, the following two core constraints are enforced:

    /// - The coefficients of $S_i$ are in the expected range.  
    /// - $U_i(\gamma) \times S_i(\gamma) =Ct_{0,i}(\gamma)$

    /// ##### Assignment
    /// * Assign evaluations to the circuit: `ai(gamma)`, `ct0_i(gamma)` for each $i$-th CRT basis.
    /// * Assign `cyclo(gamma)` to the circuit
    /// * Assign scalars `q_i` and `k0_i` to the circuit for each $i$-th CRT basis
    /// * Expose `ai(gamma)`, `ct0_i(gamma)`, `q_i`, `k0_i` as public inputs for each $i$-th CRT basis
    /// * Expose `cyclo(gamma)` as public input

    /// Since these polynomials and scalars are known to the verifier, the evaluation at $\gamma$ doesn't need to be constrained inside the circuit. Instead, this can be safely be performed (and verified) outside the circuit.

    /// ##### Range Check

    /// The coefficients of the private polynomials from each $i$-th matrix $S_i$ are checked to be in the correct range.
    /// * Range check polynomials `s`, `e`, `k1`. This has to be done only once as these polynomial are common to each $S_i$ matrix.
    /// * Range check polynomials `r1i`, `r2i` for each `Si` matrix

    /// Since negative coefficients `-z` are assigned as `p - z` to the circuit, this might result in very large coefficients. Performing the range check on such large coefficients requires large lookup tables. To avoid this, the coefficients (both negative and positive) are shifted by a constant to make them positive and then perform the range check.

    /// ##### Correct Encryption Constraint

    /// It is needed to prove that $U_i(\gamma) \times S_i(\gamma) =Ct_{0,i}(\gamma)$. This can be rewritten as `ct0i = ct0i_hat + r1i * qi + r2i * cyclo` for each $i$-th CRT basis. Where `ct0i_hat = ai * s + e + k1 * k0i`.

    /// This constrained is enforced by proving that `LHS(gamma) = RHS(gamma)`. According to the Schwartz-Zippel lemma, if this relation between polynomial when evaluated at a random point holds true, then then the polynomials are identical with high probability.
    /// * Constrain the evaluation of the polynomials `s`, `e`, `k1` at $\gamma$. This has to be done only once as these polynomial are common to each $S_i$ matrix.
    /// * Constrain the evaluation of the polynomials `r1i`, `r2i` at $\gamma$ for each `Si` matrix
    /// * Constrain that `LHS(gamma) = RHS(gamma)` for each i-th CRT basis
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

        // ASSIGNMENT

        let (ctx_gate, ctx_rlc) = builder.rlc_ctx_pair();
        let gamma = *rlc.gamma();

        let mut ais_at_gamma_assigned = vec![];
        let mut ct0is_at_gamma_assigned = vec![];
        let mut qis_assigned = vec![];
        let mut k0is_assigned = vec![];

        for z in 0..ct0is.len() {
            let ai = Poly::<F>::new(ais[z].clone());
            let ai_at_gamma = ai.eval(gamma);
            let ai_at_gamma_assigned = ctx_gate.load_witness(ai_at_gamma);
            ais_at_gamma_assigned.push(ai_at_gamma_assigned);

            let ct0i = Poly::<F>::new(ct0is[z].clone());
            let ct0i_at_gamma = ct0i.eval(gamma);
            let ct0i_at_gamma_assigned = ctx_gate.load_witness(ct0i_at_gamma);
            ct0is_at_gamma_assigned.push(ct0i_at_gamma_assigned);

            let qi_val = F::from_str_vartime(&qis[z]).unwrap();
            let qi_assigned = ctx_gate.load_witness(qi_val);
            qis_assigned.push(qi_assigned);

            let k0i_val = F::from_str_vartime(&k0is[z]).unwrap();
            let k0i_assigned = ctx_gate.load_witness(k0i_val);
            k0is_assigned.push(k0i_assigned);
        }

        // cyclo poly is equal to x^N + 1
        let cyclo_at_gamma = gamma.pow_vartime([N as u64]) + F::from(1);
        let cyclo_at_gamma_assigned = ctx_gate.load_witness(cyclo_at_gamma);

        // TODO: expose ais_at_gamma_assigned, ct0is_at_gamma_assigned, cyclo_at_gamma_assigned, qis_assigned, k0is_assigned as public inputs

        // RANGE CHECK
        s_assigned.range_check(ctx_gate, range, S_BOUND);
        e_assigned.range_check(ctx_gate, range, E_BOUND);
        k1_assigned.range_check(ctx_gate, range, K1_BOUND);

        for z in 0..ct0is.len() {
            r2is_assigned[z].range_check(ctx_gate, range, R2_BOUNDS[z]);
            r1is_assigned[z].range_check(ctx_gate, range, R1_BOUNDS[z]);
        }

        // CORRECT ENCRYPTION CONSTRAINT

        let s_at_gamma = s_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let e_at_gamma = e_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let k1_at_gamma = k1_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);

        let gate = range.gate();

        // For each `i` Prove that LHS(gamma) = RHS(gamma)
        // LHS = ct0i(gamma)
        // RHS = ai(gamma) * s(gamma) + e(gamma) + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)
        for z in 0..ct0is.len() {
            let r1i_at_gamma = r1is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);
            let r2i_at_gamma = r2is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);

            // rhs = ai(gamma) * s(gamma) + e(gamma)
            let rhs = gate.mul_add(ctx_gate, ais_at_gamma_assigned[z], s_at_gamma, e_at_gamma);

            // rhs = rhs + k1(gamma) * k0i
            let rhs = gate.mul_add(ctx_gate, k1_at_gamma, k0is_assigned[z], rhs);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = gate.mul_add(ctx_gate, r1i_at_gamma, qis_assigned[z], rhs);

            // rhs = rhs + r2i(gamma) * cyclo(gamma)
            let rhs = gate.mul_add(ctx_gate, r2i_at_gamma, cyclo_at_gamma_assigned, rhs);
            let lhs = ct0is_at_gamma_assigned[z];

            // LHS(gamma) = RHS(gamma)
            let res = gate.is_equal(ctx_gate, lhs, rhs);
            gate.assert_is_const(ctx_gate, &res, &F::from(1));
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
                    location: FailureLocation::OutsideRegion { row: 104462 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 104472 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475585 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475605 }
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
                        offset: 1
                    }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475585 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475605 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475613 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475633 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475641 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475661 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475669 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475689 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475697 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475717 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475725 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475745 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475753 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475773 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475781 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475801 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475809 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475829 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475837 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475857 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475865 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475885 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475893 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475913 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475921 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475941 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475949 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475969 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 1475977 }
                },
            ])
        );
    }
}
