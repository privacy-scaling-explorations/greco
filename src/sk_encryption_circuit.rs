use crate::constants::sk_enc_constants_4096_2x55_65537::{
    E_BOUND, K0IS, K1_BOUND, N, QIS, R1_BOUNDS, R2_BOUNDS, S_BOUND,
};
use crate::poly::{Poly, PolyAssigned};
use axiom_eth::rlc::{
    chip::RlcChip,
    circuit::{builder::RlcCircuitBuilder, instructions::RlcCircuitInstructions, RlcCircuitParams},
};
use halo2_base::{
    gates::{circuit::BaseCircuitParams, GateInstructions, RangeChip, RangeInstructions},
    utils::ScalarField,
    QuantumCell::Constant,
};
use serde::Deserialize;

/// Helper function to define the parameters of the RlcCircuit. This is a non-optimized configuration that makes use of a single advice column. Use this for testing purposes only.
pub fn test_params() -> RlcCircuitParams {
    RlcCircuitParams {
        base: BaseCircuitParams {
            k: 21,
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
/// * `k1`: scaled message polynomial.
/// * `r2is`: list of r2i polynomials for each i-th CRT basis .
/// * `r1is`: list of r1i polynomials for each CRT i-th CRT basis.
/// * `ais`: list of ai polynomials for each CRT i-th CRT basis.
/// * `ct0is`: list of ct0i (first component of the ciphertext cti) polynomials for each CRT i-th CRT basis.
#[derive(Deserialize, Clone)]
pub struct BfvSkEncryptionCircuit {
    s: Vec<String>,
    e: Vec<String>,
    k1: Vec<String>,
    r2is: Vec<Vec<String>>,
    r1is: Vec<Vec<String>>,
    ais: Vec<Vec<String>>,
    ct0is: Vec<Vec<String>>,
}

/// Payload returned by the first phase of the circuit to be reused in the second phase
pub struct Payload<F: ScalarField> {
    s_assigned: PolyAssigned<F>,
    e_assigned: PolyAssigned<F>,
    k1_assigned: PolyAssigned<F>,
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
    /// * polynomials `r1i`, `r2i` are assigned to the witness table for each $S_i$ matrix

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
            k1_assigned,
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
    /// * Assign evaluations to the circuit: `ai(gamma)`, `ct0i(gamma)` for each $U_i$ matrix
    /// * Assign `cyclo(gamma)` to the circuit. This has to be done only once as the cyclotomic polynomial is common to each $U_i$ matrix
    /// * Expose `ai(gamma)`, `ct0i(gamma)` for each $U_i$ matrix
    /// * Expose `cyclo(gamma)` as public input

    /// Since the polynomials `cyclo`, `ai` and `ct0i` are known to the verifier, the evaluation at $\gamma$ doesn't need to be constrained inside the circuit. Instead, this can be safely be performed (and verified) outside the circuit.

    /// ##### Range Check

    /// The coefficients of the private polynomials from each $i$-th matrix $S_i$ are checked to be in the correct range
    /// * Range check polynomials `s`, `e`, `k1`. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * Range check polynomials `r1i`, `r2i` for each $S_i$ matrix

    /// Since negative coefficients `-z` are assigned as `p - z` to the circuit, this might result in very large coefficients. Performing the range check on such large coefficients requires large lookup tables. To avoid this, the coefficients (both negative and positive) are shifted by a constant to make them positive and then perform the range check.

    /// ##### Evaluation at $\gamma$ Constraint

    /// Contrary to the polynomials `cyclo`, `ai` and `ct0i`, the polynomials belonging to each $S_i$ matrix  are not known by the verifier. Therefore, their evaluation at $\gamma$ must be constrained inside the circuit.

    /// * Constrain the evaluation of the polynomials `s`, `e`, `k1` at $\gamma$. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * Constrain the evaluation of the polynomials `r1i`, `r2i` at $\gamma$ for each $S_i$ matrix

    /// ##### Correct Encryption Constraint

    /// It is needed to prove that $U_i(\gamma) \times S_i(\gamma) =Ct_{0,i}(\gamma)$. This can be rewritten as `ct0i = ct0i_hat + r1i * qi + r2i * cyclo`, where `ct0i_hat = ai * s + e + k1 * k0i`.

    /// This constrained is enforced by proving that `LHS(gamma) = RHS(gamma)`. According to the Schwartz-Zippel lemma, if this relation between polynomial when evaluated at a random point holds true, then then the polynomials are identical with high probability. Note that `qi` and `k0i` (for each $U_i$ matrix) are constants to the circuit encoded during key generation.
    /// * Constrain that `ct0i(gamma) = ai(gamma) * s(gamma) + e(gamma) + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)` for each $i$-th CRT basis
    fn virtual_assign_phase1(
        builder: &mut RlcCircuitBuilder<F>,
        range: &RangeChip<F>,
        rlc: &RlcChip<F>,
        payload: Self::FirstPhasePayload,
    ) {
        let Payload {
            s_assigned,
            e_assigned,
            k1_assigned,
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
        let mut qi_constants = vec![];
        let mut k0i_constants = vec![];

        for z in 0..ct0is.len() {
            let ai = Poly::<F>::new(ais[z].clone());
            let ai_at_gamma = ai.eval(gamma);
            let ai_at_gamma_assigned = ctx_gate.load_witness(ai_at_gamma);
            ais_at_gamma_assigned.push(ai_at_gamma_assigned);

            let ct0i = Poly::<F>::new(ct0is[z].clone());
            let ct0i_at_gamma = ct0i.eval(gamma);
            let ct0i_at_gamma_assigned = ctx_gate.load_witness(ct0i_at_gamma);
            ct0is_at_gamma_assigned.push(ct0i_at_gamma_assigned);

            let qi_constant = Constant(F::from_str_vartime(QIS[z]).unwrap());
            qi_constants.push(qi_constant);

            let k0i_constant = Constant(F::from_str_vartime(K0IS[z]).unwrap());
            k0i_constants.push(k0i_constant);
        }

        // cyclo poly is equal to x^N + 1
        let cyclo_at_gamma = gamma.pow_vartime([N as u64]) + F::from(1);
        let cyclo_at_gamma_assigned = ctx_gate.load_witness(cyclo_at_gamma);

        // RANGE CHECK
        s_assigned.range_check(ctx_gate, range, S_BOUND);
        e_assigned.range_check(ctx_gate, range, E_BOUND);
        k1_assigned.range_check(ctx_gate, range, K1_BOUND);

        for z in 0..ct0is.len() {
            r2is_assigned[z].range_check(ctx_gate, range, R2_BOUNDS[z]);
            r1is_assigned[z].range_check(ctx_gate, range, R1_BOUNDS[z]);
        }

        // EVALUATION AT GAMMA CONSTRAINT

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

            // CORRECT ENCRYPTION CONSTRAINT

            // rhs = ai(gamma) * s(gamma) + e(gamma)
            let rhs = gate.mul_add(ctx_gate, ais_at_gamma_assigned[z], s_at_gamma, e_at_gamma);

            // rhs = rhs + k1(gamma) * k0i
            let rhs = gate.mul_add(ctx_gate, k1_at_gamma, k0i_constants[z], rhs);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = gate.mul_add(ctx_gate, r1i_at_gamma, qi_constants[z], rhs);

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

    use super::test_params;
    use crate::{
        constants::sk_enc_constants_4096_2x55_65537::R1_BOUNDS,
        sk_encryption_circuit::BfvSkEncryptionCircuit,
    };
    use axiom_eth::rlc::{circuit::builder::RlcCircuitBuilder, utils::executor::RlcExecutor};
    use halo2_base::{
        gates::circuit::CircuitBuilderStage,
        halo2_proofs::{
            dev::{FailureLocation, MockProver, VerifyFailure},
            halo2curves::bn256::Fr,
            plonk::{keygen_pk, keygen_vk, Any, SecondPhase},
        },
        utils::{
            fs::gen_srs,
            testing::{check_proof, gen_proof},
        },
    };
    use std::{fs::File, io::Read};

    #[cfg(feature = "bench")]
    use axiom_eth::halo2curves::bn256::Bn256;
    #[cfg(feature = "bench")]
    use halo2_base::halo2_proofs::poly::kzg::commitment::ParamsKZG;
    #[cfg(feature = "bench")]
    use prettytable::{row, Table};

    #[test]
    pub fn test_sk_enc_valid() {
        // 1. Define the inputs of the circuit
        let file_path = "src/data/sk_enc_4096_2x55_65537.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let sk_enc_circuit = serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

        // 2. Build the circuit for MockProver using the test parameters
        let rlc_circuit_params = test_params();
        let mut mock_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0)
                .use_params(rlc_circuit_params.clone());
        mock_builder.base.set_lookup_bits(8);

        let rlc_circuit = RlcExecutor::new(mock_builder, sk_enc_circuit);

        // 3. Run the mock prover. The circuit should be satisfied
        MockProver::run(
            rlc_circuit_params.base.k.try_into().unwrap(),
            &rlc_circuit,
            vec![],
        )
        .unwrap()
        .assert_satisfied();
    }

    #[test]
    pub fn test_sk_enc_full_prover() {
        // 1. Define the inputs of the circuit.
        // Since we are going to use this circuit instance for key gen, we can use an input file in which all the coefficients are set to 0
        let file_path_zeroes = "src/data/sk_enc_4096_2x55_65537_zeroes.json";
        let mut file = File::open(file_path_zeroes).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let empty_sk_enc_circuit = serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

        // 2. Generate (unsafe) trusted setup parameters
        // Here we are setting a small k for optimization purposes
        let k = 14;
        let kzg_params = gen_srs(k as u32);

        // 3. Build the circuit for key generation,
        let mut key_gen_builder =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Keygen, 0).use_k(k);
        key_gen_builder.base.set_lookup_bits(k - 1); // lookup bits set to `k-1` as suggested [here](https://docs.axiom.xyz/protocol/zero-knowledge-proofs/getting-started-with-halo2#technical-detail-how-to-choose-lookup_bits)

        let rlc_circuit = RlcExecutor::new(key_gen_builder, empty_sk_enc_circuit.clone());

        // The parameters are auto configured by halo2 lib to fit all the columns into the `k`-sized table
        let rlc_circuit_params = rlc_circuit.0.calculate_params(Some(9));

        // 4. Generate the verification key and the proving key
        let vk = keygen_vk(&kzg_params, &rlc_circuit).unwrap();
        let pk = keygen_pk(&kzg_params, vk, &rlc_circuit).unwrap();
        let break_points = rlc_circuit.0.builder.borrow().break_points();
        drop(rlc_circuit);

        // 5. Generate the proof, here we pass the actual inputs
        let mut proof_gen_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Prover, 0)
                .use_params(rlc_circuit_params);
        proof_gen_builder.base.set_lookup_bits(k - 1);

        let file_path = "src/data/sk_enc_4096_2x55_65537.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let sk_enc_circuit = serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

        let rlc_circuit = RlcExecutor::new(proof_gen_builder, sk_enc_circuit.clone());

        rlc_circuit
            .0
            .builder
            .borrow_mut()
            .set_break_points(break_points);
        let proof = gen_proof(&kzg_params, &pk, rlc_circuit);

        // 6. Verify the proof
        check_proof(&kzg_params, pk.get_vk(), &proof, true);
    }

    #[test]
    pub fn test_sk_enc_invalid_range() {
        // 1. Define the inputs of the circuit
        let file_path = "src/data/sk_enc_4096_2x55_65537.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let mut sk_enc_circuit = serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

        // 2. Invalidate the circuit by setting the value of a coefficient of the polynomial `r1is[0]` to be out of range
        let out_of_range_coeff = R1_BOUNDS[0] + 1;
        sk_enc_circuit.r1is[0][0] = out_of_range_coeff.to_string();

        // 3. Build the circuit for MockProver
        let rlc_circuit_params = test_params();
        let mut mock_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0)
                .use_params(rlc_circuit_params.clone());
        mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let rlc_circuit = RlcExecutor::new(mock_builder, sk_enc_circuit);

        // 4. Run the mock prover
        let invalid_mock_prover = MockProver::run(
            rlc_circuit_params.base.k.try_into().unwrap(),
            &rlc_circuit,
            vec![],
        )
        .unwrap();

        // 5. Assert that the circuit is not satisfied
        // In particular, it should fail the range check enforced in the second phase for the first coefficient of r1is[0] and the equality check in the second phase for the 0-th basis
        assert_eq!(
            invalid_mock_prover.verify(),
            Err(vec![
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 393180 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 393190 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 905111 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 905131 }
                },
            ])
        );
    }

    #[test]
    pub fn test_sk_enc_invalid_polys() {
        // 1. Define the inputs of the circuit
        let file_path = "src/data/sk_enc_4096_2x55_65537.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let mut sk_enc_circuit = serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

        // 2. Invalidate the circuit by setting a different `s` polynomial
        let invalid_s = vec!["1".to_string(); 1024];
        sk_enc_circuit.s = invalid_s;

        // 3. Build the circuit for MockProver
        let rlc_circuit_params = test_params();
        let mut mock_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0)
                .use_params(rlc_circuit_params.clone());
        mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let rlc_circuit = RlcExecutor::new(mock_builder, sk_enc_circuit);

        // 4. Run the mock prover
        let invalid_mock_prover = MockProver::run(
            rlc_circuit_params.base.k.try_into().unwrap(),
            &rlc_circuit,
            vec![],
        )
        .unwrap();

        // 5. Assert that the circuit is not satisfied
        // In particular, it should fail the equality check (LHS=RHS) in the second phase for each i-th CRT basis
        assert_eq!(
            invalid_mock_prover.verify(),
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
                    location: FailureLocation::OutsideRegion { row: 871319 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 871339 }
                },
                VerifyFailure::Permutation {
                    column: (Any::advice_in(SecondPhase), 1).into(),
                    location: FailureLocation::OutsideRegion { row: 871347 }
                },
            ])
        );
    }

    #[test]
    #[cfg(feature = "bench")]
    pub fn bench_sk_enc_full_prover() {
        let file_path = "src/data/sk_enc_4096_2x55_65537";

        pub struct Config {
            kzg_params: ParamsKZG<Bn256>,
            k: usize,
        }

        // Generate unsafe parameters for different values of k
        let mut configs = vec![];
        for k in 12..=18 {
            let kzg_params = gen_srs(k as u32);
            let config = Config { kzg_params, k };
            configs.push(config)
        }

        // Prepare a table to display results
        let mut table = Table::new();
        table.add_row(row![
            "K",
            "VK Generation Time",
            "PK Generation Time",
            "Proof Generation Time",
            "Proof Verification Time"
        ]);

        for config in &configs {
            println!("Running bench for k={}", config.k);
            // 1. Define the inputs of the circuit.
            // Since we are going to use this circuit instance for key gen, we can use an input file in which all the coefficients are set to 0
            let file_path_zeroes = format!("{}_zeroes.json", file_path);
            let mut file = File::open(file_path_zeroes).unwrap();
            let mut data = String::new();
            file.read_to_string(&mut data).unwrap();
            let empty_sk_enc_circuit =
                serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

            // 2. Build the circuit for key generation,
            let mut key_gen_builder =
                RlcCircuitBuilder::from_stage(CircuitBuilderStage::Keygen, 0).use_k(config.k);
            key_gen_builder.base.set_lookup_bits(config.k - 1); // lookup bits set to `k-1` as suggested [here](https://docs.axiom.xyz/protocol/zero-knowledge-proofs/getting-started-with-halo2#technical-detail-how-to-choose-lookup_bits)

            let rlc_circuit = RlcExecutor::new(key_gen_builder, empty_sk_enc_circuit.clone());

            // The parameters are auto configured by halo2 lib to fit all the columns into the `k`-sized table
            let rlc_circuit_params = rlc_circuit.0.calculate_params(Some(9));

            // 3. Generate the verification key and the proving key
            let timer = std::time::Instant::now();
            let vk = keygen_vk(&config.kzg_params, &rlc_circuit).unwrap();
            let vk_gen_time = timer.elapsed();
            let timer = std::time::Instant::now();
            let pk = keygen_pk(&config.kzg_params, vk, &rlc_circuit).unwrap();
            let pk_gen_time = timer.elapsed();
            let break_points = rlc_circuit.0.builder.borrow().break_points();
            drop(rlc_circuit);

            // 4. Generate the proof, here we pass the actual inputs
            let mut proof_gen_builder: RlcCircuitBuilder<Fr> =
                RlcCircuitBuilder::from_stage(CircuitBuilderStage::Prover, 0)
                    .use_params(rlc_circuit_params);
            proof_gen_builder.base.set_lookup_bits(config.k - 1);

            let file_path = format!("{}.json", file_path);
            let mut file = File::open(file_path).unwrap();
            let mut data = String::new();
            file.read_to_string(&mut data).unwrap();
            let sk_enc_circuit = serde_json::from_str::<BfvSkEncryptionCircuit>(&data).unwrap();

            let rlc_circuit = RlcExecutor::new(proof_gen_builder, sk_enc_circuit.clone());

            rlc_circuit
                .0
                .builder
                .borrow_mut()
                .set_break_points(break_points);
            let timer = std::time::Instant::now();
            let proof = gen_proof(&config.kzg_params, &pk, rlc_circuit);
            let proof_gen_time = timer.elapsed();

            // 5. Verify the proof
            let timer = std::time::Instant::now();
            check_proof(&config.kzg_params, pk.get_vk(), &proof, true);
            let proof_verification_time = timer.elapsed();

            table.add_row(row![
                config.k,
                format!("{:?}", vk_gen_time),
                format!("{:?}", pk_gen_time),
                format!("{:?}", proof_gen_time),
                format!("{:?}", proof_verification_time),
            ]);
        }
        println!("bfv params: {:?}", file_path);
        table.printstd();
    }
}
