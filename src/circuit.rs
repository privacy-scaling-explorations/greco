use crate::utils::evaluate_poly;
use axiom_eth::rlc::{
    chip::RlcChip,
    circuit::{builder::RlcCircuitBuilder, instructions::RlcCircuitInstructions, RlcCircuitParams},
};
use halo2_base::{
    gates::{circuit::BaseCircuitParams, RangeChip, RangeInstructions},
    utils::ScalarField,
    AssignedValue,
    QuantumCell::{Constant, Existing},
};

const K: usize = 16;

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

/// DummyCircuit is a simple circuit that takes two polynomials `a` and `b` and computes their product `c = a * b`.
///
/// * `a` polynomial coefficients where coefficients
/// * `b` polynomial coefficients where coefficients
/// * `c` polynomial coefficients such that `c = a * b`
pub struct DummyCircuit<F: ScalarField> {
    a: Vec<F>,
    b: Vec<F>,
    c: Vec<F>,
}

/// Payload returned by the first phase of the circuit to be reused in the second phase
pub struct Payload<F: ScalarField> {
    a: Vec<F>,
    b: Vec<F>,
    c: Vec<F>,
    a_assigned: Vec<AssignedValue<F>>,
    b_assigned: Vec<AssignedValue<F>>,
    c_assigned: Vec<AssignedValue<F>>,
}

impl<F: ScalarField> RlcCircuitInstructions<F> for DummyCircuit<F> {
    type FirstPhasePayload = Payload<F>;

    // Phase 0
    // * Assign the polynomials to the circuit
    // * Range check the coefficients of the polynomials
    fn virtual_assign_phase0(
        &self,
        builder: &mut RlcCircuitBuilder<F>,
        range: &RangeChip<F>,
    ) -> Self::FirstPhasePayload {
        let ctx = builder.base.main(0);
        let a_assigned = ctx.assign_witnesses(self.a.clone());
        let b_assigned = ctx.assign_witnesses(self.b.clone());
        let c_assigned = ctx.assign_witnesses(self.c.clone());

        // a coefficients should be max 4 bits
        for a_coeff in a_assigned.iter() {
            range.range_check(ctx, *a_coeff, 4)
        }

        // b coefficients should be max 4 bits
        for b_coeff in b_assigned.iter() {
            range.range_check(ctx, *b_coeff, 4)
        }

        // when multiplying two polynomials of same degree `n`
        // the coefficient corresponding to `x^n` is the sum of `n+1` products between coefficients of `a` and `b`
        // the max value of the product of two 4 bits coefficients is 8 bits
        // each sum could add a furher bit to the resulting coefficient of c
        let c_max_bits = 8 + (self.a.len() + 1);

        for c_coeff in c_assigned.iter() {
            range.range_check(ctx, *c_coeff, c_max_bits)
        }

        Payload {
            a: self.a.clone(),
            b: self.b.clone(),
            c: self.c.clone(),
            a_assigned,
            b_assigned,
            c_assigned,
        }
    }

    // Phase 1
    // * Evaluate each polynomial at the point `gamma`
    // * Assert that `a(gamma) * b(gamma) = c(gamma)`
    fn virtual_assign_phase1(
        builder: &mut RlcCircuitBuilder<F>,
        _: &RangeChip<F>,
        rlc: &RlcChip<F>,
        payload: Self::FirstPhasePayload,
    ) {
        let Payload {
            a,
            b,
            c,
            a_assigned,
            b_assigned,
            c_assigned,
        } = payload;

        let (ctx_gate, ctx_rlc) = builder.rlc_ctx_pair();
        let gamma = *rlc.gamma();

        let rlc_trace_a = rlc.compute_rlc_fixed_len(ctx_rlc, a_assigned);
        let rlc_trace_b = rlc.compute_rlc_fixed_len(ctx_rlc, b_assigned);
        let rlc_trace_c = rlc.compute_rlc_fixed_len(ctx_rlc, c_assigned);

        // enforce gate a(gamma) * b(gamma) - c(gamma) = 0
        ctx_gate.assign_region(
            [
                Constant(F::from(0)),
                Existing(rlc_trace_a.rlc_val),
                Existing(rlc_trace_b.rlc_val),
                Existing(rlc_trace_c.rlc_val),
            ],
            [0],
        );

        // Assert the correctness of the RLC
        let expected_poly_a_eval_at_gamma = evaluate_poly(&a, gamma);
        let expected_poly_b_eval_at_gamma = evaluate_poly(&b, gamma);
        let expected_poly_c_eval_at_gamma = evaluate_poly(&c, gamma);

        assert_eq!(rlc_trace_a.rlc_val.value(), &expected_poly_a_eval_at_gamma);
        assert_eq!(rlc_trace_b.rlc_val.value(), &expected_poly_b_eval_at_gamma);
        assert_eq!(rlc_trace_c.rlc_val.value(), &expected_poly_c_eval_at_gamma);
    }
}

#[cfg(test)]
mod test {

    use super::{test_params, K};
    use crate::circuit::DummyCircuit;
    use crate::utils::mul;
    use axiom_eth::rlc::{circuit::builder::RlcCircuitBuilder, utils::executor::RlcExecutor};
    use halo2_base::{
        gates::circuit::CircuitBuilderStage,
        halo2_proofs::{
            dev::{FailureLocation, MockProver, VerifyFailure},
            halo2curves::bn256::Fr,
        },
    };
    use rand::Rng;

    #[test]
    pub fn test_dummy_circuit_valid() {
        let mut rng = rand::thread_rng();

        // 1. Define the input polynomials
        // The polynomials are defined in coefficients form starting from the highest degree term.
        let a = vec![
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
        ];
        let b = vec![
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
        ];
        let c = mul(a.clone(), b.clone());

        // 2. Build the circuit
        let params = test_params();
        let mut builder =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0).use_params(params);
        builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let circuit = RlcExecutor::new(builder, DummyCircuit { a, b, c });

        // 3. Run the mock prover
        MockProver::run(K as u32, &circuit, vec![])
            .unwrap()
            .assert_satisfied();
    }

    #[test]
    pub fn test_dummy_circuit_invalid_range() {
        let mut rng = rand::thread_rng();

        // 1. Define the input polynomials
        // The polynomials are defined in coefficients form starting from the highest degree term.
        // A coefficient of `b` is greater than 4 bits and will cause the circuit to fail
        let a = vec![
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
        ];
        let out_of_range_b = vec![
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(16), // 5 bits
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
        ];
        let c = mul(a.clone(), out_of_range_b.clone());

        // 2. Build the circuit
        let params = test_params();
        let mut builder =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0).use_params(params);
        builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let circuit = RlcExecutor::new(
            builder,
            DummyCircuit {
                a,
                b: out_of_range_b,
                c,
            },
        );

        // 3. Run the mock prover
        let invalid_prover = MockProver::run(K as u32, &circuit, vec![]).unwrap();

        // 4. Assert that the circuit is not satisfied
        assert_eq!(
            invalid_prover.verify(),
            Err(vec![VerifyFailure::Lookup {
                name: "lookup".to_string(),
                lookup_index: 0,
                location: FailureLocation::OutsideRegion { row: 38 }
            },])
        );
    }

    #[test]
    pub fn test_dummy_circuit_invalid_poly_mult() {
        let mut rng = rand::thread_rng();

        // 1. Define the input polynomials
        // The polynomials are defined in coefficients form starting from the highest degree term.
        let a = vec![
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
        ];
        let b = vec![
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
            Fr::from(rng.gen_range(0..16) as u64),
        ];
        let mut invalid_c = mul(a.clone(), b.clone());
        invalid_c[0] = Fr::from(rng.gen_range(0..16) as u64); // Change the first coefficient of c

        // 2. Build the circuit
        let params = test_params();
        let mut builder =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0).use_params(params);
        builder.base.set_lookup_bits(8); // Set the lookup bits to 8

        let circuit = RlcExecutor::new(builder, DummyCircuit { a, b, c: invalid_c });

        // 3. Run the mock prover
        let invalid_prover = MockProver::run(K as u32, &circuit, vec![]).unwrap();

        // 4. Assert that the circuit is not satisfied
        // Uncomment it and the test should fail. TODO: Add better handling of the error
        // invalid_prover.assert_satisfied();
    }
}
