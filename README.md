# Greco

Circuit for proving the correct encryption under BFV fully homomorphic encryption scheme. Note that this can be also generalized to any RLWE-based FHE scheme. Based on https://hackmd.io/@gaussian/r1W98Kqqa.

### Generate Parameters

To generate the parameters for the secret key proof of encryption run the following command:

```bash
python3 scripts/circuit_sk.py -n 1024 -qis '[                                                   
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
    1152921504580894721
]' -t 65537  -output ./src/data/sk_enc_input.json
```

Where `-n` is the degree of the polynomial, `-qis` is the list of moduli qis such that qis[i] is the modulus of the i-th CRT basis of the modulus q of the ciphertext space, `-t` is the plaintext modulus and `-output` is path to the output file. The value of `𝜎` for the gaussian distribution is set to 3.2 by default.

You can modify these parameters to fit your needs.

As a result, the file `./src/data/sk_enc_input.json` will be generated including the input to the circuit. On top of that, the terminal will print the constants that need to be added to the halo2 circuit.

```rust
pub const N: usize = 1024;
pub const E_BOUND: u64 = 19;
pub const S_BOUND: u64 = 1;
pub const R1_BOUNDS: [u64; 15] = [1321, 12139, 1692, 1530, 19009, 17587, 3417, 15539, 24450, 19013, 24041, 5934, 31437, 16662, 15909];
pub const R2_BOUNDS: [u64; 15] = [576460752303292416, 576460752299360256, 576460752298508288, 576460752297984000, 576460752297820160, 576460752296706048, 576460752296411136, 576460752296214528, 576460752294969344, 576460752293265408, 576460752292773888, 576460752291823616, 576460752290938880, 576460752290709504, 576460752290447360];
pub const K1_BOUND: u64 = 32768;
```

These constants are then to be added to `src/constants/sk_enc.rs` file. Note that we define them as constants and not as inputs because these should be known at compile time.

### Circuit

```
cargo build
cargo test --release -- --nocapture
```

The halo2 circuit is based on a fork of `axiom-eth` that implements two minor changes:

- `RlcCircuit` and `RlcExecutor` are included into a utils mod such that they can be consumed outside of the crate 
- The `RlcCircuitInstructions` are modified to enable equality constraints on instance column in Second Phase