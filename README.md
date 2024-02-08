# Greco

Circuit for proving the correct encryption under BFV fully homomorphic encryption scheme. Note that this can be also generalized to any RLWE-based FHE scheme. Based on https://hackmd.io/@gaussian/r1W98Kqqa.

### Generate Parameters

To generate the parameters for the secret key proof of encryption circuit run the following command:

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
]' -t 65537  -output_input ./src/data/sk_enc_input.json -output_constants ./src/constants/sk_enc.rs
```

Where `-n` is the degree of the polynomial, `-qis` is the list of moduli qis such that qis[i] is the modulus of the i-th CRT basis of the modulus q of the ciphertext space, `-t` is the plaintext modulus and `-output_input` is path to the json file containing the circuit input and `-output_constants` is the path to the rust file containing the circuit generic constants. The value of `𝜎` for the gaussian distribution is set to 3.2 by default.

You can modify these parameters to fit your needs.

As a result:
- A file `./src/data/sk_enc_input.json` is generated including the input to the circuit that can be used for testing. It includes a random secret key, a random plaintext message and the corresponding ciphertext encrypted under the secret key.
- A file `./src/constants/sk_enc.rs` is generated including the generic constants for the circuit. Note that we separate them from the input because these should be known at compile time.

### Circuit

```
cargo build
cargo test --release -- --nocapture
```

The halo2 circuit is based on a fork of `axiom-eth` that implements two minor changes:

- `RlcCircuit` and `RlcExecutor` are included into a utils mod such that they can be consumed outside of the crate 
- The `RlcCircuitInstructions` are modified to enable equality constraints on instance column in Second Phase

Further testing, incorporating the whole flow of generating random parameters and random input and generating a proof can be run with:

```
python3 scripts/test.py 20
```

Where `20` is number of times the test should be run. Any error is added to the `scripts/error_log.txt` file.