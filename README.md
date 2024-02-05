# halo2-fhe

Circuit for proving the correct encryption under BFV fully homomorphic encryption scheme. Note that this can be also generalized to any RLWE-based FHE scheme. Based on https://hackmd.io/@gaussian/r1W98Kqqa.

### POC

The repository contains a python proof of concet for the secret key proof of encryption and for public key proof of encryption circuits

```python
cd python_poc
python3 circuit_sk.py
python3 circuit_pk.py
```

### Halo2

```
cargo build
cargo test --release -- --nocapture
```