hyperloglog
===========

A [HyperLogLog](https://static.googleusercontent.com/media/research.google.com/en/us/pubs/archive/40671.pdf) implementation in Rust, with bias correction.

Installation: use [Cargo](http://crates.io):

```toml
[dependencies]
hyperloglog = "0"
```

Usage:

```rust
let mut hll = HyperLogLog::new(error_rate);
hll.insert(&"test1");
hll.insert(&"test2");
let card_estimation = hll.len();

let mut hll2 = HyperLogLog::new_from_template(&hll);
hll2.insert(&"test3");

hll.merge(&hll2);
```

You can save the intermediate HLL object if you enable the `with_serde` feature:

```toml
[dependencies]
hyperloglog = { version = "0", features = "with_serde" }
```

Then you can serialize/deserialize the HyperLogLog struct:

```rust
let hll_ser: String = hll.serialize();
let hll_de: HyperLogLog::<String>::deserialize(&hll_ser);
```
