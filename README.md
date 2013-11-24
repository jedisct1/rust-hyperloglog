hyperloglog
===========

A [HyperLogLog](https://static.googleusercontent.com/media/research.google.com/en/us/pubs/archive/40671.pdf) implementation in Rust, with bias correction.

Installation:

```bash
rustpkg install github.com/jedisct1/rust-hyperloglog
```

Usage:

```rust
let mut hll = HyperLogLog::new(error_rate);
hll.add(~"test1");
hll.add(~"test2");
let card_estimation = hll.card();
```
