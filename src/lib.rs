// (C)opyleft 2013-2021 Frank Denis

//! HyperLogLog implementation for Rust
#![crate_name = "hyperloglog"]
#![warn(non_camel_case_types, non_upper_case_globals, unused_qualifications)]
#![allow(non_snake_case)]
#![allow(clippy::unreadable_literal)]

mod weights;
use weights::{BIAS_DATA, RAW_ESTIMATE_DATA, THRESHOLD_DATA};
use std::hash::{Hash, Hasher};
use siphasher::sip::SipHasher13;

/// A HyperLogLog counter
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "mem_dbg", derive(mem_dbg::MemDbg, mem_dbg::MemSize))]
pub struct HyperLogLog {
    alpha: f64,
    p: u8,
    number_of_registers: usize,
    registers: Vec<u8>,
    sip: SipHasher13,
}

impl HyperLogLog {
    /// Create a new `HyperLogLog` counter with the given error rate and seed.
    #[must_use]
    pub fn new_deterministic(error_rate: f64, seed: u128) -> Self {
        let key0 = (seed >> 64) as u64;
        let key1 = seed as u64;
        assert!(error_rate > 0.0 && error_rate < 1.0);
        let p = (f64::log2(1.04 / error_rate) * 2.0).ceil() as u8;
        assert!(p <= 18);
        assert!(p >= 4);
        let alpha = Self::get_alpha(p);
        let number_of_registers = 1usize << p;
        HyperLogLog {
            alpha,
            p,
            number_of_registers,
            registers: vec![0; number_of_registers],
            sip: SipHasher13::new_with_keys(key0, key1),
        }
    }

    /// Create a new `HyperLogLog` counter with the given error rate and a random
    /// seed.
    #[must_use]
    pub fn new(error_rate: f64) -> Self {
        let seed: u128 = rand::random();
        Self::new_deterministic(error_rate, seed)
    }

    /// Create a new `HyperLogLog` counter with the same parameters as an
    /// existing one.
    #[must_use]
    pub fn new_from_template(hll: &HyperLogLog) -> Self {
        HyperLogLog {
            alpha: hll.alpha,
            p: hll.p,
            number_of_registers: hll.number_of_registers,
            registers: vec![0; hll.number_of_registers],
            sip: hll.sip,
        }
    }

    /// Insert a new value into the `HyperLogLog` counter.
    pub fn insert<V: Hash>(&mut self, value: &V) {
        let mut sip = self.sip;
        value.hash(&mut sip);
        let x = sip.finish();
        self.insert_by_hash_value(x);
    }

    /// Insert a new u64 value into the `HyperLogLog` counter.
    pub fn insert_by_hash_value(&mut self, x: u64) {
        let j = x as usize & (self.number_of_registers - 1);
        let w = x >> self.p;
        let rho = Self::get_rho(w, 64 - self.p);
        let mjr = &mut self.registers[j];
        if rho > *mjr {
            *mjr = rho;
        }
    }

    /// Return the cardinality of the `HyperLogLog` counter.
    #[must_use]
    pub fn len(&self) -> f64 {
        let number_of_zero_registers = bytecount::count(&self.registers, 0);
        if number_of_zero_registers > 0 {
            let estimate = self.number_of_registers as f64 * (self.number_of_registers as f64 / number_of_zero_registers as f64).ln();
            if estimate <= Self::get_threshold(self.p) {
                return estimate
            }
        }
        self.ep()
    }

    /// Return `true` if the `HyperLogLog` counter is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0.0
    }

    /// Merge another `HyperLogLog` counter into the current one.
    pub fn merge(&mut self, src: &HyperLogLog) {
        assert!(src.p == self.p);
        assert!(src.number_of_registers == self.number_of_registers);
        let mut sip1 = src.sip;
        let mut sip2 = self.sip;
        42.hash(&mut sip1);
        42.hash(&mut sip2);
        assert_eq!(sip1.finish(), sip2.finish(), "The two SipHasher do not seem to have the same seed - Use new_deterministic instead of new to create the HyperLogLog.");
        for i in 0..self.number_of_registers {
            let (src_mir, mir) = (src.registers[i], &mut self.registers[i]);
            if src_mir > *mir {
                *mir = src_mir;
            }
        }
    }

    /// Wipe the `HyperLogLog` counter.
    pub fn clear(&mut self) {
        self.registers.fill(0);
    }

    fn get_threshold(p: u8) -> f64 {
        THRESHOLD_DATA[p as usize - 4]
    }

    pub fn precision(&self) -> u8 {
        self.p
    }

    fn get_alpha(p: u8) -> f64 {
        assert!(p >= 4);
        assert!(p <= 18);
        match p {
            4 => 0.673,
            5 => 0.697,
            6 => 0.709,
            _ => 0.7213 / (1.0 + 1.079 / (1usize << (p as usize)) as f64),
        }
    }

    fn bit_length(x: u64) -> u8 {
        (64 - x.leading_zeros()) as u8
    }

    fn get_rho(w: u64, max_width: u8) -> u8 {
        let rho = max_width - Self::bit_length(w) + 1;
        assert!(rho > 0);
        rho
    }

    fn estimate_bias(estimate: f64, p: u8) -> f64 {
        let bias_vector = BIAS_DATA[(p - 4) as usize];
        let estimate_vector = RAW_ESTIMATE_DATA[(p - 4) as usize];

        // Since the estimates are sorted, we can use a partition point to find the nearest neighbors
        let partition_point = estimate_vector.partition_point(|&x| x < estimate);

        let mut min = if partition_point > 6 {
            partition_point - 6
        } else {
            0
        };
        let mut max = core::cmp::min(partition_point + 6, estimate_vector.len());

        while max - min != 6 {
            let (min_val, max_val) = (estimate_vector[min], estimate_vector[max - 1]);
            if 2.0 * estimate - min_val > max_val {
                min += 1;
            } else {
                max -= 1;
            }
        }

        (min..max).map(|i| bias_vector[i]).sum::<f64>() / 6.0
    }

    fn ep(&self) -> f64 {
        let sum: f64 = self.registers.iter().map(|&x| 2.0f64.powi(-(x as i32))).sum();
        let estimate = self.alpha * (self.number_of_registers * self.number_of_registers) as f64 / sum;
        if estimate <= (5 * self.number_of_registers) as f64 {
            estimate - Self::estimate_bias(estimate, self.p)
        } else {
            estimate
        }
    }
}

#[cfg(feature = "serde")]
#[test]
fn hyperloglog_serialize() {
    let hll = HyperLogLog::new(0.00408);
    let bytes = bincode::serialize(&hll).unwrap();
    let _: HyperLogLog = bincode::deserialize(&bytes).unwrap();
}

#[test]
fn hyperloglog_test_simple() {
    let mut hll = HyperLogLog::new(0.00408);
    let keys = ["test1", "test2", "test3", "test2", "test2", "test2"];
    for k in &keys {
        hll.insert(k);
    }
    assert!((hll.len().round() - 3.0).abs() < std::f64::EPSILON);
    assert!(!hll.is_empty());
    hll.clear();
    assert!(hll.is_empty());
    assert!(hll.len() == 0.0);
}

#[test]
fn hyperloglog_test_merge() {
    let mut hll = HyperLogLog::new(0.00408);
    let keys = ["test1", "test2", "test3", "test2", "test2", "test2"];
    for k in &keys {
        hll.insert(k);
    }
    assert!((hll.len().round() - 3.0).abs() < std::f64::EPSILON);

    let mut hll2 = HyperLogLog::new_from_template(&hll);
    let keys2 = ["test3", "test4", "test4", "test4", "test4", "test1"];
    for k in &keys2 {
        hll2.insert(k);
    }
    assert!((hll2.len().round() - 3.0).abs() < std::f64::EPSILON);

    hll.merge(&hll2);
    assert!((hll.len().round() - 4.0).abs() < std::f64::EPSILON);
}

