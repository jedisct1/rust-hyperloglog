// (C)opyleft 2013-2019 Frank Denis

//! HyperLogLog implementation for Rust
#![crate_name = "hyperloglog"]
#![warn(non_camel_case_types, non_upper_case_globals, unused_qualifications)]
#![allow(non_snake_case)]
#![allow(clippy::unreadable_literal)]

mod data;
use rand;

use siphasher::sip::SipHasher13;
use std::cmp::Ordering::{Equal, Greater, Less};
use std::hash::{Hash, Hasher};
use std::iter::repeat;
use std::marker::PhantomData;

#[cfg_attr(feature = "with_serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HyperLogLog<V> {
    alpha: f64,
    p: u8,
    m: usize,
    M: Vec<u8>,
    key0: u64,
    key1: u64,
    #[cfg_attr(feature = "with_serde", serde(skip))]
    sip: SipHasher13,
    v_phantom: PhantomData<V>,
}

impl<V> HyperLogLog<V>
where
    V: Hash,
{
    pub fn new(error_rate: f64) -> Self {
        let key0 = rand::random();
        let key1 = rand::random();
        Self::new_from_keys(error_rate, key0, key1)
    }

    pub fn new_from_keys(error_rate: f64, key0: u64, key1: u64) -> Self {
        assert!(error_rate > 0.0 && error_rate < 1.0);
        let sr = 1.04 / error_rate;
        let p = f64::ln(sr * sr).ceil() as u8;
        assert!(p <= 64);
        let alpha = Self::get_alpha(p);
        let m = 1usize << p;

        HyperLogLog {
            alpha,
            p,
            m,
            M: repeat(0u8).take(m).collect(),
            key0,
            key1,
            sip: SipHasher13::new_with_keys(key0, key1),
            v_phantom: PhantomData,
        }
    }

    pub fn new_from_template(hll: &HyperLogLog<V>) -> Self {
        HyperLogLog {
            alpha: hll.alpha,
            p: hll.p,
            m: hll.m,
            key0: hll.key0,
            key1: hll.key1,
            M: repeat(0u8).take(hll.m).collect(),
            sip: hll.sip,
            v_phantom: PhantomData,
        }
    }

    pub fn insert(&mut self, value: &V) {
        let sip = &mut self.sip.clone();
        value.hash(sip);
        let x = sip.finish();
        let j = x as usize & (self.m - 1);
        let w = x >> self.p;
        let rho = Self::get_rho(w, 64 - self.p);
        let mjr = &mut self.M[j];
        if rho > *mjr {
            *mjr = rho;
        }
    }

    pub fn len(&self) -> f64 {
        let V = Self::vec_count_zero(&self.M);
        if V > 0 {
            let H = self.m as f64 * (self.m as f64 / V as f64).ln();
            if H <= Self::get_threshold(self.p) {
                H
            } else {
                self.ep()
            }
        } else {
            self.ep()
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0.0
    }

    pub fn merge(&mut self, src: &HyperLogLog<V>) {
        assert!(src.p == self.p);
        assert!(src.m == self.m);
        let sip1 = &mut src.sip.clone();
        let sip2 = &mut self.sip.clone();
        42.hash(sip1);
        42.hash(sip2);
        assert!(sip1.finish() == sip2.finish());
        for i in 0..self.m {
            let (src_mir, mir) = (src.M[i], &mut self.M[i]);
            if src_mir > *mir {
                *mir = src_mir;
            }
        }
    }

    pub fn clear(&mut self) {
        self.M.iter_mut().all(|x| {
            *x = 0;
            true
        });
    }

    fn get_threshold(p: u8) -> f64 {
        data::THRESHOLD_DATA[p as usize]
    }

    fn get_alpha(p: u8) -> f64 {
        assert!(p >= 4 && p <= 16);
        match p {
            4 => 0.673,
            5 => 0.697,
            6 => 0.709,
            _ => 0.7213 / (1.0 + 1.079 / (1usize << (p as usize)) as f64),
        }
    }

    fn bit_length(x: u64) -> u8 {
        let mut bits: u8 = 0;
        let mut xm = x;
        while xm != 0 {
            bits += 1;
            xm >>= 1;
        }
        bits
    }

    fn get_rho(w: u64, max_width: u8) -> u8 {
        let rho = max_width - Self::bit_length(w) + 1;
        assert!(rho > 0);
        rho
    }

    fn vec_count_zero(v: &[u8]) -> usize {
        bytecount::count(v, 0)
    }

    fn estimate_bias(E: f64, p: u8) -> f64 {
        let bias_vector = data::BIAS_DATA[(p - 4) as usize];
        let nearest_neighbors =
            Self::get_nearest_neighbors(E, data::RAW_ESTIMATE_DATA[(p - 4) as usize]);
        let sum = nearest_neighbors
            .iter()
            .fold(0.0, |acc, &neighbor| acc + bias_vector[neighbor]);
        sum / nearest_neighbors.len() as f64
    }

    fn get_nearest_neighbors(E: f64, estimate_vector: &[f64]) -> Vec<usize> {
        let ev_len = estimate_vector.len();
        let mut r: Vec<(f64, usize)> = repeat((0.0f64, 0usize)).take(ev_len).collect();
        for i in 0..ev_len {
            let dr = E - estimate_vector[i];
            r[i] = (dr * dr, i);
        }
        r.sort_by(|a, b| {
            if a < b {
                Less
            } else if a > b {
                Greater
            } else {
                Equal
            }
        });
        r.truncate(6);
        r.iter()
            .map(|&ez| match ez {
                (_, b) => b,
            })
            .collect()
    }

    fn ep(&self) -> f64 {
        let sum = self
            .M
            .iter()
            .fold(0.0, |acc, &x| acc + 2.0f64.powi(-(x as i32)));
        let E = self.alpha * (self.m * self.m) as f64 / sum;
        if E <= (5 * self.m) as f64 {
            E - Self::estimate_bias(E, self.p)
        } else {
            E
        }
    }
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

#[cfg(feature = "with_serde")]
#[test]
fn hyperloglog_test_serialize() {
    let mut hll = HyperLogLog::new(0.00408);
    let keys = ["test1", "test2", "test3", "test2", "test2", "test2"];
    for k in &keys {
        hll.insert(k);
    }

    let hll_str: String = serde_json::to_string(&hll).unwrap();

    let hll_de = match serde_json::from_str::<HyperLogLog<String>>(&hll_str) {
        Ok(hll) => hll,
        Err(_) => panic!("Failed to deserialize"),
    };

    assert!((hll.len() - hll_de.len()).abs() < std::f64::EPSILON);
}

#[test]
fn hyperloglog_test_merge_with_keys() {
    let key0 = rand::random();
    let key1 = rand::random();

    let mut hll = HyperLogLog::new_from_keys(0.00408, key0, key1);
    let keys = ["test1", "test2", "test3", "test2", "test2", "test2"];
    for k in &keys {
        hll.insert(k);
    }
    assert!((hll.len().round() - 3.0).abs() < std::f64::EPSILON);

    let mut hll2 = HyperLogLog::new_from_keys(0.00408, key0, key1);
    let keys2 = ["test3", "test4", "test4", "test4", "test4", "test1"];
    for k in &keys2 {
        hll2.insert(k);
    }
    assert!((hll2.len().round() - 3.0).abs() < std::f64::EPSILON);

    hll.merge(&hll2);
    assert!((hll.len().round() - 4.0).abs() < std::f64::EPSILON);
}
