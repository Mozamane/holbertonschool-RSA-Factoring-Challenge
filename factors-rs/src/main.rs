// (C) 2019 Fernando Gonzalez-Morales
// That's right.
#![allow(irrefutable_let_patterns)]
use std::{env, fs};
use ramp::Int;
use rayon::prelude::*;


fn main() {
    let natural: Vec<Int> = [1,2, 3,4, 5,6, 7,8,9,10, 11,12, 13,14,15,16, 17,18, 19,20,21,22, 23,24,25,26,27,28, 29,30,
                            31,32,33,34,35,36, 37,38,39,40, 41,42, 43,44,45,46, 47,48,49,50,51,52, 53,54,55,56,57,58,59,60, 61,62,63,64,65,66, 67,68,69,70,
                            71,72, 73,74,75,76,77,78, 79,80, 81, 82, 83,84,85,86,87,88, 89,90,91,92,93,94,95,96, 97].into_iter()
                                                   .map(|x| Int::from(*x))
                                                   .collect();
    let filename = env::args().nth(1).unwrap();
    let contents = fs::read_to_string(filename).unwrap();
    let numbers: Vec<Int> = contents.par_lines()
                            .filter_map(|line| line.parse().ok())
                            .collect();
    let (small, large): (Vec<_>, Vec<_>) = numbers.into_par_iter()
                                                  .partition(|n| *n < 10000);
    large.par_iter().for_each(|n| {
        let (p, q) = pollard_brent(n);
        println!("{}={}*{}", n, p, q);
    });
    small.par_iter().for_each(|n| {
        let (p, q) = look_up(n, &primes);
        println!("{}={}*{}", n, p, q);
    });
}

/* https://stackoverflow.com/a/2274520/9221785 */
fn look_up(n: &Int, ps: &Vec<Int>) -> (Int, Int) {
    for p in ps.iter() {
		if let (q, r) = n.divmod(p) {
            if r == 0 { return (p.clone(), q) };
        }
    }
    (Int::one(), n.clone())
}

/* Refs:
 * +http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.117.1230&rep=rep1&type=pdf
 * +https://comeoncodeon.wordpress.com/2010/09/18/pollard-rho-brent-integer-factorization/
 */
fn pollard_brent(num: &Int) -> (Int, Int) {
    let (mut x, mut y) = (Int::from(2), Int::from(2));
    let (mut d, mut q) = (Int::one(), Int::one());
    let (mut ys, mut r) = (Int::zero(), 1);
    const M: i32 = 71;
    let f = |i: &Int| (i.pow_mod(&Int::from(2), &num) + M) % num;
    while d == 1 {
        x = f(&y);
        for _ in 0..r {
            y = f(&y);
        }
        let mut k = 0_i32;
        while k < r && d == 1 {
            for _ in 0..(M.min(r - k)) {
                y = f(&y);
                q = q * (&x - &y).abs() % num;
            }
            ys = y.clone();
            d = q.gcd(num);
            k += M;
        }
        r *= 2;
    }
    if d == *num {
        loop {
            ys = f(&ys);
            d = ((&x - &ys).abs()).gcd(num);
            if d > 1 {
                break;
            }
        }
    }
    (d.clone(), num/d)
}
