#![no_main]
use libfuzzer_sys::fuzz_target;
extern crate hcompress;
use crate::hcompress::*;
use crate::hcompress::lib2::*;

fuzz_target!(|data: Data| {
    let mut d = data.d.iter().map(|&x| x as i32).collect::<Vec<i32>>();
    let mut input: Vec<i32> = d.clone();

    let y = 1;
    let x = input.len();
    let scale = 0; //lossless

    if x < 10 { return }

    let mut compressed: Vec<u8> = Vec::with_capacity(x*y*100);
    let res = fits_hcompress(&mut input, y, x, scale, &mut compressed);

    let mut uncompressed: Vec<i32> = vec![0;x*y];
    let res = fits_hdecompress(&compressed, 0, &mut uncompressed);

    assert_eq!(d, uncompressed);
});
