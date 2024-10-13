#![no_main]
use libfuzzer_sys::fuzz_target;
extern crate hcompress;
use crate::hcompress::*;

fuzz_target!(|data: Data16| {
    let mut d = data.d.iter().map(|&x| x as i32).collect::<Vec<i32>>();
    let mut input: Vec<i32> = d.clone();

    let x = 1;
    let y = input.len();
    let scale = 0; //lossless

    if y < 10 { return }

    let mut compressed: Vec<u8> = Vec::with_capacity(x*y*100);
    let mut encoder = hcompress::write::HCEncoder::new();
    let res = encoder.write(&mut input, y, x, scale, &mut compressed);

    let mut uncompressed: Vec<i32> = vec![0;x*y];
    let mut decoder = hcompress::read::HCDecoder::new();
    let res = decoder.read(&compressed, 0, &mut uncompressed);

    assert_eq!(d, uncompressed);
});
