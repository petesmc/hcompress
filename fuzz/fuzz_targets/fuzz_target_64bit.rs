#![no_main]
use libfuzzer_sys::fuzz_target;
use bytemuck::cast_slice_mut;
extern crate hcompress;
use crate::hcompress::*;

fuzz_target!(|data: Data32| {
    let mut d = data.d.iter().map(|&x| x as i64).collect::<Vec<i64>>();
    let mut input: Vec<i64> = d.clone();

    let y = 1;
    let x = input.len();
    let scale = 0; //lossless

    if x < 10 { return }

    let mut compressed: Vec<u8> = Vec::with_capacity(x*y*100);
    let mut encoder = hcompress::write::HCEncoder::new();

    let res = encoder.write64(&mut input, y, x, scale, &mut compressed);

    let mut uncompressed: Vec<i32> = vec![0;x*y*2];
    let mut decoder = hcompress::read::HCDecoder::new();
    let res = decoder.read64(&compressed, 0, cast_slice_mut(&mut uncompressed));

    assert_eq!(data.d, uncompressed[..(x*y)]);
});
