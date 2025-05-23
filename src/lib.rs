#[cfg(test)]
#[macro_use]
extern crate quickcheck;

pub mod read;
pub mod write;

pub(crate) const CODE_MAGIC: [u8; 2] = [0xDD, 0x99];

#[derive(Clone, Debug)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct Data16 {
    pub d: Vec<i16>,
}

#[derive(Clone, Debug)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct Data32 {
    pub d: Vec<i32>,
}

#[cfg(test)]
mod tests {

    use bytemuck::cast_slice_mut;
    use quickcheck::TestResult;

    use super::*;

    quickcheck! {
        fn encode_decode_32bit(input_data: Vec<i16>) -> TestResult {
            let d = input_data.iter().map(|&x| x as i32).collect::<Vec<i32>>();
            let mut input: Vec<i32> = d.clone();

            let y = 1;
            let x = input.len();
            let scale = 0; //lossless

            if x < 10 {
                return TestResult::discard();
            }

            let mut compressed: Vec<u8> = Vec::with_capacity(x * y * 100);
            let mut encoder = crate::write::HCEncoder::new(&mut compressed);

            let _ = encoder.write(&mut input, y, x, scale);

            let mut uncompressed: Vec<i32> = vec![0; x * y];
            let mut decoder = crate::read::HCDecoder::new();
            let _ = decoder.read(&compressed, 0, &mut uncompressed);

            if d == uncompressed {
                return TestResult::passed();
            } else {
                return TestResult::failed();
            }
        }


        fn encode_decode_64bit(input_data: Vec<i32>) -> TestResult {
            let d = input_data.iter().map(|&x| x as i64).collect::<Vec<i64>>();
            let mut input: Vec<i64> = d.clone();

            let y = 1;
            let x = input.len();
            let scale = 0; //lossless

            if x < 2 {
                return TestResult::discard();
            }

            let mut compressed: Vec<u8> = Vec::with_capacity(x * y * 100);
            let mut encoder = crate::write::HCEncoder::new(&mut compressed);

            let _ = encoder.write64(&mut input, y, x, scale);

            let mut uncompressed: Vec<i32> = vec![0; x * y * 2];
            let mut decoder = crate::read::HCDecoder::new();
            let _ = decoder.read64(&compressed, 0, cast_slice_mut(&mut uncompressed));

            if input_data == uncompressed[..(x*y)] {
                return TestResult::passed();
            } else {
                return TestResult::failed();
            }
        }
    }
}
