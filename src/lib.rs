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
    pub d: Vec<i32>
}