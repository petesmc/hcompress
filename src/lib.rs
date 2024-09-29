pub mod read;
pub mod write;

#[derive(Clone, Debug)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct Data {
    pub d: Vec<i16>,
    //pub bs: u8,
}
