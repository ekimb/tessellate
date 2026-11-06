#[derive(Clone, Copy, Debug)]
pub struct Anchor {
    pub q: u32,
    pub r: u32,
    pub len: u32,
}

pub struct Chain {
    pub rec: usize,
    pub minus: bool,
    pub anchors: Vec<Anchor>,
}
