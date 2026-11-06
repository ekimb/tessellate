mod block;
pub mod budgets;
mod codec;
mod edit;
mod indexed;
pub mod fastx;
pub mod packed;
mod reference;
pub mod runs;
mod sequence;
pub mod store;
mod tile;
pub mod tq;
mod view;

pub use packed::{Packed, PackedView};
pub(crate) use reference::ReferenceCache;
pub use sequence::parse;
pub use tile::{AMB, Tile};
pub use view::View;
pub const DEFAULT_RATE: usize = 32;

pub const CHUNK: usize = 1 << 18;
const _: () = assert!(CHUNK % 256 == 0);
pub const EMPTY: u8 = 3;
