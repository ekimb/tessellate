use super::packed::DecodedBlock;
use std::cell::RefCell;
use std::collections::{HashMap, VecDeque};
use std::rc::Rc;
use std::sync::atomic::{AtomicU64, Ordering};

type Key = (u64, usize, usize);
#[derive(Default)]
struct Boundaries {
    blocks: HashMap<Key, Rc<Vec<usize>>>,
    order: VecDeque<Key>,
}
thread_local! { static CACHE: RefCell<Boundaries> = RefCell::new(Boundaries::default()); }

const MAX_BLOCKS: usize = 1024;
const MAX_BLOCK_BYTES: usize = 8 << 20;

struct DecodedCache {
    blocks: HashMap<Key, Rc<DecodedBlock>>,
    order: VecDeque<Key>,
    bytes: usize,
    max_blocks: usize,
    max_bytes: usize,
}
impl DecodedCache {
    fn new(max_blocks: usize, max_bytes: usize) -> Self {
        Self {
            blocks: HashMap::new(),
            order: VecDeque::new(),
            bytes: 0,
            max_blocks,
            max_bytes,
        }
    }
    fn insert(&mut self, key: Key, block: Rc<DecodedBlock>) {
        let size = block.size_in_bytes();
        if self.max_blocks == 0 || size > self.max_bytes || self.blocks.contains_key(&key) {
            return;
        }
        while self.blocks.len() >= self.max_blocks || self.bytes > self.max_bytes - size {
            let old = self
                .order
                .pop_front()
                .expect("nonempty decoded block cache");
            self.bytes -= self.blocks.remove(&old).unwrap().size_in_bytes();
        }
        self.bytes += size;
        self.order.push_back(key);
        self.blocks.insert(key, block);
    }
}
thread_local! {
    static DECODED: RefCell<DecodedCache> =
        RefCell::new(DecodedCache::new(MAX_BLOCKS, MAX_BLOCK_BYTES));
}
#[derive(Clone, Copy)]
pub(crate) struct ReferenceCache(u64);
impl Default for ReferenceCache {
    fn default() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        Self(NEXT.fetch_add(1, Ordering::Relaxed))
    }
}

impl ReferenceCache {
    pub(super) fn decoded(&self, rec: usize, block: usize) -> Option<Rc<DecodedBlock>> {
        DECODED.with(|c| c.borrow().blocks.get(&(self.0, rec, block)).cloned())
    }

    pub(super) fn insert_decoded(&self, rec: usize, block: usize, decoded: Rc<DecodedBlock>) {
        DECODED.with(|c| c.borrow_mut().insert((self.0, rec, block), decoded));
    }

    pub(super) fn boundaries(&self, rec: usize, block: usize) -> Option<Rc<Vec<usize>>> {
        CACHE.with(|c| c.borrow().blocks.get(&(self.0, rec, block)).cloned())
    }

    pub(super) fn insert_boundaries(&self, rec: usize, block: usize, bounds: Vec<usize>) {
        let key = (self.0, rec, block);
        CACHE.with(|c| {
            let mut c = c.borrow_mut();
            if c.blocks.len() >= MAX_BLOCKS {
                if let Some(old) = c.order.pop_front() {
                    c.blocks.remove(&old);
                }
            }
            c.order.push_back(key);
            c.blocks.insert(key, Rc::new(bounds));
        });
    }
}
