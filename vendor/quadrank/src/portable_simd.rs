#[derive(Clone, Copy)]
#[repr(C, align(32))]
pub struct __m256i(pub [u8; 32]);

#[derive(Clone, Copy)]
#[repr(C, align(16))]
pub struct __m128i(pub [u32; 4]);

#[inline(always)]
pub unsafe fn _mm256_shuffle_epi8(a: __m256i, b: __m256i) -> __m256i {
    let mut r = [0u8; 32];
    for lane in 0..2 {
        let base = 16 * lane;
        for j in 0..16 {
            let idx = b.0[base + j];
            r[base + j] = if idx & 0x80 != 0 {
                0
            } else {
                a.0[base + (idx & 0x0f) as usize]
            };
        }
    }
    __m256i(r)
}

#[inline(always)]
pub unsafe fn _mm256_sad_epu8(a: __m256i, b: __m256i) -> __m256i {
    let mut r = [0u8; 32];
    for k in 0..4 {
        let mut sum: u64 = 0;
        for j in 0..8 {
            let i = 8 * k + j;
            sum += a.0[i].abs_diff(b.0[i]) as u64;
        }
        r[8 * k..8 * k + 8].copy_from_slice(&sum.to_le_bytes());
    }
    __m256i(r)
}

#[inline(always)]
fn lanes64(x: __m256i) -> [u64; 4] {
    core::array::from_fn(|k| u64::from_le_bytes(x.0[8 * k..8 * k + 8].try_into().unwrap()))
}

#[inline(always)]
fn from_lanes64(l: [u64; 4]) -> __m256i {
    let mut r = [0u8; 32];
    for k in 0..4 {
        r[8 * k..8 * k + 8].copy_from_slice(&l[k].to_le_bytes());
    }
    __m256i(r)
}

#[inline(always)]
pub unsafe fn _mm256_unpacklo_epi64(a: __m256i, b: __m256i) -> __m256i {
    let (a, b) = (lanes64(a), lanes64(b));
    from_lanes64([a[0], b[0], a[2], b[2]])
}

#[inline(always)]
pub unsafe fn _mm256_unpackhi_epi64(a: __m256i, b: __m256i) -> __m256i {
    let (a, b) = (lanes64(a), lanes64(b));
    from_lanes64([a[1], b[1], a[3], b[3]])
}

#[inline(always)]
pub unsafe fn _mm_sign_epi32(a: __m128i, b: __m128i) -> __m128i {
    __m128i(core::array::from_fn(|i| {
        let s = b.0[i] as i32;
        if s < 0 {
            a.0[i].wrapping_neg()
        } else if s == 0 {
            0
        } else {
            a.0[i]
        }
    }))
}
