use rsomics_bamio::raw::RawRecord;

use crate::lookup::{
    CIGAR_DEL, CIGAR_DIFF, CIGAR_EQUAL, CIGAR_INS, CIGAR_MATCH, CIGAR_REF_SKIP, CIGAR_SOFT_CLIP,
};

/// Sentinel slot values for deletions and ref-skips.
pub(crate) const SLOT_DEL: (u8, u8) = (16, 0);
pub(crate) const SLOT_SKIP: (u8, u8) = (17, 0);

/// One active read in the pileup buffer.
///
/// CIGAR is flattened into `slots` indexed by reference offset from `beg`:
/// each slot is `(base4, qual)`. O(1) column lookup with no per-column allocation.
///
/// `base4` values: 0–15 = seq_nt16; 16 = deletion; 17 = ref-skip.
pub(crate) struct ActiveRead {
    /// Reference start (inclusive, 0-based).
    pub(crate) beg: i64,
    /// Reference end (exclusive, 0-based).
    pub(crate) end: i64,
    /// tid of the read.
    pub(crate) tid: i32,
    /// Per-reference-position `(base4, qual)`, indexed by `(ref_pos - beg)`.
    pub(crate) slots: Vec<(u8, u8)>,
}

impl ActiveRead {
    /// Build from a raw record; `slots` is a recycled Vec from the pool.
    pub(crate) fn new(rec: RawRecord, mut slots: Vec<(u8, u8)>) -> Self {
        let beg = i64::from(rec.alignment_start());
        let tid = rec.reference_sequence_id();
        let qscores: &[u8] = rec.quality_scores();

        let span = {
            let mut s: i64 = 0;
            for (op, l) in rec.cigar_ops() {
                match op {
                    CIGAR_MATCH | CIGAR_DEL | CIGAR_REF_SKIP | CIGAR_EQUAL | CIGAR_DIFF => {
                        s += i64::from(l);
                    }
                    _ => {}
                }
            }
            s as usize
        };
        slots.clear();
        slots.reserve(span.saturating_sub(slots.capacity()));

        let mut qpos: usize = 0;
        for (op, l) in rec.cigar_ops() {
            let l = l as usize;
            match op {
                CIGAR_MATCH | CIGAR_EQUAL | CIGAR_DIFF => {
                    for i in 0..l {
                        let qi = qpos + i;
                        let b4 = rec.seq_nibble(qi);
                        let q = if qi < qscores.len() { qscores[qi] } else { 0 };
                        slots.push((b4, q));
                    }
                    qpos += l;
                }
                CIGAR_DEL => {
                    slots.extend(std::iter::repeat_n(SLOT_DEL, l));
                }
                CIGAR_REF_SKIP => {
                    slots.extend(std::iter::repeat_n(SLOT_SKIP, l));
                }
                CIGAR_INS | CIGAR_SOFT_CLIP => {
                    qpos += l;
                }
                _ => {}
            }
        }

        let end = beg + slots.len() as i64;
        ActiveRead {
            beg,
            end,
            tid,
            slots,
        }
    }

    /// Return the slot Vec to the pool.
    pub(crate) fn retire(self, pool: &mut Vec<Vec<(u8, u8)>>) {
        if pool.len() < 128 {
            pool.push(self.slots);
        }
    }

    /// O(1) lookup: `(base4, qual, is_del, is_refskip)` at reference `pos`.
    #[inline]
    pub(crate) fn at(&self, pos: i64) -> (u8, u8, bool, bool) {
        let off = (pos - self.beg) as usize;
        if off >= self.slots.len() {
            return (0, 0, true, false);
        }
        let (b4, q) = self.slots[off];
        match b4 {
            16 => (16, 0, true, false), // deletion
            17 => (16, 0, true, true),  // ref-skip
            b => (b, q, false, false),
        }
    }
}
