use std::io::Write;
use std::num::NonZero;
use std::path::Path;

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::call::simple_call;
use crate::lookup::{FLAG_DUP, FLAG_QCFAIL, FLAG_SECONDARY, FLAG_UNMAP};
use crate::output::write_fasta;
use crate::pileup::ActiveRead;

/// Default exclude-flags matching samtools consensus (UNMAP|SECONDARY|QCFAIL|DUP).
pub const DEFAULT_EXCL_FLAGS: u16 = FLAG_UNMAP | FLAG_SECONDARY | FLAG_QCFAIL | FLAG_DUP;

#[derive(Debug, Clone)]
pub struct ConsensusOpts {
    /// Weight each base by its quality score (`--use-qual`).
    pub use_qual: bool,
    /// Minimum base quality to count a base (`--min-BQ`).
    pub min_qual: u8,
    /// Minimum depth to emit a non-N call (`-d`/`--min-depth`).
    pub min_depth: u32,
    /// Minimum fraction of total score the called allele(s) must reach (`-c`/`--call-fract`).
    pub call_fract: f64,
    /// Minimum ratio of second-best to best score to call a het (`-H`/`--het-fract`).
    pub het_fract: f64,
    /// Enable IUPAC ambiguity codes (`-A`/`--ambig`).
    pub ambig: bool,
    /// Emit deletion (`*`) bases (`--show-del yes`).
    pub show_del: bool,
    /// FASTA line wrap length (`-l`/`--line-len`).
    pub line_len: usize,
    /// Reads with any of these FLAG bits are excluded.
    pub excl_flags: u16,
    /// Include only reads with any of these FLAG bits set (0 = no filter).
    pub incl_flags: u16,
    /// Minimum mapping quality (`--min-MQ`).
    pub min_mapq: u8,
}

impl Default for ConsensusOpts {
    fn default() -> Self {
        Self {
            use_qual: false,
            min_qual: 0,
            min_depth: 1,
            call_fract: 0.75,
            het_fract: 0.5,
            ambig: false,
            show_del: false,
            line_len: 70,
            excl_flags: DEFAULT_EXCL_FLAGS,
            incl_flags: 0,
            min_mapq: 0,
        }
    }
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct ConsensusStats {
    pub sequences: u64,
    pub positions: u64,
}

/// Compute simple-mode FASTA consensus for every covered contig.
///
/// Algorithm: `bam_consensus.c::calculate_consensus_simple` + `basic_fasta`,
/// tag 1.23.1. Single-threaded decode + CIGAR walk; `workers` is accepted but
/// ignored — bgzf MultithreadedReader interacts poorly with header-then-raw-read
/// sequences across platforms.
pub fn consensus(
    bam_path: &Path,
    out: &mut dyn Write,
    opts: &ConsensusOpts,
    workers: NonZero<usize>,
) -> Result<ConsensusStats> {
    let mut stats = ConsensusStats::default();

    let _ = workers;
    let st = NonZero::<usize>::new(1).unwrap();
    let mut reader = rsomics_bamio::open_with_workers(bam_path, st)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();

    let reader_mut = reader.get_mut();
    // Two pre-allocated record buffers swapped to avoid per-record heap alloc.
    let mut read_buf = RawRecord::default();
    let mut next_rec: Option<RawRecord> = None;

    let mut active: Vec<ActiveRead> = Vec::with_capacity(512);
    // Minimum end among active reads — avoids O(n) retain when cursor hasn't
    // reached the nearest read end.
    let mut active_min_end: i64 = i64::MAX;
    let mut col_buf: Vec<(u8, u8)> = Vec::with_capacity(512);
    // Free list of recycled slot Vecs — avoids per-read malloc for CIGAR slots.
    // Pre-filled so the initial burst of reads doesn't hit the global allocator.
    let mut slot_pool: Vec<Vec<(u8, u8)>> = (0..64).map(|_| Vec::with_capacity(300)).collect();

    let mut cur_tid: i32 = -1;
    let mut cur_seq: Vec<u8> = Vec::new();
    let mut last_pos: i64 = -1;

    let mut cursor_tid: i32 = -1;
    let mut cursor_pos: i64 = 0;

    macro_rules! read_next {
        () => {{
            let n = rsomics_bamio::raw::read_record(reader_mut, &mut read_buf)?;
            if n > 0 {
                // Swap: read_buf reuses next_rec's allocation; next_rec gets the fresh record.
                let mut tmp = next_rec.take().unwrap_or_default();
                std::mem::swap(&mut read_buf, &mut tmp);
                next_rec = Some(tmp);
                true
            } else {
                next_rec = None;
                false
            }
        }};
    }

    read_next!();

    loop {
        #[allow(clippy::while_let_loop)]
        loop {
            let Some(ref nr) = next_rec else { break };
            let flag = nr.flags();
            if nr.reference_sequence_id() < 0 || flag & FLAG_UNMAP != 0 {
                read_next!();
                continue;
            }
            if opts.excl_flags != 0 && flag & opts.excl_flags != 0 {
                read_next!();
                continue;
            }
            if opts.incl_flags != 0 && flag & opts.incl_flags == 0 {
                read_next!();
                continue;
            }
            if nr.mapping_quality() < opts.min_mapq {
                read_next!();
                continue;
            }

            let rtid = nr.reference_sequence_id();
            let rbeg = i64::from(nr.alignment_start());

            if cursor_tid < 0 {
                cursor_tid = rtid;
                cursor_pos = rbeg;
            }

            if rtid > cursor_tid || (rtid == cursor_tid && rbeg > cursor_pos) {
                break;
            }

            let cur_nr = next_rec.take().unwrap();
            let recycled = slot_pool.pop().unwrap_or_default();
            let ar = ActiveRead::new(cur_nr, recycled);
            if ar.end < active_min_end {
                active_min_end = ar.end;
            }
            active.push(ar);

            read_next!();
        }

        if cursor_pos >= active_min_end {
            let mut i = 0;
            while i < active.len() {
                if active[i].tid != cursor_tid || active[i].end <= cursor_pos {
                    let retired = active.swap_remove(i);
                    retired.retire(&mut slot_pool);
                } else {
                    i += 1;
                }
            }
            active_min_end = active.iter().map(|r| r.end).fold(i64::MAX, i64::min);
        }

        if active.is_empty() {
            active_min_end = i64::MAX;
            let Some(ref nr) = next_rec else {
                break;
            };
            let rtid = nr.reference_sequence_id();
            let rbeg = i64::from(nr.alignment_start());
            if rtid != cursor_tid && cur_tid >= 0 && !cur_seq.is_empty() {
                let name = ref_names.get(cur_tid as usize).ok_or_else(|| {
                    RsomicsError::InvalidInput(format!("tid {cur_tid} not in header"))
                })?;
                write_fasta(out, name, &cur_seq, opts.line_len).map_err(RsomicsError::Io)?;
                stats.sequences += 1;
                cur_seq.clear();
                cur_tid = -1;
                last_pos = -1;
            }
            cursor_tid = rtid;
            cursor_pos = rbeg;
            continue;
        }

        let emit_tid = cursor_tid;
        let emit_pos = cursor_pos;

        if emit_tid != cur_tid {
            if cur_tid >= 0 && !cur_seq.is_empty() {
                let name = ref_names.get(cur_tid as usize).ok_or_else(|| {
                    RsomicsError::InvalidInput(format!("tid {cur_tid} not in header"))
                })?;
                write_fasta(out, name, &cur_seq, opts.line_len).map_err(RsomicsError::Io)?;
                stats.sequences += 1;
                cur_seq.clear();
            }
            cur_tid = emit_tid;
            last_pos = -1;
        }

        col_buf.clear();
        for ar in &mut active {
            let (b4, q, is_del, is_refskip) = ar.at(emit_pos);
            if is_refskip {
                continue;
            }
            col_buf.push((if is_del { 16 } else { b4 }, q));
        }

        let (cb, _cq) = simple_call(&col_buf, opts);

        if cb == b'*' && !opts.show_del {
            last_pos = emit_pos;
            cursor_pos += 1;
            continue;
        }

        // Fill positions between last covered and this one with N.
        if !cur_seq.is_empty() && emit_pos > last_pos + 1 {
            let gap = (emit_pos - last_pos - 1) as usize;
            cur_seq.extend(std::iter::repeat_n(b'N', gap));
            stats.positions += gap as u64;
        }

        cur_seq.push(cb);
        stats.positions += 1;
        last_pos = emit_pos;
        cursor_pos += 1;
    }

    if cur_tid >= 0 && !cur_seq.is_empty() {
        let name = ref_names
            .get(cur_tid as usize)
            .ok_or_else(|| RsomicsError::InvalidInput(format!("tid {cur_tid} not in header")))?;
        write_fasta(out, name, &cur_seq, opts.line_len).map_err(RsomicsError::Io)?;
        stats.sequences += 1;
    }

    Ok(stats)
}
