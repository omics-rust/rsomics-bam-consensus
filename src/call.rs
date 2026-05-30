use crate::driver::ConsensusOpts;
use crate::lookup::{N_BASES, SEQI2A, SEQI2C, SEQI2G, SEQI2T};

/// `bam_consensus.c::calculate_consensus_simple` — returns ASCII base call and
/// quality (0–100 as a percentage of total score).
///
/// `reads`: `(base4, qual)` pairs; `base4` is 4-bit seq_nt16, gap/del encoded as 16.
pub fn simple_call(reads: &[(u8, u8)], opts: &ConsensusOpts) -> (u8, u8) {
    // IUPAC het table: `const char *het = "NACMGRSVTWYHKDBN*ac?g???t???????";`
    // from bam_consensus.c `calculate_consensus_simple`.
    const HET: &[u8; 32] = b"NACMGRSVTWYHKDBN*ac?g???t???????";

    let mut score = [0u64; N_BASES]; // A C G T *
    let mut tot_depth: u32 = 0;

    for &(b4, q) in reads {
        if q < opts.min_qual {
            continue;
        }
        let w: u64 = if opts.use_qual { q as u64 } else { 1 };
        if b4 < 16 {
            let b = b4 as usize;
            let qa = SEQI2A[b] as u64 * w;
            let qc = SEQI2C[b] as u64 * w;
            let qg = SEQI2G[b] as u64 * w;
            let qt = SEQI2T[b] as u64 * w;
            if qa > 0 {
                score[0] += qa;
            }
            if qc > 0 {
                score[1] += qc;
            }
            if qg > 0 {
                score[2] += qg;
            }
            if qt > 0 {
                score[3] += qt;
            }
        } else {
            score[4] += 8 * w;
        }
        tot_depth += 1;
    }

    let tscore: u64 = score.iter().sum();

    // Find best (call1) and second-best (call2).
    // Slot i → seqi bit = 1 << i  (0→A=1, 1→C=2, 2→G=4, 3→T=8, 4→*=16).
    let mut call1: usize = 15; // seqi 15 = N
    let mut call2: usize = 15;
    let mut score1: u64 = 0;
    let mut score2: u64 = 0;

    for (i, &s) in score.iter().enumerate() {
        if s > score1 {
            score2 = score1;
            call2 = call1;
            score1 = s;
            call1 = 1 << i;
        } else if s > score2 {
            score2 = s;
            call2 = 1 << i;
        }
    }

    let mut used_base = call1;
    let mut used_score = score1;

    if opts.ambig && score1 > 0 && (score2 as f64) >= opts.het_fract * score1 as f64 {
        used_base |= call2;
        used_score += score2;
    }

    // Depth / fraction gate.
    if tot_depth < opts.min_depth
        || (tscore > 0 && (used_score as f64) < opts.call_fract * tscore as f64)
    {
        used_base = if call1 == 16 { 16 } else { 0 };
    }

    let cb = HET[used_base.min(31)];
    let cq: u8 = if used_base != 0 && tscore > 0 {
        (100 * used_score / tscore).min(100) as u8
    } else {
        0
    };

    (cb, cq)
}
