//! FASTA consensus from a coordinate-sorted BAM — simple-mode port of
//! `samtools consensus -m simple`.
//!
//! Per reference position the engine accumulates a weighted score for each of
//! A, C, G, T, and * (deletion); weight per read is `qual` (or 1 when
//! `use_qual = false`) times a compatibility factor from the seqi2A/C/G/T
//! tables. The highest-scoring allele is the call; if the second-highest
//! scores ≥ `het_fract × score1` and `ambig` is on, the two are OR-ed into an
//! IUPAC code. If total depth < `min_depth` or the call's fraction of total
//! score < `call_fract`, the call is N (or * for a gap).
//!
//! Reference source: `samtools/bam_consensus.c` (MIT), tag 1.23.1.
//!
//! Pileup engine: custom lightweight walker over coordinate-sorted BAM records,
//! modelled on `samtools/bam_plbuf.c` / `consensus_pileup.c`. Avoids per-column
//! `Vec` allocations and HashMap overlap tracking; single-threaded decode + CIGAR walk.

mod call;
mod driver;
mod lookup;
mod output;
mod pileup;

pub use call::simple_call;
pub use driver::{ConsensusOpts, ConsensusStats, DEFAULT_EXCL_FLAGS, consensus};
