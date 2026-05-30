// seqi2{A,C,G,T}: index is 4-bit seq_nt16 (A=1 C=2 G=4 T=8).
// Pure bases carry weight 8; 2-base ambiguity 4; 3-base 2; N 1.
// Source: bam_consensus.c static arrays, lines ~1908–1911 in 1.23.1.
//                              * A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
pub(crate) const SEQI2A: [u8; 16] = [0, 8, 0, 4, 0, 4, 0, 2, 0, 4, 0, 2, 0, 2, 0, 1];
pub(crate) const SEQI2C: [u8; 16] = [0, 0, 8, 4, 0, 0, 4, 2, 0, 0, 4, 2, 0, 0, 2, 1];
pub(crate) const SEQI2G: [u8; 16] = [0, 0, 0, 0, 8, 4, 4, 1, 0, 0, 0, 0, 4, 2, 2, 1];
pub(crate) const SEQI2T: [u8; 16] = [0, 0, 0, 0, 0, 0, 0, 0, 8, 4, 4, 2, 8, 2, 2, 1];

/// Score/depth slots: 0=A 1=C 2=G 3=T 4=* (gap).
pub(crate) const N_BASES: usize = 5;

// SAM FLAG bits (SAMv1 §1.4).
pub(crate) const FLAG_UNMAP: u16 = 0x4;
pub(crate) const FLAG_SECONDARY: u16 = 0x100;
pub(crate) const FLAG_QCFAIL: u16 = 0x200;
pub(crate) const FLAG_DUP: u16 = 0x400;

// CIGAR op codes (BAM packed encoding, low nibble).
pub(crate) const CIGAR_MATCH: u8 = 0;
pub(crate) const CIGAR_INS: u8 = 1;
pub(crate) const CIGAR_DEL: u8 = 2;
pub(crate) const CIGAR_REF_SKIP: u8 = 3;
pub(crate) const CIGAR_SOFT_CLIP: u8 = 4;
pub(crate) const CIGAR_EQUAL: u8 = 7;
pub(crate) const CIGAR_DIFF: u8 = 8;
