use std::io::Write;

/// Write a FASTA record to `out`, wrapping at `line_len`.
pub(crate) fn write_fasta(
    out: &mut dyn Write,
    name: &str,
    seq: &[u8],
    line_len: usize,
) -> std::io::Result<()> {
    if seq.is_empty() {
        return Ok(());
    }
    writeln!(out, ">{name}")?;
    for chunk in seq.chunks(line_len) {
        out.write_all(chunk)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}
