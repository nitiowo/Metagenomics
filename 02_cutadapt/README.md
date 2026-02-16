# Primer demultiplexing and trimming

Cutadapt was used in 3 steps:
1. `02a_reverse_demux/` — separate by reverse primer (HCO2198 vs SSUR)
2. `02b_forward_demux/` — split HCO2198 reads by forward primer (LCO1490 vs mICOintF)
3. `02c_rescue_demux/` — check unmatched reads for anything missed in step 1

- `demux_submit_array.job` — SGE job array that runs all three steps (one job per sample per primer)
- `post_demux_qc.job` — SGE job to perform QC on final demux output