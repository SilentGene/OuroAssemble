# OuroAssemble — detailed usage

## 1. Input

One PacBio HiFi metagenome per run: FASTQ or FASTA, optionally gzipped. Short reads
and other long-read types are not supported.

## 2. Running

```bash
ouroassemble run -r reads.fq.gz -o results -t 32 [-c config.yaml] [options] [-- snakemake args]
```

`ouroassemble run` writes the effective user config to `results/ouroassemble.config.yaml`
and calls Snakemake with `--software-deployment-method conda`. The workflow can also be
run directly:

```bash
snakemake -s $(python -c "import ouroassemble;print(ouroassemble.SNAKEFILE)") \
    -d results --configfile my.yaml --sdm conda --cores 32
```

### Memory

* `-m/--memory GB` (config `max_memory_gb`): upper limit for every job. Each step uses
  `min(<step>.mem_mb, limit)`; locally the summed memory of concurrent jobs also stays
  below the limit (`--resources mem_mb=...`). With a cluster profile it caps each job's
  request (e.g. `mqsub -m`), which is a hard limit there: a job exceeding it is killed.
* Per step: `metamdbg`, `hifiasm_meta`, `metaflye`, `hicanu`, `mapping`, `checkm2`,
  `gtdbtk`, `genomad` each take `mem_mb` in the config. HiCanu's `maxMemory` follows it.
* hifiasm-meta needs > 16 GB even for small inputs (default Bloom filter `-f 37`).
* CPUs: each step's `threads` is a fixed number, except HiCanu's default `threads: all`,
  which takes every core of a local run (`-t`). HiCanu then starts once the other
  assemblers have freed their cores and runs alone, which suits its long, well
  parallelised overlap stage. With `--profile` it requests `hicanu.cluster_threads`
  (default 32) instead. Any step can use `threads: all`, or set HiCanu back to a number.

### On a cluster

Pass a Snakemake profile with `--profile`. Each rule declares `threads`,
`resources.mem_mb` and `resources.runtime` (minutes) taken from the config, so they can
be tuned per step:

```yaml
metaflye: {threads: 48, mem_mb: 200000, runtime: 2880}
gtdbtk:   {threads: 32, mem_mb: 160000}
```

### Resuming an interrupted run

If a run stops part-way (walltime reached, node failure, Ctrl-C), run the same command
again with `--resume`, optionally with more threads or memory:

```bash
ouroassemble run -r reads.fq.gz -o out -c my.yaml -t 36 -m 200 --resume
```

* Snakemake's stale lock is removed, and only missing files are made: finished steps are
  kept even if the config or OuroAssemble changed in between (`--rerun-triggers mtime`).
  Make sure no other run is using the same output directory.
* The assemblers continue from their own intermediate files in `02_assembly/<assembler>/`:

  | Assembler | Resumes from |
  |---|---|
  | HiCanu | its last finished task (e.g. individual overlap jobs) |
  | metaFlye | its last finished stage (`--resume`) |
  | hifiasm-meta | its saved read correction and overlaps (`*.bin`); if these are unusable, it starts over |
  | metaMDBG | cannot resume (starts over; it is the fastest assembler) |

  A signature of the reads and assembler options is kept in
  `02_assembly/<assembler>/.ouroassemble_resume`. If the reads or options change, that
  directory is emptied and the assembler starts from scratch instead, so stale
  intermediate files are never reused. Threads and memory are not part of the signature.
* To chain a resume job after a cluster job that may hit its walltime, make it depend on
  the first job ending in any state (PBS: `-W depend=afterany:<jobid>`).

## 3. Steps and files

| Step | Directory | Files |
|---|---|---|
| read filtering | `01_reads/` | `reads.filtered.fq.gz`, `read_stats.tsv` |
| assembly | `02_assembly/<assembler>/` | native assembler output |
| standardisation | `03_contigs/<assembler>/` | `contigs.fa` (IDs `<assembler>__<id>`), `contigs.tsv`, `circular_candidates.fa` |
| trim + rotate | `04_candidates/<assembler>/` | `trimmed.fa`, `rotated.fa`, `validation_ref.fa`, `trim_rotate.tsv`, `junctions.bed` |
| junction validation | `05_validation/<assembler>/` | `reads.bam`, `junction_eval.tsv`, `suspicious_sites.bed`, `igv_batch.txt`, `igv_sites.bed`, (`igv_report.html`) |
| classification | `06_classify/` | `candidates.tsv`, `genomes/`, `checkm2/`, `gtdbtk/`, `genomad/`, `classification.tsv` |
| dereplication | `07_derep/` | `skani_edges.tsv`, `clusters.tsv` |
| merged assembly | `08_merge/` | `pool.tsv`, `pool.fa`, `pool.paf`, `selection.tsv`, `merged_assembly.fa`, `reads.bam`, `mapping_stats.tsv` (+ `05_validation/<assembler>/contig_scan.tsv`) |
| report | `final/` | `report.html`, `circular_genomes.tsv`, FASTAs |

### 3.1 Circular markers

| Assembler | File | Marker |
|---|---|---|
| metaMDBG | `contigs.fasta.gz` | header `circular=yes` |
| hifiasm-meta | `asm.p_ctg.gfa` | `S`-line ID ends with `c`; `asm.rescue.fa` with `hifiasm_meta.include_rescue: true` (rows get `source=rescue`) |
| metaFlye | `assembly_info.txt` | column 4 `circ.` = `Y` (column 5 `repeat` is kept) |
| HiCanu | `asm.contigs.fasta` | header `suggestCircular=yes` |

### 3.2 Overlap trimming and rotation

Some assemblers (notably HiCanu) keep a copy of the contig start at its end. The longest
exact suffix = prefix match is searched first (`trim_rotate.overlap_search` bp), then an
approximate one with minimap2 (`asm5`, identity ≥ `min_identity`). An overlap ≥
`min_overlap` bp is removed once. The trimmed contig `A+B` is rotated to `B+A`; the
former end→start seam is at `junction_pos` (0-based) in `rotated.fa`.

### 3.3 Junction validation

Reads are mapped with `minimap2 -ax map-hifi -Y` to `validation_ref.fa` = every contig of
the same assembly, with candidates replaced by their rotated versions, so reads that
belong elsewhere receive low MAPQ instead of being forced onto the candidate.

`junction_eval.tsv` columns:

| Column | Meaning |
|---|---|
| `spanning_reads` | reads with MAPQ ≥ `min_mapq` covering `J ± flank` in one alignment, no indel ≥ `min_clip` near J |
| `depth_at_junction` / `unique_depth_at_junction` | all / MAPQ ≥ `min_mapq` alignments covering J |
| `span_ratio` | spanning / unique depth |
| `clipped_at_junction` | reads clipped (≥ `min_clip`) within ± 2×`clip_window` of J |
| `low_mapq_frac` | 1 − unique depth / depth at J (high → repeat) |
| `n_clip_sites` | other positions where ≥ `min_clip_reads` and ≥ `clip_frac` × depth reads are clipped |
| `n_low_cov_regions` | windows (`cov_window`) below `low_cov_frac` × median depth |
| `status`, `reason` | see below |

| Status | Rule |
|---|---|
| PASS | spanning ≥ `min_spanning`, span ratio ≥ `min_span_ratio`, no suspicious sites |
| PENDING | depth at J < `min_depth` (`low_coverage_at_junction`); low-MAPQ fraction ≥ `repeat_low_mapq_frac` (`repeat_at_junction`); or junction fine but suspicious sites elsewhere (`suspicious_sites_elsewhere`) |
| FAIL | otherwise (`junction_not_spanned`) |

Manual inspection: `igv -b 05_validation/<assembler>/igv_batch.txt` loads the reference,
alignments and sites and saves a snapshot of every junction and suspicious site. Set
`junction.igv_reports: true` for an HTML igv-reports page per assembler.

### 3.4 Classification

Candidates with status PASS (and PENDING if `classify.classify_pending`) are *selected*.

* ≥ `min_chrom_len`: each written to its own FASTA in `06_classify/genomes/` and run
  through CheckM2 and GTDB-Tk. `chromosome` = completeness ≥ `min_completeness` and
  contamination ≤ `max_contamination`; otherwise `large_circular`. The genome size of
  the GTDB-Tk closest reference is read from the GTDB metadata files and reported as
  `ref_genome_length` / `length_ratio`.
* < `min_chrom_len` (or all with `genomad_on_all`): geNomad → `plasmid`, `virus` or
  `other_circular`.

### 3.5 Dereplication

Selected sequences ≥ `derep.small_len` are compared with `skani triangle`; pairs
involving shorter sequences with minimap2 against the doubled (circular) target. Links
need ANI ≥ `min_ani` and aligned fraction ≥ `min_af` in both directions; clusters are
single-linkage. The representative is chosen by: junction status → CheckM2
(completeness − 5 × contamination) → spanning reads → fewer suspicious sites → length
closest to the cluster median → `assembler_priority`.

### 3.6 Merged assembly

1. Every contig of every assembly gets a read-support scan (`contig_scan.tsv`: clip
   pile-ups and low-coverage windows, same parameters as the junction check).
2. Pool (`pool.tsv`), ranked by tier → suspicious sites → length → `assembler_priority`:
   * tier 0: dereplication representatives with junction PASS (overlap-trimmed sequence);
   * tier 1: representatives with PENDING status;
   * tier 2: linear contigs ≥ `merge.min_contig_len`, plus circular candidates whose
     junction FAILed (`category=circular_unverified`, if `merge.include_failed_circular`).
   Non-representative members of a circle cluster are not pooled.
3. `minimap2 -x asm5 -c -X` all-vs-all of the pool. Going down the ranking, a contig is
   dropped when ≥ `merge.max_covered_frac` of its length is covered by already kept
   contigs with alignment blocks ≥ `merge.min_block` bp at ≥ `merge.min_identity`.
   `selection.tsv` records `covered_frac` and `redundant_with` for every pooled contig.
4. Reads are mapped to the merged assembly; `% reads mapped` is reported for it and
   for every single assembly (from the validation alignments).

Contigs are selected, never cut, joined or polished: a region present only as part of a
larger contig that was dropped is not rescued, and partially overlapping contigs
(< `max_covered_frac`) are both kept. Disable the stage with `merge: {enabled: false}`.

## 4. Configuration reference

See the annotated defaults: `ouroassemble init-config -o config.yaml`.

## 5. Limitations

* Rotation does not move the sequence start to *dnaA* / a replication origin.
* Collinearity against the closest reference is not computed (only genome length).
* PENDING candidates are flagged, not resolved; inspect them in IGV.
* A circle that an assembler collapsed from a tandem repeat (e.g. two near-identical
  copies → one) can still PASS, because reads shorter than the true molecule wrap around
  the collapsed version consistently. In the example data metaMDBG collapses the 7.3 kb
  *Buchnera* pTrp plasmid (two trpEG copies) to 3.6 kb; hifiasm-meta keeps both copies.
  Compare lengths across assemblers in `final/circular_genomes.tsv` (same cluster or
  not) when copy number matters. OuroAssemble never trims overlaps longer than
  `trim_rotate.max_overlap_frac` × length itself, so it does not cause such collapses.
