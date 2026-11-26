# 06 — classify validated circular candidates
#      large (>= min_chrom_len): one FASTA per contig -> CheckM2 + GTDB-Tk
#      small: geNomad (plasmid / virus)


rule prepare_classification:
    input:
        contigs=expand("03_contigs/{asm}/contigs.tsv", asm=ASSEMBLERS),
        trim=expand("04_candidates/{asm}/trim_rotate.tsv", asm=ASSEMBLERS),
        junction=expand("05_validation/{asm}/junction_eval.tsv", asm=ASSEMBLERS),
        trimmed=expand("04_candidates/{asm}/trimmed.fa", asm=ASSEMBLERS),
    output:
        tsv="06_classify/candidates.tsv",
        genomes=directory("06_classify/genomes"),
        genomad="06_classify/genomad_input.fa",
        selected="06_classify/selected.fa",
        derep_large="06_classify/derep_large.fa",
    params:
        script=SCRIPTS / "prepare_classification.py",
        small_len=config["derep"]["small_len"],
        min_chrom_len=CLS["min_chrom_len"],
        pending="--classify-pending" if CLS["classify_pending"] else "",
        genomad_all="--genomad-on-all" if CLS["genomad_on_all"] else "",
    resources:
        mem_mb=cap_mem(8000),
        runtime=60,
    log:
        "logs/prepare_classification.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} --contigs {input.contigs} --trim {input.trim} "
        "--junction {input.junction} --trimmed {input.trimmed} "
        "--min-chrom-len {params.min_chrom_len} {params.pending} {params.genomad_all} "
        "--out-tsv {output.tsv} --out-genomes {output.genomes} "
        "--out-genomad {output.genomad} --out-selected {output.selected} "
        "--small-len {params.small_len} --out-derep-large {output.derep_large} > {log} 2>&1"


rule checkm2:
    input:
        "06_classify/genomes",
    output:
        "06_classify/checkm2/quality_report.tsv",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        db=DB["checkm2"],
    threads: job("checkm2")["threads"]
    resources:
        mem_mb=cap_mem(job("checkm2")["mem_mb"]),
        runtime=job("checkm2")["runtime"],
    log:
        "logs/checkm2.log",
    benchmark:
        "benchmarks/checkm2.tsv"
    conda:
        "../envs/checkm2.yaml"
    shell:
        "if ls {input}/*.fa > /dev/null 2>&1; then "
        "  checkm2 predict --input {input} --extension fa --output-directory {params.outdir} "
        "  --threads {threads} --database_path {params.db} --force > {log} 2>&1; "
        "else "
        "  mkdir -p {params.outdir}; "
        "  printf 'Name\\tCompleteness\\tContamination\\n' > {output}; "
        "  echo 'no genomes to assess' > {log}; "
        "fi"


rule gtdbtk:
    input:
        "06_classify/genomes",
    output:
        "06_classify/gtdbtk/summary.tsv",
    params:
        outdir=lambda w, output: os.path.join(os.path.dirname(output[0]), "out"),
        db=DB["gtdbtk"],
        ani="--skip_ani_screen" if config["gtdbtk"]["skip_ani_screen"] else "",
        pplacer=config["gtdbtk"]["pplacer_cpus"],
    threads: job("gtdbtk")["threads"]
    resources:
        mem_mb=cap_mem(job("gtdbtk")["mem_mb"]),
        runtime=job("gtdbtk")["runtime"],
    log:
        "logs/gtdbtk.log",
    benchmark:
        "benchmarks/gtdbtk.tsv"
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        "if ls {input}/*.fa > /dev/null 2>&1; then "
        "  export GTDBTK_DATA_PATH={params.db}; "
        "  gtdbtk classify_wf --genome_dir {input} --extension fa "
        "  --out_dir {params.outdir} --cpus {threads} --pplacer_cpus {params.pplacer} "
        "  {params.ani} --force > {log} 2>&1; "
        "fi; "
        "files=$(ls {params.outdir}/gtdbtk.*.summary.tsv 2>/dev/null || true); "
        "if [ -n \"$files\" ]; then "
        "  awk 'FNR==1 && NR!=1 {{next}} {{print}}' $files > {output}; "
        "else "
        "  printf 'user_genome\\tclassification\\tclosest_genome_reference\\n' > {output}; "
        "fi"


rule genomad:
    input:
        "06_classify/genomad_input.fa",
    output:
        plasmid="06_classify/genomad/plasmid_summary.tsv",
        virus="06_classify/genomad/virus_summary.tsv",
    params:
        outdir=lambda w, output: os.path.join(os.path.dirname(output.plasmid), "out"),
        db=DB["genomad"],
        extra=config["genomad"]["extra"],
    threads: job("genomad")["threads"]
    resources:
        mem_mb=cap_mem(job("genomad")["mem_mb"]),
        runtime=job("genomad")["runtime"],
    log:
        "logs/genomad.log",
    benchmark:
        "benchmarks/genomad.tsv"
    conda:
        "../envs/genomad.yaml"
    shell:
        "if [ -s {input} ]; then "
        "  genomad end-to-end --cleanup --threads {threads} {params.extra} "
        "  {input} {params.outdir} {params.db} > {log} 2>&1; "
        "  cp {params.outdir}/genomad_input_summary/genomad_input_plasmid_summary.tsv {output.plasmid}; "
        "  cp {params.outdir}/genomad_input_summary/genomad_input_virus_summary.tsv {output.virus}; "
        "else "
        "  printf 'seq_name\\tplasmid_score\\n' > {output.plasmid}; "
        "  printf 'seq_name\\tvirus_score\\n' > {output.virus}; "
        "  echo 'no sequences to screen' > {log}; "
        "fi"


def classification_inputs(wildcards):
    inp = {"candidates": "06_classify/candidates.tsv"}
    if not CLS["skip_checkm2"]:
        inp["checkm2"] = "06_classify/checkm2/quality_report.tsv"
    if not CLS["skip_gtdbtk"]:
        inp["gtdbtk"] = "06_classify/gtdbtk/summary.tsv"
    if not CLS["skip_genomad"]:
        inp["plasmid"] = "06_classify/genomad/plasmid_summary.tsv"
        inp["virus"] = "06_classify/genomad/virus_summary.tsv"
    return inp


def classification_opts(wildcards, input):
    opts = []
    for key in ("checkm2", "gtdbtk", "plasmid", "virus"):
        if key in input.keys():
            opts.append(f"--{key} {input[key]}")
    if DB["gtdb_metadata"] and not CLS["skip_gtdbtk"]:
        opts.append("--gtdb-metadata " + " ".join(DB["gtdb_metadata"]))
    return " ".join(opts)


rule summarise_classification:
    input:
        unpack(classification_inputs),
    output:
        "06_classify/classification.tsv",
    params:
        script=SCRIPTS / "summarise_classification.py",
        opts=classification_opts,
        min_chrom_len=CLS["min_chrom_len"],
        comp=CLS["min_completeness"],
        cont=CLS["max_contamination"],
    resources:
        mem_mb=cap_mem(8000),
        runtime=120,
    log:
        "logs/summarise_classification.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} --candidates {input.candidates} {params.opts} "
        "--min-chrom-len {params.min_chrom_len} --min-completeness {params.comp} "
        "--max-contamination {params.cont} --out {output} > {log} 2>&1"
