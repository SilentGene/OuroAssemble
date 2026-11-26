# 03/04 — standardise assemblies, extract circular candidates,
#          trim duplicated end overlaps and rotate (A+B -> B+A)


def standardise_args(wildcards):
    """assembler-specific inputs for scripts/standardise_assembly.py"""
    asm = wildcards.asm
    if asm == "metamdbg":
        return {"fasta": "02_assembly/metamdbg/contigs.fasta.gz"}
    if asm == "hifiasm_meta":
        return {"gfa": "02_assembly/hifiasm_meta/asm.p_ctg.gfa"}
    if asm == "metaflye":
        return {
            "fasta": "02_assembly/metaflye/assembly.fasta",
            "info": "02_assembly/metaflye/assembly_info.txt",
        }
    if asm == "hicanu":
        return {"fasta": "02_assembly/hicanu/asm.contigs.fasta"}
    raise ValueError(asm)


def standardise_opts(wildcards, input):
    opts = []
    for key in ("fasta", "gfa", "info"):
        if key in input.keys():
            opts.append(f"--{key} {input[key]}")
    if wildcards.asm == "hifiasm_meta" and config["hifiasm_meta"]["include_rescue"]:
        opts.append("--rescue 02_assembly/hifiasm_meta/asm.rescue.fa")
    return " ".join(opts)


rule standardise_assembly:
    input:
        unpack(standardise_args),
    output:
        fasta="03_contigs/{asm}/contigs.fa",
        tsv="03_contigs/{asm}/contigs.tsv",
        candidates="03_contigs/{asm}/circular_candidates.fa",
    params:
        opts=standardise_opts,
        min_len=config["candidates"]["min_circular_len"],
        script=SCRIPTS / "standardise_assembly.py",
    resources:
        mem_mb=cap_mem(8000),
        runtime=120,
    log:
        "logs/standardise_assembly/{asm}.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} --assembler {wildcards.asm} {params.opts} "
        "--min-circular-len {params.min_len} "
        "--out-fasta {output.fasta} --out-tsv {output.tsv} "
        "--out-candidates {output.candidates} > {log} 2>&1"


rule trim_rotate:
    input:
        contigs="03_contigs/{asm}/contigs.fa",
        candidates="03_contigs/{asm}/circular_candidates.fa",
    output:
        trimmed="04_candidates/{asm}/trimmed.fa",
        rotated="04_candidates/{asm}/rotated.fa",
        ref="04_candidates/{asm}/validation_ref.fa",
        tsv="04_candidates/{asm}/trim_rotate.tsv",
        bed="04_candidates/{asm}/junctions.bed",
    params:
        script=SCRIPTS / "trim_rotate.py",
        search=config["trim_rotate"]["overlap_search"],
        min_overlap=config["trim_rotate"]["min_overlap"],
        min_identity=config["trim_rotate"]["min_identity"],
        max_frac=config["trim_rotate"]["max_overlap_frac"],
    resources:
        mem_mb=cap_mem(8000),
        runtime=120,
    log:
        "logs/trim_rotate/{asm}.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} --contigs {input.contigs} --candidates {input.candidates} "
        "--overlap-search {params.search} --min-overlap {params.min_overlap} "
        "--min-identity {params.min_identity} --max-overlap-frac {params.max_frac} "
        "--out-trimmed {output.trimmed} --out-rotated {output.rotated} "
        "--out-ref {output.ref} --out-tsv {output.tsv} --out-bed {output.bed} > {log} 2>&1"
