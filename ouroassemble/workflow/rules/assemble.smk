# 02 — assembly with four HiFi metagenome assemblers (each optional)


rule assemble_metamdbg:
    input:
        FILTERED_READS,
    output:
        "02_assembly/metamdbg/contigs.fasta.gz",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        extra=config["metamdbg"]["extra"],
    threads: job("metamdbg")["threads"]
    resources:
        mem_mb=cap_mem(job("metamdbg")["mem_mb"]),
        runtime=job("metamdbg")["runtime"],
    log:
        "logs/assemble_metamdbg.log",
    benchmark:
        "benchmarks/assemble_metamdbg.tsv"
    conda:
        "../envs/metamdbg.yaml"
    shell:
        "metaMDBG asm --out-dir {params.outdir} --in-hifi {input} "
        "--threads {threads} {params.extra} > {log} 2>&1"


rule assemble_hifiasm_meta:
    input:
        FILTERED_READS,
    output:
        "02_assembly/hifiasm_meta/asm.p_ctg.gfa",
    params:
        prefix=lambda w, output: output[0].removesuffix(".p_ctg.gfa"),
        extra=config["hifiasm_meta"]["extra"],
    threads: job("hifiasm_meta")["threads"]
    resources:
        mem_mb=cap_mem(job("hifiasm_meta")["mem_mb"]),
        runtime=job("hifiasm_meta")["runtime"],
    log:
        "logs/assemble_hifiasm_meta.log",
    benchmark:
        "benchmarks/assemble_hifiasm_meta.tsv"
    conda:
        "../envs/hifiasm_meta.yaml"
    shell:
        "hifiasm_meta -t {threads} -o {params.prefix} {params.extra} {input} > {log} 2>&1"


rule assemble_metaflye:
    input:
        FILTERED_READS,
    output:
        fasta="02_assembly/metaflye/assembly.fasta",
        info="02_assembly/metaflye/assembly_info.txt",
    params:
        outdir=lambda w, output: os.path.dirname(output.fasta),
        extra=config["metaflye"]["extra"],
    threads: job("metaflye")["threads"]
    resources:
        mem_mb=cap_mem(job("metaflye")["mem_mb"]),
        runtime=job("metaflye")["runtime"],
    log:
        "logs/assemble_metaflye.log",
    benchmark:
        "benchmarks/assemble_metaflye.tsv"
    conda:
        "../envs/flye.yaml"
    shell:
        "flye --pacbio-hifi {input} --meta --out-dir {params.outdir} "
        "--threads {threads} {params.extra} > {log} 2>&1"


rule assemble_hicanu:
    input:
        FILTERED_READS,
    output:
        "02_assembly/hicanu/asm.contigs.fasta",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        genome_size=config["hicanu"]["genome_size"],
        params=config["hicanu"]["params"],
        extra=config["hicanu"]["extra"],
        mem_gb=lambda w, resources: max(1, int(resources.mem_mb / 1024)),
    threads: job("hicanu")["threads"]
    resources:
        mem_mb=cap_mem(job("hicanu")["mem_mb"]),
        runtime=job("hicanu")["runtime"],
    log:
        "logs/assemble_hicanu.log",
    benchmark:
        "benchmarks/assemble_hicanu.tsv"
    conda:
        "../envs/canu.yaml"
    shell:
        "canu -p asm -d {params.outdir} genomeSize={params.genome_size} "
        "-pacbio-hifi {input} {params.params} "
        "useGrid=false maxThreads={threads} maxMemory={params.mem_gb}g "
        "{params.extra} > {log} 2>&1"
