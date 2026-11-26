# 02 — assembly with four HiFi metagenome assemblers (each optional)
#
# Resuming: HiCanu, metaFlye and hifiasm-meta continue from the files left in their work
# directory when a job is interrupted (e.g. killed at walltime) and simply run again.
# scripts/resume_guard.sh keeps a signature of reads + options in that directory and
# empties it instead when they changed, so stale state is never reused. Snakemake only
# deletes the declared output files, never these work directories. metaMDBG (1.4) cannot
# resume and starts over.

RESUME_GUARD = SCRIPTS / "resume_guard.sh"


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
        outdir=lambda w, output: os.path.dirname(output[0]),
        extra=config["hifiasm_meta"]["extra"],
        guard=RESUME_GUARD,
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
        # hifiasm-meta reuses its *.bin files (read correction, overlaps) by itself; if
        # they are unusable (e.g. truncated by the kill), redo everything with -i.
        """
        mode=$(bash {params.guard} {params.outdir} {input} {params.extra:q})
        echo "[ouroassemble] $(date): $mode run" >> {log}
        if ! hifiasm_meta -t {threads} -o {params.prefix} {params.extra} {input} >> {log} 2>&1; then
            [ "$mode" = resume ] || exit 1
            echo "[ouroassemble] $(date): resume failed, starting over (-i)" >> {log}
            hifiasm_meta -i -t {threads} -o {params.prefix} {params.extra} {input} >> {log} 2>&1
        fi
        """


rule assemble_metaflye:
    input:
        FILTERED_READS,
    output:
        fasta="02_assembly/metaflye/assembly.fasta",
        info="02_assembly/metaflye/assembly_info.txt",
    params:
        outdir=lambda w, output: os.path.dirname(output.fasta),
        extra=config["metaflye"]["extra"],
        guard=RESUME_GUARD,
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
        # Flye resumes from its last completed stage (params.json marks a started run)
        """
        mode=$(bash {params.guard} {params.outdir} {input} {params.extra:q})
        resume=""
        if [ "$mode" = resume ] && [ -s {params.outdir}/params.json ]; then resume="--resume"; fi
        echo "[ouroassemble] $(date): $mode run $resume" >> {log}
        flye --pacbio-hifi {input} --meta --out-dir {params.outdir} \
            --threads {threads} {params.extra} $resume >> {log} 2>&1
        """


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
        guard=RESUME_GUARD,
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
        # Canu resumes by itself when run again in the same directory (finished overlap
        # jobs, stores, ... are kept); threads/memory may differ from the first attempt.
        """
        mode=$(bash {params.guard} {params.outdir} {input} \
            {params.genome_size:q} {params.params:q} {params.extra:q})
        echo "[ouroassemble] $(date): $mode run" >> {log}
        canu -p asm -d {params.outdir} genomeSize={params.genome_size} \
            -pacbio-hifi {input} {params.params} \
            useGrid=false maxThreads={threads} maxMemory={params.mem_gb}g \
            {params.extra} >> {log} 2>&1
        """
