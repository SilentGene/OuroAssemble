# OuroAssemble — database download workflow (`ouroassemble download-db`).
# Run with --directory <db dir>; every target `<name>.done` contains the
# absolute path to put into the config's `databases:` section.

GTDBTK_URL = config.get(
    "gtdbtk_url",
    "https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/"
    "gtdbtk_package/full_package/gtdbtk_data.tar.gz",
)


rule checkm2_db:
    output:
        "checkm2.done",
    log:
        "logs/checkm2_db.log",
    conda:
        "envs/checkm2.yaml"
    shell:
        "checkm2 database --download --path checkm2 > {log} 2>&1 && "
        "realpath $(find checkm2 -name '*.dmnd' | head -n 1) > {output}"


rule gtdbtk_db:
    output:
        "gtdbtk.done",
    params:
        url=GTDBTK_URL,
    threads: 4
    log:
        "logs/gtdbtk_db.log",
    conda:
        "envs/download.yaml"
    shell:
        "mkdir -p gtdbtk && "
        "curl -L --fail --retry 3 -o gtdbtk/gtdbtk_data.tar.gz {params.url} > {log} 2>&1 && "
        "tar -I pigz -xf gtdbtk/gtdbtk_data.tar.gz -C gtdbtk >> {log} 2>&1 && "
        "rm gtdbtk/gtdbtk_data.tar.gz && "
        "realpath $(ls -d gtdbtk/release* | sort | tail -n 1) > {output}"


rule genomad_db:
    output:
        "genomad.done",
    log:
        "logs/genomad_db.log",
    conda:
        "envs/genomad.yaml"
    shell:
        "genomad download-database genomad > {log} 2>&1 && "
        "realpath genomad/genomad_db > {output}"
