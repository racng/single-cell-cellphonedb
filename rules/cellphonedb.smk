rule prep_counts:
    input:
        adata=config['adata']
    output:
        counts=temp("cellphonedb/{sample}/counts.txt"),
        meta=temp("cellphonedb/{sample}/meta.txt")
        # counts="cellphonedb/{sample}/counts.txt",
        # meta="cellphonedb/{sample}/meta.txt"
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/prep_counts.py"

rule run_cellphoendb:
    input:
        counts='cellphonedb/{sample}/counts.txt',
        meta='cellphonedb/{sample}/meta.txt'
    output:
        "cellphonedb/{sample}/deconvoluted.txt",
        "cellphonedb/{sample}/pvalues.txt",
        "cellphonedb/{sample}/means.txt",
        "cellphonedb/{sample}/significant_means.txt",
    log:
        "cellphonedb/{sample}/cellphonedb.statistical_analysis.log"
    params:
        outdir="cellphonedb/{sample}",
        iterations=config['cellphonedb']['args']['iterations'],
        threshold=config['cellphonedb']['args']['threshold']
    threads: config['cellphonedb']['args']['threads']
    conda:
        "../envs/cellphonedb.yaml"
    shell:
        "cellphonedb method statistical_analysis {input.meta} {input.counts} "
        "--output-path {params.outdir} "
        "--threads={threads} "
        "--iterations={params.iterations} "
        "--threshold={params.threshold} &> {log}"
        

        