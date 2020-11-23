checkpoint samples:
    input: 
        adata=config['adata']
    output:
        directory("{outdir}/samples")
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/read_adata.py"

rule prep_counts:
    input:
        adata=config['adata']
    output:
        counts=temp("{outdir}/cellphonedb/{sample}/counts.txt"),
        meta=temp("{outdir}/cellphonedb/{sample}/meta.txt")
        # counts="cellphonedb/{sample}/counts.txt",
        # meta="cellphonedb/{sample}/meta.txt"
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/prep_counts.py"

rule run_cellphoendb:
    input:
        counts='{outdir}/cellphonedb/{sample}/counts.txt',
        meta='{outdir}/cellphonedb/{sample}/meta.txt'
    output:
        "{outdir}/cellphonedb/{sample}/deconvoluted.txt",
        "{outdir}/cellphonedb/{sample}/pvalues.txt",
        "{outdir}/cellphonedb/{sample}/means.txt",
        "{outdir}/cellphonedb/{sample}/significant_means.txt",
    log:
        "{outdir}/cellphonedb/{sample}/cellphonedb.statistical_analysis.log"
    params:
        outdir="{outdir}/cellphonedb/{sample}",
        gene_id_type=config['cellphonedb']['args']['counts-data'],
        iterations=config['cellphonedb']['args']['iterations'],
        threshold=config['cellphonedb']['args']['threshold']
    threads: config['cellphonedb']['args']['threads']
    conda:
        "../envs/cellphonedb.yaml"
    shell:
        "cellphonedb method statistical_analysis {input.meta} {input.counts} "
        "--counts-data {params.gene_id_type} "
        "--output-path {params.outdir} "
        "--threads={threads} "
        "--iterations={params.iterations} "
        "--threshold={params.threshold} &> {log}"
        

def aggregate_results(wildcards):
    checkpoint_output = checkpoints.samples.get(**wildcards).output[0]
    return expand("{outdir}/cellphonedb/{sample}/means.txt", 
        outdir=config['outdir'],
        sample=glob_wildcards(os.path.join(checkpoint_output, "{sample}.txt")).sample)

rule aggregate:
    input:
        aggregate_results
    output:
        touch("{outdir}/finished.txt")
