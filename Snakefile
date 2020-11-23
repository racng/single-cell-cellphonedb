configfile: "config.yaml"

include: "rules/cellphonedb.smk"    

rule all:
    input:
        expand("{outdir}/finished.txt", outdir=config['outdir'])        
