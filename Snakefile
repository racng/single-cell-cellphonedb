# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


#configfile: "config.yaml"
#report: "report/workflow.rst"

# include: "rules/init.smk"
workdir: config['workdir']
include: "rules/cellphonedb.smk"    

rule all:
    input:
        expand("cellphonedb/{sample}/means.txt", sample=config['sample']['groups'])
        
