import scanpy as sc
import os

# read adata
adata = sc.read(snakemake.input['adata'])
# get all available samples in samples key
samples_cat = adata.obs[snakemake.config['sample']['key']].unique()

# use all samples if none specified
if snakemake.config['sample']['groups']==[]:
    samples = samples_cat
else:
    # check specified samples exist
    samples = snakemake.config['sample']['groups']
    missing = set(samples) - set(samples_cat)
    if len(missing)>0:
        raise ValueError("Samples not found: ", missing)

# write samples
path = f"{snakemake.config['outdir']}/samples"
if not os.path.isdir(path):
    os.system("mkdir -p " + path)
for sample in samples:
    with open("{}/{}.txt".format(path, sample), "w") as f:
        f.write(sample)
