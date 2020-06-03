import scanpy as sc
import os

# read adata
adata = sc.read(snakemake.input['adata'])
# get all available samples in samples key
samples_cat = adata.obs[snakemake.config['sample']['key']].astype('category').cat.categories
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
if not os.path.isdir("samples"):
    os.system("mkdir samples")
for sample in samples:
    with open("samples/{}.txt".format(sample), "w") as f:
        f.write(sample)
