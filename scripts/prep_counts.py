import scanpy as sc
import pandas as pd

def subset_adata(
    adata,
    celltype_key,
    celltype_groups=None, 
    sample_key=None,
    sample=None
    ):
    adata = adata.copy()
    # remove nan
    adata = adata[~adata.obs[celltype_key].isnull()]
    # remove string converted nan values 
    adata = adata[adata.obs[celltype_key]!='nan']
    
    # include only certain celltypes if specified
    if celltype_groups is not None and len(celltype_groups)>0:
        adata = adata[adata.obs[celltype_key].isin(celltype_groups)]

    # get sample subset 
    if sample_key is not None:
        adata = adata[adata.obs[sample_key]==sample]
        if adata.shape[0] == 0:
            raise ValueError('No cells for {} in {}'.format(sample, sample_key))
    
    return adata

def write_counts(
    adata, 
    gene_key, 
    path, 
    ):
    
    # make count matrix
    df_expr_matrix = pd.DataFrame(adata.X.T.toarray())
    # set cell IDs as columns
    df_expr_matrix.columns = adata.obs.index
    # set adata var index or a var column (ex. EnsembleIDs, gene names) as index
    if gene_key=='index':
        df_expr_matrix.set_index(adata.var_names, inplace=True)
    else:
        df_expr_matrix.set_index(adata.var[gene_key], inplace=True)
    # write count matrix
    df_expr_matrix.to_csv(path, sep='\t')

def write_meta(
    adata,
    celltype_key,
    path
    ):
    # get cell IDs and annotation
    df_meta = pd.DataFrame({
        'Cell': list(adata.obs.index), 
        'cell_type': list(adata.obs[celltype_key])
        })
    df_meta.to_csv(path, sep='\t', index=False)
    


celltype_key = snakemake.config['celltype']['key']
celltype_groups= snakemake.config['celltype']['groups']
gene_key = snakemake.config['gene']['key']
sample_key = snakemake.config['sample']['key']
sample = snakemake.wildcards['sample']

path_counts = snakemake.output['counts']
path_meta = snakemake.output['meta']

adata = sc.read(snakemake.input['adata'])
subset = subset_adata(adata, celltype_key, celltype_groups, sample_key, sample)
write_counts(subset, gene_key, path_counts)
write_meta(subset, celltype_key, path_meta)
