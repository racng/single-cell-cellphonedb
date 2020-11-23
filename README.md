# Snakemake workflow: single-cell-cellphonedb

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)

Snakemake workflow for running [cellphonedb](https://github.com/Teichlab/cellphonedb) on a combined annotated dataset (AnnData object) saved as h5ad file.
Workflow structure is modelled after cookiecutter snakemake [template](https://github.com/snakemake-workflows/cookiecutter-snakemake-workflow).  

## Authors

* Rachel Ng (@racng)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing `config.yaml`.
- `outdir`: Output directory (default: `output`)
- `adata`: Path to `h5ad` file
- `gene`: Key in adata.var that contains gene ids. Put `index` if using adata.var_names
- `celltype:key`: Key in adata.obs that contains cell type labels
- `celltype:groups`: Cell types to use (Put `[]` if using all)
- `sample:key`: Key in adata.obs that contains sample labels
- `sample:groups`: Samples to analyze (Put `[]` if using all)
- `cellphonedb:args:`: Additional arguments for the cellphonedb cli (counts-data, iterations, threshold, threads)
### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [mamba](https://github.com/mamba-org/mamba):

    conda install -c conda-forge mamba
    mamba create -c conda-forge -c bioconda -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally using `$N` cores (or 'all' cores)

    snakemake --use-conda --cores $N

Use alternate config files using the `--configfile` option.

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for further details.

### Step 5: Investigate results

Cellphonedb outputs are organized by sample subfolders under the specified `outdir`.




