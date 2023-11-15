# [DownClim - Downscale Climate Projections](https://sylvainschmitt.github.io/DownClim/)
Sylvain Schmitt -
Nov 13, 2023

- [Installation](#installation)
- [Credentials](#credentials)
- [Usage](#usage)
- [Workflow](#workflow)
- [Data](#data)

[`snakemake` &
`singularity`](https://github.com/sylvainschmitt/snakemake_singularity)
workflow to downscale climate projections

**Description.**

![Workflow.](dag/dag.svg)

# Installation

This workflow is built on:

- [x] Python ≥3.5
- [x] Snakemake ≥5.24.1
- [x] Singularity ≥3.7.3

Once installed simply clone the workflow:

``` bash
git clone git@github.com:sylvainschmitt/DownClim.git
cd DownClim
```

# Credentials

*Detail how to get the different credentials.*

# Usage

``` bash
module load bioinfo/Snakemake/7.20.0 # for test on nod depending on your HPC
snakemake -np # dry run
snakemake --dag | dot -Tsvg > dag/dag.svg # dag
snakemake -j 1 --resources mem_mb=10000 # local run (test)
sbatch job_muse.sh # HPC run with slurm
```

# Workflow

### [rule](https://github.com/sylvainschmitt/smkTemplate/blob/main/rules/rule.smk)

- Script:
  [`rule.R`](https://github.com/sylvainschmitt/smkTemplate/blob/main/scripts/rule.R)

Dummy example of a rule based on an R script.

# Data

*Detail the different data sources.*
