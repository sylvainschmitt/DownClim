# DownClim - Downscale Climate Projections
Sylvain Schmitt -
Dec 5, 2023

- [Installation](#installation)
- [Credentials](#credentials)
- [Usage](#usage)
- [Workflow](#workflow)
  - [Country](#country)
  - [CHELSA](#chelsa)
  - [CORDEX](#cordex)
  - [Downscaling](#downscaling)
- [Data](#data)
- [Results](#results)

[`snakemake` &
`singularity`](https://github.com/sylvainschmitt/snakemake_singularity)
workflow to downscale climate projections

**Description.**

![Workflow.](dag/dag.svg)

# Installation

This workflow is built on:

- [x] Python ≥3.5
- [x] Snakemake ≥5.24.1
- [x] Conda ≥ *to precise*

Once installed simply clone the workflow:

``` bash
git clone git@github.com:sylvainschmitt/DownClim.git
cd DownClim
```

# Credentials

Data are retrieve from the [Institut Pierre-Simon Laplace
node](https://esgf-node.ipsl.upmc.fr/search/cordex-ipsl/). You need
first to [create an
account](https://esgf.github.io/esgf-user-support/user_guide.html#create-an-account)
on this page ([create
account](https://esgf-node.ipsl.upmc.fr/user/add/?next=http://esgf-node.ipsl.upmc.fr/search/cordex-ipsl/)
link at corner right).

Then you’ll need to register credentials locally to use the workflow. A
[help
page](https://esgf.github.io/esgf-user-support/user_guide.html?highlight=credentials%20pem#access-data-with-the-command-line-via-opendap)
is available. In linux, you need `myproxy-logon` installed and to run
this command line (with your user name):

``` bash
myproxy-logon -sesgf-node.ipsl.upmc.fr -l {user_name} -b -T -t 72 -o ~/.esg/credentials.pem
```

To run the workflow on a cluster, you can simply copy your local
credentials to the cluster. For instance:

``` bash
cp ~/.esg/credentials.pem /my_cluster/
```

# Usage

``` bash
module load bioinfo/Snakemake/7.20.0 # for test on nod depending on your HPC
snakemake -np # dry run
snakemake --dag | dot -Tsvg > dag/dag.svg # dag
snakemake -j 1 --resources mem_mb=10000 # local run (test)
sbatch job_muse.sh # HPC run with slurm
```

# Workflow

## Country

### [get_country](https://github.com/sylvainschmitt/DownClim/blob/xarray/rules/get_country.py)

- Script:
  [`get_country.py`](https://github.com/sylvainschmitt/DownClim/blob/xarray/scripts/get_country.py)
- Environment:
  [`gadm.yml`](https://github.com/sylvainschmitt/DownClim/blob/xarray/envs/gadm.yml)

Python script to get country limits with GADM and to define sampling
points in the land for evaluation.

## CHELSA

### [get_chelsa](https://github.com/sylvainschmitt/DownClim/blob/xarray/rules/get_chelsa.py)

- Script:
  [`get_chelsa.py`](https://github.com/sylvainschmitt/DownClim/blob/xarray/scripts/get_chelsa.py)
- Environment:
  [`xarray.yml`](https://github.com/sylvainschmitt/DownClim/blob/xarray/envs/xarray.yml)

Python script to download, crop, adjust CHELSA monthly variables.

### [summarise_chelsa](https://github.com/sylvainschmitt/DownClim/blob/xarray/rules/summarise_chelsa.py)

- Script:
  [`summarise_chelsa.py`](https://github.com/sylvainschmitt/DownClim/blob/xarray/scripts/summarise_chelsa.py)
- Environment: *to be defined*

Python script to summarise CHELSA monthly variables on a defined period.

## CORDEX

### [get_cordex](https://github.com/sylvainschmitt/DownClim/blob/xarray/rules/get_cordex.py)

- Script:
  [`get_cordex.py`](https://github.com/sylvainschmitt/DownClim/blob/xarray/scripts/get_cordex.py)
- Environment:
  [`xarray.yml`](https://github.com/sylvainschmitt/DownClim/blob/xarray/envs/xarray.yml)

Python script to download, crop, reproject, and adjust CORDEX monthly
variables.

### [summarise_cordex](https://github.com/sylvainschmitt/DownClim/blob/xarray/rules/summarise_cordex.py)

- Script:
  [`summarise_cordex.py`](https://github.com/sylvainschmitt/DownClim/blob/xarray/scripts/summarise_cordex.py)
- Environment: *to be defined*

Python script to summarise CORDEX monthly variables on a defined period.

### [get_anomalies](https://github.com/sylvainschmitt/DownClim/blob/xarray/rules/get_anomalies.py)

- Script:
  [`get_anomalies.py`](https://github.com/sylvainschmitt/DownClim/blob/xarray/scripts/get_anomalies.py)
- Environment: *to be defined*

Python script to compute CORDEX monthly anomalies between to defined
periods.

## Downscaling

### [downscale](https://github.com/sylvainschmitt/DownClim/blob/xarray/rules/downscale.py)

- Script:
  [`downscale.py`](https://github.com/sylvainschmitt/DownClim/blob/xarray/scripts/downscale.py)
- Environment: *to be defined*

Python script to compute CORDEX downscaled monthly values with bias
correction on CHELSA.

### [evaluate](https://github.com/sylvainschmitt/DownClim/blob/xarray/rules/evaluate.py)

- Script:
  [`evaluate.py`](https://github.com/sylvainschmitt/DownClim/blob/xarray/scripts/evaluate.py)
- Environment: *to be defined*

Python script to evaluate CORDEX downscaling versus raw projection on
the evluation period of CHELSA.

# Data

[**CORDEX**](https://cordex.org/)**: Coordinated Regional Climate
Downscaling Experiment**

*The CORDEX vision is to advance and coordinate the science and
application of regional climate downscaling through global
partnerships.*

[**CHELSA V2.1.1**](https://chelsa-climate.org/)**: Climatologies at
high resolution for the earth’s land surface areas**

*CHELSA (Climatologies at high resolution for the earth’s land surface
areas) is a very high resolution (30 arc sec, ~1km) global downscaled
climate data set currently hosted by the Swiss Federal Institute for
Forest, Snow and Landscape Research WSL. It is built to provide free
access to high resolution climate data for research and application, and
is constantly updated and refined.*

# Results
