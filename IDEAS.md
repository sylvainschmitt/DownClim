# DownClim - Downscale Climate Projections
S Schmitt, G Vieilledent , T Arsouze, A Mauri -
Nov 13, 2023

- [Introduction](#introduction)
- [Data](#data)
  - [Baseline](#baseline)
  - [Climate projections](#climate-projections)
- [Workflow](#workflow)
- [Examples](#examples)
- [Assesment](#assesment)
- [Next steps](#next-steps)
- [Poeple](#poeple)
- [References](#references)

# Introduction

The purpose of `DownClim` is to offer a tool for regional and national
climate projections including the mechanistic ‘dynamic’ downscaling of
the CORDEX initiative. `DownClim` is opposed to the direct statistical
downscaling of global climate projections found in WorldClim and CHELSA.
`DownClim` is directly inspired of AFRICLIM (Platts, Omeny, and Marchant
2014). The approach is justified by an improvement in regional
projections of CORDEX compared to CMIP (insert few references), although
it can increase uncertainty and sometimes be less reliable (insert
indian monsoon reference). The tool is an automated `singularity`
workflow easily reproducible and scalable associated to `singularity`
images for enhance reproducibility and portability.

As a first step we will develop the pipeline for three study countries
(French Guiana, Ivory Coast and New Caledonia) over three continents,
with a single baseline (CHELSA) and a reduced number of global (GCM) and
regional (RCM) climate models for a few scenarios (RCP or SSP). Behind
the choice of countries we have different hypotheses. (1) we expect New
Caledonia regional topography to be important for projections. (2) we
expect Ivory Coast regional climates to be important for projections.
(3) we assume French Guiana to be a control with mainly lowland and
homogeneous climate. To demonstrate the aim of the tool, we want to
validate and/or assess the projections for each country using (1)
projections on observed year (e.g. 2005\>2020, even if the differences
between scenarios will not be significant), and (2) the difference with
forest climatic extent with this product versus classical products.

# Data

## Baseline

*WorldClim, Chelsa, ERA5-Land, CRU, TAMSAT, CHIRPS, …*

We will use CHELSA (version, reference) for the example with
spatial_resolution and temporal_resolution. Ideally the tool should be
built for different baselines (next steps).

## Climate projections

Depending on temporal resolution of the baseline and the sutdy area, we
should list available data on CORDEX regarding:

- global climate models (GCM)
- regional climate models (RCM)
- temporal resolution (AFRICLIM averages over 20 years but we could use
  finer)
- spatial resolution (a priori the finest)
- scenarios (RCP or SSP)

# Workflow

Ideas on the fly:

- snakemake
- singularity
- a bit of everything
- GRASS GIS (more than R terra) for raster operation and interpolation
- R
- R ncdf4 (or more appropriate NetCDF tools)
- rmarkdown for automatic reporting of the pipeline?

# Examples

- French Guiana, Ivory Coast, New Caledonia
- Chelsa v2.1 - 1km - annual?
- All RCM x GCM available at 1km - annual?

# Assesment

- predicted vs observed 2005-2023
- tropical moist forest potential distribution based on CHELSA future vs
  DownClim projections

# Next steps

*To develop once first steps done.*

# Poeple

- Sylvain Schmitt (sylvain.schmitt@cirad.fr)
- Ghislain Vieilledent (ghislain.vieilledent@cirad.fr)
- Thomas Arsouze (thomas.arsouze@cirad.fr)
- Achille Mauri (mauri.achille@gmail.com)

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-platts2014" class="csl-entry">

Platts, Philip J., Peter A. Omeny, and Rob Marchant. 2014. “AFRICLIM:
High-Resolution Climate Projections for Ecological Applications in
Africa.” *African Journal of Ecology* 53 (1): 103–8.
<https://doi.org/10.1111/aje.12180>.

</div>

</div>
