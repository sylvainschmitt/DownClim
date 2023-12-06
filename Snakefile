configfile: "config/config_dag.yml"

import pandas as pd
proj = pd.read_table(config["projections"]).set_index("projection", drop=False)

# rules #
rule all:
   input:
      # downsaled projection
      # expand("results/projection/downscaled/{proj}_{baseline}_{period_proj}_{period_base}_{ds_method}.nc",
      #        proj=proj.projection,
      #        baseline=config["baseline"],
      #        period_base=config["hist_years"],
      #        period_proj=config["proj_years"],
      #        ds_method=config["ds_method"]),
      # evaluation
      expand("results/evaluation/{proj}_{baseline}_{period_eval}_{period_base}_{base_eval}_{ds_method}.tsv",
             proj=proj.projection,
             baseline=config["baseline"],
             period_base=config["hist_years"],
             period_eval=config["eval_years"],
             base_eval=config["base_eval"],
             ds_method=config["ds_method"]),
      # ensemble
      expand("results/projection/ensemble/{baseline}_{period_proj}_{period_base}_{ds_method}_{ens_method}.nc",
             baseline=config["baseline"],
             period_base=config["hist_years"],
             period_proj=config["proj_years"],
             base_eval=config["base_eval"],
             ds_method=config["ds_method"],
             ens_method=config["ens_method"])

## area ##
include: "rules/get_area.py"

## baseline ##
include: "rules/get_chelsa2.py"
include: "rules/get_worldclim2.py"
include: "rules/get_cru4.py"

## projection ##
include: "rules/get_proj.py"

## downscaling ##
include: "rules/downscale_bc.py"
include: "rules/downscale_qt.py"

## evaluation ##
include: "rules/evaluate_bc.py"
include: "rules/evaluate_qt.py"

## ensemble ##
include: "rules/ensemble_sma.py"
include: "rules/ensemble_bma.py"
