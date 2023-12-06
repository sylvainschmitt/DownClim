configfile: "config/config_test.yml"

import pandas as pd
exp = pd.read_table(config["projections"]).set_index("exp", drop=False)

# rules #
rule all:
   input:
      expand("results/cordex/downscaled/{exp}_{period_proj}_{period_base}.nc",
             exp=exp.exp,
             period_base=config["hist_years"],
             period_proj=config["proj_years"]),
      expand("results/evaluation/{exp}_{period_eval}_{period_base}.tsv",
             exp=exp.exp,
             period_base=config["hist_years"],
             period_eval=config["eval_years"])

## area ##
include: "rules/get_area.py"

## baselin ##
include: "rules/get_chelsa.py"

## projection ##
include: "rules/get_cordex.py"

## downscaling ##
include: "rules/downscale_bc.py"

## evaluation ##
include: "rules/evaluate_bc.py"
