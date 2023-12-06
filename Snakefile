configfile: "config/config_test.yml"
experiment = "config/experiment_test.tsv"

import pandas as pd
exp = pd.read_table(experiment).set_index("exp", drop=False)

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

## country ##
include: "rules/get_country.py"

## chelsa ##
include: "rules/get_chelsa.py"
include: "rules/summarise_chelsa.py"

## cordex ##
include: "rules/get_cordex.py"
include: "rules/summarise_cordex.py"
include: "rules/get_anomalies.py"

## downscaling ##
include: "rules/downscale.py"
include: "rules/evaluate.py"
# include: "rules/eval_all.py"
# include: "rules/eval_tab.py"
