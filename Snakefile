configfile: "config/config_test.yml"
experiment = "config/experiment_test.tsv"

import pandas as pd
exp = pd.read_table(experiment).set_index("exp", drop=False)

# rules #

rule all:
   input:
      expand("results/downscaled/proj_{var}_{exp}_{month}.tif",
              exp=exp.exp,
              month=config["months"],
              var=config["variables"]),
      "results/evaluation/summary.tsv"

## country ##
include: "rules/get_bb.py"
include: "rules/sampling_pts.py"

## chelsa ##
include: "rules/retrieve_chelsa.py"
include: "rules/crop_chelsa.py"
include: "rules/summarise_chelsa_hist.py"
include: "rules/summarise_chelsa_eval.py"

## cordex ##
include: "rules/get_script_cordex.py"
include: "rules/retrieve_cordex.py"
include: "rules/merge_cordex.py"
include: "rules/project_cordex.py"
include: "rules/crop_cordex.py"
include: "rules/summarise_cordex_hist.py"
include: "rules/summarise_cordex_eval.py"
include: "rules/summarise_cordex_proj.py"
include: "rules/cordex_anomalies.py"

## downscaling ##
include: "rules/downscale.py"
include: "rules/evaluate.py"
include: "rules/eval_all.py"
include: "rules/eval_tab.py"
