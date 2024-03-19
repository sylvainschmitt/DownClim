configfile: "config/config_ex.yml"

import pandas as pd
proj = pd.read_table(config["projections"])
dom = pd.read_table(config["domains"])
areas = config["area"]
dom = dom.query("area == @areas")
proj_dom = proj.merge(dom, on='domain', how='left')
proj_dom["id"] = proj_dom.area + "_" + proj_dom.project + "_" \
                 + proj_dom.domain + "_" + proj_dom.institute \
                 + "_" + proj_dom.model + "_" +  proj_dom.experiment + \
                 "_" +  proj_dom.ensemble + "_" + proj_dom.rcm + \
                 "_" + proj_dom.downscaling

# rules #
rule all:
   input:
      # projections
      expand("results/projection/downscaled/{proj}_{baseline}_{aggregation}_{period_proj}_{period_base}_{ds_method}.nc",
             proj=proj_dom.id,
             baseline=config["baseline"],
             aggregation=config["aggregation"],
             period_base=config["hist_years"],
             period_proj=config["proj_years"],
             ds_method=config["ds_method"]),
      # evaluations
      "results/evaluation/histograms.tsv",
      "results/evaluation/evaluations.tsv"

## area ##
include: "rules/get_area.py"

## baseline ##
include: "rules/get_chelsa2.py"
include: "rules/aggregate_base.py"

## projection ##
include: "rules/get_cordex.py"
include: "rules/get_cmip6.py"
include: "rules/aggregate_proj.py"

## downscaling ##
include: "rules/downscale_bc.py"

## evaluation ##
include: "rules/hist_base.py"
include: "rules/hist_proj_raw.py"
include: "rules/hist_proj_ds.py"
include: "rules/merge_hist.py"
include: "rules/eval_proj_raw.py"
include: "rules/eval_proj_ds.py"
include: "rules/merge_eval.py"
