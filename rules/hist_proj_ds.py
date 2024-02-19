rule hist_proj_ds:
    input:
        "results/projection/downscaled/{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_eval}_{ds_method}.nc",
        "results/countries/{area}.shp"
    output:
        "results/evaluation/hist/ds/{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_eval}_{ds_method}.tsv"
    log:
        "results/logs/hist_ds_{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_eval}_{ds_method}.log"
    benchmark:
        "results/benchmarks/hist_ds_{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_eval}_{ds_method}.benchmark.txt"
    threads: 1
    params:
      area="{area}",
      origin="{project}",
      type="downscaled", # other is raw
      domain="{domain}",
      institute="{institute}",
      model="{model}",
      experiment="{experiment}",
      ensemble="{ensemble}",
      rcm="{rcm}",
      downscaling="{downscaling}",
      base="{base}",
      aggregation="{aggregation}",
      period_proj="{period_proj}",
      period_eval="{period_eval}",
      ds_method="{ds_method}"
    script:
      "../scripts/eval_hist.py"
