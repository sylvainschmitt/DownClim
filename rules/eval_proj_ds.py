rule eval_proj_ds:
    input:
        "results/projection/downscaled/{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_eval}_{ds_method}.nc",
        "results/{base}/means/{area}_{base_eval}_{aggregation}_{period_eval}.nc",
        "results/countries/{area}.shp"
    output:
        "results/evaluation/eval/ds/{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_eval}_{ds_method}_{base_eval}.tsv"
    log:
        "results/logs/eval_ds_{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_eval}_{ds_method}_{base_eval}.log"
    benchmark:
        "results/benchmarks/eval_ds_{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_eval}_{ds_method}_{base_eval}.benchmark.txt"
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
      ds_method="{ds_method}",
      base_eval="{base_eval}"
    script:
      "../scripts/eval_metrics.py"
