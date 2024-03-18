rule eval_proj_raw:
    input:
        "results/projection/means/{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_eval}.nc",
        "results/{base}/means/{area}_{base_eval}_{aggregation}_{period_eval}.nc",
        "results/countries/{area}.shp"
    output:
        "results/evaluation/eval/raw/{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_eval}_{base_eval}.tsv"
    log:
        "results/logs/eval_raw_{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_eval}_{base_eval}.log"
    benchmark:
        "results/benchmarks/eval_raw_{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_eval}_{base_eval}.benchmark.txt"
    threads: 1
    params:
      area="{area}",
      origin="{project}",
      type="raw",
      domain="{domain}",
      institute="{institute}",
      model="{model}",
      experiment="{experiment}",
      ensemble="{ensemble}",
      rcm="{rcm}",
      downscaling="{downscaling}",
      base="{base}",
      aggregation="{aggregation}",
      period_proj="none",
      period_eval="{period_eval}",
      ds_method="none",
      base_eval="{base_eval}"
    script:
      "../scripts/eval_metrics.py"
