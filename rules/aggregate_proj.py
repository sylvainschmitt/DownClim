rule aggregate_proj:
    input:
        "results/projection/raw/{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}"
    output:
        "results/projection/means/{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period}.nc"
    log:
        "results/logs/aggregate_proj_{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period}.log"
    benchmark:
        "results/benchmarks/aggregate_proj_{area}_{project}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        period="{period}",
        aggregation="{aggregation}",
        area="{area}"
    script:
      "../scripts/aggregate_proj.py"
