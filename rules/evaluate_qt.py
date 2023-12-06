rule evaluate_qt:
    input:
        proj="results/projection/raw/{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}.nc",
        base_eval="results/{base_eval}/raw/{area}_{base_eval}.nc",
        ds="results/projection/downscaled/{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{baseline}_{period_eval}_{period_base}_qt.nc",
        pts="results/countries/{area}_pts.shp"
    output:
        "results/evaluation/{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{baseline}_{period_eval}_{period_base}_{base_eval}_qt.tsv"
    log:
        "results/logs/evaluate_qt_{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{baseline}_{period_eval}_{period_base}_{base_eval}.log"
    benchmark:
        "results/benchmarks/evaluate_qt_{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{baseline}_{period_eval}_{period_base}_{base_eval}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        period_base="{period_base}",
        period_eval="{period_eval}"
    script:
      "../scripts/evaluate_qt.py"
