rule downscale_qt:
    input:
        proj="results/projection/raw/{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}.nc",
        base="results/{baseline}/raw/{area}_{baseline}.nc"
    output:
        "results/projection/downscaled/{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{baseline}_{period_proj}_{period_base}_qt.nc"
    log:
        "results/logs/downscale_qt_{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{baseline}_{period_proj}_{period_base}.log"
    benchmark:
        "results/benchmarks/downscale_qt_{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{baseline}_{period_proj}_{period_base}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        period_base="{period_base}",
        period_proj="{period_proj}"
    script:
      "../scripts/downscale_qt.py"
