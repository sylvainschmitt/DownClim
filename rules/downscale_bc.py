rule downscale_bc:
    input:
        proj_proj="results/projection/means/{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{aggregation}_{period_proj}.nc",
        proj_base="results/projection/means/{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{aggregation}_{period_base}.nc",
        base_base="results/{base}/means/{area}_{base}_{aggregation}_{period_base}.nc"
    output:
        "results/projection/downscaled/{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_base}_bc.nc"
    log:
        "results/logs/downscale_bc_{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_base}.log"
    benchmark:
        "results/benchmarks/downscale_bc_{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}_{aggregation}_{period_proj}_{period_base}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        period_base="{period_base}",
        period_proj="{period_proj}"
    script:
      "../scripts/downscale_bc.py"
