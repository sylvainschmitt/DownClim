rule ensemble_sma:
    input:
        expand("results/projection/downscaled/{proj}_{baseline}_{period_proj}_{period_base}_{ds_method}.nc",
               proj=proj.projection, allow_missing=True)
    output:
        "results/projection/ensemble/{baseline}_{period_proj}_{period_base}_{ds_method}_sma.nc"
    log:
        "results/logs/ensemble_sma_{baseline}_{period_proj}_{period_base}_{ds_method}.log"
    benchmark:
        "results/benchmarks/ensemble_sma_{baseline}_{period_proj}_{period_base}_{ds_method}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        period_base="{period_base}",
        period_proj="{period_proj}"
    script:
      "../scripts/ensemble_sma.py"
