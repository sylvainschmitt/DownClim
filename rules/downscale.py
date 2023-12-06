rule downscale:
    input:
        "results/cordex/anomaly/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_proj}_{period_base}.nc",
        "results/chelsa/mean/{country}_CHELSA_{period_base}_V.2.1.tif"
    output:
        "results/cordex/downscaled/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_proj}_{period_base}.nc"
    log:
        "results/logs/downscaled_{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_proj}_{period_base}.log"
    benchmark:
        "results/benchmarks/downscaled_{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_proj}_{period_base}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        period_base="{period_base}",
        period_proj="{period_proj}"
    script:
      "../scripts/downscale.py"
