rule evaluate:
    input:
        "results/cordex/downscaled/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_eval}_{period_base}.nc",
        "results/cordex/mean/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_eval}.nc",
        "results/chelsa/mean/{country}_CHELSA_{period_eval}_V.2.1.tif",
        "results/countries/{country}_pts.shp"
    output:
        "results/evaluation/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_eval}_{period_base}.tsv"
    log:
        "results/logs/evaluate_{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_eval}_{period_base}.log"
    benchmark:
        "results/benchmarks/evaluate_{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_eval}_{period_base}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        period_base="{period_base}",
        period_proj="{period_eval}"
    script:
      "../scripts/evaluate.py"
