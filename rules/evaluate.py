rule evaluate:
    input:
        "results/downscaled/eval_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}_{month}.tif",
        "results/cordex/eval/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc",
        "results/chelsa/eval/{var}_{country}_{month}.tif",
        "results/countries/{country}_sampling.tsv"
    output:
        "results/evaluation/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}_{month}.tsv"
    log:
        "results/logs/evaluate_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}_{month}.log"
    benchmark:
        "results/benchmarks/evaluate_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}_{month}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        var="{var}",
        month="{month}",
        country="{country}"
    script:
      "../scripts/evaluate_rast.R"
