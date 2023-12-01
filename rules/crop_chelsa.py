rule crop_chelsa:
    input:
        "results/chelsa/raw/{var}_{month}_{year}.tif",
        "results/countries/{country}_gadm"
    output:
        temp("results/chelsa/cropped/{var}_{country}_{month}_{year}.tif")
    log:
        "results/logs/crop_chelsa_{var}_{country}_{month}_{year}.log"
    benchmark:
        "results/benchmarks/crop_chelsa_{var}_{country}_{month}_{year}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        country="{country}"
    script:
      "../scripts/crop_chelsa.R"
