rule summarise_chelsa_hist:
    input:
        expand("results/chelsa/cropped/{var}_{country}_{month}_{year}.tif",
                year=config["hist_years"], allow_missing = True)
    output:
        "results/chelsa/hist/{var}_{country}_{month}.tif"
    log:
        "results/logs/summarise_chelsa_hist_{var}_{country}_{month}.log"
    benchmark:
        "results/benchmarks/summarise_chelsa_hist_{var}_{country}_{month}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    script:
      "../scripts/summarise_mean_rast.R"
