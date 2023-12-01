rule summarise_chelsa_eval:
    input:
        expand("results/chelsa/cropped/{var}_{country}_{month}_{year}.tif",
                year=config["eval_years"], allow_missing = True)
    output:
        "results/chelsa/eval/{var}_{country}_{month}.tif"
    log:
        "results/logs/summarise_chelsa_eval_{var}_{country}_{month}.log"
    benchmark:
        "results/benchmarks/summarise_chelsa_eval_{var}_{country}_{month}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    script:
      "../scripts/summarise_mean_rast.R"
