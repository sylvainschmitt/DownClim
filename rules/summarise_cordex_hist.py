rule summarise_cordex_hist:
    input:
        "results/cordex/cropped/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    output:
        "results/cordex/hist/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    log:
        "results/logs/summarise_cordex_hist_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/summarise_cordex_hist_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        years=config["hist_years"]
    script:
      "../scripts/summarise_mean_nc.R"
