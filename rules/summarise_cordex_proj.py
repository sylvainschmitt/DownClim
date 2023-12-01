rule summarise_cordex_proj:
    input:
        "results/cordex/cropped/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    output:
        "results/cordex/proj/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    log:
        "results/logs/summarise_cordex_proj_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/summarise_cordex_proj_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        years=config["proj_years"]
    script:
      "../scripts/summarise_mean_nc.R"
