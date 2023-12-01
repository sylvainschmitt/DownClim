rule summarise_cordex_eval:
    input:
        "results/cordex/cropped/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    output:
        "results/cordex/eval/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    log:
        "results/logs/summarise_cordex_eval_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/summarise_cordex_eval_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        years=config["eval_years"]
    script:
      "../scripts/summarise_mean_nc.R"
