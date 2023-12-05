rule summarise_cordex:
    input:
        expand("results/cordex/raw/{country}_{domain}_{gcm}_{rcm}_{rcp}_{years}.nc",
               years=config["cordex_years"], allow_missing = True)
    output:
        "results/cordex/mean/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period}.nc"
    log:
        "results/logs/summarise_cordex_{period}_{country}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/summarise_cordex_{period}_{country}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        period="{period}"
    script:
      "../scripts/summarise_cordex.py"
