rule get_script_cordex:
    output:
        "results/cordex/wget/{gcm}_{rcm}_{rcp}.sh"
    log:
        "results/logs/get_script_cordex_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/get_script_cordex_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        gcm="{gcm}",
        rcm="{rcm}",
        rcp="{rcp}",
        domain=config["domain"]
    script:
      "../scripts/get_script_cordex.py"
