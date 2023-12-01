rule get_script_cordex:
    output:
        "results/cordex/wget/{var}_{domain}_{gcm}_{rcm}_{rcp}.sh"
    log:
        "results/logs/get_script_cordex_{var}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/get_script_cordex_{var}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        domain="{domain}",
        gcm="{gcm}",
        rcm="{rcm}",
        rcp="{rcp}",
        var="{var}"
    script:
      "../scripts/get_script_cordex.py"
