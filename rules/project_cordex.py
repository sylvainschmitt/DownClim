rule project_cordex:
    input:
        "results/cordex/merged/{var}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    output:
        temp("results/cordex/projected/{var}_{domain}_{gcm}_{rcm}_{rcp}.nc")
    log:
        "results/logs/project_cordex_{var}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/project_cordex_{var}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    shell:
      "cdo remapbil,global_0.22 {input[0]} {output[0]}"
