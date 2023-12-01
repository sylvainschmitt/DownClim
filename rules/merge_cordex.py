rule merge_cordex:
    input:
        "results/cordex/raw/{var}_{domain}_{gcm}_{rcm}_historical",
        "results/cordex/raw/{var}_{domain}_{gcm}_{rcm}_{rcp}"
    output:
        temp("results/cordex/merged/{var}_{domain}_{gcm}_{rcm}_{rcp}_historical.nc"),
        temp("results/cordex/merged/{var}_{domain}_{gcm}_{rcm}_{rcp}_only.nc"),
        temp("results/cordex/merged/{var}_{domain}_{gcm}_{rcm}_{rcp}.nc")
    log:
        "results/logs/merge_cordex_{var}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/merge_cordex_{var}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    shell:
      "cdo mergetime {input[0]}/*.nc {output[0]} ; "
      "cdo mergetime {input[1]}/*.nc {output[1]} ; "
      "cdo mergetime {output[0]} {output[1]} {output[2]}"
