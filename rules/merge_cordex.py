rule merge_cordex:
    input:
        "results/cordex/raw/{gcm}_{rcm}_{rcp}"
    output:
        "results/cordex/merged/{gcm}_{rcm}_{rcp}.nc"
    log:
        "results/logs/merge_cordex_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/merge_cordex_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    shell:
      "cdo mergetime {input[0]}/*.nc {output[0]}"
