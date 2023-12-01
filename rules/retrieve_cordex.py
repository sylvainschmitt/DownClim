rule retrieve_cordex:
    input:
        "results/cordex/wget/{var}_{domain}_{gcm}_{rcm}_{rcp}.sh"
    output:
        directory("results/cordex/raw/{var}_{domain}_{gcm}_{rcm}_{rcp}")
    log:
        "results/logs/retrieve_cordex_{var}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/retrieve_cordex_{var}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        script="{var}_{domain}_{gcm}_{rcm}_{rcp}.sh"
    shell:
      # "module load java ; " # only for muse, need to add a condition
      # "ESG_HOME='.esg/' ; " # only for muse, need to add a condition
      "mkdir -p {output} ; "
      "cp {input} {output} ; "
      "cd {output} ; "
      "bash {params.script} -u ; "
      "bash {params.script}"
