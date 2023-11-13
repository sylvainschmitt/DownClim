rule rule:
    output:
        directory("results/{key}")
    log:
        "results/logs/rule_{key}.log"
    benchmark:
        "results/benchmarks/rule_{key}.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-troll/releases/download/0.0.1/sylvainschmitt-singularity-troll.latest.sif" # dummy
    threads: 1
    script:
        "../scripts/rule.R"
        