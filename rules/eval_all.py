rule eval_all:
    input:
        expand("results/evaluation/{var}_{exp}_{month}.tsv",
              exp=exp.exp,
              month=config["months"],
              var=config["variables"])
    output:
        "results/evaluation/all.tsv"
    log:
        "results/logs/eval_all.log"
    benchmark:
        "results/benchmarks/eval_all.benchmark.txt"
    threads: 1
    script:
      "../scripts/eval_all.R"
