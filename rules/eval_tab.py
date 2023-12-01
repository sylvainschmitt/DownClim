rule eval_tab:
    input:
        "results/evaluation/all.tsv"
    output:
        "results/evaluation/summary.tsv"
    log:
        "results/logs/eval_tab.log"
    benchmark:
        "results/benchmarks/eval_tab.benchmark.txt"
    threads: 1
    script:
      "../scripts/eval_tab.R"
