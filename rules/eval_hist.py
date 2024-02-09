rule eval_hist:
    input:
        "results/projection/downscaled/{data}.nc"
    output:
        "results/evaluation/hist/{data}.nc"
    log:
        "results/logs/eval_hist_{data}.log"
    benchmark:
        "results/benchmarks/eval_hist_{data}.benchmark.txt"
    threads: 1
    script:
      "../scripts/eval_hist.py"
