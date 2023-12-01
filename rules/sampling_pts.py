rule sampling_pts:
    input:
        "results/countries/{country}_gadm"
    output:
        "results/countries/{country}_sampling.tsv"
    log:
        "results/logs/sampling_pts_{country}.log"
    benchmark:
        "results/benchmarks/sampling_pts_{country}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        country="{country}",
        log10_eval_pts=config["log10_eval_pts"]
    script:
      "../scripts/sampling_pts.R"
