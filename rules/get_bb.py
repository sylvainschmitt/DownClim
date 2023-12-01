rule get_bb:
    output:
        "results/countries/{country}_bb.txt",
        directory("results/countries/{country}_gadm"),
        "results/countries/{country}.png"
    log:
        "results/logs/get_bb_{country}.log"
    benchmark:
        "results/benchmarks/get_bb_{country}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        country="{country}"
    script:
      "../scripts/get_bb.R"
