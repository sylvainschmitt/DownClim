rule summarise_chelsa:
    input:
        "results/chelsa/raw/{country}_CHELSA_1981-2010_V.2.1.nc"
    output:
        "results/chelsa/mean/{country}_CHELSA_{period}_V.2.1.tif"
    log:
        "results/logs/summarise_chelsa_{country}_{period}.log"
    benchmark:
        "results/benchmarks/summarise_chelsa_{country}_{period}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        period="{period}"
    script:
      "../scripts/summarise_chelsa.py"
