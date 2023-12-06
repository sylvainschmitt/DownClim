rule get_country:
    output:
        "results/countries/{country}.shp",
        "results/countries/{country}.png",
        "results/countries/{country}_pts.shp",
        "results/countries/{country}_pts.png"
    log:
        "results/logs/get_country_{country}.log"
    benchmark:
        "results/benchmarks/get_country_{country}.benchmark.txt"
    threads: 1
    conda:
        "../envs/gadm.yml"
    resources:
        mem_mb=1000
    params:
        country="{country}",
        log10_eval_pts=config["log10_eval_pts"]
    script:
      "../scripts/get_country.py"
