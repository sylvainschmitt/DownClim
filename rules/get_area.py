rule get_area:
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
    resources:
        mem_mb=1000
    conda:
        "../envs/gadm.yml"
    params:
        country="{country}",
        log10_eval_pts=config["log10_eval_pts"]
    script:
      "../scripts/get_area.py"
