rule get_area:
    output:
        "results/countries/{area}.shp",
        "results/countries/{area}.png"
    log:
        "results/logs/get_country_{area}.log"
    benchmark:
        "results/benchmarks/get_country_{area}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/gadm.yml"
    params:
        area="{area}",
        log10_eval_pts=config["log10_eval_pts"]
    script:
      "../scripts/get_area.py"
