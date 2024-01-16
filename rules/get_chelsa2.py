rule get_chelsa2:
    input:
        expand("results/countries/{area}.shp", area=config["area"])
    output:
        expand("results/chelsa2/raw/{area}_chelsa2.nc", area=config["area"])
    log:
        "results/logs/get_chelsa2.log"
    benchmark:
        "results/benchmarks/get_chelsa2.benchmark.txt"
    threads: 6
    resources:
        mem_mb=10000
    conda:
        "../envs/xarray.yml"
    params:
        area=config["area"],
        variables=config["variables"],
        time_frequency=config["time_frequency"],
        base_years=config["base_years"],
        tmp="results/chelsa2/raw/tmp/"
    script:
      "../scripts/get_chelsa2.py"
