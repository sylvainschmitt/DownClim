rule get_chelsa2:
    input:
        "results/countries/{area}.shp"
    output:
        "results/chelsa2/raw/{area}_chelsa2.nc"
    log:
        "results/logs/get_chelsa2_{area}.log"
    benchmark:
        "results/benchmarks/get_chelsa2_{area}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        area="{area}",
        variables=config["variables"],
        time_frequency=config["time_frequency"],
        base_years=config["base_years"]
    shell:
      "../scripts/get_chelsa2.py"
