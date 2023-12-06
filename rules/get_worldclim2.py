rule get_worldclim2:
    input:
        "results/countries/{area}.shp"
    output:
        "results/worldclim2/raw/{area}_worldclim2.nc"
    log:
        "results/logs/get_worldclim2_{area}.log"
    benchmark:
        "results/benchmarks/get_worldclim2_{area}.benchmark.txt"
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
      "../scripts/get_worldclim2.py"
