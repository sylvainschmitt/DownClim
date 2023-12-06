rule get_cru4:
    input:
        "results/countries/{area}.shp"
    output:
        "results/cru4/raw/{area}_cru4.nc"
    log:
        "results/logs/get_cru4_{area}.log"
    benchmark:
        "results/benchmarks/get_cru4_{area}.benchmark.txt"
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
      "../scripts/get_cru4.py"
