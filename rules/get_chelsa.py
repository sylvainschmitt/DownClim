rule get_chelsa:
    input:
        "results/countries/{country}.shp"
    output:
        "results/chelsa/raw/{country}_CHELSA_1981-2010_V.2.1.nc"
    log:
        "results/logs/get_chelsa_{country}.log"
    benchmark:
        "results/benchmarks/get_chelsa_{country}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        country="{country}",
        variables=config["variables"]
    shell:
      "../scripts/get_chelsa.py"
