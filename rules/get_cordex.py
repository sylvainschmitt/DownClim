rule get_cordex:
    input:
        "results/countries/{country}.shp"
    output:
        "results/cordex/raw/{country}_{domain}_{gcm}_{rcm}_{rcp}_{years}.nc"
    log:
        "results/logs/get_cordex_{country}_{domain}_{gcm}_{rcm}_{rcp}_{years}.log"
    benchmark:
        "results/benchmarks/get_cordex_{country}_{domain}_{gcm}_{rcm}_{rcp}_{years}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        country="{country}",
        domain="{domain}",
        gcm="{gcm}",
        rcm="{rcm}",
        rcp="{rcp}",
        years="{years}",
        variables=config["variables"]
    shell:
      "../scripts/get_cordex.py"
