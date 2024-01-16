def list_areas(wildcards):
    return dom[dom.domain == wildcards.domain].area

def list_cordex_in(wildcards):
    return expand("results/countries/{area}.shp", 
                   area=dom[dom.domain == wildcards.domain].area)

rule get_cordex:
    input:
        list_cordex_in
    output:
        directory("results/projection/raw/CORDEX_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}")
    log:
        "results/logs/get_CORDEX_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}.log"
    benchmark:
        "results/benchmarks/get_cordex_CORDEX_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}.benchmark.txt"
    threads: 10
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        area=list_areas,
        activity="{activity}",
        domain="{domain}",
        institute="{institute}",
        model="{model}",
        experiment="{experiment}",
        ensemble="{ensemble}",
        rcm="{rcm}",
        downscaling="{downscaling}",
        variables=config["variables"],
        time_frequency=config["time_frequency"],
        esgf_credential=config["esgf_credential"]
    script:
      "../scripts/get_cordex.py"
