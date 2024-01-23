def list_areas(wildcards):
    return dom[dom.domain == wildcards.domain].area

def list_cordex_in(wildcards):
    return expand("results/{base}/raw/{area}_{base}.nc", 
                   area=dom[dom.domain == wildcards.domain].area, allow_missing=True)

rule get_cordex:
    input:
        list_cordex_in
    output:
        directory("results/projection/raw/CORDEX_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}/")
    log:
        "results/logs/get_CORDEX_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}.log"
    benchmark:
        "results/benchmarks/get_cordex_CORDEX_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}_{base}.benchmark.txt"
    threads: 10
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        area=list_areas,
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
