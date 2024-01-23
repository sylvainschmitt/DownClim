rule get_cmip6:
    input:
        expand("results/chelsa2/raw/{area}_chelsa2.nc", area=config["area"])
    output:
        directory("results/projection/raw/CMIP6_world_{institute}_{model}_{experiment}_{ensemble}_none_none_{base}/")
    log:
        "results/logs/get_cmip6_world_{institute}_{model}_{experiment}_{ensemble}_none_none_{base}.log"
    benchmark:
        "results/benchmarks/get_cmip6_world_{institute}_{model}_{experiment}_{ensemble}_none_none_{base}.benchmark.txt"
    threads: 10
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        area=config["area"],
        institute="{institute}",
        model="{model}",
        experiment="{experiment}",
        ensemble="{ensemble}",
        variables=config["variables"],
        time_frequency=config["time_frequency"],
        esgf_credential=config["esgf_credential"]
    script:
      "../scripts/get_cmip6.py"
