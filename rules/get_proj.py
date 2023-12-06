rule get_proj:
    input:
        "results/countries/{area}.shp"
    output:
        "results/projection/raw/{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}.nc"
    log:
        "results/logs/get_cordex_{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}.log"
    benchmark:
        "results/benchmarks/get_cordex_{area}_{project}_{activity}_{domain}_{institute}_{model}_{experiment}_{ensemble}_{rcm}_{downscaling}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        area="{area}",
        project="{project}",
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
        proj_years=config["esgf_years"]
    shell:
      "../scripts/get_proj.py"
