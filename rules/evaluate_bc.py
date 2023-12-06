rule evaluate_bc:
    input:
        proj=expand("results/cordex/raw/{country}_{domain}_{gcm}_{rcm}_{rcp}_{years}.nc",
                     years=config["cordex_years"], allow_missing = True),
        base="results/chelsa/raw/{country}_CHELSA_1981-2010_V.2.1.nc",
        ds="results/cordex/downscaled/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_eval}_{period_base}.nc",
        pts="results/countries/{country}_pts.shp"
    output:
        "results/evaluation/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_eval}_{period_base}.tsv"
    log:
        "results/logs/evaluate_bc_{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_eval}_{period_base}.log"
    benchmark:
        "results/benchmarks/evaluate_bc_{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_eval}_{period_base}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        period_base="{period_base}",
        period_eval="{period_eval}"
    script:
      "../scripts/evaluate_bc.py"
