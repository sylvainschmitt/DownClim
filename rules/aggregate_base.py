rule aggregate_base:
    input:
        "results/{base}/raw/{area}_{base}.nc"
    output:
        "results/{base}/means/{area}_{base}_{aggregation}_{period}.nc"
    log:
        "results/logs/aggregate_{base}_{area}_{aggregation}_{period}.log"
    benchmark:
        "results/benchmarks/aggregate_{base}_{area}_{aggregation}_{period}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    params:
        period="{period}",
        aggregation="{aggregation}"
    script:
      "../scripts/aggregate_base.py"
