rule get_anomalies:
    input:
        "results/cordex/mean/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_base}.nc",
        "results/cordex/mean/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_proj}.nc"
    output:
        "results/cordex/anomaly/{country}_{domain}_{gcm}_{rcm}_{rcp}_{period_proj}_{period_base}.nc"
    log:
        "results/logs/get_anomalies_{period_proj}_{period_base}_{country}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/get_anomalies_{period_proj}_{period_base}_{country}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        period_base="{period_base}",
        period_proj="{period_proj}"
    script:
      "../scripts/ge_anomalies.py"
