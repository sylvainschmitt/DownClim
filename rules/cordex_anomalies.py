rule cordex_anomalies:
    input:
        "results/cordex/{type}/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc",
        "results/cordex/hist/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    output:
        "results/cordex/{type}_bias/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    log:
        "results/logs/cordex_anomalies_{type}_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/cordex_anomalies_{type}_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        var="{var}"
    script:
      "../scripts/compute_anomalies_nc.R"
