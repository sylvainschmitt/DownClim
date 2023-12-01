rule downscale:
    input:
        "results/cordex/{type}_bias/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc",
        "results/chelsa/hist/{var}_{country}_{month}.tif"
    output:
        "results/downscaled/{type}_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}_{month}.tif"
    log:
        "results/logs/downscaled_{type}_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}_{month}.log"
    benchmark:
        "results/benchmarks/downscaled_{type}_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}_{month}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        var="{var}",
        month="{month}"
    script:
      "../scripts/downscale_rast.R"
