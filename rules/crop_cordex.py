rule crop_cordex:
    input:
        "results/countries/{country}_bb.txt",
        "results/cordex/projected/{var}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    output:
        "results/cordex/cropped/{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.nc"
    log:
        "results/logs/crop_cordex_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/crop_cordex_{var}_{country}_{domain}_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    shell:
      "source {input[0]} ; "
      "cdo -sellonlatbox,$xmin,$xmax,$ymin,$ymax {input[1]} {output[0]}"
