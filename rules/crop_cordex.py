rule crop_cordex:
    input:
        "results/cordex/merged/{gcm}_{rcm}_{rcp}.nc"
    output:
        "results/cordex/cropped/{gcm}_{rcm}_{rcp}.nc"
    log:
        "results/logs/crop_cordex_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/crop_cordex_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    params:
        lonmin=config["lonmin"],
        lonmax=config["lonmax"],
        latmin=config["latmin"],
        latmax=config["latmax"]
    resources:
        mem_mb=1000
    shell:
      "cdo -sellonlatbox,{params.lonmin},{params.lonmax},{params.latmin},{params.latmax} {input[0]} {output[0]}"
