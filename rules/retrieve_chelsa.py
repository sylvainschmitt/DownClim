rule retrieve_chelsa:
    output:
        "results/chelsa/raw/{var}_{month}_{year}.tif"
    log:
        "results/logs/retrieve_chelsa_{var}_{month}_{year}.log"
    benchmark:
        "results/benchmarks/retrieve_chelsa_{var}_{month}_{year}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        var="{var}",
        month="{month}",
        year="{year}"
    shell:
      "base_url='https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly' ; "
      "url=$base_url'/{params.var}/CHELSA_{params.var}_{params.month}_{params.year}_V.2.1.tif' ; "
      "wget $url -O {output[0]}"
