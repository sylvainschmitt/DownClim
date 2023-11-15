# IPSL: https://esgf-node.ipsl.upmc.fr/search/cordex-ipsl/
# CDS: https://cds.climate.copernicus.eu/cdsapp#!/dataset/projections-cordex-domains-single-levels?tab=overview
# 

rule retrieve_proj:
   input:
        "wget/{gcm}_{rcm}_{rcp}.sh"
    output:
        directory("results/cordex/{gcm}_{rcm}_{rcp}") # to later make as temporary
    log:
        "results/logs/retrieve_proj_{gcm}_{rcm}_{rcp}.log"
    benchmark:
        "results/benchmarks/retrieve_proj_{gcm}_{rcm}_{rcp}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
        script="{gcm}_{rcm}_{rcp}.sh"
    shell:
        "module load java ; " # only for muse, need to add a condition
        "ESG_HOME='.esg/' ; " # only for muse, need to add a condition
        "mkdir -p {output} ; "
        "cp {input} {output} ; "
        "cd {output} ; "
        "bash {params.script} -u ; "
        "bash {params.script}"
