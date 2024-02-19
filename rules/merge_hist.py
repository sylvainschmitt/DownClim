rule merge_hist:
    input:
      ds = expand("results/evaluation/hist/ds/{proj}_{baseline}_{aggregation}_{period_eval}_{period_base}_{ds_method}.tsv",
                      proj=proj_dom.id,
                      baseline=config["baseline"],
                      aggregation=config["aggregation"],
                      period_base=config["hist_years"],
                      period_eval=config["eval_years"],
                      base_eval=config["base_eval"],
                      ds_method=config["ds_method"]),
      raw = expand("results/evaluation/hist/raw/{proj}_{baseline}_{aggregation}_{period_eval}.tsv",
                      proj=proj_dom.id,
                      baseline=config["baseline"],
                      aggregation=config["aggregation"],
                      period_eval=config["eval_years"]),
      base = expand("results/evaluation/hist/base/{area}_{base}_{aggregation}_{period_eval}.tsv",
                      area=config["area"],
                      base=config["baseline"],
                      aggregation=config["aggregation"],
                      period_eval=config["eval_years"])
    output:
      "results/evaluation/histograms.tsv"
    log:
      "results/logs/merge_hist.log"
    benchmark:
      "results/benchmarks/merge_hist.benchmark.txt"
    script:
      "../scripts/eval_hist.py"
