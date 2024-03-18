rule merge_eval:
    input:
      ds = expand("results/evaluation/eval/ds/{proj}_{baseline}_{aggregation}_{period_eval}_{period_base}_{ds_method}_{base_eval}.tsv",
                      proj=proj_dom.id,
                      baseline=config["baseline"],
                      aggregation=config["aggregation"],
                      period_base=config["hist_years"],
                      period_eval=config["eval_years"],
                      ds_method=config["ds_method"],
                      base_eval=config["base_eval"]),
      raw = expand("results/evaluation/eval/raw/{proj}_{baseline}_{aggregation}_{period_eval}_{base_eval}.tsv",
                      proj=proj_dom.id,
                      baseline=config["baseline"],
                      aggregation=config["aggregation"],
                      period_eval=config["eval_years"],
                      base_eval=config["base_eval"])
    output:
      "results/evaluation/evaluations.tsv"
    log:
      "results/logs/merge_eval.log"
    benchmark:
      "results/benchmarks/merge_eval.benchmark.txt"
    script:
      "../scripts/merge_eval.py"
