configfile: "config/config.yml"

rule all:
   input:
         expand("results/{key}",
                 key=config["key"])
                
# Rules #
include: "rules/rule.smk"
