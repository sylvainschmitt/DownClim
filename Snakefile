configfile: "config/config_test.yml"

rule all:
   input:
         expand("results/cordex/{gcm}_{rcm}_{rcp}",
                 gcm=config["gcm"],
                 rcm=config["rcm"],
                 rcp=config["rcp"])
                
# Rules #
include: "rules/retrieve_proj.py"
