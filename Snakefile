configfile: "config/config_test.yml"

rule all:
   input:
         expand("results/cordex/cropped/{gcm}_{rcm}_{rcp}.nc",
                 gcm=config["gcm"],
                 rcm=config["rcm"],
                 rcp=config["rcp"])
                
# Rules #

## CORDEX ##
include: "rules/get_script_cordex.py"
include: "rules/retrieve_cordex.py"
include: "rules/merge_cordex.py"
include: "rules/crop_cordex.py"

