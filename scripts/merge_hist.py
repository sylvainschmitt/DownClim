# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
ds =  snakemake.input.ds
raw =  snakemake.input.raw
base =  snakemake.input.base
out_file = snakemake.output[0]

# libs
import pandas as pd 

# code
a = ds + raw + base
# remove 0 !
pd.concat(a).to_csv(out_file, sep="\t", index=False)
