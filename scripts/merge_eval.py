# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
ds =  snakemake.input.ds
raw =  snakemake.input.raw
out_file = snakemake.output[0]

# libs
import pandas as pd 

# code
ds_all = pd.concat(list(map(lambda t: pd.read_csv(t, sep="\t"), ds)))
raw_all = pd.concat(list(map(lambda t: pd.read_csv(t, sep="\t"), raw)))
all = pd.concat([ds_all, raw_all])
all.to_csv(out_file, sep="\t", index=False)
