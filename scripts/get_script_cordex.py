# esgf-pylient https://esgf-pyclient.readthedocs.io/en/latest/notebooks/demo/subset-cmip5.html

# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
domain = snakemake.params.domain
gcm = snakemake.params.gcm
rcm = snakemake.params.rcm
rcp = snakemake.params.rcp
var = snakemake.params.var
out_file = snakemake.output[0]

# test
# domain = "AUS-22"
# gcm = "MOHC-HadGEM2-ES"
# rcm = "CCLM5-0-15"
# rcp = "historical"
# var = "tas"

# code
from pyesgf.search import SearchConnection

server = 'https://esgf-node.ipsl.upmc.fr/esg-search/'
conn = SearchConnection(server, distrib=True)

ctx = conn.new_context(
    facets='*',
    project = 'CORDEX',
    domain = domain,
    driving_model = gcm,
    rcm_name= rcm,
    time_frequency = 'mon',
    experiment = rcp,
    variable = var
    )
if(ctx.hit_count == 0):
  raise SystemExit('The query has no results')
# ctx.search()[0].dataset_id
f = open(out_file, "w")
f.write(ctx.get_download_script())
