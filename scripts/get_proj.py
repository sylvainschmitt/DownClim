# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area_file = snakemake.intput[0]
domain = snakemake.params.domain
area = snakemake.params.area
project = snakemake.params.project
activity = snakemake.params.activity
domain = snakemake.params.domain
institute = snakemake.params.institute
model = snakemake.params.model
experiment = snakemake.params.experiment
ensemble = snakemake.params.ensemble
rcm = snakemake.params.rcm
downscaling = snakemake.params.downscaling
variables = snakemake.params.variables
time_frequency = snakemake.params.time_frequency
proj_years = snakemake.params.proj_years
nc_file = snakemake.output[0]

# test
raise Exception("Currently in development!")

# code
from pyesgf.search import SearchConnection
# esgf-pylient https://esgf-pyclient.readthedocs.io/en/latest/notebooks/demo/subset-cmip5.html

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
f.write(ctx.get_download_script()) # downloading wget script must be replaced by direct reading of data online
