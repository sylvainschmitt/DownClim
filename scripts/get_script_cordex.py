# esgf-pylient https://esgf-pyclient.readthedocs.io/en/latest/notebooks/demo/subset-cmip5.html

log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

from pyesgf.search import SearchConnection

server = 'https://esgf-node.ipsl.upmc.fr/esg-search/'
conn = SearchConnection(server, distrib=True)

ctx = conn.new_context(
    project = 'CORDEX',
    domain = snakemake.params.domain,
    driving_model = snakemake.params.gcm,
    rcm_name= snakemake.params.rcm,
    time_frequency = 'mon',
    variable = 'tas',
    experiment = snakemake.params.rcp
    )
if(ctx.hit_count == 0):
  raise SystemExit('The query has no results')
# ctx.search()[0].dataset_id
f = open( snakemake.output[0], "w")
f.write(ctx.get_download_script())
