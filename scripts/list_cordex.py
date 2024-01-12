# test
project = "CORDEX"
activity = "none"
domains = ["SAM-22", "AFR-22", "AUS-22"]
experiments = ["rcp26", "rcp85"]
variables = ["tas", "tasmin", "tasmax", "pr"]
time_frequency = "mon"
file = "test.tsv"

# libs
from pyesgf.search import SearchConnection
import pandas as pd

# list
server = 'https://esgf-node.ipsl.upmc.fr/esg-search/'
conn = SearchConnection(server, distrib=True)
# list of keys: https://esgf-node.llnl.gov/esg-search/search?project=CORDEX&facets=*&limit=0
ctx = conn.new_context(
   facets='*',
   project = project,
   domain = domains,
   experiment = experiments,
   time_frequency = time_frequency,
   variable = variables
  )  
results = ctx.search()
datasets = [res.dataset_id for res in results]
pd.DataFrame(dict(dataset=datasets)).to_csv(file, sep="\t")
