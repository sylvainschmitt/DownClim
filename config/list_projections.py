# libs
from pyesgf.search import SearchConnection
import gcsfs
import pandas as pd

# cordex
project = "CORDEX"
activity = "none"
domains = ["SAM-22", "AFR-22", "AUS-22"]
experiments = ["rcp26", "rcp85"]
variables = ["tas", "tasmin", "tasmax", "pr"]
time_frequency = "mon"
server = 'https://esgf-node.ipsl.upmc.fr/esg-search/'
conn = SearchConnection(server, distrib=True)
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
df = pd.DataFrame(dict(dataset=datasets))
df[["dataset", "datanode"]]  = df['dataset'].str.split('|', expand=True)
df[['project', 'product', 'domain', 'institute', 'model', 'experiment',
    'ensemble', 'rcm', 'downscaling', 'time_frequency', 'variable', 'version']] = df['dataset'].str.split('.', expand=True)
df.project = df.project.str.upper()
cordex = df[['project', 'domain', 'institute', 'model', 'experiment', 'ensemble', 'rcm', 'downscaling']].drop_duplicates()

# cmip6
experiments = ["ssp126", "ssp585"]
variables =  ["tas", "tasmin", "tasmax", "pr"]
gcs = gcsfs.GCSFileSystem(token='anon')
df = pd.read_csv('https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv')
df_ta = df.query("activity_id == 'ScenarioMIP' & table_id == 'Amon' & variable_id == @variables & experiment_id == @experiments")
df_ta = df_ta.rename(columns = {"institution_id": "institute", "source_id": "model", "experiment_id": "experiment", "member_id": "ensemble"})
df_ta.insert(0, "project", "CMIP6")
df_ta.insert(0, "domain", "world")
df_ta.insert(0, "rcm", "none")
df_ta.insert(0, "downscaling", "none")
cmip6 = df_ta[['project', 'domain', 'institute', 'model', 'experiment', 'ensemble', 'rcm', 'downscaling']].drop_duplicates().groupby(["institute", "model", 'experiment',]).tail(1)

# save
pd.concat([cordex, cmip6]).to_csv("config/projections_all.tsv", sep="\t", index=False)
