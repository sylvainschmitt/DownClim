# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area = snakemake.params.area
log10_eval_pts = snakemake.params.log10_eval_pts
area_file = snakemake.output[0]
area_fig = snakemake.output[1]

# test
# area = "Côte-d'ivoire"
# log10_eval_pts = 4

# libs
import pygadm
import matplotlib.pyplot as plt
import re

# code
area = re.sub("-", " ", area)
code = pygadm.AdmNames(area).GID_0[0]
gdf = pygadm.AdmItems(admin = code)

# country
gdf.plot()
plt.savefig(area_fig)
gdf.to_file(area_file)
