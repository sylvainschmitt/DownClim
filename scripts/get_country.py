# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
country = snakemake.params.country
log10_eval_pts = snakemake.params.log10_eval_pts
country_file = snakemake.output[0]
country_fig = snakemake.output[1]
pts_file = snakemake.output[2]
pts_fig = snakemake.output[3]

# test
# country = "New-Caledonia"
# log10_eval_pts = 4

# libs
# from gadm import GADMDownloader
import pygadm
import matplotlib.pyplot as plt
import re

# code
country = re.sub("-", " ", country)
# downloader = GADMDownloader(version="4.0")
# gdf = downloader.get_shape_data_by_country_name(country_name=country, ad_level=0)
code = pygadm.AdmNames(country).GID_0[0]
gdf = pygadm.AdmItems(admin = code)
pts = gdf.sample_points(pow(10, log10_eval_pts))

# country
gdf.plot()
plt.savefig(country_fig)
gdf.to_file(country_file)

# points
pts.plot()
plt.savefig(pts_fig)
pts.to_file(pts_file)
