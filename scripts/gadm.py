from gadm import GADMDownloader
downloader = GADMDownloader(version="4.0")
country_name = "New Caledonia"
ad_level = 0
gdf = downloader.get_shape_data_by_country_name(country_name=country_name, ad_level=ad_level)
