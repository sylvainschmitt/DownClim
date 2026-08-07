from __future__ import annotations

from collections.abc import Hashable, Iterable
from enum import Enum
from pathlib import Path
from typing import Any

import geopandas as gpd
import xarray as xr
import xesmf as xe

from .aoi import get_aoi_informations
from .dataset.simulations import PeriodKind, SimulationCatalog, SimulationFile
from .dataset.utils import (
    Aggregation,
    DataProduct,
    _check_output_dir,
    check_input_dir,
    climatology_filename,
    get_regridder,
)
from .logging_config import get_logger

logger = get_logger(__name__)


class DownscaleMethod(Enum):
    """Class to define the downscaling methods available."""

    BIAS_CORRECTION = "bias_correction"
    QUANTILE_MAPPING = "quantile_mapping"
    DYNAMICAL = "dynamical"

    @classmethod
    def _missing_(cls, value: Any) -> None:
        msg = f"Unknown or not implemented downscaling method '{value}'. Right now only 'bias_correction' is implemented."
        raise ValueError(msg)


def bias_correction(
    baseline: xr.Dataset, historical: xr.Dataset, projection: xr.Dataset
) -> xr.Dataset:
    """Bias correction of the projections. All datasets must be on the same grid.

    Args:
        baseline (xr.Dataset): Baseline product dataset.
        historical (xr.Dataset): Historical period dataset.
        projection (xr.Dataset): Projection period dataset.

    Returns:
        xr.Dataset: Bias corrected projections.
    """

    # Compute anomalies between projection and historical
    anomalies = projection - historical
    # Add anomalies to the baseline
    projection = baseline + anomalies

    if "pr" in projection:
        rel_anomalies_pr = anomalies["pr"] / (historical["pr"] + 1)
        projection["pr"] = baseline["pr"] * (1 + rel_anomalies_pr)

    return projection


def _check_populate_simulations(
    simulations: list[str] | None, aoi_n: str, input_dir: str, dataproduct: DataProduct
) -> list[str]:
    """Check and populate simulations for a specific AOI and data product.

    Args:
        simulations (list[str] | None): List of simulations to downscale.
        aoi_n (str): AOI name.
        input_dir (str): Input directory where the simulation files are located.
        dataproduct (DataProduct): Data product.

    Returns:
        list[str]: List of populated simulations for the AOI and data product.
    """
    if simulations is None:
        simulations = [
            str(p)
            for p in Path(f"{input_dir}/{dataproduct.product_name}").glob(
                f"{aoi_n}_{dataproduct.product_name}*.nc"
            )
        ]
        msg = f"{dataproduct.product_name.upper()} simulations to downscale not provided. Using all files found in {input_dir}/{dataproduct.product_name}."
        logger.warning(msg)
    if simulations == []:
        msg = f"No {dataproduct.product_name.upper()} simulations to downscale found."
        logger.warning(msg)
    return simulations


def _populate_simulations(
    cmip6_simulations_to_downscale: list[str] | None,
    cordex_simulations_to_downscale: list[str] | None,
    aoi_n: str,
    input_dir: str,
) -> tuple[list[str], list[str]]:
    """
    Populate CMIP6 and CORDEX simulations to downscale for a given AOI.

    Args:
        cmip6_simulations_to_downscale (list[str] | None): List of CMIP6 simulations to downscale. Defaults to None,
        cordex_simulations_to_downscale (list[str] | None): List of CORDEX simulations to downscale. Defaults to None,
        aoi_n (str): AOI name.
        input_dir (str): Input directory where the simulation files are located.

    Returns:
        tuple[list[str], list[str]]: Tuple of populated CMIP6 and CORDEX simulations to downscale for the AOI.
    """

    logger.info("   Checking simulations to downscale for AOI: %s", aoi_n)
    if not cmip6_simulations_to_downscale:
        cmip6_simulations_to_downscale = _check_populate_simulations(
            cmip6_simulations_to_downscale, aoi_n, input_dir, DataProduct.CMIP6
        )
    if not cordex_simulations_to_downscale:
        cordex_simulations_to_downscale = _check_populate_simulations(
            cordex_simulations_to_downscale, aoi_n, input_dir, DataProduct.CORDEX
        )
    return cmip6_simulations_to_downscale, cordex_simulations_to_downscale


def _get_downscaling_grid_file(
    downscaling_grid_file: list[str] | None,
    aoi_n: str,
    input_dir: str,
    baseline_product: DataProduct,
    i: int,
) -> Path:
    """
    Get the downscaling grid path for a given AOI.

    Args:
        downscaling_grid_file (list[str] | None): Path to the grid file on which to downscale. One for each aoi. Defaults to None, meaning
        the grid will be extracted from the baseline product for each aoi.
        aoi_n (str): AOI name.
        input_dir (str): Input directory where the simulation files are located.
        baseline_product (DataProduct): Baseline product to use.
        i (int): Index of the grid file to retrieve in the list.

    Returns:
        Path: Path to the downscaling grid file.
    """
    logger.info("       Checking downscaling grid file...")
    if downscaling_grid_file is None:
        downscaling_grid_file_aoi = f"{input_dir}/{baseline_product.product_name}/{baseline_product.product_name}_{aoi_n}_grid.nc"
        msg = f"Downscaling grid file not provided. Using default grid file {downscaling_grid_file_aoi} which is extracted from {baseline_product.product_name}"
        logger.warning(msg)
    else:
        downscaling_grid_file_aoi = downscaling_grid_file[i]
    if not Path(downscaling_grid_file_aoi).is_file():
        msg = f"Downscaling grid file {downscaling_grid_file_aoi} not found. Please provide a valid downscaling grid file."
        logger.error(msg)
        raise FileNotFoundError(msg)
    return Path(downscaling_grid_file_aoi)


def _get_baseline_historical_data(
    baseline_product: DataProduct,
    aoi_n: str,
    input_dir: str,
    aggregation: Aggregation,
    historical_period: tuple[int, int],
) -> xr.Dataset:
    """
    Get the baseline historical data for a given AOI.

    Args:
        baseline_product (DataProduct): Baseline product to use.
        aoi_n (str): AOI name.
        input_dir (str): Input directory where the baseline files are located.
        aggregation (Aggregation): Aggregation method to use.
        historical_period (tuple[int, int]): Historical period (start, end).

    Returns:
        xr.Dataset: Baseline historical data.
    """
    logger.info("       Checking baseline historical data...")
    baseline_file = climatology_filename(
        f"{input_dir}/{baseline_product.product_name}",
        aoi_n,
        baseline_product,
        aggregation,
        historical_period,
    )
    if not Path(baseline_file).is_file():
        msg = f"""Baseline historical data not found: {baseline_product.product_name} for {aoi_n} should be located in {baseline_file}.
        Please download it first by using the `downclim.downclim.DownClimContext.download_data` method."""
        raise FileNotFoundError(msg)
    return xr.open_dataset(baseline_file)


def _regrid_baseline_data(
    ds_baseline: xr.Dataset,
    downscaling_grid: xr.Dataset,
    baseline_product: DataProduct,
    downscaling_grid_file_aoi: str,
    aoi_n: str,
    input_dir: str,
    output_dir: str,
) -> xr.Dataset:
    """
    Regrid baseline data on the downscaling grid. Nothing is done is baseline grid is the downscaling grid.

    Args:
        ds_baseline (xr.Dataset): Baseline data.
        downscaling_grid (xr.Dataset): Downscaling grid that will receive baseline data.
        baseline_product (DataProduct): Path to the baseline grid file.
        downscaling_grid_file_aoi (str): Path to the downscaling grid file.
        aoi_n (str): AOI name.
        input_dir (str): Input directory where the baseline files are located.
        output_dir (str): Output directory where the downscaled files will be saved.

    Returns:
        xr.Dataset: Regridded baseline data.
    """
    baseline_grid_file = f"{input_dir}/{baseline_product.product_name}/{baseline_product.product_name}_{aoi_n}_grid.nc"
    baseline_grid = xr.open_dataset(baseline_grid_file)
    if baseline_grid.equals(downscaling_grid):
        logger.info("       Baseline grid matches downscaling grid")
        return ds_baseline
    logger.info("       Regridding baseline data on the downscaling grid.")
    regridder = get_regridder(
        ds_baseline,
        downscaling_grid,
        baseline_grid_file,
        downscaling_grid_file_aoi,
        f"{output_dir}/..",
    )
    return regridder(ds_baseline, keep_attrs=True)  # type: ignore[no-any-return]


def _netcdf_encoding(ds: xr.Dataset) -> dict[Hashable, dict[str, Any]]:
    """Compressed float32 encoding for all the data variables of a dataset.

    zlib is lossless; float32 is the standard precision for stored climate data.
    """
    return {
        var: {"zlib": True, "complevel": 5, "dtype": "float32"} for var in ds.data_vars
    }


def _regrid_historical_data(
    historical_file: str | Path,
    historical_period: tuple[int, int],
    downscaling_grid: xr.Dataset,
    downscaling_grid_file: Path,
    aoi_n: str,
    output_dir: str,
    data_product: str,
) -> xr.Dataset:
    """
    Regrid historical data onto the downscaling grid and save to netcdf

    Args:
        historical_file (str | Path): Path to the historical file.
        historical_period (tuple[int, int]): Historical period.
        downscaling_grid (xr.Dataset): Downscaling grid.
        downscaling_grid_file (Path): Path to the downscaling grid file.
        aoi_n (str): AOI name.
        output_dir (str): Output directory.
        data_product (str): Name of the data product.

    Returns:
        xr.Dataset: Regridded historical data.
    """
    # Open historical datasets and interpolate the data onto downscaling grid
    logger.info(
        "       Regridding historical data %s, period %s, for AOI: %s.",
        historical_file,
        historical_period,
        aoi_n,
    )
    historical_regridded_file = f"{output_dir}/{data_product}/{Path(historical_file).stem}-{Path(downscaling_grid_file).stem}.nc"
    if Path(historical_regridded_file).is_file():
        logger.warning(
            "        Regridded historical dataset for %s already exists: %s. No action taken.",
            historical_file,
            historical_regridded_file,
        )
        ds_historical_regridded = xr.open_dataset(historical_regridded_file)
    else:
        logger.info(
            "       Regridding historical dataset for %s: %s.",
            historical_file,
            historical_regridded_file,
        )
        ds_historical = xr.open_dataset(historical_file)
        regridder = xe.Regridder(ds_historical, downscaling_grid, "bilinear")
        ds_historical_regridded = regridder(ds_historical, keep_attrs=True)
        ds_historical_regridded.to_netcdf(
            historical_regridded_file,
            encoding=_netcdf_encoding(ds_historical_regridded),
        )
    return ds_historical_regridded


def _downscale_period(
    period: str,
    simulations_product: DataProduct,
    baseline_product: DataProduct,
    downscaling_grid_file: str,
    files_to_downscale: list[SimulationFile],
    aoi_n: str,
    ds_baseline_regridded: xr.Dataset,
    ds_historical_regridded: xr.Dataset,
    regridder: xe.Regridder,
    output_dir: str,
    method: DownscaleMethod,
) -> None:
    """
    Downscales CMIP6 / CORDEX files found for a dedicated set of requirements and a dedicated period.

    Args:
        period (str): Period to downscale.
        simulations_product (DataProduct): type of simulations to downscale, CMIP6 or CORDEX
        baseline_product (DataProduct): baseline product
        downscaling_grid_file (str): path to the downscaling grid
        files_to_downscale (list[SimulationFile]): Simulations to downscale.
        aoi_n (str): AOI name.
        ds_baseline_regridded (xr.Dataset): Baseline product on downscaling grid.
        ds_historical_regridded (xr.Dataset): Regridded historical data onto downscaling grid.
        regridder (xe.Regridder): Regridder to use from original grid to downscaling grid.
        output_dir (str): Output directory.
        method (DownscaleMethod): Downscaling method to use.

    """
    logger.info(
        "       Regridding dataset for %s, period: %s, for AOI: %s.",
        [simulation_file.path for simulation_file in files_to_downscale],
        period,
        aoi_n,
    )
    for simulation_file in files_to_downscale:
        downscaled_file = f"{output_dir}/{simulations_product.product_name}/{Path(simulation_file.path).stem}-downscaled-{baseline_product.product_name}_baseline-{Path(downscaling_grid_file).stem}.nc"

        if Path(downscaled_file).is_file():
            logger.warning(
                "        Downscaled dataset for %s already exists: %s. No action taken.",
                simulation_file.path,
                downscaled_file,
            )
            continue
        ds_to_downscale = xr.open_dataset(simulation_file.path)
        ds_to_downscale_regridded = regridder(ds_to_downscale, keep_attrs=True)

        # Downscale
        logger.info(
            "       Downscaling dataset %s, period: %s, for AOI: %s using method: %s.",
            simulation_file.path,
            period,
            aoi_n,
            method.value,
        )
        if method == DownscaleMethod.BIAS_CORRECTION:
            ds_downscaled = bias_correction(
                ds_baseline_regridded,
                ds_historical_regridded,
                ds_to_downscale_regridded,
            )

        # Save downscaled dataset
        logger.info("       Saving downscaled dataset into: %s", downscaled_file)
        ds_downscaled.attrs.update(simulation_file.to_attrs())
        ds_downscaled.to_netcdf(
            downscaled_file, encoding=_netcdf_encoding(ds_downscaled)
        )


def run_downscaling(
    aoi: list[gpd.GeoDataFrame],
    historical_period: tuple[int, int],
    evaluation_period: tuple[int, int],
    projection_period: tuple[int, int],
    baseline_product: DataProduct,
    cmip6_simulations_to_downscale: list[str] | None = None,
    cordex_simulations_to_downscale: list[str] | None = None,
    downscaling_grid_file: list[str] | None = None,
    periods_to_downscale: Iterable[str] | None = None,
    aggregation: Aggregation = Aggregation.MONTHLY_MEAN,  # type: ignore[assignment]
    method: DownscaleMethod = DownscaleMethod.BIAS_CORRECTION,
    input_dir: str | None = None,
    output_dir: str | None = None,
) -> None:
    """
    Run the downscaling process.

    Args:
        aoi (list[gpd.GeoDataFrame]): List of areas of interest.
        historical_period (tuple[int, int]): Baseline period (start, end).
        evaluation_period (tuple[int, int]): Evaluation period (start, end).
        projection_period (tuple[int, int]): Projection period (start, end).
        baseline_product (DataProduct): Baseline product to use.
        cmip6_simulations_to_downscale (list[str] | None): List of CMIP6 simulations to downscale. Defaults to None,
        which means all available CMIP6 simulations in "<input_dir>".
        cordex_simulations_to_downscale (list[str] | None): List of CORDEX simulations to downscale. Defaults to None,
        which means all available CORDEX simulations in "<input_dir>".
        downscaling_grid_file (list[str] | None): Path to the grid file on which to downscale. One for each aoi. Defaults to None, meaning
        the grid will be extracted from the baseline product for each aoi.
        periods_to_downscale (list[str] | None): List of periods to downscale. Can be any combination of ['evaluation', 'projection']. Defaults to None, meaning all periods will be downscaled.
        aggregation (Aggregation, optional): Aggregation method to use. Defaults to Aggregation.MONTHLY_MEAN.
        method (DownscaleMethod, optional): Downscaling method to use. Defaults to DownscaleMethod.BIAS_CORRECTION.
    input_dir (str, optional): Input directory for the data. Only used if "<cmip6_simulations_to_downscale>" or "<cordex_simulations_to_downscale>" are None. Defaults to "./results".
        output_dir (str, optional): Output directory for the results. Defaults to "./results/downscaled".

    Raises:
        FileNotFoundError: If a required file is not found.
        ValueError: If a required parameter is invalid.
    """

    logger.info("Starting downscaling process...")
    # Check input directory
    input_dir = check_input_dir(input_dir, "./results")

    # Create output directory
    output_dir = _check_output_dir(
        output_dir, "./results/downscaled", ["cmip6", "cordex", "../regridder"]
    )

    # Check the downscaling method
    if method != DownscaleMethod.BIAS_CORRECTION:
        msg = "Method not implemented yet, only bias_correction is available."
        raise ValueError(msg)

    # Get AOIs information
    aoi_name, _ = get_aoi_informations(aoi)

    # Define periods to downscale
    logger.info("Checking periods to downscale...")
    if periods_to_downscale is None:
        periods_to_downscale = ["evaluation", "projection"]
        logger.warning(
            "Periods to downscale not provided. Using default periods %s.",
            periods_to_downscale,
        )
    if not all(
        period in ["evaluation", "projection"] for period in periods_to_downscale
    ):
        msg = f"Invalid periods found in {periods_to_downscale}. Please provide valid periods : ['evaluation', 'projection']."
        raise ValueError(msg)

    for i, aoi_n in enumerate(aoi_name):
        # Populate simulations to downscale
        cmip6_simulations_to_downscale, cordex_simulations_to_downscale = (
            _populate_simulations(
                cmip6_simulations_to_downscale,
                cordex_simulations_to_downscale,
                aoi_n,
                input_dir,
            )
        )

        # Get the downscaling grid
        downscaling_grid_path = _get_downscaling_grid_file(
            downscaling_grid_file, aoi_n, input_dir, baseline_product, i
        )
        downscaling_grid = xr.open_dataset(downscaling_grid_path)

        # Get baseline historical data
        ds_baseline = _get_baseline_historical_data(
            baseline_product, aoi_n, input_dir, aggregation, historical_period
        )

        # Interpolate baseline data on downscaling grid (if needed)
        ds_baseline_regridded = _regrid_baseline_data(
            ds_baseline,
            downscaling_grid,
            baseline_product,
            str(downscaling_grid_path),
            aoi_n,
            input_dir,
            output_dir,
        )

        # Build the catalog of simulations to downscale
        catalog = SimulationCatalog.from_files(
            [*cmip6_simulations_to_downscale, *cordex_simulations_to_downscale],
            {
                PeriodKind.HISTORICAL: historical_period,
                PeriodKind.EVALUATION: evaluation_period,
                PeriodKind.PROJECTION: projection_period,
            },
        )
        logger.info(
            "   Simulations to downscale for AOI %s: %s",
            aoi_n,
            catalog.simulations(),
        )

        for simulation in catalog.simulations():
            for historical_file in catalog.historical_files(simulation):
                ds_historical = xr.open_dataset(historical_file.path)
                regridder = get_regridder(
                    ds_historical,
                    downscaling_grid,
                    str(historical_file.path),
                    str(downscaling_grid_path),
                    f"{output_dir}/../regridder",
                )
                ds_historical_regridded = _regrid_historical_data(
                    historical_file.path,
                    historical_period,
                    downscaling_grid,
                    downscaling_grid_path,
                    aoi_n,
                    output_dir,
                    simulation.product.product_name,
                )

                for period in periods_to_downscale:
                    # Check existing files matching the context to downscale
                    files_to_downscale = catalog.matching_files(
                        simulation, PeriodKind(period)
                    )
                    _downscale_period(
                        period,
                        simulation.product,
                        baseline_product,
                        str(downscaling_grid_path),
                        files_to_downscale,
                        aoi_n,
                        ds_baseline_regridded,
                        ds_historical_regridded,
                        regridder,
                        output_dir,
                        method,
                    )
