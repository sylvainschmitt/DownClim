from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from downclim.dataset.utils import DataProduct
from downclim.evaluation import (
    _check_populate_simulations,
    _get_evaluation_grid_file,
    compute_correlation,
    compute_mean,
    compute_rmse,
    compute_std,
)


def _make_test_datasets() -> tuple[xr.Dataset, xr.Dataset]:
    time = xr.cftime_range("2000-01", periods=3, freq="ME", calendar="standard")
    ref = xr.Dataset({"tas": (("time", "lat", "lon"), np.ones((3, 2, 2)))})
    ref = ref.assign_coords(time=time)
    ev = xr.Dataset({"tas": (("time", "lat", "lon"), np.ones((3, 2, 2)) * 2)})
    ev = ev.assign_coords(time=time)
    return ref, ev


def test_compute_rmse():
    ref, ev = _make_test_datasets()
    result = compute_rmse(ref, ev)
    assert "tas" in result
    expected = ((ref["tas"] - ev["tas"]) ** 2).mean(dim=["lat", "lon"]) ** 0.5
    np.testing.assert_array_almost_equal(result["tas"].values, expected.values)


def test_compute_mean():
    ref, ev = _make_test_datasets()
    result = compute_mean(ref, ev)
    assert "tas" in result
    expected = (ev["tas"] - ref["tas"]).mean(dim=["lat", "lon"])
    np.testing.assert_array_almost_equal(result["tas"].values, expected.values)


def test_compute_std():
    ref, ev = _make_test_datasets()
    result = compute_std(ref, ev)
    assert "tas" in result


def test_compute_correlation():
    ref, ev = _make_test_datasets()
    result = compute_correlation(ref, ev)
    assert "tas" in result


@pytest.mark.network
def test_run_evaluation():
    pass  # requires full data pipeline


def test_get_evaluation_grid_file_none():
    grid_file = _get_evaluation_grid_file(
        None, DataProduct.CHELSA, "Vanuatu", "./results/downscaled", 0
    )
    assert grid_file == "./results/downscaled/../chelsa/chelsa_Vanuatu_grid.nc"


def test_get_evaluation_grid_file_dataproduct():
    grid_file = _get_evaluation_grid_file(
        DataProduct.CHELSA, DataProduct.CHIRPS, "Vanuatu", "./results/downscaled", 0
    )
    assert grid_file == "./results/downscaled/../chelsa/chelsa_Vanuatu_grid.nc"


def test_get_evaluation_grid_file_list():
    grids = ["./grid_nc.nc", "./grid_ic.nc"]
    assert (
        _get_evaluation_grid_file(grids, DataProduct.CHELSA, "Vanuatu", "in", 0)
        == "./grid_nc.nc"
    )
    assert (
        _get_evaluation_grid_file(grids, DataProduct.CHELSA, "Vanuatu", "in", 1)
        == "./grid_ic.nc"
    )


def test_check_populate_simulations(tmp_path):
    cmip6_dir = tmp_path / "cmip6"
    cmip6_dir.mkdir()
    matching = cmip6_dir / (
        "Vanuatu_cmip6_IPSL_IPSL-CM6A-LR_ssp585_r1i1p1f1_"
        "monthly-mean_2006-01-01_2019-12-31-downscaled-chelsa_baseline-"
        "chelsa_Vanuatu_grid.nc"
    )
    matching.touch()
    other = cmip6_dir / (
        "Vanuatu_cmip6_IPSL_IPSL-CM6A-LR_ssp585_r1i1p1f1_"
        "monthly-mean_2006-01-01_2019-12-31-downscaled-chelsa_baseline-"
        "chelsa_Vanuatu_grid_v2.nc"
    )
    other.touch()
    simulations = _check_populate_simulations(
        [str(matching), str(other)],
        "Vanuatu",
        str(tmp_path),
        DataProduct.CMIP6,
        (2006, 2019),
        "chelsa_Vanuatu_grid",
    )
    assert simulations == [str(matching)]
