"""Tests for the simulations data structure."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import numpy as np
import pytest
import xarray as xr
from pydantic import ValidationError

from downclim.dataset.cordex import _get_filename_from_cordex_context
from downclim.dataset.simulations import (
    PeriodKind,
    Simulation,
    SimulationCatalog,
    SimulationFile,
)
from downclim.dataset.utils import Aggregation, DataProduct


def _write_file(
    tmp_path: Path, name: str, simulation: Simulation, period: tuple[int, int]
) -> Path:
    """Write a minimal netCDF file carrying the simulation metadata as attrs."""
    path = tmp_path / name
    ds = xr.Dataset({"tas": (("time",), np.ones(12))})
    ds.attrs.update(
        SimulationFile(simulation=simulation, period=period, path=path).to_attrs()
    )
    ds.to_netcdf(path, engine="scipy")
    return path


def _cmip6_sim(experiment: str = "historical") -> Simulation:
    return Simulation(
        product=DataProduct.CMIP6,
        aoi_n="Vanuatu",
        institute="IPSL",
        source="IPSL-CM6A-LR",
        ensemble="r1i1p1f1",
        experiment=experiment,
    )


class TestSimulation:
    def test_required_fields(self):
        sim = _cmip6_sim()
        assert sim.product is DataProduct.CMIP6
        assert sim.aoi_n == "Vanuatu"
        assert sim.aggregation is Aggregation.MONTHLY_MEAN

    def test_extra_fields_forbidden(self):
        with pytest.raises(ValidationError):
            Simulation(product=DataProduct.CMIP6, aoi_n="Vanuatu", unknown="field")  # type: ignore[call-arg]

    def test_matching_key_ignores_experiment(self):
        historical = _cmip6_sim(experiment="historical")
        projection = _cmip6_sim(experiment="ssp245")
        assert historical.matching_key() == projection.matching_key()

    def test_matching_key_differs_between_products(self):
        cordex = Simulation(
            product=DataProduct.CORDEX,
            aoi_n="Vanuatu",
            domain="EUR-11",
            driving_model="IPSL-CM6A-LR",
            rcm_name="ALADIN63",
            ensemble="r1i1p1",
            rcm_version="v1",
        )
        assert _cmip6_sim().matching_key() != cordex.matching_key()


class TestFilenameParsing:
    def test_cmip6_from_filename(self):
        sim = Simulation.from_filename(
            "results/Vanuatu_cmip6_IPSL_IPSL-CM6A-LR_historical_r1i1p1f1_monthly-mean_1980-01-01_2005-12-31.nc"
        )
        assert sim.product is DataProduct.CMIP6
        assert sim.aoi_n == "Vanuatu"
        assert sim.institute == "IPSL"
        assert sim.source == "IPSL-CM6A-LR"
        assert sim.experiment == "historical"
        assert sim.ensemble == "r1i1p1f1"

    def test_cordex_from_filename(self):
        sim = Simulation.from_filename(
            "results/Vanuatu_cordex_EUR-11_IPSL-CM6A-LR_ALADIN63_r1i1p1_v1_monthly-mean_1980-01-01_2005-12-31.nc"
        )
        assert sim.product is DataProduct.CORDEX
        assert sim.domain == "EUR-11"
        assert sim.driving_model == "IPSL-CM6A-LR"
        assert sim.rcm_name == "ALADIN63"
        assert sim.ensemble == "r1i1p1"
        assert sim.rcm_version == "v1"

    def test_from_filename_invalid(self):
        with pytest.raises(
            ValueError, match=r"(Unknown data product|Cannot parse simulation context)"
        ):
            Simulation.from_filename("not_a_valid_filename.nc")

    def test_cordex_filename_builder_roundtrip(self, tmp_path):
        """The CORDEX filename must be parseable back into the same context."""
        output_file = _get_filename_from_cordex_context(
            str(tmp_path),
            "Vanuatu",
            DataProduct.CORDEX,
            "EUR-11",
            "IPSL-CM6A-LR",
            "ALADIN63",
            "r1i1p1",
            "v1",
            Aggregation.MONTHLY_MEAN,  # type: ignore[arg-type]
            "1980-01-01",
            "2005-12-31",
        )
        simulation = Simulation.from_filename(output_file)
        assert simulation.product is DataProduct.CORDEX
        assert simulation.domain == "EUR-11"
        assert simulation.rcm_version == "v1"


class TestAttrs:
    def test_attrs_roundtrip(self):
        sim = _cmip6_sim()
        assert Simulation.from_attrs(sim.to_attrs()) == sim

    def test_from_attrs_missing_returns_none(self):
        assert Simulation.from_attrs({"title": "something"}) is None

    def test_simulation_file_from_file_uses_attrs(self, tmp_path):
        sim = _cmip6_sim(experiment="historical")
        path = _write_file(tmp_path, "sim.nc", sim, (1980, 2005))
        simulation_file = SimulationFile.from_file(path)
        assert simulation_file.simulation == sim
        assert simulation_file.period == (1980, 2005)

    def test_simulation_file_from_file_falls_back_to_filename(self, tmp_path):
        name = (
            "Vanuatu_cmip6_IPSL_IPSL-CM6A-LR_historical_r1i1p1f1_"
            "monthly-mean_1980-01-01_2005-12-31.nc"
        )
        path = _write_file(tmp_path, name, _cmip6_sim(), (1980, 2005))
        # Remove the downclim attributes to force filename parsing
        with xr.open_dataset(path, engine="scipy") as ds:
            ds.load()
            ds.attrs = {}
            ds.to_netcdf(path, engine="scipy")
        simulation_file = SimulationFile.from_file(path)
        assert simulation_file.simulation.product is DataProduct.CMIP6
        assert simulation_file.period == (1980, 2005)


class TestCatalog:
    PERIODS: ClassVar[dict[PeriodKind, tuple[int, int]]] = {
        PeriodKind.HISTORICAL: (1980, 2005),
        PeriodKind.EVALUATION: (2006, 2019),
        PeriodKind.PROJECTION: (2071, 2100),
    }

    def test_matching_across_periods(self, tmp_path):
        historical = _write_file(
            tmp_path,
            "Vanuatu_cmip6_IPSL_IPSL-CM6A-LR_historical_r1i1p1f1_"
            "monthly-mean_1980-01-01_2005-12-31.nc",
            _cmip6_sim(experiment="historical"),
            (1980, 2005),
        )
        projection = _write_file(
            tmp_path,
            "Vanuatu_cmip6_IPSL_IPSL-CM6A-LR_ssp245_r1i1p1f1_"
            "monthly-mean_2071-01-01_2100-12-31.nc",
            _cmip6_sim(experiment="ssp245"),
            (2071, 2100),
        )
        catalog = SimulationCatalog.from_files([historical, projection], self.PERIODS)

        simulations = catalog.simulations()
        assert len(simulations) == 1
        simulation = simulations[0]

        assert [f.path for f in catalog.historical_files(simulation)] == [historical]
        assert [
            f.path for f in catalog.matching_files(simulation, PeriodKind.PROJECTION)
        ] == [projection]

    def test_matching_files_missing_raises(self, tmp_path):
        historical = _write_file(
            tmp_path,
            "Vanuatu_cmip6_IPSL_IPSL-CM6A-LR_historical_r1i1p1f1_"
            "monthly-mean_1980-01-01_2005-12-31.nc",
            _cmip6_sim(),
            (1980, 2005),
        )
        catalog = SimulationCatalog.from_files([historical], self.PERIODS)
        with pytest.raises(FileNotFoundError):
            catalog.matching_files(catalog.simulations()[0], PeriodKind.PROJECTION)

    def test_distinct_simulations(self, tmp_path):
        first = _write_file(
            tmp_path,
            "Vanuatu_cmip6_IPSL_IPSL-CM6A-LR_historical_r1i1p1f1_"
            "monthly-mean_1980-01-01_2005-12-31.nc",
            _cmip6_sim(),
            (1980, 2005),
        )
        second = Simulation(
            product=DataProduct.CMIP6,
            aoi_n="Vanuatu",
            institute="IPSL",
            source="IPSL-CM6A-LR",
            ensemble="r2i1p1f1",
            experiment="historical",
        )
        second_file = _write_file(
            tmp_path,
            "Vanuatu_cmip6_IPSL_IPSL-CM6A-LR_historical_r2i1p1f1_"
            "monthly-mean_1980-01-01_2005-12-31.nc",
            second,
            (1980, 2005),
        )
        catalog = SimulationCatalog.from_files([first, second_file], self.PERIODS)
        assert len(catalog.simulations()) == 2
