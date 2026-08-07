"""Typed models describing local simulation files and their inventory.

- :class:`Simulation`: the product-agnostic identity of a simulation run.
- :class:`SimulationFile`: a single local file, binding a simulation identity,
  a period and a path.
- :class:`SimulationCatalog`: the inventory of local files used by the
  downscaling workflow.

Simulation metadata is stored in the netCDF global attributes of the files
(prefixed with ``downclim_``) and used as the primary source of information,
with downclim-style filename parsing as a fallback.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from datetime import datetime as dt
from enum import Enum
from pathlib import Path
from typing import Any

import xarray as xr
from pydantic import BaseModel, ConfigDict, ValidationError

from ..logging_config import get_logger
from .utils import Aggregation, DataProduct

logger = get_logger(__name__)

#: Prefix used for the netCDF global attributes storing a simulation context.
ATTR_PREFIX = "downclim_"

#: Fields identifying a simulation run regardless of experiment / period.
_IDENTITY_FIELDS = (
    "institute",
    "source",
    "ensemble",
    "experiment",
    "domain",
    "driving_model",
    "rcm_name",
    "rcm_version",
)

#: Filename layout (parts after the AOI prefix) for each data product.
_FIELDS_BY_PRODUCT: dict[str, tuple[str, ...]] = {
    "cmip6": (
        "aoi_n",
        "data_product",
        "institute",
        "source",
        "experiment",
        "ensemble",
        "aggregation",
        "tmin",
        "tmax",
    ),
    "cordex": (
        "aoi_n",
        "data_product",
        "domain",
        "driving_model",
        "rcm_name",
        "ensemble",
        "rcm_version",
        "aggregation",
        "tmin",
        "tmax",
    ),
}


class PeriodKind(Enum):
    """Kind of period a local simulation file belongs to."""

    HISTORICAL = "historical"
    EVALUATION = "evaluation"
    PROJECTION = "projection"


def _data_product_from_name(name: str) -> DataProduct:
    """Return the DataProduct enum member matching ``name`` (its product_name)."""
    for product in DataProduct:
        if product.product_name == name:
            return product
    msg = f"Unknown data product '{name}'."
    raise ValueError(msg)


def _year(value: Any) -> int:
    """Extract the year from a string such as '2005' or '2005-12-31'."""
    text = str(value)
    for fmt in ("%Y-%m-%d", "%Y"):
        try:
            return dt.strptime(text, fmt).year
        except ValueError:
            continue
    msg = f"Cannot extract a year from '{value}'."
    raise ValueError(msg)


def _parse_filename(filename: str | Path) -> tuple[dict[str, str], DataProduct]:
    """Split a downclim-style filename into a context mapping and its data product.

    Two layouts are supported, identified by the data product token:
    - cmip6:  {aoi}_{product}_{institute}_{source}_{experiment}_{ensemble}_{aggregation}_{tmin}_{tmax}
    - cordex: {aoi}_{product}_{domain}_{driving_model}_{rcm_name}_{ensemble}_{rcm_version}_{aggregation}_{tmin}_{tmax}

    Returns:
        tuple[dict[str, str], DataProduct]: filename context and data product.
    """
    parts = Path(filename).name.split(".nc")[0].split("_")
    if len(parts) < 2:
        msg = f"Cannot parse simulation context from filename '{filename}'."
        raise ValueError(msg)
    product = _data_product_from_name(parts[1])
    fields = _FIELDS_BY_PRODUCT[product.product_name]
    if len(parts) != len(fields):
        msg = (
            f"Cannot parse simulation context from filename '{filename}': "
            f"expected {len(fields)} parts, got {len(parts)}."
        )
        raise ValueError(msg)
    return dict(zip(fields, parts, strict=True)), product


def _read_attrs(path: Path) -> dict[str, Any]:
    """Read the global attributes of a netCDF file without loading its data."""
    for engine in (None, "scipy"):
        try:
            with xr.open_dataset(path, decode_cf=False, engine=engine) as ds:
                return dict(ds.attrs)
        except (OSError, ValueError, RuntimeWarning):
            logger.debug(
                "Could not read attributes of %s with engine %r.", path, engine
            )
    return {}


class Simulation(BaseModel):
    """Product-agnostic identity of a simulation run.

    Fields that do not apply to a data product (e.g. ``domain`` for a CMIP6
    simulation) are ``None``.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    product: DataProduct
    aoi_n: str
    institute: str | None = None
    source: str | None = None
    ensemble: str | None = None
    experiment: str | None = None
    domain: str | None = None
    driving_model: str | None = None
    rcm_name: str | None = None
    rcm_version: str | None = None
    aggregation: Aggregation = Aggregation.MONTHLY_MEAN  # type: ignore[assignment]

    def matching_key(self) -> tuple[str | None, ...]:
        """Fields identifying a simulation run regardless of experiment / period."""
        identity_fields = (
            ("source", "ensemble")
            if self.product is DataProduct.CMIP6
            else ("domain", "driving_model", "rcm_name", "ensemble", "rcm_version")
        )
        return (
            self.product.name,
            self.aoi_n,
            self.aggregation.value,
            *(getattr(self, field) for field in identity_fields),
        )

    def brief(self) -> str:
        """Short human-readable description of the simulation."""
        key = self.matching_key()
        identity = " / ".join(str(value) for value in key[3:])
        return f"{key[0]} simulation ({identity})"

    def to_attrs(self) -> dict[str, str]:
        """Serialize the simulation into netCDF global attributes."""
        attrs: dict[str, str] = {}
        for name, value in self.model_dump(exclude_none=True).items():
            if name == "product":
                attrs[f"{ATTR_PREFIX}{name}"] = self.product.product_name
            elif name == "aggregation":
                attrs[f"{ATTR_PREFIX}{name}"] = self.aggregation.value
            else:
                attrs[f"{ATTR_PREFIX}{name}"] = str(value)
        return attrs

    @classmethod
    def from_attrs(cls, attrs: Mapping[str, Any]) -> Simulation | None:
        """Reconstruct a Simulation from netCDF global attributes.

        Returns None if no ``downclim_`` attributes are present.
        """
        if not any(key.startswith(ATTR_PREFIX) for key in attrs):
            return None
        values: dict[str, Any] = {}
        for name in cls.model_fields:
            key = f"{ATTR_PREFIX}{name}"
            if key not in attrs:
                continue
            raw = attrs[key]
            if name == "product":
                values[name] = _data_product_from_name(str(raw))
            elif name == "aggregation":
                values[name] = Aggregation(str(raw))
            else:
                values[name] = raw
        return cls(**values)

    @classmethod
    def from_context(cls, context: Mapping[str, str]) -> Simulation:
        """Build a Simulation from a filename context mapping."""
        values: dict[str, Any] = {
            "product": _data_product_from_name(context["data_product"]),
            "aoi_n": context["aoi_n"],
            "aggregation": Aggregation(context["aggregation"]),
        }
        for field in _IDENTITY_FIELDS:
            if field in context:
                values[field] = context[field]
        return cls(**values)

    @classmethod
    def from_filename(cls, filename: str | Path) -> Simulation:
        """Parse a Simulation from a downclim-style netcdf filename."""
        context, _ = _parse_filename(filename)
        return cls.from_context(context)


class SimulationFile(BaseModel):
    """A single local simulation file: identity + period + path."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    simulation: Simulation
    period: tuple[int, int]
    path: Path

    def to_attrs(self) -> dict[str, str]:
        """Serialize the simulation and its period into netCDF global attributes."""
        attrs = self.simulation.to_attrs()
        attrs[f"{ATTR_PREFIX}tmin"] = str(self.period[0])
        attrs[f"{ATTR_PREFIX}tmax"] = str(self.period[1])
        return attrs

    @classmethod
    def from_file(cls, path: str | Path) -> SimulationFile:
        """Build a SimulationFile from a netCDF file.

        The netCDF global attributes are used as the primary metadata source,
        with downclim-style filename parsing as a fallback.
        """
        path = Path(path)
        attrs = _read_attrs(path)
        simulation = None
        try:
            simulation = Simulation.from_attrs(attrs)
        except (ValueError, ValidationError):
            logger.warning(
                "Invalid downclim attributes in %s, falling back to filename parsing.",
                path,
            )
        if simulation is not None and f"{ATTR_PREFIX}tmin" in attrs:
            period = (
                _year(attrs[f"{ATTR_PREFIX}tmin"]),
                _year(attrs[f"{ATTR_PREFIX}tmax"]),
            )
            return cls(simulation=simulation, period=period, path=path)
        context, _ = _parse_filename(path)
        simulation = Simulation.from_context(context)
        period = (_year(context["tmin"]), _year(context["tmax"]))
        return cls(simulation=simulation, period=period, path=path)


class SimulationCatalog:
    """Local file inventory, grouped by data product and period.

    The catalog is built once per area of interest and holds the files to
    downscale. Simulations are matched by :meth:`Simulation.matching_key`, so
    that historical, evaluation and projection files of the same simulation run
    are associated regardless of their experiment.
    """

    def __init__(
        self,
        files: Iterable[SimulationFile],
        periods: Mapping[PeriodKind, tuple[int, int]],
    ) -> None:
        self._files = list(files)
        self._periods = dict(periods)

    @classmethod
    def from_files(
        cls,
        paths: Iterable[str | Path],
        periods: Mapping[PeriodKind, tuple[int, int]],
    ) -> SimulationCatalog:
        """Build a catalog from the paths of local simulation files."""
        return cls((SimulationFile.from_file(path) for path in paths), periods)

    @property
    def files(self) -> list[SimulationFile]:
        """All the files of the catalog."""
        return list(self._files)

    def simulations(self) -> list[Simulation]:
        """One representative simulation per unique simulation key."""
        by_key: dict[tuple[str | None, ...], Simulation] = {}
        for file in self._files:
            by_key.setdefault(file.simulation.matching_key(), file.simulation)
        return list(by_key.values())

    def historical_files(self, simulation: Simulation) -> list[SimulationFile]:
        """Files of the historical period matching a simulation."""
        return self._files_for(simulation, PeriodKind.HISTORICAL)

    def matching_files(
        self, simulation: Simulation, period: PeriodKind
    ) -> list[SimulationFile]:
        """Files of a given period matching a simulation.

        Raises:
            FileNotFoundError: if no file matches the simulation for this period.
        """
        files = self._files_for(simulation, period)
        if len(files) == 0:
            start, end = self._periods[period]
            msg = (
                f"No matching files found for {simulation.brief()} in the "
                f"{period.value} period ({start}-{end})."
            )
            raise FileNotFoundError(msg)
        return files

    def _files_for(
        self, simulation: Simulation, period: PeriodKind
    ) -> list[SimulationFile]:
        wanted = self._periods[period]
        key = simulation.matching_key()
        return [
            file
            for file in self._files
            if file.simulation.matching_key() == key and file.period == wanted
        ]
