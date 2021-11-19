import os
from dataclasses import dataclass

import pandas as pd
import pooch
from scipy.io import loadmat
from xarray import DataArray, Dataset

from .utils import add_attributes_and_rename_variables, add_cf_attributes


@dataclass
class Standardizer:
    """Standardize raw data"""

    raw_data_path: str = None

    @property
    def raw_pooch(self):
        pooch_obj = pooch.create(
            path=self.raw_data_path or pooch.os_cache("NORVAL"),
            base_url="https://gws-access.jasmin.ac.uk/public/jmmp/NORVAL/",
            registry=None,
        )
        pooch_obj.load_registry(
            os.path.join(os.path.dirname(__file__), "registry_raw.txt")
        )
        return pooch_obj

    @property
    def kogur(self) -> Dataset:
        """Standardized Kögur dataset"""

        # Open mat file
        filename = self.raw_pooch.fetch("Kogur/all_gridded.mat")
        mat = loadmat(filename, squeeze_me=True, mat_dtype=True)

        # Dimensions
        coords = {
            name: DataArray(mat[name], dims=name) for name in ("tvec", "dvec", "xvec")
        }

        # Variables
        variables = {
            name: DataArray(mat[name], dims=tuple(coords), coords=coords)
            for name in ("ugrid", "vgrid", "ws", "Sfinal1", "Tfinal1", "PDfinal1", "wc")
        }

        # Initialize dataset
        ds = Dataset(variables, coords=coords, attrs={"featureType": "timeSeries"})
        ds = ds.rename_dims(xvec="station")
        ds["station"] = ds["station"]

        # Convert matlab time
        ds["tvec"] = pd.to_datetime(ds["tvec"].values - 719529, unit="D")
        ds["tvec"] = ds["tvec"].dt.round("H")

        # Add coordinates
        moorings = Dataset(
            {var: DataArray(mat[var], dims="dist") for var in ("lon", "lat", "dist")}
        )
        moorings = moorings.rename(dist="xvec")
        interp = moorings.interp(xvec=ds["xvec"])
        for coord in ("lon", "lat"):
            interp[coord][-1] = moorings[coord][-1]
            ds = ds.assign_coords({coord: interp[coord]})
        ds["xvec"] = ds["xvec"] * 1.0e3

        # Manually add CF attributes
        attrs_dict = dict(
            tvec={"standard_name": "time", "long_name": "time vector for gridded data"},
            dvec={
                "standard_name": "depth",
                "long_name": "depth of gridded product",
                "positive": "down",
            },
            xvec={"long_name": "distance vector for gridded product", "units": "m"},
            station={"long_name": "station id"},
            lon={"standard_name": "longitude"},
            lat={"standard_name": "latitude"},
            ugrid={"standard_name": "sea_water_x_velocity", "long_name": "u velocity"},
            vgrid={"standard_name": "sea_water_y_velocity", "long_name": "v velocity"},
            ws={"long_name": "through section velocity", "units": "m s-1"},
            Sfinal1={"standard_name": "sea_water_salinity", "long_name": "Salinity"},
            Tfinal1={
                "standard_name": "sea_water_potential_temperature",
                "long_name": "Potential Temperature",
            },
            PDfinal1={
                "standard_name": "sea_water_sigma_theta",
                "long_name": "Potential Density",
            },
            wc={"long_name": "cross section velocity", "units": "m s-1"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs_dict)

        # Rename
        ds = ds.rename(
            xvec="dist", ws="through_section_velocity", wc="cross_section_velocity"
        )

        return add_cf_attributes(ds)

    def _interpolate_latrabjarg_bathymetry(self, dist_coord: DataArray) -> Dataset:
        """Interpolate Látrabjarg bathymetry using distances from sill"""

        # Open mat file
        filename = self.raw_pooch.fetch("Latrabjarg/LatBat.mat")
        mat = loadmat(filename, squeeze_me=True, mat_dtype=True)

        # Create dataset
        ds = Dataset(
            {
                name: DataArray(data, dims="regdist")
                for name, data in mat.items()
                if not name.startswith("_")
            }
        )

        # Remove duplicates
        ds = ds.groupby("regdist").mean()

        # Center at the sill
        ds["regdist"] = ds["regdist"] - ds["regdist"][ds["regbat"].argmax()]

        # Manually add CF attributes
        attrs_dict = dict(
            reglat={"standard_name": "latitude"},
            reglon={"standard_name": "longitude"},
            regbat={"standard_name": "sea_floor_depth_below_geoid"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs_dict)

        # Interpolate
        ds = ds.rename(regdist=dist_coord.name)
        ds = ds.interp(**{dist_coord.name: dist_coord.values})
        ds = ds.rename_dims({dist_coord.name: "station"})
        return ds

    def _initialize_latrabjarg(self, mat) -> Dataset:
        """Initialize Látrabjarg dataset"""

        # Dimensions
        coords = {
            name: DataArray(
                mat[name][0 if name == "X" else ..., 0 if name == "Y" else ...],
                dims=name,
            )
            for name in ("Y", "X")
        }

        # Variables
        variables = {
            name: DataArray(mat[name], dims=tuple(coords), coords=coords)
            for name in mat
            if (name not in coords and not name.startswith("_") and name != "PCENT")
        }

        # Initialize dataset
        ds = Dataset(variables, coords=coords, attrs={"featureType": "timeSeries"})
        ds = ds.rename_dims(X="station")
        ds["station"] = ds["station"]

        # Center at the sill
        notnull_count = ds["OrtVel" if "OrtVel" in ds else "SIG"].notnull().sum("Y")
        sill_station = (
            ds["station"].where(notnull_count == notnull_count.max()).mean().round()
        )
        ds["X"] = ds["X"] - ds["X"][int(sill_station)]

        # Add coordinates
        ds = ds.merge(self._interpolate_latrabjarg_bathymetry(ds["X"]))
        ds = ds.set_coords(["longitude", "latitude"])
        ds["X"] = ds["X"] * 1.0e3

        return ds

    @property
    def latrabjarg_climatology(self) -> Dataset:
        """Standardized Látrabjarg climatology dataset"""

        # Open mat file
        filename = self.raw_pooch.fetch("Latrabjarg/meanFields_MastrpoleEtAl2017.mat")
        mat = loadmat(filename, squeeze_me=True, mat_dtype=True)

        # Create dataset
        ds = self._initialize_latrabjarg(mat)

        # Manually add CF attributes
        attrs_dict = dict(
            Y={
                "standard_name": "depth",
                "positive": "down",
            },
            X={"long_name": "distance from sill", "units": "m"},
            station={"long_name": "station id"},
            num_data={"long_name": "occupations"},
            SIG={
                "standard_name": "sea_water_sigma_theta",
            },
            THE={
                "standard_name": "sea_water_potential_temperature",
            },
            SAL={"standard_name": "sea_water_salinity"},
            NSQ={"standard_name": "square_of_brunt_vaisala_frequency_in_sea_water"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs_dict)

        # Rename
        ds = ds.rename(X="dist")

        return add_cf_attributes(ds)

    @property
    def latrabjarg_survey(self) -> Dataset:
        """Standardized Látrabjarg survey dataset"""

        # Open mat file
        filename = self.raw_pooch.fetch("Latrabjarg/Velocities_VageEtAl2011.mat")
        mat = loadmat(filename, squeeze_me=True, mat_dtype=True)

        # Create dataset
        ds = self._initialize_latrabjarg(mat)

        # Manually add CF attributes
        attrs_dict = dict(
            Y={
                "standard_name": "depth",
                "positive": "down",
            },
            X={"long_name": "distance from sill", "units": "m"},
            station={"long_name": "station id"},
            OrtVel={"long_name": "through section velocity", "units": "m s-1"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs_dict)

        # Rename
        ds = ds.rename(X="dist", OrtVel="through_section_velocity")

        return add_cf_attributes(ds)
