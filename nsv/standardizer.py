import os
from dataclasses import dataclass

import pandas as pd
import pooch
import xarray as xr
from scipy.io import loadmat
from xarray import DataArray, Dataset

from .utils import add_attributes_and_rename_variables, add_cf_attributes, compute_pt0


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
        ds["xvec"] = ds["xvec"]

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
                "units": "degree_C",
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
        ds["X"] = ds["X"]

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
            X={"long_name": "distance from sill", "units": "km"},
            station={"long_name": "station id"},
            num_data={"long_name": "occupations"},
            SIG={
                "standard_name": "sea_water_sigma_theta",
            },
            THE={
                "standard_name": "sea_water_potential_temperature",
            },
            SAL={"standard_name": "sea_water_practical_salinity"},
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
        ds["station"] = ds["station"]

        # Manually add CF attributes
        attrs_dict = dict(
            Y={
                "standard_name": "depth",
                "positive": "down",
            },
            X={"long_name": "distance from sill", "units": "km"},
            station={"long_name": "station id"},
            OrtVel={"long_name": "through section velocity", "units": "m s-1"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs_dict)

        # Rename
        ds = ds.rename(X="dist", OrtVel="through_section_velocity")

        return add_cf_attributes(ds)

    @property
    def fim_1m(self) -> Dataset:
        """Standardized FIM in 1m depth bins"""

        return self._fim(1)

    @property
    def fim_25m(self) -> Dataset:
        """Standardized FIM in 25m depth bins"""

        return self._fim(25)

    def _fim(self, resolution) -> Dataset:
        """Standardized FIM dataset"""

        # Open NetCDF file
        filename = self.raw_pooch.fetch(f"FIM/FIM_CTD_{resolution}m.nc")
        ds = xr.open_mfdataset(filename)
        ds.attrs = {"featureType": "timeSeries"}

        # Drop variables
        ds = ds.drop(
            (
                var
                for var in ds.variables
                if var.split("_")[0] in ["stn", "ctd", "cruise"]
            )
        )

        # Manually add CF attributes
        attrs_dict = dict(
            pressure={"standard_name": "depth", "positive": "down", "units": "m"},
            station={"long_name": "station id"},
            lat={"standard_name": "latitude"},
            lon={"standard_name": "longitude"},
            bathy={"standard_name": "sea_floor_depth_below_geoid"},
            temp={
                "standard_name": "sea_water_temperature",
            },
            sal={"standard_name": "sea_water_practical_salinity"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs_dict)

        # Compute potential temperature
        ds["potential_temperature"] = compute_pt0(ds)

        return add_cf_attributes(ds)

    @property
    def osnap(self) -> Dataset:
        """Standardized OSNAP dataset"""

        # Open NetCDF file
        filenames = [
            self.raw_pooch.fetch(f"OSNAP/OSNAP_{file}_201408_201604_2018.nc")
            for file in ("Gridded_TS", "Transports")
        ]
        ds = xr.open_mfdataset(filenames)
        ds.attrs = {"featureType": "timeSeries"}

        # Rename dimensions
        ds = ds.rename_dims(LONGITUDE="station", LATITUDE="station")
        ds["station"] = ds["station"]

        # Manually add CF attributes
        attrs = dict(
            station={"long_name": "station id"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs)

        # Compute potential temperature
        ds["potential_temperature"] = compute_pt0(ds)

        return add_cf_attributes(ds)

    @property
    def ovide(self) -> Dataset:
        """Standardized OVIDE dataset"""

        # Open NetCDF file
        filename = self.raw_pooch.fetch("OVIDE/Daniault_2016_climatology.nc")
        ds = xr.open_mfdataset(filename)
        ds.attrs = {"featureType": "timeSeries"}

        # Rename dimensions
        ds = ds.rename_dims(Dist="station", DistTr="mid_station")
        ds = ds.rename(Distance="distance", Dist_Transport="mid_distance")
        ds = ds.drop("Year")
        ds["station"] = ds["station"]
        ds["mid_station"] = ds["mid_station"] + 0.5

        # STATION
        # Manually add CF attributes
        attrs = dict(
            Depth={"standard_name": "depth"},
            station={"long_name": "station id"},
            Longitude={"standard_name": "longitude"},
            Latitude={"standard_name": "latitude"},
            Bottom={"standard_name": "sea_floor_depth_below_geoid"},
            SAL={"standard_name": "sea_water_practical_salinity"},
            TPOT={"standard_name": "sea_water_potential_temperature"},
            SIG0={"standard_name": "sea_water_sigma_theta"},
        )
        ds_station = ds.drop_dims("mid_station")
        ds_station = add_attributes_and_rename_variables(ds_station, attrs)

        # Sigmas
        pressures = [1, 2, 4]
        sig = xr.concat(
            [ds_station[f"SIG{i}"] for i in pressures], "reference_pressure"
        )
        sig["reference_pressure"] = DataArray(
            [i * 1.0e7 for i in pressures],
            dims="reference_pressure",
            attrs={"standard_name": "reference_pressure"},
        )
        sig.attrs["standard_name"] = "sea_water_sigma_theta"
        ds_station = ds_station.drop((f"SIG{i}" for i in pressures))
        ds_station["SIG"] = sig

        # MID STATION
        # Manually add CF attributes
        attrs = dict(
            Depth={"standard_name": "depth"},
            mid_station={"long_name": "mid station id"},
            Lon_Transport={"standard_name": "longitude"},
            Lat_Transport={"standard_name": "latitude"},
            Bottom_Transport={"standard_name": "sea_floor_depth_below_geoid"},
        )
        ds_mid = ds.drop_dims("station")
        ds_mid = add_attributes_and_rename_variables(ds_mid, attrs)
        ds_mid = ds_mid.rename(
            {
                var: f"mid_{var}"
                for var, da in ds_mid.data_vars.items()
                if "mid_station" in da.dims
            }
        )

        return add_cf_attributes(ds_station.merge(ds_mid))
