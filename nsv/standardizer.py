import os
from dataclasses import dataclass
from string import ascii_uppercase

import pandas as pd
import pooch
import xarray as xr
from scipy.io import loadmat
from xarray import DataArray, Dataset

from .utils import (
    add_attributes_and_rename_variables,
    compute_pt0,
    dms2d,
    final_cleanup_before_returning,
)


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
    @final_cleanup_before_returning
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
        ds = Dataset(variables, coords=coords)
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
            xvec="distance", ws="through_section_velocity", wc="cross_section_velocity"
        )

        return ds

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
        ds = Dataset(variables, coords=coords)
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
    @final_cleanup_before_returning
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
        ds = ds.rename(X="distance")

        return ds

    @property
    @final_cleanup_before_returning
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
        ds = ds.rename(X="distance", OrtVel="through_section_velocity")

        return ds

    @final_cleanup_before_returning
    def fim(self, resolution: int) -> Dataset:
        """
        Standardized FIM dataset

        Args:
            sec (int): resolution {1, 25}

        Returns:
            Dataset: Standardized dataset
        """

        # Check input
        ok_res = [1, 25]
        if resolution not in ok_res:
            raise ValueError(
                f"{resolution} is not available. Available sections: {ok_res!r}"
            )

        # Open NetCDF file
        filename = self.raw_pooch.fetch(f"FIM/FIM_CTD_{resolution}m.nc")
        ds = xr.open_mfdataset(filename)

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

        return ds

    @property
    @final_cleanup_before_returning
    def osnap(self) -> Dataset:
        """Standardized OSNAP dataset"""

        # Open mat file
        filename = self.raw_pooch.fetch(
            "OSNAP/OSNAP_Gridded_TSV_201408_201805_2021.mat"
        )
        mat = loadmat(filename, squeeze_me=True, mat_dtype=True)

        # Dimensions
        coords = {
            name: DataArray(mat["osnap"][name.upper()].item(), dims=name)
            for name in ("lon", "depth", "time")
        }
        coords2 = {
            name: DataArray(mat["osnap"][name.upper()].item(), dims=name)
            for name in ("lon", "depth")
        }
        coords1 = {
            name: DataArray(mat["osnap"][name.upper()].item(), dims=name)
            for name in ("lon",)
        }

        # Variables
        variables = {
            name: DataArray(
                mat["osnap"][name.upper()].item(), dims=tuple(coords), coords=coords
            )
            for name in ("theta", "psal", "velo")
        }
        variables["area"] = DataArray(
            mat["osnap"]["AREA"].item(), dims=tuple(coords2), coords=coords2
        )
        variables["longitude"] = DataArray(
            mat["osnap"]["LON"].item(), dims=tuple(coords1), coords=coords1
        )
        variables["latitude"] = DataArray(
            mat["osnap"]["LAT"].item(), dims=tuple(coords1), coords=coords1
        )
        # Initialize dataset
        ds = Dataset(variables, coords=coords)
        ds = ds.rename_dims(lon="station")
        ds["station"] = ds["station"]
        ds = ds.drop_vars("lon")

        # Add coordinates
        ds = ds.set_coords(["longitude", "latitude"])

        # Manually add CF attributes
        attrs_dict = dict(
            time={"standard_name": "time", "long_name": "time vector for gridded data"},
            depth={
                "standard_name": "depth",
                "long_name": "depth of gridded product",
                "positive": "down",
            },
            station={"long_name": "station id"},
            longitude={"standard_name": "longitude"},
            latitude={"standard_name": "latitude"},
            velo={"long_name": "velocity normal to the section", "units": "m s-1"},
            psal={"standard_name": "sea_water_practical_salinity"},
            theta={
                "standard_name": "sea_water_potential_temperature",
                "long_name": "Potential Temperature",
                "units": "degree_C",
            },
        )
        ds = add_attributes_and_rename_variables(ds, attrs_dict)

        return ds

    @property
    @final_cleanup_before_returning
    def ovide(self) -> Dataset:
        """Standardized OVIDE dataset"""

        # Open NetCDF file
        filename = self.raw_pooch.fetch("OVIDE/Daniault_2016_climatology.nc")
        ds = xr.open_mfdataset(filename)

        # Rename dimensions
        ds = ds.rename_dims(Dist="station", DistTr="mid_station")
        ds = ds.rename(Distance="distance", Dist_Transport="mid_distance")
        ds = ds.drop("Year")
        ds["station"] = ds["station"]
        ds["mid_station"] = ds["mid_station"] + 0.5

        # STATION
        # Manually add CF attributes
        attrs = dict(
            Depth={"standard_name": "depth", "positive": "down"},
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
            Depth={"standard_name": "depth", "positive": "down"},
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

        return ds_station.merge(ds_mid)

    @property
    @final_cleanup_before_returning
    def ho2000(self) -> Dataset:
        """Standardized Hansen & Osterhus 2000 dataset"""

        # Read variables
        dataframes = {}
        for var in ("Tem", "Sal"):
            filename = self.raw_pooch.fetch(f"Hansen_Osterhus_2000/{var}_data.csv")
            df = pd.read_csv(filename, index_col="Depth")
            dataframes[var] = df.apply(pd.to_numeric, errors="coerce")
        ds = Dataset(dataframes).rename(dim_1="Stations")

        # Assign coordinates
        coords_fname = self.raw_pooch.fetch("Hansen_Osterhus_2000/Stations_coord.csv")
        coords_df = pd.read_csv(coords_fname).set_index("Stations")
        coords_ds = Dataset(coords_df).drop("Stations")
        ds = ds.assign_coords(
            {
                pref: dms2d(*[coords_ds[f"{pref} {suff}"] for suff in ("deg", "min")])
                for pref in ("LON", "LAT")
            }
        )

        # Rename dimensions
        ds = ds.rename(Stations="station")

        # Manually add CF attributes
        attrs = dict(
            Depth={
                "standard_name": "depth",
                "positive": "down",
            },
            station={"long_name": "station id"},
            LON={"standard_name": "longitude"},
            LAT={"standard_name": "latitude"},
            Sal={"standard_name": "sea_water_practical_salinity"},
            Tem={"standard_name": "sea_water_temperature"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs)

        # Compute potential temperature
        ds["potential_temperature"] = compute_pt0(ds)

        return ds

    @final_cleanup_before_returning
    def m82_1(self, sec_id) -> Dataset:
        """
        Standardized Meteor cruise M82/1

        Args:
        sec (int): Section ID {1, 2, 3, 4, 5, 6, 7, 8, 9}

        Returns:
            Dataset: Standardized dataset
        """

        # Check input
        sec_dict = {
            1: list(range(322, 333)) + list(range(334, 337)),
            2: range(337, 345),
            3: range(347, 363),
            4: list(range(363, 379)),
            5: range(379, 387),
            6: (list(range(387, 396)) + list(range(401, 407))),
            7: range(407, 427),
            8: list(range(437, 445)) + list(range(430, 426)),
            9: list(range(451, 457)),
        }
        if sec_id not in sec_dict:
            raise ValueError(
                f"{sec_id} is not available. Available sections: {list(sec_dict)!r}"
            )

        # Open CSV file
        filename = self.raw_pooch.fetch("Quadfasel_2018/M82-1_CTD.tab")
        df = pd.read_csv(filename, sep="\t", header=176, parse_dates=["Date/Time"])

        # Extract section
        events = [f"M82/1_{i}" for i in sec_dict[sec_id]]
        df = df[df["Event"].isin(events)]

        # Transform to Dataset (round depths and average duplicates)
        df["Depth water [m]"] = df["Depth water [m]"].round()
        indexes = ["Event", "Depth water [m]"]
        df = (
            df.set_index(indexes)
            .groupby(indexes)
            .mean(numeric_only=False)
            .sort_values(indexes)
        )
        ds = df.to_xarray()

        # Rename dimensions and assign coordinates
        ds = ds.assign_coords(
            {
                coord: ds[coord].mean("Depth water [m]")
                for coord in ["Date/Time", "Longitude", "Latitude", "Elevation [m]"]
            }
        )
        ds = ds.reset_coords("Elevation [m]").rename(Event="station")

        # Manually add CF attributes
        ds = ds.rename(
            {"DO [ml/l]": "DO", "O2 [µmol/l]": "O2", "Transmission [%]": "Transmission"}
        )
        ds["Elevation [m]"] *= -1
        attrs = {
            "Depth water [m]": {
                "standard_name": "depth",
                "positive": "down",
            },
            "Date/Time": {"standard_name": "time"},
            "station": {"long_name": "station id"},
            "Longitude": {"standard_name": "longitude"},
            "Latitude": {"standard_name": "latitude"},
            "Temp [°C]": {"standard_name": "sea_water_temperature"},
            "Sal": {"standard_name": "sea_water_practical_salinity"},
            "Press [dbar]": {"standard_name": "sea_water_pressure"},
            "DO": {"long_name": "oxygen dissolved", "units": "ml/l"},
            "O2": {"long_name": "oxygen", "units": "µmol/l"},
            "Transmission": {"long_name": "transmission of light", "units": "%"},
            "Elevation [m]": {"standard_name": "sea_floor_depth_below_geoid"},
        }
        ds = add_attributes_and_rename_variables(ds, attrs)

        # Compute potential temperature
        ds["potential_temperature"] = compute_pt0(ds)

        return ds

    @property
    @final_cleanup_before_returning
    def eel(self) -> Dataset:
        """Standardized EEL dataset"""

        # Open CSV file
        data_name = self.raw_pooch.fetch("EEL/csv_ctdgrid/EELCTDandLADCP_3Dfield.csv")
        coords_name = self.raw_pooch.fetch("EEL/csv_ctdgrid/EELCTDandLADCP_refpos.csv")
        pd_args = {"sep": ",", "index_col": None, "header": 0}
        df_coords = pd.read_csv(coords_name, **pd_args)
        df_data = pd.read_csv(data_name, **pd_args)

        # Transform to Dataset
        indexes = ["Refdist", "Year", "Depth"]
        ds_coords = df_coords.set_index("Refdist").sort_values("Refdist").to_xarray()
        ds_data = df_data.set_index(indexes).sort_values(indexes).to_xarray()
        ds = xr.merge([ds_coords, ds_data])

        # Rename dimensions and assign coordinates
        ds = ds.rename(Staname="station", Refdist="distance")
        ds["station"] = ds["station"].isel(Year=0, Depth=0)
        ds = ds.swap_dims(distance="station")
        ds["Year"] = pd.to_datetime(ds["Year"].astype(str))

        # Manually add CF attributes
        attrs = dict(
            distance={"long_name": "distance", "units": "km"},
            station={"long_name": "station id"},
            Year={"standard_name": "time"},
            Depth={"standard_name": "depth", "positive": "down"},
            CruiseID={"long_name": "cruise id"},
            LonSta={"standard_name": "longitude"},
            LatSta={"standard_name": "latitude"},
            DepthSta={"standard_name": "sea_floor_depth_below_geoid"},
            PTMP={"standard_name": "sea_water_potential_temperature"},
            PSAL={"standard_name": "sea_water_practical_salinity"},
            Sigma0={"standard_name": "sea_water_sigma_theta"},
            Vrel={"long_name": "relative velocity", "units": "m s-1"},
            Vladcp={"long_name": "LADCP through section velocity", "units": "m s-1"},
            Vabs={"long_name": "absolute velocity", "units": "m s-1"},
            Vladcpalong={"long_name": "LADCP cross section velocity", "units": "m s-1"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs)

        return ds

    @final_cleanup_before_returning
    def kn203_2(self, sec_id: str) -> Dataset:
        """
        Standardized Knorr cruise KN203-2

        Args:
            sec_id (str): Section ID {"A", "B", "C", "D", "E"}

        Returns:
            Dataset: Standardized dataset
        """

        # Check input
        sec_id = sec_id.upper()
        ok_sec = list(ascii_uppercase[:5])
        if sec_id not in ok_sec:
            raise ValueError(
                f"{sec_id} is not available. Available sections: {ok_sec!r}"
            )

        # Open CSV file
        fname = self.raw_pooch.fetch(
            f"Semper_2020/CTD/Section{sec_id}_NEIceland_CTD.tab"
        )
        with open(fname, "r") as f:
            for skiprow, row in enumerate(f):
                if row.startswith("*"):
                    break
        df = pd.read_csv(fname, sep="\t", header=skiprow + 1, parse_dates=["Date/Time"])

        # Transform to Dataset
        indexes = ["Event", "Press [dbar]"]
        ds = df.set_index(indexes).sort_values(indexes).to_xarray()

        # Rename dimensions and assign coordinates
        ds = ds.assign_coords(
            {
                coord: ds[coord].mean("Press [dbar]")
                for coord in ["Date/Time", "Longitude", "Latitude"]
            }
        )
        ds = ds.rename(Event="station")

        # Manually add CF attributes
        attrs = {
            "Date/Time": {"standard_name": "time"},
            "Press [dbar]": {"standard_name": "depth", "positive": "down"},
            "station": {"long_name": "station id"},
            "Longitude": {"standard_name": "longitude"},
            "Latitude": {"standard_name": "latitude"},
            "Sal": {"standard_name": "sea_water_practical_salinity"},
            "Temp [°C]": {"standard_name": "sea_water_temperature"},
        }
        ds = add_attributes_and_rename_variables(ds, attrs)

        # Compute potential temperature
        ds["potential_temperature"] = compute_pt0(ds)

        return ds

    @final_cleanup_before_returning
    def pos503(self, sec_id) -> Dataset:
        """
        Standardized Meteor cruise POS503

        Args:
        sec (int): Section ID {1, 2, 3, 4, 5}

        Returns:
            Dataset: Standardized dataset
        """

        # Check input
        sec_dict = {
            1: list(range(423, 435)),
            2: list(range(435, 447)),
            3: list(range(447, 456)),
            4: list(range(456, 464)),
            5: list(range(464, 497)),
        }
        if sec_id not in sec_dict:
            raise ValueError(
                f"{sec_id} is not available. Available sections: {list(sec_dict)!r}"
            )

        # Open CSV file
        filename = self.raw_pooch.fetch("Quadfasel_2018/POS503_CTD.tab")
        df = pd.read_csv(filename, sep="\t", header=156, parse_dates=["Date/Time"])

        # Extract section
        events = [f"POS503_{i}-1" for i in sec_dict[sec_id]]
        df = df[df["Event"].isin(events)]

        # Transform to Dataset (round depths and average duplicates)
        df["Depth water [m]"] = df["Depth water [m]"].round()
        indexes = ["Event", "Depth water [m]"]
        df = (
            df.set_index(indexes)
            .groupby(indexes)
            .mean(numeric_only=False)
            .sort_values(indexes)
        )
        ds = df.to_xarray()

        # Rename dimensions and assign coordinates
        ds = ds.assign_coords(
            {
                coord: ds[coord].mean("Depth water [m]")
                for coord in ["Date/Time", "Longitude", "Latitude", "Elevation [m]"]
            }
        )
        ds = ds.reset_coords("Elevation [m]").rename(Event="station")

        # Manually add CF attributes
        ds = ds.rename({"DO [ml/l]": "DO", "O2 [µmol/l]": "O2"})
        ds["Elevation [m]"] *= -1
        attrs = {
            "Depth water [m]": {
                "standard_name": "depth",
                "positive": "down",
            },
            "Date/Time": {"standard_name": "time"},
            "station": {"long_name": "station id"},
            "Longitude": {"standard_name": "longitude"},
            "Latitude": {"standard_name": "latitude"},
            "Temp [°C]": {"standard_name": "sea_water_temperature"},
            "Sal": {"standard_name": "sea_water_practical_salinity"},
            "Press [dbar]": {"standard_name": "sea_water_pressure"},
            "DO": {"long_name": "oxygen dissolved", "units": "ml/l"},
            "O2": {"long_name": "oxygen", "units": "µmol/l"},
            "Elevation [m]": {"standard_name": "sea_floor_depth_below_geoid"},
        }
        ds = add_attributes_and_rename_variables(ds, attrs)

        # Compute potential temperature
        ds["potential_temperature"] = compute_pt0(ds)

        return ds
