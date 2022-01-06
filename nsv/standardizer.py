import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pooch
import xarray as xr
from gsw import z_from_p
from scipy.io import loadmat
from xarray import DataArray, Dataset

from .utils import (
    add_attributes_and_rename_variables,
    add_cf_attributes,
    compute_pt0,
    dms2d,
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

    @property
    def HO2000(self) -> Dataset:
        """Standardized Hansen & Osterhus 2000 dataset"""

        Tfilename = self.raw_pooch.fetch("Hansen_Osterhus_2000/Tem_data.csv")
        Sfilename = self.raw_pooch.fetch("Hansen_Osterhus_2000/Sal_data.csv")
        Cfilename = self.raw_pooch.fetch("Hansen_Osterhus_2000/Stations_coord.csv")

        df_T = pd.read_csv(Tfilename)
        df_S = pd.read_csv(Sfilename)
        df_C = pd.read_csv(Cfilename)
        dfT = df_T.loc[:, "S1":"S15"]
        dfS = df_S.loc[:, "S1":"S15"]
        dfD = df_T.loc[:, "Depth"]
        dfC = df_C.loc[:, "LON deg":"LAT min"]
        dfT = dfT.apply(pd.to_numeric, errors="coerce")
        dfS = dfS.apply(pd.to_numeric, errors="coerce")
        dfD = dfD.apply(pd.to_numeric, errors="coerce")
        dfC = dfC.apply(pd.to_numeric, errors="coerce")

        # Converting Lat and Lon
        COORD = dfC.to_numpy()
        lat = []
        lon = []
        for n in range(COORD.shape[0]):
            lon.append(dms2d(COORD[n, 0], COORD[n, 1], 0.0))
            lat.append(dms2d(COORD[n, 2], COORD[n, 3], 0.0))

        # Initialize dataset
        ds = Dataset(
            data_vars=dict(
                TEM=(["Depth", "station"], dfT.to_numpy()),
                SAL=(["Depth", "station"], dfS.to_numpy()),
                Longitude=(
                    [
                        "station",
                    ],
                    lon,
                ),
                Latitude=(
                    [
                        "station",
                    ],
                    lat,
                ),
            ),
            coords=dict(
                station=(["station"], range(1, 16)),
                Depth=(["Depth"], dfD.to_numpy()),
            ),
            attrs=dict(description="Hansen & Osterhus 2000 dataset"),
        )

        # Manually add CF attributes
        attrs = dict(
            Depth={"standard_name": "depth"},
            station={"long_name": "station id"},
            Longitude={"standard_name": "longitude"},
            Latitude={"standard_name": "latitude"},
            SAL={"standard_name": "sea_water_practical_salinity"},
            TEM={"standard_name": "sea_water_temperature"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs)

        # Compute potential temperature
        ds["potential_temperature"] = compute_pt0(ds)

        return add_cf_attributes(ds)

    @property
    def Q2018_sec1(self) -> Dataset:
        """Standardized Section 1 of Quadfasel et al 2018 dataset"""

        return self._Q2018(0)

    @property
    def Q2018_sec2(self) -> Dataset:
        """Standardized Section 2 of Quadfasel et al 2018 dataset"""

        return self._Q2018(1)

    @property
    def Q2018_sec3(self) -> Dataset:
        """Standardized Section 3 of Quadfasel et al 2018 dataset"""

        return self._Q2018(2)

    @property
    def Q2018_sec4(self) -> Dataset:
        """Standardized Section 4 of Quadfasel et al 2018 dataset"""

        return self._Q2018(3)

    @property
    def Q2018_sec5(self) -> Dataset:
        """Standardized Section 5 of Quadfasel et al 2018 dataset"""

        return self._Q2018(4)

    @property
    def Q2018_sec6(self) -> Dataset:
        """Standardized Section 6 of Quadfasel et al 2018 dataset"""

        return self._Q2018(5)

    @property
    def Q2018_sec7(self) -> Dataset:
        """Standardized Section 7 of Quadfasel et al 2018 dataset"""

        return self._Q2018(6)

    @property
    def Q2018_sec8(self) -> Dataset:
        """Standardized Section 8 of Quadfasel et al 2018 dataset"""

        return self._Q2018(7)

    @property
    def Q2018_sec9(self) -> Dataset:
        """Standardized Section 9 of Quadfasel et al 2018 dataset"""

        return self._Q2018(8)

    def _Q2018(self, section) -> Dataset:
        """Standardized Quadfasel et al 2018 dataset"""

        filename = self.raw_pooch.fetch("Quadfasel_2018/M82-1_CTD.tab")
        df = pd.read_csv(filename, sep="\t", header=176)

        # Creating Sections
        sec1 = list(range(322, 333)) + list(range(334, 337))
        sec2 = range(337, 345)
        sec3 = range(347, 363)
        sec4 = list(range(363, 379))[::-1]
        sec5 = range(379, 387)
        sec6 = (list(range(387, 396)) + list(range(401, 407)))[::-1]
        sec7 = range(407, 427)
        sec8 = list(range(437, 445)) + list(range(430, 426, -1))
        sec9 = list(range(451, 457))[::-1]

        SECS = [sec1, sec2, sec3, sec4, sec5, sec6, sec7, sec8, sec9]
        sec = SECS[section]

        # Finding stations of the chosen section
        Stations = []
        for st in range(len(sec)):
            s = "M82/1_" + str(sec[st])
            Stations.append(s)
            if st == 0:
                df_s = df[df["Event"] == s]
            else:
                df_s = pd.concat([df_s, df[df["Event"] == s]])

        # Creating variables array
        Depths = np.unique(np.around(df_s["Depth water [m]"], decimals=1))[::1]
        nx = len(Stations)
        nk = len(Depths)
        Tem = np.ones((nk, nx)) * -99999.9
        Sal = np.ones((nk, nx)) * -99999.9
        Time = []
        Lon = []
        Lat = []
        Baty = []

        for s in range(len(Stations)):
            DF = df_s[df_s["Event"] == Stations[s]]
            Time.append(np.unique(DF["Date/Time"])[0])
            Lon.append(np.unique(DF["Longitude"])[0])
            Lat.append(np.unique(DF["Latitude"])[0])
            Baty.append(-np.unique(DF["Elevation [m]"])[0])
            for z in range(len(Depths)):
                depth_z = DF["Depth water [m]"]
                df_z = DF[np.around(depth_z, decimals=1) == Depths[z]]
                if len(df_z) > 0:
                    if len(df_z) > 1:
                        df_T = df_z["Temp [°C]"]
                        df_S = df_z["Sal"]
                        Tem[z, s] = df_T[df_T.index[0]]
                        Sal[z, s] = df_S[df_S.index[0]]
                    else:
                        Tem[z, s] = df_z["Temp [°C]"]
                        Sal[z, s] = df_z["Sal"]
                elif Depths[z] <= depth_z[depth_z.index[-1]]:
                    Tem[z, s] = np.nan
                    Sal[z, s] = np.nan

        # Creating Dataset
        ds = Dataset(
            data_vars=dict(
                TEM=(["Depth", "station"], Tem),
                SAL=(["Depth", "station"], Sal),
            ),
            coords=dict(
                station=(["station"], sec),
                Depth=(["Depth"], Depths),
            ),
            attrs=dict(
                description="Section "
                + str(section + 1)
                + " of Quadfasel et al 2018 dataset"
            ),
        )

        # Dropping Nan and masking land
        ds = ds.dropna("Depth")
        ds = ds.where(ds > -99999.9, np.nan)

        # Adding Time, Topography, Lon and Lat variables
        daTime = DataArray(Time, coords=[sec], dims=["station"])
        daBaty = DataArray(Baty, coords=[sec], dims=["station"])
        daLon = DataArray(Lon, coords=[sec], dims=["station"])
        daLat = DataArray(Lat, coords=[sec], dims=["station"])

        ds["Time"] = daTime
        ds["Topography"] = daBaty
        ds["Longitude"] = daLon
        ds["Latitude"] = daLat

        # Manually add CF attributes
        attrs = dict(
            Depth={"standard_name": "depth"},
            station={"long_name": "station id"},
            Longitude={"standard_name": "longitude"},
            Latitude={"standard_name": "latitude"},
            SAL={"standard_name": "sea_water_practical_salinity"},
            TEM={"standard_name": "sea_water_temperature"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs)

        # Compute potential temperature
        ds["potential_temperature"] = compute_pt0(ds)

        return add_cf_attributes(ds)

    @property
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
        ds.attrs = {"featureType": "timeSeries"}

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

        return add_cf_attributes(ds)

    @property
    def S2020_secA(self) -> Dataset:
        """Standardized Section A of Semper et al 2020 dataset"""

        return self._S2020("A")

    @property
    def S2020_secB(self) -> Dataset:
        """Standardized Section B of Semper et al 2020 dataset"""

        return self._S2020("B")

    @property
    def S2020_secC(self) -> Dataset:
        """Standardized Section C of Semper et al 2020 dataset"""

        return self._S2020("C")

    @property
    def S2020_secD(self) -> Dataset:
        """Standardized Section D of Semper et al 2020 dataset"""

        return self._S2020("D")

    @property
    def S2020_secE(self) -> Dataset:
        """Standardized Section E of Semper et al 2020 dataset"""

        return self._S2020("E")

    def _S2020(self, sec) -> Dataset:
        """Standardized Semper et al 2020 dataset"""

        CTDfile = self.raw_pooch.fetch(
            "Semper_2020/CTD/Section" + sec + "_NEIceland_CTD.tab"
        )
        skiprow = 0
        a = open(CTDfile)
        for line in a:
            if line[0] == "*":
                skiprow += 1
                break
            skiprow += 1
        dfCTD = pd.read_csv(CTDfile, sep="\t", header=skiprow)

        # Creating variables array
        Stations = np.unique(dfCTD["Event"])
        Depths = -z_from_p(
            dfCTD["Press [dbar]"],
            dfCTD["Latitude"],
            dfCTD["Press [dbar]"] * 0.0,
            dfCTD["Press [dbar]"] * 0.0,
        )
        Depths = np.unique(np.around(Depths, decimals=1))[::1]
        nx = len(Stations)
        nk = len(Depths)
        Tem = np.ones((nk, nx)) * -99999.9
        Sal = np.ones((nk, nx)) * -99999.9
        Time = []
        Lon = []
        Lat = []

        for s in range(len(Stations)):
            DF = dfCTD[dfCTD["Event"] == Stations[s]]
            Time.append(np.unique(DF["Date/Time"])[0])
            Lon.append(np.unique(DF["Longitude"])[0])
            Lat.append(np.unique(DF["Latitude"])[0])
            for z in range(len(Depths)):
                depth_z = -z_from_p(
                    DF["Press [dbar]"],
                    DF["Latitude"],
                    DF["Press [dbar]"] * 0.0,
                    DF["Press [dbar]"] * 0.0,
                )

                df_z = DF[np.around(depth_z, decimals=1) == Depths[z]]
                if len(df_z) > 0:
                    if len(df_z) > 1:
                        df_T = df_z["Temp [°C]"]
                        df_S = df_z["Sal"]
                        Tem[z, s] = df_T[df_T.index[0]]
                        Sal[z, s] = df_S[df_S.index[0]]
                    else:
                        Tem[z, s] = df_z["Temp [°C]"]
                        Sal[z, s] = df_z["Sal"]
                elif Depths[z] <= depth_z[depth_z.index[-1]]:
                    Tem[z, s] = np.nan
                    Sal[z, s] = np.nan

        # Creating Dataset
        ds = Dataset(
            data_vars=dict(
                TEM=(["Depth", "station"], Tem),
                SAL=(["Depth", "station"], Sal),
            ),
            coords=dict(
                station=(["station"], range(nx)),
                Depth=(["Depth"], Depths),
            ),
            attrs=dict(description="Section " + sec + " of Semper et al 2020 dataset"),
        )

        # Dropping Nan and masking land
        ds = ds.dropna("Depth")
        ds = ds.where(ds > -99999.9, np.nan)

        # Adding Time, Topography, Lon and Lat variables
        daTime = DataArray(Time, coords=[range(nx)], dims=["station"])
        daLon = DataArray(Lon, coords=[range(nx)], dims=["station"])
        daLat = DataArray(Lat, coords=[range(nx)], dims=["station"])

        ds["Time"] = daTime
        ds["Longitude"] = daLon
        ds["Latitude"] = daLat

        # Manually add CF attributes
        attrs = dict(
            Depth={"standard_name": "depth"},
            station={"long_name": "station id"},
            Longitude={"standard_name": "longitude"},
            Latitude={"standard_name": "latitude"},
            SAL={"standard_name": "sea_water_practical_salinity"},
            TEM={"standard_name": "sea_water_temperature"},
        )
        ds = add_attributes_and_rename_variables(ds, attrs)

        # Compute potential temperature
        ds["potential_temperature"] = compute_pt0(ds)

        return add_cf_attributes(ds)
