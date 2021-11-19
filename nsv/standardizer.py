import os
from dataclasses import dataclass

import cf_xarray  # noqa: F401
import pandas as pd
import pooch
from scipy.io import loadmat
from xarray import DataArray, Dataset


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
        """Standardized kogur dataset"""

        # Open mat file
        filename = self.raw_pooch.fetch("Kogur/all_gridded.mat")
        mat = loadmat(filename, squeeze_me=True, mat_dtype=True)

        # Dimensions
        coord_names = (
            "tvec",
            "dvec",
            "xvec",
        )
        coords = {name: DataArray(mat[name], dims=name) for name in coord_names}

        # Variables
        var_names = (
            "ugrid",
            "vgrid",
            "ws",
            "Sfinal1",
            "Tfinal1",
            "PDfinal1",
            "wc",
        )
        variables = {
            name: DataArray(mat[name], dims=coord_names, coords=coords)
            for name in var_names
        }

        # Initialize dataset
        ds = Dataset(variables, coords=coords, attrs={"featureType": "timeSeries"})
        ds = ds.rename_dims(xvec="station")
        ds["station"] = ds["station"]

        # Convert matlab time
        ds["tvec"] = pd.to_datetime(ds["tvec"].values - 719529, unit="D")
        ds["tvec"] = ds["tvec"].dt.round("H")

        # Add coordinates
        mooring = Dataset(
            {var: DataArray(mat[var], dims="dist") for var in ("lon", "lat", "dist")}
        )
        mooring = mooring.rename(dist="xvec")
        interp = mooring.interp(xvec=ds["xvec"])
        for coord in ("lon", "lat"):
            interp[coord][-1] = mooring[coord][-1]
            ds = ds.assign_coords({coord: interp[coord]})

        # Manually add CF attributes
        attrs_dict = dict(
            tvec={"standard_name": "time", "long_name": "time vector for gridded data"},
            dvec={
                "standard_name": "depth",
                "long_name": "depth of gridded product",
                "positive": "down",
            },
            xvec={"long_name": "distance vector for gridded product"},
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
        for var, attrs in attrs_dict.items():
            ds[var].attrs = attrs
            if "standard_name" in attrs:
                ds = ds.rename({var: attrs["standard_name"].replace("sea_water_", "")})

        # Automagically add CF attributes
        ds = ds.rename(xvec="dist")
        ds = ds.cf.guess_coord_axis().cf.add_canonical_attributes()

        return ds
