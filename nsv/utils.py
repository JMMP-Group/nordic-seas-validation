import cf_xarray  # noqa: F401
from xarray import Dataset


def add_cf_attributes(ds: Dataset) -> Dataset:
    return ds.cf.guess_coord_axis().cf.add_canonical_attributes()


def add_attributes_and_rename_variables(ds: Dataset, attrs_dict: dict) -> Dataset:
    for var, attrs in attrs_dict.items():
        ds[var].attrs = {**ds[var].attrs, **attrs}
        if "standard_name" in attrs:
            ds = ds.rename(
                {
                    var: attrs["standard_name"]
                    .replace("sea_water_", "")
                    .replace("in_sea_water", "")
                }
            )

    return ds
