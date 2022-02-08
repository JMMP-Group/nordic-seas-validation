from dataclasses import dataclass

import cf_xarray  # noqa: F401
import numpy as np
import xoak  # noqa: F401
from xarray import DataArray, Dataset


@dataclass
class SectionFinder:
    """
    Parameters
    ----------
    ds_domain: Dataset
        domain_cfg dataset
    """

    ds_domain: Dataset

    def __post_init__(self):
        self._grids = {}

    @property
    def grids(self) -> dict:
        """Dictionary mapping each grid to a dataset with its coordinates"""

        if not self._grids:
            for grid in ("u", "v", "t", "f"):
                ds = Dataset(
                    coords={
                        "lon": self.ds_domain[f"glam{grid}"],
                        "lat": self.ds_domain[f"gphi{grid}"],
                    }
                )
                ds = ds.squeeze(drop=True)
                for key, value in ds.sizes.items():
                    ds[f"{key}_index"] = DataArray(range(value), dims=key)
                self._grids[grid] = ds.cf.guess_coord_axis()

        return self._grids

    def nearest_neighbor(self, lons, lats, grid: str) -> Dataset:
        """
        Given the coordinates defining a section, find the nearest points
        on a model grid.

        Args:
            lons (1D array-like): Longitudes defining a section
            lats (1D array-like): Latitudes defining a section
            grid (string): Model grid `{"u", "v", "t", "f"}`

        Returns:
            Dataset
        """

        if not self.grids[grid].xoak.index:
            self._grids[grid].xoak.set_index(("lat", "lon"), "sklearn_geo_balltree")

        return self.grids[grid].xoak.sel(lat=DataArray(lats), lon=DataArray(lons))

    nearest_neighbour = nearest_neighbor

    def zigzag_section(self, lons, lats, grid: str) -> Dataset:
        """
        Given the coordinates defining a section, find the correspoinding zigzag section
        on a model grid.

        Args:
            lons (1D array-like): Longitudes defining a section
            lats (1D array-like): Latitudes defining a section
            grid (string): Model grid `{"u", "v", "t", "f"}`

        Returns:
            Dataset
        """

        def diff_and_inds_where_insert(ix, iy):

            dx, dy = (np.diff(ii) for ii in (ix, iy))
            inds = np.argwhere(np.abs(dx) + np.abs(dy) > 1).squeeze()

            return dx, dy, inds

        # Extract indexes
        ds = self.nearest_neighbor(lons, lats, grid)
        ix, iy = (ds[f"{i}_index"].data for i in ("x", "y"))

        # Remove duplicates
        mask = np.abs(np.diff(ix)) + np.abs(np.diff(iy)) == 0
        ix, iy = (np.delete(ii, np.argwhere(mask)) for ii in (ix, iy))

        # Initialize variables
        dx, dy, inds = diff_and_inds_where_insert(ix, iy)

        # Stop when inds is empty
        while inds.size:

            # Extract diffs
            dx, dy = (di[inds] for di in (dx, dy))

            # Insert new values
            # Single diagonal step: move along y direction first
            mask = np.abs(dx * dy) == 1
            ix = np.insert(ix, inds + 1, ix[inds] + (dx / 2).astype(int))
            iy = np.insert(
                iy, inds + 1, iy[inds] + np.where(mask, dy, (dy / 2).astype(int))
            )

            # Prepare for next iteration
            dx, dy, inds = diff_and_inds_where_insert(ix, iy)

        # Find new lat/lon
        ds = self.grids[grid].isel(
            x=DataArray(ix, dims=ds.dims), y=DataArray(iy, dims=ds.dims)
        )

        return self.nearest_neighbor(ds["lon"], ds["lat"], grid)

    def velocity_points_along_zigzag_section(self, lons, lats) -> dict:
        """
        Given the coordinates defining a section, find the corrisponding velocity points
        along a zigzag section (f-grid). Useful to compute accurate volume fluxes.

        Args:
            lons (1D array-like): Longitudes defining a section
            lats (1D array-like): Latitudes defining a section

        Returns:
            dict: Dictionary mapping u/v grids to their coordinates and indexes.
        """

        # ZigZag path along f
        ds = self.zigzag_section(lons, lats, "f")

        # Find dimension name
        dim = list(ds.dims)[0]
        ds[dim] = ds[dim]

        # Compute diff and find max index of each pair
        ds_diff = ds.diff(dim)
        ds_roll = ds.rolling({dim: 2}).max().dropna(dim)

        # Fill dictionary
        ds_dict = {}
        for grid in ("u", "v"):

            # Apply mask
            # Meridional: u - Zonal: v
            ds = ds_roll.where(
                ds_diff[f"{'y' if grid == 'u' else 'x'}_index"], drop=True
            )
            da_dim = ds[dim]

            if not ds.sizes[dim]:
                # Empty: either zonal or meridional
                continue

            # Find new lat/lon
            ds = self.grids[grid].isel(
                x=DataArray(ds["x_index"].astype(int), dims=dim),
                y=DataArray(ds["y_index"].astype(int), dims=dim),
            )
            ds = self.nearest_neighbor(ds["lon"], ds["lat"], grid)

            # Assign coordinate - useful to concatenate u and v after extraction
            ds_dict[grid] = ds.assign_coords({dim: da_dim - 1})

        return ds_dict
