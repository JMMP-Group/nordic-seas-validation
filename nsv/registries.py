import os

import pooch

RAW = pooch.create(
    path=pooch.os_cache("NORVAL"),
    base_url="https://gws-access.jasmin.ac.uk/public/jmmp/NORVAL/",
    registry=None,
)
RAW.load_registry(os.path.join(os.path.dirname(__file__), "registry_raw.txt"))
