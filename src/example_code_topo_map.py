# -*- coding: utf-8 -*-

"""
    Author: Andrew Lindstrom
    Date: 2025-03-12
    Purpose: isochrones for public transit in Portland, OR to create topo-style map
"""


import os
import logging
import json

import geopandas as gpd
import pandas as pd

from isochrone import Isochrone

def main():
    """Read downloaded centroids of Portland, OR census tracts and determine a set of isochrones for them
    """
    os.makedirs(os.path.join(os.path.dirname(os.path.dirname(__file__)),"data","logs"),exist_ok=True)
    logging.basicConfig(
        level=logging.DEBUG,
        filename=os.path.join(os.path.dirname(os.path.dirname(__file__)),"data","logs",".log")
    )
    gtfs_data_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "data",
        "trimet"
    )
    with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "mapbox.env"), "r") as f:
        pk_dict = json.load(f)
        mapbox_pk = pk_dict["mapbox_pk"]
    ic = Isochrone(gtfs_path=gtfs_data_path, mapbox_pk=mapbox_pk)

    # pioneer square location
    lat = 45.518889
    lng = -122.679227

    lines = [5, 10, 15, 20, 25, 30, 35, 40, 45]

    all_data = []
    for line in lines:
        print(f"Fetching data for iso line {line}")
        iso_line = ic.transit_isochrone(
            lat = lat,
            lng = lng,
            time = line,
            transfers = 1
        )
        all_data.append(
            {
                "iso_line": line,
                "geometry": iso_line.iloc[0]["geometry"]
            }
        )
    
    df = pd.DataFrame(all_data)
    gdf = gpd.GeoDataFrame(data=df, geometry="geometry", crs=iso_line.crs)

    gdf.to_file(
        os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "data",
            "portland_iso_lines.geojson"
        )
    )


if __name__ == "__main__":
    main()