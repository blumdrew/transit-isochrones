# -*- coding: utf-8 -*-

"""
    Author: Andrew Lindstrom
    Date: 2025-03-11
    Purpose: hillsboro jobs
"""

import os
import logging
import json

import geopandas as gpd
import pandas as pd

from isochrone import Isochrone

def main():
    with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "mapbox.env"), "r") as f:
        pk_dict = json.load(f)
        mapbox_pk = pk_dict["mapbox_pk"]
    lat = 45.521384234830755
    lng = -122.9853161854452

    ic = Isochrone(
        gtfs_path="/Users/andrewmurp/Documents/python/transit-isochrones/data/trimet",
        mapbox_pk=mapbox_pk
    )

    hb = ic.transit_isochrone(lat=lat, lng=lng, transfers=1, time=45, time_of_day="08:00:00")

    hb.to_file("/Users/andrewmurp/Documents/python/transit-isochrones/data/hb.shp")

main()