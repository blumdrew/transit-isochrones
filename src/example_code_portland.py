# -*- coding: utf-8 -*-

"""
    Author: Andrew Lindstrom
    Date: 2025-03-11
    Purpose: Get isochrones for a set of places in Portland, Oregon
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
    ct_data_pth = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "data",
        "pdx_tracts",
        "pdx_tracts.shp"
    )
    ct_data = gpd.read_file(ct_data_pth)

    gtfs_data_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "data",
        "trimet"
    )
    with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "mapbox.env"), "r") as f:
        pk_dict = json.load(f)
        mapbox_pk = pk_dict["mapbox_pk"]
    ic = Isochrone(gtfs_path=gtfs_data_path, mapbox_pk=mapbox_pk)
    # constants
    walk_speed = 3
    wd = []
    td = []
    dd = []
    for data_idx in ct_data.index:
        data_pt = ct_data.loc[data_idx]
        lat = data_pt["geometry"].y
        lng = data_pt["geometry"].x
        print(f"Processing point {data_idx} out of {ct_data.index.max()}")
        # 20 minute walk isochrone
        distance = 10 *(int(walk_speed * (20 / 60) * 1609.344) // 10)
        output_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "data",
            "pdx_tracts_output",
            f"tract_{data_pt['NAME'].replace('.','_')}"
        )
        print(f"Fetching walk data...")
        try:
            walk_data = gpd.read_file(os.path.join(output_path, "walk_isochrone.geojson"))
        except Exception:
            walk_data = ic.walking_isochrone(
                distance=distance,
                lat=lat,
                lng=lng
            )
        print(f"Fetching transit data...")
        try:
            transit_data = gpd.read_file(os.path.join(output_path, "transit_isochrone.geojson"))
        except Exception:
            try:
                transit_data = ic.transit_isochrone(
                    lat=lat,
                    lng=lng,
                    time=25,
                    transfers=0
                )
            except Exception:
                print(f"Failed to fetch data for {lat, lng}")
        print(f"Fetching drive data...")
        drive_data = ic.driving_isochrone(
            time = 15,
            lat=lat,
            lng=lng
        )
        # add some bits that will show us what it is
        walk_data["NAME"] = data_pt["NAME"]
        walk_data["GEOID"] = data_pt["GEOID"]
        transit_data["NAME"] = data_pt["NAME"]
        transit_data["GEOID"] = data_pt["GEOID"]
        drive_data["NAME"] = data_pt["NAME"]
        drive_data["GEOID"] = data_pt["GEOID"]
        
        os.makedirs(output_path, exist_ok=True)
        walk_data.to_file(os.path.join(output_path, "walk_isochrone.geojson"))
        transit_data.to_file(os.path.join(output_path, "transit_isochrone.geojson"))
        drive_data.to_file(os.path.join(output_path, "drive_isochrone.geojson"))
        wd.append(walk_data)
        td.append(transit_data)
        dd.append(drive_data)
        
    w = pd.concat(wd)
    t = pd.concat(td)
    d = pd.concat(dd)

    w = gpd.GeoDataFrame(data=w, geometry="geometry")
    t = gpd.GeoDataFrame(data=t, geometry="geometry")
    d = gpd.GeoDataFrame(data=d, geometry="geometry")

    w.to_file(os.path.join(os.path.dirname(output_path), "full_walk_data.geojson"))
    t.to_file(os.path.join(os.path.dirname(output_path), "full_transit_data.geojson"))
    d.to_file(os.path.join(os.path.dirname(output_path), "full_drive_data.geojson"))
        

main()