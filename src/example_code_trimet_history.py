# -*- coding: utf-8 -*-

"""
    Author: Andrew Lindstrom
    Date: 2025-03-18
    Purpose: historical trimet data comparisons
"""

import os
import json
import geopandas as gpd
import pandas as pd
import logging

from isochrone import Isochrone

def main():
    """Create isochrones for each regional center in Metro's 2040 Plan
    for 2013, 2019, and 2025 to explore how access to jobs on transit has changed over the years
    """
    logging.basicConfig(
        level=logging.INFO,
        filename=os.path.join(os.path.dirname(os.path.dirname(__file__)),"data","logs","trimet_history.log")
    )
    base_path = os.path.dirname(os.path.dirname(__file__))
    
    lat_lngs = {
        "downtown":(45.51888, -122.679311),
        "gresham":(45.502643, -122.427054),
        "gateway":(45.530643,-122.563609),
        "clackamas":(45.435535,-122.568143),
        "oregon_city":(45.359993,-122.604555),
        "washington_square":(45.452435,-122.778777),
        "beaverton":(45.491451,-122.801632),
        "tansabourne":(45.53474,-122.877708),
        "hillsboro":(45.52141,-122.985272)
    }
    # - if you are reusing this code, make sure to point it to where you have this data
    gtfs_paths = [
        os.path.join(base_path, "data", "trimet_2025"), # 2025 
        os.path.join(base_path, "data", "trimet_2019"), # 2019
        os.path.join(base_path, "data", "trimet_2013") # 2013
    ]
    # mapbox pk needs to be in the right spot :)
    with open(os.path.join(base_path, "mapbox.env"), "r") as f:
        pk_dict = json.load(f)
        mapbox_pk = pk_dict["mapbox_pk"]
    # some constants
    walk_speed = 3
    travel_time = 45
    time_of_day = "16:00:00"
    transfers = 1

    all_dfs = []
    for ltln_key in lat_lngs.keys():
        lat, lng = lat_lngs.get(ltln_key)
        ll_data = []
        for gtfs_path in gtfs_paths:
            year = gtfs_path[-4:]
            ic = Isochrone(gtfs_path=gtfs_path, mapbox_pk=mapbox_pk)
            gdf = ic.transit_isochrone(
                lat=lat,
                lng=lng,
                time=travel_time,
                transfers=transfers,
                time_of_day=time_of_day
            )
            gdf["center"] = ltln_key
            gdf["lat"] = lat
            gdf["lon"] = lng
            gdf["year"] = year
            gdf["travel_time"] = travel_time
            gdf.drop(
                ["type","stop_id","distance_to_start","elapsed_time","remaining_time","remaining_distance"],
                axis=1,
                inplace=True
            )
            gdf = gdf[["center","lat","lon","year","travel_time","geometry"]]
            ll_data.append(gdf)
            all_dfs.append(gdf)
        gdf_ll = gpd.GeoDataFrame(data=pd.concat(ll_data), geometry="geometry", crs=gdf.crs)
        print(gdf_ll)
        gdf_ll.to_file(os.path.join(base_path,"data",f"test_trimet_{ltln_key}.geojson"))
        

main()