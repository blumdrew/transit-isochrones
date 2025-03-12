# -*- coding: utf-8 -*-

"""
    Author: Andrew Lindstrom
    Date: 2025-03-03
    Purpose: Wrapper for Mapbox Isochrone API
"""

import os
import json
import requests
import logging

import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import geopandas as gpd
from shapely import Point, box

class Isochrone(object):
    """Wrapper for mapbox isochrone API + GTFS data to generate transit isochrones
    Also supports walking and driving isochrones
    """

    def __init__(
        self,
        gtfs_path: os.PathLike,
        mapbox_pk: str
    ) -> None:
        if os.path.isdir(gtfs_path):
            self.gtfs_path = gtfs_path
        else:
            raise FileNotFoundError("GTFS path must point to existing path")
        self.mapbox_pk = mapbox_pk
        self.mapbox_cache_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "data",
            "mapbox_cache"
        )
        os.makedirs(self.mapbox_cache_path, exist_ok=True)

        # gtfs specifics needed for sanity (mostly stop bounding box for api calls)
        self.stops = pd.read_csv(os.path.join(self.gtfs_path, "stops.txt"))
        try:
            self.calendar = pd.read_csv(os.path.join(self.gtfs_path, "calendar.txt"))
        except FileNotFoundError:
            self.calendar_dates = pd.read_csv(os.path.join(self.gtfs_path, "calendar_dates.txt"))
            self.calendar = self._handle_missing_calendar()
        self.stops_gdf = gpd.GeoDataFrame(
            data = self.stops, 
            geometry = [Point(xy) for xy in zip(self.stops["stop_lon"],self.stops["stop_lat"])]
        )
        self.stops_gdf.crs = "EPSG:4326"
        self.bbox = box(
            self.stops_gdf.bounds["minx"].min(),
            self.stops_gdf.bounds["miny"].min(),
            self.stops_gdf.bounds["maxx"].max(),
            self.stops_gdf.bounds["maxy"].max()
        )
        # forward declare
        self.trips = None
        self.stop_times = None

        
    ### internal methods
    # gtfs specific methods
    def _handle_missing_calendar(
        self
    ) -> pd.DataFrame:
        """Handle missing calendar by constructing it from calendar dates"""
        temp = self.calendar_dates.copy()
        temp["date_obj"] = pd.to_datetime(temp["date"], format="%Y%m%d")
        temp["monday"] = np.where(temp["date_obj"].dt.day_of_week == 0, 1, 0)
        temp["tuesday"] = np.where(temp["date_obj"].dt.day_of_week == 1, 1, 0)
        temp["wednesday"] = np.where(temp["date_obj"].dt.day_of_week == 2, 1, 0)
        temp["thursday"] = np.where(temp["date_obj"].dt.day_of_week == 3, 1, 0)
        temp["friday"] = np.where(temp["date_obj"].dt.day_of_week == 4, 1, 0)
        temp["saturday"] = np.where(temp["date_obj"].dt.day_of_week == 5, 1, 0)
        temp["sunday"] = np.where(temp["date_obj"].dt.day_of_week == 6, 1, 0)

        df = temp.groupby(
            by="service_id",
            as_index=False
        ).agg(
            {
                "monday":"max","tuesday":"max","wednesday":"max",
                "thursday":"max","friday":"max","saturday":"max","sunday":"max",
                "date":("min","max")
            }
        )
        df.columns = [
            "service_id","monday","tuesday","wednesday","thursday","friday",
            "saturday","sunday","start_date","end_date"
        ]
        return df

    def _calendar_from_day_type(
        self,
        day_type: str
    ) -> pd.DataFrame:
        # pass in day_type, return calendar df of those service_ids and s/e dates
        # handle case when calendar is 0 for all days... why do you do this TriMet
        all_cols = ["monday","tuesday","wednesday","thursday","friday","saturday","sunday"]
        if (self.calendar[all_cols] == 0).all().all():
            # we need to make self.calender from calendar dates..
            if self.calendar_dates is None or self.calendar_dates.empty:
                self.calendar_dates = pd.read_csv(os.path.join(self.path,"calendar_dates.txt"))
            
            self.calendar_dates["date"] = pd.to_datetime(
                self.calendar_dates["date"],
                format="%Y%m%d"
            )
            self.calendar_dates["day_of_week"] = self.calendar_dates["date"].dt.strftime("%A").str.lower()

            for dow in self.calendar_dates["day_of_week"].unique():
                self.calendar_dates[dow] = (self.calendar_dates["day_of_week"] == dow).astype(int)

            gf = self.calendar_dates.groupby(
                by=["service_id"],
                as_index=False
            ).agg(
                {
                    "monday":"sum",
                    "tuesday":"sum",
                    "wednesday":"sum",
                    "thursday":"sum",
                    "friday":"sum",
                    "saturday":"sum",
                    "sunday":"sum",
                    "date":("min","max"),
                }
            )
            gf.columns = [
                "service_id","monday","tuesday","wednesday",
                "thursday","friday","saturday","sunday",
                "start_date","end_date"
            ]
            gf[all_cols] = np.where(
                gf[all_cols] > 0,
                1,
                0
            )
            made_calender = True
        else:
            made_calender = False
        
        if day_type == "weekday":
            ref_cols = ["monday","tuesday","wednesday","thursday","friday"]
        elif day_type == "weekend":
            ref_cols = ["saturday","sunday"]
        else:
            ref_cols = [day_type]

        if not made_calender:
            return self.calendar[(self.calendar[ref_cols] == 1).all(axis=1)]
        else:
            return gf[(gf[ref_cols] == 1).all(axis=1)]
    
    # handle transfers in transit isochrone
    def _add_transfers(
        self,
        all_stops: pd.DataFrame,
        stop_times: pd.DataFrame,
        transfers_allowed: int = 1,
        max_transfer_walk_time: float = 2,
        walk_speed: float = 3,
        transfer_buffer_time: float = 1
    ) -> pd.DataFrame:
        """Internal method. Do not call outside of transit_isochrone construction
        """
        if transfers_allowed == 0:
            return all_stops
        ## stops in utm, euclidean buffer will have to do
        stops_utm = self.stops_gdf.to_crs(self.stops_gdf.estimate_utm_crs())
        buffer_dist = int(walk_speed * (max_transfer_walk_time / 60) * 1609.344)
        buffer_stops = gpd.GeoDataFrame(
            data = stops_utm["stop_id"],
            geometry = stops_utm.buffer(buffer_dist)
        )
        # add buffer geometry to all_stops
        buffer_stops["stop_id"] = buffer_stops["stop_id"].astype(str)
        all_stops["stop_id"] = all_stops["stop_id_stop_times"].astype(str)
        all_stops = all_stops.merge(
            buffer_stops,
            left_on="stop_id_stop_times",
            right_on="stop_id",
            how="inner",
            suffixes=["","_geom"]
        )
        # need for analysis
        stop_times["stop_id"] = stop_times["stop_id"].astype(str)
        sid_rid_did = stop_times[["stop_id","route_id","direction_id"]].drop_duplicates()

        full_data = [all_stops]
        while transfers_allowed > 0:
            transfer_data = []
            for idx in full_data[-1].index:
                row = full_data[-1].loc[idx]
                stops_to_check = stops_utm[stops_utm.within(row["geometry"])]
                row_stop_loc = row["geometry"].centroid
                if stops_to_check.empty:
                    continue
                #### TRANSFER LOGIC
                # first check - not the same route we are on already
                # no anti-direction transfers allowed
                stops_to_check["stop_id"] = stops_to_check["stop_id"].astype(str)
                sid_rid_did["stop_id"] = sid_rid_did["stop_id"].astype(str)
                stops_to_check = stops_to_check.merge(
                    sid_rid_did,
                    on="stop_id"
                )
                stops_to_check = stops_to_check[
                    stops_to_check["route_id"] != row["route_id"]
                ]
                if stops_to_check.empty:
                    logging.debug(f"Processed transfers for stop {idx} out of {full_data[-1].index.max()}: no new stops on different routes")
                    continue
                stops_to_check["distance_to_transfer_point"] = stops_to_check.distance(row_stop_loc)
                
                # reduce to nearest transfer point 
                stops_to_check_gf = stops_to_check.groupby(
                    by=["route_id","direction_id"],
                    as_index=False
                )[["distance_to_transfer_point"]].min()
                stops_to_check = stops_to_check.merge(
                    stops_to_check_gf,
                    on=["route_id","direction_id","distance_to_transfer_point"],
                    how="inner"
                )
                # merge next trip_id in stop_times...
                stops_to_check = stops_to_check.merge(
                    stop_times,
                    on=["stop_id","route_id","direction_id"],
                    how="inner",
                    suffixes=["","_stop_times"]
                )
                # filter to only transfers reachable after the bus arrives to initial point
                stops_to_check = stops_to_check[
                    stops_to_check["arrival_time_num"] > (
                        row["arrival_time_num_stop_times"]
                        + ((stops_to_check["distance_to_transfer_point"] / 1609.344) / (walk_speed / 60))
                        + transfer_buffer_time
                    )
                ]
                # filter to only the first trip after the transfer
                stops_to_check_gf = stops_to_check.groupby(
                    by=["stop_id","route_id","direction_id"],
                    as_index=False
                )[["arrival_time_num"]].min()
                stops_to_check = stops_to_check.merge(
                    stops_to_check_gf,
                    on=["stop_id","route_id","direction_id","arrival_time_num"],
                    how="inner"
                )
                stops_to_check["total_transfer_time"] = stops_to_check["arrival_time_num"] - row["arrival_time_num_stop_times"]
                stops_to_check["elapsed_time_transfer"] = stops_to_check["arrival_time_num"] - row["arrival_time_num"]
                stops_to_check["remaining_time_transfer"] = (
                    (row["elapsed_time"] + row["remaining_time"])
                    - stops_to_check["elapsed_time_transfer"]
                )
                # remove all transfers which are outside the time window
                stops_to_check = stops_to_check[stops_to_check["remaining_time_transfer"] > 0]
                if stops_to_check.empty:
                    logging.debug(f"Processed transfers for stop {idx} out of {full_data[-1].index.max()}: no transfers have enough time")
                    continue
                stops_to_check = stops_to_check[[
                    "stop_id","route_id","direction_id","trip_id","stop_lat","stop_lon",
                    "arrival_time","arrival_time_num","total_transfer_time","elapsed_time_transfer",
                    "remaining_time_transfer","geometry"
                ]]
                # stop_id of transfer to new route
                # route_id of route to be transferred to
                # direction_id of route to be transferred to
                # trip_id of trip to be transferred to
                # stop_lat/stop_lon of stop_id
                # arrival_time/arrival_time_num of trip leaving stop
                # elapsed_time_transfer: how much time has elapsed at start of the transfer trip
                # remaining_time_transfer: how much time remains at the start of the transfer trip
                #### END TRANSFER LOGIC
                #### TRAVEL ALONG TRANSFER LOGIC
                transfer_trips = stops_to_check.merge(
                    stop_times,
                    on="trip_id",
                    how="inner",
                    suffixes=["","_stop_times"]
                )
                # filter to only stops on or after transfer
                transfer_trips = transfer_trips[
                    transfer_trips["arrival_time_num_stop_times"] >= transfer_trips["arrival_time_num"]
                ]
                transfer_trips["elapsed_time"] = (
                    transfer_trips["arrival_time_num_stop_times"] # time transfered bus arrives at stop
                    - row["arrival_time_num"] # time transfered from bus left initial stop
                    + row["time_to_start"] # time from initial point to first bus
                )
                transfer_trips["remaining_time"] = (row["elapsed_time"] + row["remaining_time"]) - transfer_trips["elapsed_time"]
                transfer_trips = transfer_trips[transfer_trips["remaining_time"] > 0]
                
                if transfer_trips.empty:
                    logging.debug(f"Processed transfers for stop {idx} out of {full_data[-1].index.max()}: no legitimate transfers found")
                
                transfer_trips["remaining_distance"] = (
                    walk_speed * (transfer_trips["remaining_time"] / 60) * 1609.344
                ).astype(int).round(2)
                transfer_trips["remaining_distance"] = (transfer_trips["remaining_distance"] / 10).round().astype(int) * 10
                #### TRAVEL ALONG TRANSFER LOGIC END
                #### MAKE TRANSFER_TRIPS LOOKS LIKE ALL_TRIPS
                # need arrival_time and arrival_time_num to reference the FIRST stop!
                transfer_trips["arrival_time"] = row["arrival_time"]
                transfer_trips["arrival_time_num"] = row["arrival_time_num"]
                transfer_trips["distance_to_start"] = row["distance_to_start"]
                transfer_trips["time_to_start"] = row["time_to_start"]
                transfer_trips = transfer_trips[[
                    "stop_id","stop_lat","stop_lon","distance_to_start","time_to_start","route_id",
                    "direction_id","trip_id","arrival_time","arrival_time_num",
                    "stop_id_stop_times","arrival_time_stop_times","arrival_time_num_stop_times",
                    "elapsed_time","remaining_time","remaining_distance","total_transfer_time"
                ]]
                transfer_trips["transfer"] = True
                transfer_trips["transfer_time"] = transfer_trips["total_transfer_time"] + row["transfer_time"]
                transfer_trips.drop("total_transfer_time",axis=1,inplace=True)
                transfer_data.append(transfer_trips)
                logging.debug(f"Processed transfers for stop {idx} out of {full_data[-1].index.max()}")
            ##
            all_transfers = pd.concat(transfer_data, ignore_index=True)
            
            # take only the fasest transfer to each potential new stop id
            trips_to_keep = all_transfers.groupby(
                by=["stop_id_stop_times"]
            )[["remaining_time"]].max()
            all_transfers = all_transfers.merge(
                trips_to_keep,
                on=["stop_id_stop_times","remaining_time"]
            )
            # add logic buffers back in for next round of transfers
            all_transfers = all_transfers.merge(
                buffer_stops,
                left_on="stop_id_stop_times",
                right_on="stop_id",
                how="inner",
                suffixes=["","_geom"]
            )
            full_data.append(all_transfers)
            transfers_allowed -= 1
        
        final_data = pd.concat(full_data)
        final_data.drop(["stop_id_geom","geometry"],axis=1,inplace=True)
        
        # reduce to only stops with maximum time remaining and minimum transfer time
        # if not done separately, stops can be lost since one transfer may have more
        # time remaining AND a longer transfer time. But that would be fine, so do
        # time remaining first, transfer time second. Transfer time reduciton
        # needs to be done to eliminate cases where you could transfer at any number of
        # stops
        final_data_gf_1 = final_data.groupby(
            by="stop_id_stop_times",
            as_index=False
        )[["remaining_time"]].max()
        final_data = final_data.merge(
            final_data_gf_1,
            on=["stop_id_stop_times","remaining_time"]
        )
        final_data_gf_2 = final_data.groupby(
            by="stop_id_stop_times",
            as_index=False
        )[["transfer_time"]].min()
        final_data = final_data.merge(
            final_data_gf_2,
            on=["stop_id_stop_times","transfer_time"]
        )
        return final_data
    
    # static methods
    @staticmethod
    def time_str_to_num(
        time_str: str
    ) -> float:
        """Input: time_str in (24h) hh:mm:ss
        Output: number of minutes since 00:00:00
        """
        try:
            h, m, s = time_str.split(":")
        except Exception:
            return None
        return (int(h) * 60) + (int(m)) + (int(s) / 60)

    def read_stop_times(
        self,
        day_of_week: str
    ) -> pd.DataFrame:
        """Wrapper to get stop times only on a specific day or type of day"""
        calendar = self._calendar_from_day_type(day_of_week)
        # read trips
        if self.trips is None:
            trips = pd.read_csv(os.path.join(self.gtfs_path, "trips.txt"))
            self.trips = trips
        else:
            trips = self.trips
        trips = trips.merge(calendar, on="service_id", how="inner", suffixes=["","_calendar"])
        # read stop_times
        if self.stop_times is None:
            stop_times = pd.read_csv(os.path.join(self.gtfs_path, "stop_times.txt"))
            self.stop_times = stop_times
        else:
            stop_times = self.stop_times
        
        stop_times = stop_times.merge(trips, on="trip_id", how="inner", suffixes=["","_trip"])
        
        return_cols = [
            "stop_id","trip_id","route_id","direction_id",
            "shape_id","service_id","start_date","end_date",
            "arrival_time",
        ]
        return stop_times[return_cols]

    def walking_isochrone(
        self,
        distance: float,
        lat: float,
        lng: float,
        cache: bool = True
    ) -> gpd.GeoDataFrame:
        """Wrapper for mapbox API call"""
        # round lat/lng to 8 decimal points
        profile = "walking"
        lat = round(lat, 8)
        lng = round(lng, 8)
        # check that point is in box
        pt = Point(lng, lat)
        if not pt.within(self.bbox):
            logging.critical(f"lat/lng point {pt} is not within bounding box of stops in gtfs")
        # check distance is legitimate
        if distance < 1:
            distance = 10
        api_call = f"https://api.mapbox.com/isochrone/v1/mapbox/{profile}/{lng},{lat}?contours_meters={distance}&polygons=true&access_token={self.mapbox_pk}"
        cache_fname = f"{profile}_{lng}_{lat}_{distance}.geojson"
        cache_path = os.path.join(self.mapbox_cache_path, cache_fname)
        if os.path.isfile(cache_path):
            logging.info(f"Read data for {(distance, lat, lng, profile)} from cache")
            return gpd.read_file(cache_path)
        res = requests.get(api_call)
        if cache:
            with open(cache_path, 'w') as f:
                f.write(res.text)
        logging.info(f"Got data from mapbox API. Call: {(distance, lat, lng, profile)}")
        return gpd.read_file(cache_path)
    
    def transit_isochrone(
        self,
        lat: float,
        lng: float,
        time: float,
        transfers: int = 0,
        time_of_day: str = "16:00:00",
        day_of_week: str = "weekday",
        time_strictness: str = "loose",
        intial_walk_time: float = None,
        walk_speed: float = 3,
        dissolve: bool = True,
        **kwargs
    ) -> gpd.GeoDataFrame:
        """Calculate transit isochrone from a given point
        Args:
            lat (float): lattitude of starting point
            lng (float): longitude of starting point
            time (float): time (in minutes) of buffer area
        (Optional):
            transfers (int): number of transfers calculated. Defaults to 0
            time_of_day (str): HH:MM:SS string of time of day to calculate. Defaults to 08:00:00. 24h
            time_strictness (str): 'strict' or 'loose'. If 'strict', calculated based on leaving starting
                point at exactly time_of_day. If 'loose', optimize for first trip after starting time.
                Defaults to loose
            initial_walk_time (float): Maximum time to walk to first stop. If nothing is passed, will use
                the full time allocation to search. Defaults to None
            walk_speed (float): assumed walking speed in mph. Defaults to 3
            dissolve (bool): If data returned should be each polygon for each stop, or dissolved into one. Defaults to True
        (Transfer related options):
            max_transfer_walk_time (float): In minutes. Determines maximum distance away a transfer will be looked for. Defaults to 2
            transfer_buffer_time (float): In minutes. Determines how much padding is needed to make a transfer in time. Defaults to 1
        """
        max_transfer_walk_time = kwargs.get("max_transfer_walk_time", 2)
        transfer_buffer_time = kwargs.get("transfer_buffer_time", 1)
        # check for point being within bounding box
        pt = Point(lng, lat)
        if not pt.within(self.bbox):
            raise ValueError(f"lat/lng point {pt} must be within bounding box of stops in gtfs")
        ## calculate initial isochrone for stop placement
        if intial_walk_time is None:
            intial_walk_time = time
        initial_walk_distance = int(walk_speed * (intial_walk_time / 60) * 1609.344)
        logging.debug("Fetching initial stop area isochrone")
        initial_isochrone = self.walking_isochrone(initial_walk_distance, lat, lng)

        ## convert to utm
        area_utm = initial_isochrone.estimate_utm_crs()
        # will need this first isochrone later
        initial_isochrone_utm = initial_isochrone.to_crs(area_utm)
        stops = self.stops_gdf.to_crs(area_utm)

        logging.debug("Determining initial starting stops")
        initial_stops = stops[stops.within(initial_isochrone_utm.loc[0]["geometry"])]
    
        ## this must be the most expensive way to convert a point to a different proj
        start_gdf = gpd.GeoDataFrame(
            data = [1,],
            geometry = [Point(xy) for xy in zip([lng,], [lat,])],
            crs = "EPSG:4326"
        )
        start_gdf.to_crs(crs=area_utm,inplace=True)
        start_utm_point = start_gdf.loc[0]["geometry"]
        
        initial_stops["distance_to_start"] = initial_stops.distance(start_utm_point)
        initial_stops["time_to_start"] = (
            (initial_stops["distance_to_start"] / 1609.344)
            / (walk_speed / 60)
        )
        
        # now, add route
        stop_route = initial_stops.drop("geometry",axis=1)
        stop_times = self.read_stop_times(day_of_week)
        # reduce stop_times to only minimum by start_date
        stop_times_gf = stop_times.groupby(
            by=["stop_id","route_id","direction_id"],
            as_index=False
        )[["start_date"]].min()
        
        stop_times = stop_times.merge(
            stop_times_gf,
            on=["stop_id","route_id","direction_id","start_date"],
            how="inner"
        )
        stops_with_route = stop_times.drop(
            ["arrival_time","start_date","end_date","trip_id","shape_id"],
            axis=1
        ).drop_duplicates()
        stops_with_route["stop_id"] = stops_with_route["stop_id"].astype(str)
        stop_route["stop_id"] = stop_route["stop_id"].astype(str)
        stop_route = stop_route.merge(
            stops_with_route[["stop_id","route_id","direction_id"]],
            how="inner"
        )
        # closest stops by route by direction
        starting_stops = stop_route.groupby(
            by=["route_id","direction_id"],
            as_index=False
        )[["distance_to_start"]].min()
        starting_stops = stop_route.merge(
            starting_stops,
            on=["route_id","direction_id","distance_to_start"],
            how="inner"
        )
        stop_times["arrival_time_num"] = stop_times["arrival_time"].apply(self.time_str_to_num)
        # only keep stops after the starting time
        stop_times_near_tod = stop_times[stop_times["arrival_time_num"] >= self.time_str_to_num(time_of_day)]
        stop_times_near_tod["stop_id"] = stop_times_near_tod["stop_id"].astype(str)
        starting_stops_trip = starting_stops.merge(
            stop_times_near_tod,
            on=["stop_id","route_id","direction_id"],
            how="inner",
            suffixes=["","stop_times"]
        )
        # only starting stop trips that can be reached in (time + walk time)
        # only do this if time_strictness is "strict", since starting before is generally okay
        # for "loose"
        if time_strictness == "strict":
            starting_stops_trip = starting_stops_trip[
                starting_stops_trip["arrival_time_num"] 
                >= (starting_stops_trip["time_to_start"] + self.time_str_to_num(time_of_day))
            ]
        # reduce to first trip
        starting_stops_trip_gf = starting_stops_trip.groupby(
            by=["stop_id","route_id","direction_id"],
            as_index=False
        )[["arrival_time_num"]].min()
        starting_stops_trip = starting_stops_trip.merge(
            starting_stops_trip_gf,
            how="inner",
            on=["stop_id","route_id","direction_id","arrival_time_num"]
        )

        # reduce to only trips within [time] of start time
        starting_stops_trip = starting_stops_trip[
            starting_stops_trip["arrival_time_num"] 
            <= (self.time_str_to_num(time_of_day) + time)
        ]

        logging.debug("Determining en route stops on each route from starting stops")
        # merge all stops on first trip
        all_stops = starting_stops_trip.merge(
            stop_times_near_tod,
            on="trip_id",
            how="inner",
            suffixes=["","_stop_times"]
        )
        
        # calculate elapsed time since leaving
        if time_strictness == "loose":
            all_stops["elapsed_time"] = (
                all_stops["arrival_time_num_stop_times"] 
                - all_stops["arrival_time_num"]
                + all_stops["time_to_start"]
            )
        elif time_strictness == "strict":
            all_stops["elapsed_time"] = (
                all_stops["arrival_time_num_stop_times"]
                - self.time_str_to_num(time_of_day)
                # + all_stops["time_to_start"]
            )
        
        # only keep stops that are reached after time begins
        if time_strictness == "loose":
            all_stops = all_stops[all_stops["elapsed_time"] >= all_stops["time_to_start"]]
        else:
            all_stops = all_stops[
                all_stops["arrival_time_num_stop_times"] >= all_stops["arrival_time_num"]
            ]

        all_stops["remaining_time"] = time - all_stops["elapsed_time"]
        all_stops["remaining_distance"] = (
            walk_speed * (all_stops["remaining_time"] / 60) * 1609.344
        ).astype(int)
        # round to nearest 10m
        all_stops["remaining_distance"] = (all_stops["remaining_distance"] / 10).round().astype(int) * 10
    
        # only keep stops that are reached before time ends
        all_stops = all_stops[all_stops["remaining_time"].round(decimals=4) >= 0]
        # reduce to only maximum remaining time per ending stop id
        trips_to_keep = all_stops.groupby(
            by=["stop_id_stop_times"]
        )[["remaining_time"]].max()
        all_stops = all_stops.merge(
            trips_to_keep,
            on=["stop_id_stop_times","remaining_time"]
        )
        ## reduce to only columns we need to make life easier in transfer land
        all_stops = all_stops[[
            "stop_id","stop_lat","stop_lon","distance_to_start","time_to_start",
            "route_id","direction_id","trip_id","arrival_time",
            "arrival_time_num","stop_id_stop_times","arrival_time_stop_times",
            "arrival_time_num_stop_times","elapsed_time","remaining_time","remaining_distance"
        ]]
        all_stops["transfer"] = False
        all_stops["transfer_time"] = 0
        
        logging.debug("Adding transfers")
        all_stops = self._add_transfers(
            all_stops = all_stops, 
            stop_times = stop_times,
            transfers_allowed = transfers,
            walk_speed = walk_speed,
            max_transfer_walk_time = max_transfer_walk_time,
            transfer_buffer_time = transfer_buffer_time
        )
        
        # merge back stop lat lng for stop_id_stop_times
        all_stops["stop_id"] = all_stops["stop_id"].astype(str)
        st = self.stops.copy()
        st["stop_id"] = st["stop_id"].astype(str)
        all_stops = all_stops.merge(
            st[["stop_id","stop_lat","stop_lon"]],
            how="inner",
            left_on="stop_id_stop_times",
            right_on="stop_id",
            suffixes=["_start",""]
        )
        # make a bunch of mapbox api calls
        gdf_output = [
            {
                "type":"start",
                "stop_id":None,
                "lat":lat,
                "lon":lng,
                "distance_to_start":None,
                "elapsed_time":None,
                "remaining_time":intial_walk_time,
                "remaining_distance":None,
                "geometry":initial_isochrone.iloc[0]["geometry"]
            }
        ]
        for idx in all_stops.index:
            logging.debug(f"Processing mapbox data for item {idx} out of {all_stops.index.max()}")
            row = all_stops.loc[idx]
            res = self.walking_isochrone(
                distance = row["remaining_distance"],
                lat = row["stop_lat"],
                lng = row["stop_lon"]
            )
            gdf_output.append(
                {
                    "type":"bus_stop",
                    "stop_id":row["stop_id"],
                    "lat":row["stop_lat"],
                    "lon":row["stop_lon"],
                    "distance_to_start":row["distance_to_start"],
                    "elapsed_time":row["elapsed_time"],
                    "remaining_time":row["remaining_time"],
                    "remaining_distance":row["remaining_distance"],
                    "geometry":res.iloc[0]["geometry"]
                }
            )

        df = pd.DataFrame(gdf_output)
        gdf = gpd.GeoDataFrame(data = df, geometry="geometry", crs=initial_isochrone.crs)
        if dissolve:
            gdf = gdf.dissolve()
        return gdf

    def driving_isochrone(
        self,
        time: float,
        lat: float,
        lng: float,
        cache: bool = True
    ) -> gpd.GeoDataFrame:
        """Wrapper for mapbox API call"""
        # round lat/lng to 8 decimal points
        profile = "driving"
        lat = round(lat, 8)
        lng = round(lng, 8)
        # check that point is in box
        pt = Point(lng, lat)
        if not pt.within(self.bbox):
            logging.critical(f"lat/lng point {pt} is not within bounding box of stops in gtfs")
        
        api_call = f"https://api.mapbox.com/isochrone/v1/mapbox/{profile}/{lng},{lat}?contours_minutes={int(time)}&polygons=true&access_token={self.mapbox_pk}"
        cache_fname = f"{profile}_{lng}_{lat}_{int(time)}.geojson"
        cache_path = os.path.join(self.mapbox_cache_path, cache_fname)
        if os.path.isfile(cache_path):
            logging.info(f"Read data for {(int(time), lat, lng, profile)} from cache")
            return gpd.read_file(cache_path)
        res = requests.get(api_call)
        if cache:
            with open(cache_path, 'w') as f:
                f.write(res.text)
        logging.info(f"Got data from mapbox API. Call: {(int(time), lat, lng, profile)}")
        return gpd.read_file(cache_path)

