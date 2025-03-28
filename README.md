# transit-isochrones

Extension of [Mapbox Isochrone API](https://docs.mapbox.com/api/navigation/isochrone/) to public transit.

## What you need to get started

A free tier Mapbox account. Up to 100,000 API hits per month are allowed on the free tier. In the mapbox.env file in the base directory, place your public key. Reminder - do not publish your public key to GitHub directly, unless you want to pay for someone elses API hits.
The GTFS files of a public transit provider of reference. Place these in /data/gtfs/{service-provider-name}. If you are looking for historical GTFS files, here are a few options:
- [Mobility Database](https://mobilitydatabase.org/)
- [Transit Land](https://www.transit.land/)

## Using the code

The provided example code should give you a place to start. isochrone.py contains the core code, with the Isochrone class being main container.
The Isochrone class requires a mapbox public key (recommended storage location: .. mapbox.env) and one or more paths pointing to GTFS data. See "What you need to get started" for more information on how to get a key.
Additionally, download any number of GTFS files to use within the isochrone calculations to pass into the Isochrone class.

Call the `transit_isochrone` function on any lat/lng point. Required arguments:

- lat: lattitude of starting point
- lng: longitude of starting point
- time: time, in minutes, of isochrone buffer. Mapbox API has a maximum of value of 60, so that is also the maximum value here.
- time_of_day: time of day to calculate travel time at - HH:MM:SS

## Example outputs

