#!/usr/bin/env python
from auvlib.data_tools import std_data, all_data
from auvlib.bathy_maps import draw_map
import sys

all_pings = all_data.all_mbes_ping.parse_folder(sys.argv[1]) # parse folder of gsf data
all_entries = all_data.all_nav_entry.parse_folder(sys.argv[1]) # parse folder of gsf data
mbes_pings = all_data.convert_matched_entries(all_pings, all_entries) # convert to std_data pings
# std_data.write_data(mbes_pings, "submaps.cereal")

# mbes_pings = std_data.mbes_ping.read_data("submaps.cereal")

d = draw_map.BathyMapImage(mbes_pings, 500, 500) # create a bathymetry height map
d.draw_height_map(mbes_pings) # draw the height map
d.draw_indices(mbes_pings, 1000)
d.draw_track(mbes_pings) # draw the track of the vehicle
d.write_image("height_map.png") # save the height map to "height_map.png"