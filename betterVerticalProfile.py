import matplotlib.pyplot as plt
import numpy as np

import csv
from geopy import distance
import math
import operator
import os
from pathlib import Path
import sys
from typing import List

LENGTH_OF_CHART = 20
SLICE_INTERVAL = 0.1  # miles
METERS_TO_MILE = 1 / 1609.344
METERS_TO_FEET = 3.28084
MILES_TO_FEET = 5280


def calc_grade(min_elevation, max_elevation) -> float:
    slice_in_feet = SLICE_INTERVAL * 5280

    return ((max_elevation - min_elevation) / slice_in_feet) * 100


class HikingPath:
    def __init__(self, name, points):
        self.name = name

        # the raw points, in lat/lon/elevation(m)
        self.points = points

        # the points for processing
        self.distances = []
        self.elevations = []
        self.max_elevation = float("-inf")
        self.min_elevation = float("inf")
        self.elevation_gain = 0
        self.elevation_loss = 0
        self.grades = []
        self.average_grade = 0

        self.calc_data()

    def get_length(self):
        if len(self.distances) == 0:
            return 0
        return self.distances[-1]

    def get_min_max_elevation(self) -> tuple:
        return int(self.min_elevation), int(self.max_elevation)

    def plot_hike(self, ax):
        ax.plot(self.distances, self.elevations)

    def calc_data(self):
        # save this data off in case we use it for multiple plots
        self.distances.clear()
        self.elevations.clear()
        sum_of_distances = 0
        for i in range(len(self.points)):
            if i == 0:
                self.distances.append(0)
            else:
                difference = self.points[i].distance_from(self.points[i - 1])
                sum_of_distances += difference
                self.distances.append(sum_of_distances)

            elevation = self.points[i].elevation * METERS_TO_FEET
            if elevation > self.max_elevation:
                self.max_elevation = elevation
                if sum_of_distances != 0:
                    self.average_grade = int((self.max_elevation - self.min_elevation) /
                                             (sum_of_distances * MILES_TO_FEET) * 100)
            if elevation < self.min_elevation:
                self.min_elevation = elevation
            self.elevations.append(self.points[i].elevation * METERS_TO_FEET)

        self.max_elevation = int(self.max_elevation)
        self.min_elevation = int(self.min_elevation)

        # Why divide into increments of 0.1 mi? It gave decently accurate numbers
        increment = len(self.elevations) / SLICE_INTERVAL
        skip = int(int(increment) / len(self.elevations))

        for elevation_index in range(skip, len(self.elevations), skip):
            # find elevation gain. Compare the maximum to the minimum within slice
            list_slice = self.elevations[elevation_index - skip:elevation_index]

            # measure gain for gain section, ignore loss as it is negligible
            if list_slice[-1] > list_slice[0]:
                max_elevation = max(list_slice)
                local_min = min(list_slice[:len(list_slice) - list_slice.index(max_elevation)])
                self.elevation_gain += int(max_elevation - local_min)
            else:
                min_elevation = min(list_slice)
                local_max = max(list_slice[:len(list_slice) - list_slice.index(min_elevation)])
                self.elevation_loss += int(local_max - min_elevation)


class PositionData:
    def __init__(self, lat, lon, elevation):
        # in degrees
        self.lat = float(lat)
        self.lon = float(lon)

        # in meters
        self.elevation = float(elevation)

    def distance_from(self, other):
        this = (self.lat, self.lon)
        other_loc = (other.lat, other.lon)
        two_d_distance = distance.distance(this, other_loc).miles
        return math.sqrt(math.pow(two_d_distance, 2) +
                         math.pow(math.fabs((self.elevation - other.elevation) * METERS_TO_MILE), 2))


def get_all_valid_files(path) -> List:
    ret_val = []
    for file in os.listdir(path):
        if file.endswith(".csv"):
            ret_val.append(file)
    return ret_val


def main():
    if len(sys.argv) != 2:
        print("--------------------------------- HELP ---------------------------")
        print("Run the program with the directory you want to run as the argument")
        print("The directory should have well formatted csv files in them")
        return

    path = Path(sys.argv[1])
    if not path.exists() or not path.is_dir() or len(os.listdir(path=str(path))) == 0:
        print("The path, " + str(path) + ", is not valid")
        return
    else:
        print("Using directory " + str(path))

    list_of_files = get_all_valid_files(str(path))

    if len(list_of_files) == 0:
        print("There are no valid files")

    font = {'size': 10}

    plt.rc('font', **font)

    trail_list = []

    for file in list_of_files:
        with open(str(path) + "/" + file, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            pos_list = []
            for row in reader:
                pos_list.append(PositionData(row['Latitude'], row['Longitude'], row['Elevation']))
            trail_list.append(HikingPath(file[:-4], pos_list))

    save_name = str(path) + "/"
    min_elevation = float("inf")
    max_elevation = float("-inf")
    max_length = 0
    for trail in trail_list:
        # find longest trail, min elevation, max elevation

        min_elevation = min(trail.min_elevation, min_elevation)
        max_elevation = max(trail.max_elevation, max_elevation)
        max_length = max(trail.get_length(), max_length)

        save_name = save_name + trail.name.replace(" ", "") + "vs"
    save_name = save_name[:-2]
    max_length_feet = max_length * MILES_TO_FEET
    inches_per_feet = LENGTH_OF_CHART / max_length_feet
    y_axis_size = (max_elevation * inches_per_feet) * 1.25
    plt.rcParams['toolbar'] = 'None'
    f, ax = plt.subplots(1, 1, figsize=(LENGTH_OF_CHART, y_axis_size))
    max_elevation_diff = max_elevation
    if min_elevation < 0:
        max_elevation_diff += min_elevation
    ratio = max_elevation_diff / max_length_feet
    ax.set_ylim(top=(max_elevation * 1.15))
    ax.set_xlim(0, max_length)

    # find the tick lengths. We want 5 ticks for vertical, and a tick every half mile for horizontal
    y_ticks = np.arange(0, max_elevation + 10, max_elevation / 5)
    plt.yticks(y_ticks)
    x_ticks = np.arange(0, max_length, 0.5)
    plt.xticks(x_ticks)
    trail_list.sort(key=operator.attrgetter('max_elevation'), reverse=True)
    legend_names = []
    for trail in trail_list:
        trail.plot_hike(ax)
        legend_names.append(trail.name)
    ax.legend(legend_names)
    ax.set_aspect(1.0 / ax.get_data_ratio() * ratio, adjustable='box-forced')
    plt.savefig(save_name + ".png")
    plt.show()
    return


if __name__ == '__main__':
    main()
