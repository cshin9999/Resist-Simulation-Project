import numpy as np
import random
from scipy import spatial
from scipy.stats import gamma
import time
import pandas as pd
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# ---------------------------------THIS IS TO FIND INTENSITY OVER DISTANCE---------------------------------------------
# combining intensity profile of both 28nm and 32nm and halfing it to get the intensity profile for 30nm;
def intensity_30nm_profile():
    # all the data from excel, the columns have 'nan's however because intensity28 has less distance data points
    df = pd.read_csv('30nm.csv')

    distance32 = df.values[:, 0] * 60
    intensity32 = df.values[:, 1]
    distance28 = df.values[:, 2] * 60
    intensity28 = df.values[:, 3]

    # gettings rid of nans
    distance28 = distance28[~np.isnan(distance28)]
    intensity28 = intensity28[~np.isnan(intensity28)]

    # interpolation functions based on intensity data
    interpolation_32nm = interp1d(distance32, intensity32, assume_sorted=False)
    interpolation_28nm = interp1d(distance28, intensity28, assume_sorted=False)

    # making empty array to be filled based on distance in 28nm; plug 28nm distance values into interpolated 32nm function
    intensity32_interpolated = np.zeros(len(distance28))

    # getting intensity of 32nm( y data) based on interpolation at each distance of 28nm points(the x data)
    for count, item in enumerate(distance28):
        intensity32_interpolated[count] = interpolation_32nm(item)

    # also getting interpolation function of 30nm
    intensity30 = (intensity28 + intensity32_interpolated) / 2
    global interpolation_30nm
    interpolation_30nm = interp1d(distance28, intensity30, assume_sorted=False)

    # print(1/interpolation_30nm(15)) # 17.0236 multiplier is there to make sure everything within 15nm radius gets clearing dose

    # shows the intensity profile of 28nm, 30nm, and 32nm
    """
    plt.scatter(distance28, intensity32_interpolated, c = 'red')
    plt.plot(distance28, intensity30, 'green')
    plt.scatter(distance28, intensity28, c = 'blue')
    plt.show()
    """

    # shows the intensity profile of 30nm
    """
    plt.plot(distance28, adjusted_interpolation_30nm(distance28), 'green')
    plt.show()
    """


def pag_and_base():
    for i in range(len(coordinate_list)):
        # 12.5% chance for
        if random.random() < acid_concentration:
            coordinate_list[i][3] = 1

        # 3.33% chance for
        if random.random() < base_concentration:
            coordinate_list[i][4] = 1


def primary_electrons():
    # all coordinates in coordinate_list where z coordinate = 0 (meaning top of photoresist)
    surface_xy_coordinates = coordinate_list[np.in1d(coordinate_list[:, 2], 0)]
    # number of primary electrons in the column
    column_xy_and_primary_electrons = surface_xy_coordinates[:, :3]

    total_number_of_photons = []
    for i in range(len(column_xy_and_primary_electrons)):
        # x and y coordinate for each array value in the column_xy_and_primary_electrons
        x = column_xy_and_primary_electrons[i][0]
        y = column_xy_and_primary_electrons[i][1]
        # try is there because over 30nm radius range, the intensity function has no data to go off on, interp1d will be out of bounds
        try:
            distance = np.sqrt((x ** 2) + (y ** 2))

            # the function tells percentage of photons in demicals dropped in area
            fraction_of_photons = interpolation_30nm(distance)

            # number of photons per nm^2 calculated to be 3.87 at center
            # 17.0236384626 is so that at edge of the 15nm radius, the photoresist blocks gets clearing dose
            photons_in_xy_column = photons_per_nm_squared * fraction_of_photons * 17.0236384626

            # number of photons absorbed, will round to nearest integer; Beer Lambert law  I=I0 exp(-Alpha *t);
            # it will be 1 - I0 exp(-Alpha *t) since we want to know how much is absorbed instead of how much brightness at the end
            expected_photons_absorbed_in_xy_column = photons_in_xy_column * (1 - np.exp(-4.9 * 0.03))

            # photons absorbed in the column through poisson, the input is the lambda
            actual_photons_absorbed_in_xy_column = round(np.random.poisson(expected_photons_absorbed_in_xy_column))

            # print(expected_photons_absorbed_in_xy_column)
            # stores the actual photons absorbed in the xy column in the array
            column_xy_and_primary_electrons[i][2] = actual_photons_absorbed_in_xy_column

            # this part is for printing number of primary electrons in the radius 15
            if distance <= 15:
                total_number_of_photons.append(column_xy_and_primary_electrons[i][2])

        except:
            pass

    for i in range(len(column_xy_and_primary_electrons)):
        # electrons in the column is column_xy_and_primary_electrons[i][2]
        # searches coordinate list for the blocks that are in the column
        # filter searches matching x coordiantes and matching y coordinates

        blocks_in_column = np.where(np.in1d(coordinate_list[:, 0], column_xy_and_primary_electrons[i][0]) &
                                    np.in1d(coordinate_list[:, 1], column_xy_and_primary_electrons[i][1]))[0]
        count = 0
        # puts electrons in the columns depending on how many were said to have been absorbed in the column
        while count < column_xy_and_primary_electrons[i][2]:
            choosen_block = random.choice(blocks_in_column)
            coordinate_list[choosen_block, 5] += 1
            count += 1


# random coordinates relative to distance of a coordinate, also has random angle/direction relative to the initial coordinate
def random_spherical(x_coordinate, y_coordinate, z_coordinate, radius, npoints):
    theta = 2 * np.pi * np.random.rand(npoints)
    phi = np.arccos(2 * np.random.rand(npoints) - 1)
    x = x_coordinate + radius * np.cos(theta) * np.sin(phi)
    y = y_coordinate + radius * np.sin(theta) * np.sin(phi)
    z = z_coordinate + radius * np.cos(phi)
    return np.vstack([x, y, z]).T

def secondary_electrons():
    # all coordinates that contain at least 1 primary electron
    primary_electron_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 5], [0])]

    # organizes the list of coordinates_containing secondary electrons
    secondary_electron_location_list = []

    for i in range(len(primary_electron_coordinates)):
        # we're assuming the electron travels at a distance probability distribution of a gamma with alpha 10 and beta 0.45, but this did not curve fit well into the original data
        # secondary electrons are generated; Also primary electron loses energy while moving about
        secondary_electron_location = random_spherical(primary_electron_coordinates[i][0],
                                                       primary_electron_coordinates[i][1],
                                                       primary_electron_coordinates[i][2],
                                                       gamma.rvs(10, scale=0.45),
                                                       sum([np.random.poisson(secondary_electrons_per_primary) for j
                                                            in
                                                            range(primary_electron_coordinates[i][5])]) +
                                                       primary_electron_coordinates[i][5])

        secondary_electron_location_list.append(secondary_electron_location)

    secondary_electron_location_list = np.concatenate(np.asarray(secondary_electron_location_list))

    # finding closest integer coordinate in the meshgrid for the coordinate containing secondary electron
    meshgrid_coordinates = np.mgrid[-60:60, -60:60, -30:30].reshape(3, -1).T
    tree = spatial.cKDTree(meshgrid_coordinates)
    I = tree.query(secondary_electron_location_list)
    secondary_electron_containing_coordinates = meshgrid_coordinates[I[1], :]

    # secondary electrons going lower than photoresist gets absorbed by substrate and leaves the system
    # secondary electrons going higher than photoresist flies off and leaves the system
    secondary_electron_containing_coordinates = secondary_electron_containing_coordinates[
        (0 >= secondary_electron_containing_coordinates[:, 2]) & (
            secondary_electron_containing_coordinates[:, 2] >= -29)]

    # the coordinate going over the 61x61 bounds in the xy comes back on the other side until it travels the full distance
    for i in range(len(secondary_electron_containing_coordinates)):

        if secondary_electron_containing_coordinates[i][0] >= 30:
            secondary_electron_containing_coordinates[i][0] = secondary_electron_containing_coordinates[i][0] - 30
        elif secondary_electron_containing_coordinates[i][0] <= -30:
            secondary_electron_containing_coordinates[i][0] = secondary_electron_containing_coordinates[i][0] + 30

        elif secondary_electron_containing_coordinates[i][1] >= 30:
            secondary_electron_containing_coordinates[i][1] = secondary_electron_containing_coordinates[i][1] - 30
        elif secondary_electron_containing_coordinates[i][1] <= -30:
            secondary_electron_containing_coordinates[i][1] = secondary_electron_containing_coordinates[i][1] + 30

    for i in range(len(coordinate_list)):
        list_of_electrons_in_block = secondary_electron_containing_coordinates[
            np.in1d(secondary_electron_containing_coordinates[:, 0], coordinate_list[i][0]) &
            np.in1d(secondary_electron_containing_coordinates[:, 1], coordinate_list[i][1]) &
            np.in1d(secondary_electron_containing_coordinates[:, 2], coordinate_list[i][2])]
        coordinate_list[i][6] += len(list_of_electrons_in_block)

def acid_no_neighbors():
    for i in range(len(coordinate_list)):
        if coordinate_list[i][3] > 0 and coordinate_list[i][6] > 0:
            coordinate_list[i][3] -= 1
            coordinate_list[i][6] -= 1
            coordinate_list[i][7] += 1

def acid_neighbors():
    coordinate_list_xyz = coordinate_list[:, 0:3]

    tree = spatial.cKDTree(coordinate_list_xyz, leafsize=100)

    for i in range(len(coordinate_list_xyz)):
        # idx is a list of indexes in the coordinate_list that are neighbors since both coordinate_list_xyz and
        # coordinate_list share same order in terms of indexes per xyz coordinate
        idx = tree.query_ball_point(coordinate_list_xyz[i], neighbor_radius)

        # list of indexes that are neighbors
        # print(idx)

        # prints xyz coordinates of neighbors
        # print(coordinate_list_xyz[idx])

        # prints the neighbor array
        # print(np.take(coordinate_list, idx, axis=0))

        # checks if voxel has secondary electrons--------------------------------------------------------------------
        # gets all indexes of neighbors that have at least 1 PAG
        while coordinate_list[i][6] > 0:
            filtered_idx = []
            for j in range(len(idx)):
                if coordinate_list[idx[j]][3] > 0:
                    filtered_idx.append(idx[j])
            # all indexes of neighbors that have at least 1 PAG
            # print(filtered_idx)
            if len(filtered_idx) > 0:
                choosen_idx = random.choice(filtered_idx)
                coordinate_list[i][6] -= 1
                coordinate_list[choosen_idx][7] += 1
                coordinate_list[choosen_idx][3] -= 1
            else:
                break

def developer_spread():
    coordinate_list_xyz = coordinate_list[:, 0:3]

    tree = spatial.cKDTree(coordinate_list_xyz, leafsize=100)

    for i in range(len(coordinate_list_xyz)):

        # idx is a list of indexes in the coordinate_list that are neighbors since both coordinate_list_xyz and
        # coordinate_list share same order in terms of indexes per xyz coordinate
        idx = tree.query_ball_point(coordinate_list_xyz[i], 1.0)

        # list of indexes that are neighbors
        # print(idx)

        # prints xyz coordinates of neighbors
        # print(coordinate_list_xyz[idx])

        # prints the neighbor array
        # print(np.take(coordinate_list, idx, axis=0))

        # checks if voxel has developer
        # gets all indexes of face neighbors and mark them as touching developer
        if coordinate_list[i][9] > 0:
            coordinate_list[idx, 8] = 1


def clear():
    coordinate_list_xyz = coordinate_list[:, 0:3]

    tree = spatial.cKDTree(coordinate_list_xyz, leafsize=100)

    for i in range(len(coordinate_list_xyz)):

        if coordinate_list[i][8] == 0:
            continue
        # idx is a list of indexes in the coordinate_list that are neighbors since both coordinate_list_xyz and
        # coordinate_list share same order in terms of indexes per xyz coordinate
        idx = tree.query_ball_point(coordinate_list_xyz[i], 1.8)

        # list of indexes that are neighbors
        # print(idx)

        # prints xyz coordinates of neighbors
        # print(coordinate_list_xyz[idx])

        # prints the neighbor array
        # print(np.take(coordinate_list, idx, axis=0))

        # checks if voxel is touching developer and more than 9 acids in the diagonal vicinity
        # if so, it now is cleared and has developer
        nearby_acid_sum = np.sum(coordinate_list[idx, 7])
        nearby_base_sum = np.sum(coordinate_list[idx, 4])
        nearby_developer_sum = np.sum(coordinate_list[idx, 9])

        if coordinate_list[i][8] > 0 and (
                    ((0.38 * nearby_acid_sum) - (0.38 * 0.6 * nearby_base_sum)) + (
                            nearby_developer_sum / 13.5)) >= 1.0:
            coordinate_list[i][8] = 0
            coordinate_list[i][9] = 1

def plot_3dfigure(data, color):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    ax.scatter(data[:, 0], data[:, 1], data[:, 2], c=color, s=20, linewidths=None)

    plt.show()

if __name__ == "__main__":
    # ------------------------------------------------------------------------------------------------------------------
    # CONSTANTS
    neighbor_radius = 1.0
    secondary_electrons_per_primary = 1.75
    acid_concentration = 0.2
    base_concentration = 0.033
    photons_per_nm_squared = 3.87
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # PHOTORESIST AND DEVELOPER LAYER CREATION

    # 61 by 61 by (30 + 1) # the +1 is the uppermost layer that is all developer and no photoresist
    # makes a 2d array with a list of xyz coordinates
    coordinates = np.mgrid[-30:31, -30:31, -29:2].reshape(3, -1).T

    # makes a 2d array with all the parameter values; int32 so it's all integer instead of floats
    blank_parameters = np.zeros(shape=(len(coordinates), 7), dtype=np.int32)

    # merges the coordinates and parameters array into a 2d array
    # 3 values for coordinates, 7 values for each parameter
    # order:x,y,z, PAG, Base, primary electron, secondary electron, acid, developer touching, developer
    coordinate_list = np.column_stack((coordinates, blank_parameters))

    developer_layer = coordinate_list[coordinate_list[:, 2] > 0, :]
    developer_layer[:, 9] = developer_layer[:, 9] + 1
    coordinate_list = coordinate_list[coordinate_list[:, 2] <= 0, :]
    # ------------------------------------------------------------------------------------------------------------------


    intensity_30nm_profile()

    pag_and_base()

    primary_electrons()

    secondary_electrons()

    acid_no_neighbors()

    acid_neighbors()

    acid_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 7], [0])]

    # merges the developer layer and wafer layers. Then turns the arrays into the original order
    coordinate_list = np.concatenate((coordinate_list, developer_layer), axis=0)
    proper_order_index = np.lexsort((coordinate_list[:, 2], coordinate_list[:, 1], coordinate_list[:, 0]))
    coordinate_list = coordinate_list[proper_order_index]

    previous_coordinate_list = None
    while (False == np.array_equal(previous_coordinate_list, coordinate_list)):
        previous_coordinate_list = np.copy(coordinate_list)
        developer_spread()
        clear()

        # print(previous_coordinate_list)

    developer_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 9], [0])]
    #removes the uppermost developer layer so only photoresist layers are present
    developer_coordinates = developer_coordinates[developer_coordinates[:, 2] < 1, :]



    """
    pag_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 3], [0])]
    data1 = np.asarray(pag_coordinates)
    plot_3dfigure(data1, 'r')
    
    base_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 4], [0])]
    data2 = np.asarray(base_coordinates)
    plot_3dfigure(data2, 'b')
    
    primary_electron_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 5], [0])]
    data3 = np.asarray(primary_electron_coordinates)
    plot_3dfigure(data3, 'y')
    
    secondary_electron_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 6], [0])]
    data4 = np.asarray(secondary_electron_coordinates)
    plot_3dfigure(data4, 'k')
    
    acid_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 7], [0])]
    data5 = np.asarray(acid_coordinates)
    plot_3dfigure(data5, 'g')
    """

    data6 = np.asarray(developer_coordinates)
    plot_3dfigure(data6, 'm')

    #3721 comes from 61 * 61
    print("Volume Cleared")
    print(np.sum(coordinate_list[:, 9]) - 3721)
