import numpy as np
import random
from scipy import spatial
from scipy.stats import gamma
import time

value = None
while True:
    value = input("Type: '0' for no neighbor check, '1' for face neighbor check, '2' for diagonal neighbor check:")
    try:
       value = int(value)
    except ValueError:
       print('Valid input, please')
       continue
    if value in {0,1,2}:
       break
    else:
       print('Invalid range, please: 0, 1, or 2')

neighbor_radius = None
if value == 1:
    neighbor_radius = 1.0
if value == 2:
    neighbor_radius = 1.8


start_time = time.clock()

# 60x60x30
# makes a 2d array with a list of xyz coordinates
coordinates = np.indices((60, 60, 30)).T.reshape(-1, 3)
""" alternatively: coordinates = np.mgrid[0:60, 0:60, 0:30].reshape(3, -1).T"""

# makes a 2d array with all the parameter values; int32 so it's all integer instead of floats
blank_parameters = np.zeros(shape=(len(coordinates), 6), dtype=np.int32)

# merges the coordinates and parameters array into a 2d array
# 3 values for coordinates, 6 values for each parameter
# order:x,y,z, PAG, Base, primary electron, secondary electron, acid, cleared
coordinate_list = np.column_stack((coordinates, blank_parameters))

def pag_and_base():
    for i in range(len(coordinate_list)):
        # 12.5% chance for
        if random.random() < 0.125:
            coordinate_list[i][3] = 1

        # 4% chance for
        if random.random() < 0.04:
            coordinate_list[i][4] = 1


# all coordinates where z coordinate = 0
surface_xy_coordinates = coordinate_list[np.in1d(coordinate_list[:, 2], 0)]
# x and y coordinates and primary electron parameter
column_xy_and_electrons = surface_xy_coordinates[:, :3]


def primary_electrons_in_column():
    for i in range(len(column_xy_and_electrons)):
        # x and y coordinate for each array value in the column_xy_and_electrons
        x = column_xy_and_electrons[i][0]
        y = column_xy_and_electrons[i][1]
        # center for the gaussian is at (30,30) in x,y
        # the function tells percentage of photons in demicals dropped in area
        fraction_of_photons = np.exp(-(((x - 30) ** 2) / 400 + ((y - 30) ** 2) / 400))
        # number of photons per nm^2 calculated to be 13.57 at center
        photons_in_xy_column = 13.57 * fraction_of_photons
        # number of photons absorbed, will round to nearest integer; Beer Lambert law  I=I0 exp(-Alpha *t);
        # it will be 1 - I0 exp(-Alpha *t) since we want to know how much is absorbed instead of how much brightness at the end
        expected_photons_absorbed_in_xy_column = photons_in_xy_column * (1 - np.exp(-5 * 0.03))
        # photons absorbed in the column through poisson, the input is the lamda
        actual_photons_absorbed_in_xy_column = round(np.random.poisson(expected_photons_absorbed_in_xy_column))
        # print(expected_photons_absorbed_in_xy_column)
        # stores the actual photons absorbed in the xy column in the array
        column_xy_and_electrons[i][2] = actual_photons_absorbed_in_xy_column


"""
#use this function to find the actual photons absorbed at the top of the center column
print(column_xy_and_electrons[1830])
"""


# list of electron values
def primary_electrons():
    for i in range(len(column_xy_and_electrons)):
        # electrons in the column is column_xy_and_electrons[i][2]
        # searches coordinate list for the blocks that are in the column
        # filter searches matching x coordiantes and matching y coordinates

        blocks_in_column = np.where(np.in1d(coordinate_list[:, 0], column_xy_and_electrons[i][0]) &
                                    np.in1d(coordinate_list[:, 1], column_xy_and_electrons[i][1]))[0]
        count = 0
        # puts electrons in the columns depending on how many were said to have been absorbed in the column
        while count < column_xy_and_electrons[i][2]:
            choosen_block = random.choice(blocks_in_column)
            coordinate_list[choosen_block, 5] += 1
            count += 1


pag_and_base()
print("done")
primary_electrons_in_column()
print("done")
primary_electrons()
print("done")

pag_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 3], [0])]
base_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 4], [0])]
primary_electron_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 5], [0])]

#random coordinates relative to distance of a coordinate, also has random angle/direction relative to the initial coordinate
def random_spherical(x_coordinate, y_coordinate, z_coordinate, radius, npoints):
    theta = 2 * np.pi * np.random.rand(npoints)
    phi = np.arccos(2 * np.random.rand(npoints) - 1)
    x = x_coordinate + radius * np.cos(theta) * np.sin(phi)
    y = y_coordinate + radius * np.sin(theta) * np.sin(phi)
    z = z_coordinate + radius * np.cos(phi)
    return np.vstack([x,y,z]).T


def secondary_electrons():
    #organizes the list of coordinates_containing secondary electrons
    rough_secondary_electron_coordinate_list = []

    #print(primary_electron_coordinates)
    for i in range(len(primary_electron_coordinates)):
        #we're assuming the electron travels at a distance gamma with alpha 10 and beta 0.45, but this did not curve fit well into the original data
        #we're also assuming the number of electrons generated is based on the poisson distribution of 5
        a = random_spherical(primary_electron_coordinates[i][0], primary_electron_coordinates[i][1],
                             primary_electron_coordinates[i][2], gamma.rvs(10, scale=0.45),
                             sum([np.random.poisson(5) for j in range(primary_electron_coordinates[i][5])]) + primary_electron_coordinates[i][5])
        # note replace 5 in the 5 + primary_electron_coordinates[i][5] to the poisson thing and the radius from 3 to poisson
        rough_secondary_electron_coordinate_list.append(a)

    rough_secondary_electron_coordinate_list = np.concatenate(np.asarray(rough_secondary_electron_coordinate_list))
    #print(rough_secondary_electron_coordinate_list)
    #print(len(rough_secondary_electron_coordinate_list))


    #finding closest integer coordinate in the meshgrid for the coordinate containing secondary electron
    meshgrid_coordinates = np.mgrid[-60:60, -60:60, -30:30].reshape(3, -1).T
    tree = spatial.cKDTree(meshgrid_coordinates)
    I = tree.query(rough_secondary_electron_coordinate_list)
    secondary_electron_containing_coordinates = meshgrid_coordinates[I[1], :]
    #deletes coordinates that have a z value below 0 or 30 or above

    #print(type(secondary_electron_containing_coordinates))
    #print(secondary_electron_containing_coordinates)

    secondary_electron_containing_coordinates = secondary_electron_containing_coordinates[
        (30 > secondary_electron_containing_coordinates[:, 2]) & (secondary_electron_containing_coordinates[:, 2] >= 0)]

    #print(secondary_electron_containing_coordinates)

    #the coordinate going over the 60x60 bounds in the xy comes back on the other side until it travels the full distance
    for i in range(len(secondary_electron_containing_coordinates)):

        if secondary_electron_containing_coordinates[i][0] >= 60:
            secondary_electron_containing_coordinates[i][0] = secondary_electron_containing_coordinates[i][0] - 60
        if secondary_electron_containing_coordinates[i][0] < 0:
            secondary_electron_containing_coordinates[i][0] = secondary_electron_containing_coordinates[i][0] + 60

        if secondary_electron_containing_coordinates[i][1] >= 60:
            secondary_electron_containing_coordinates[i][1] = secondary_electron_containing_coordinates[i][1] - 60
        if secondary_electron_containing_coordinates[i][1] < 0:
            secondary_electron_containing_coordinates[i][1] = secondary_electron_containing_coordinates[i][1] + 60

    #print(secondary_electron_containing_coordinates)

    for i in range(len(coordinate_list)):
        list_of_electrons_in_block = secondary_electron_containing_coordinates[
            np.in1d(secondary_electron_containing_coordinates[:, 0], coordinate_list[i][0]) &
            np.in1d(secondary_electron_containing_coordinates[:, 1], coordinate_list[i][1]) &
            np.in1d(secondary_electron_containing_coordinates[:, 2], coordinate_list[i][2])]
        coordinate_list[i][6] += len(list_of_electrons_in_block)
        #print(len(list_of_electrons_in_block))

secondary_electrons()

secondary_electron_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 6], [0])]
#print(secondary_electron_coordinates)

def acid_no_neighbors():
    for i in range(len(coordinate_list)):
        if coordinate_list[i][3] > 0 and coordinate_list[i][6] > 0:
            coordinate_list[i][3] -= 1
            coordinate_list[i][6] -= 1
            coordinate_list[i][7] += 1

def acid_face_neighbors():
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
        #print(np.take(coordinate_list, idx, axis=0))

        # gets all indexes of neighbors that have at least 1 PAG
        while coordinate_list[i][6] > 0:
            filtered_idx = []
            for j in range(len(idx)):
                if coordinate_list[idx[j]][3] > 0:
                    filtered_idx.append(idx[j])
            # all indexes of neighbors that have at least 1 PAG
            #print(filtered_idx)
            if len(filtered_idx) > 0:
                choosen_idx = random.choice(filtered_idx)
                coordinate_list[i][6] -= 1
                coordinate_list[choosen_idx][7] += 1
                coordinate_list[choosen_idx][3] -= 1
            else:
                break


acid_no_neighbors()
if value != 0:
    acid_face_neighbors()

#total number of acids in all coordinates combined
print(np.sum(coordinate_list[:, 7]))

acid_coordinates = coordinate_list[~np.in1d(coordinate_list[:, 7], [0])]

print(time.clock() - start_time, "seconds")

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def plot_3dfigure(data, color):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    ax.scatter(data[:,0],data[:,1],data[:,2], c=color, s=20, linewidths=None)

    plt.show()


data1 = np.asarray(pag_coordinates)
plot_3dfigure(data1, 'r')
data2 = np.asarray(base_coordinates)
plot_3dfigure(data2, 'b')
data3 = np.asarray(primary_electron_coordinates)
plot_3dfigure(data3, 'y')
data4 = np.asarray(secondary_electron_coordinates)
plot_3dfigure(data4, 'k')
data5 = np.asarray(acid_coordinates)
plot_3dfigure(data5, 'g')


#np.savetxt('something.csv', coordinate_list)

#np.save('coordinate_array_list', coordinate_list)
#np.save('coordinate_array_2', primary_electron_coordinates)
#a = np.load('outfile_name.npy')