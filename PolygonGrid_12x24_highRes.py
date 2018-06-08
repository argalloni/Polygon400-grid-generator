#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 2018

@author: Alessandro Galloni

This script generates a 12x24 grid of stimuli (228 patterns) to use as input 
for the Mightex Polygon400. 

Each pattern contains one stimulus spot, in the form of an array of 0s with 1s at the spot location.
Different spot locations are separated by one line of 0s. 
The size of the array thus determines the grid resolution and the relative sizes of the
light spots and spaces.

After generating an ordered array of stim patterns, the order of the patterns 
is randomly shuffled and the pattern with the greatest separation between 
subsequent patterns is chosen. This array is the exported to a .txt file that
can be read by the Polygon400.

"""
import numpy as np
import random as rd
import math
#import scipy as sp
#import scipy.spatial as spt
#import matplotlib
#import matplotlib.pyplot as plt
import time

np.set_printoptions(threshold=np.nan)

# ---------------------------------------------------------------
"""                   Grid parameters                         """
# ---------------------------------------------------------------

output_file = 'PolygonGrid_test'
stim_rows = 24
stim_columns = 12
grid_resolution = 3     # relative size of spots compared to the blank spaces between spots

all_stims = 1   # include pattern with all spots: Yes = 1, No = 0
separation = 1      #0=ordered, 1=maxseparated
start_pattern = 137    #start with stim in the middle of the grid (0-indexed)     
min_distance = 5
max_distance = 40 
iterations = 1000000


# ---------------------------------------------------------------
"""                        Functions                           """
# ---------------------------------------------------------------

#Calculate 1-indexed row and column of stim from 1-indexed pattern number
def ordered_gridposition(pattern_number): 
    row_position_1indexed = math.ceil(pattern_number/stim_columns)
    row_counter = row_position_1indexed-1
    column_position_1indexed = (pattern_number-row_counter*stim_columns)
    return(row_position_1indexed,column_position_1indexed)
    
#Calculate distance between two stim patterns 
def grid_distance(grid1,grid2):  
    index1 = np.where(grid1 == 1)
    index2 = np.where(grid2 == 1)
    rowdist = index1[0][0] - index2[0][0]
    coldist = index1[1][0] - index2[1][0]
    return(np.sqrt((rowdist)**2+(coldist)**2))

#Calcualte distance between two stim patterns, with the two ordered pattern numbers as input
def position_distance(pattern_number1,pattern_number2):   
    position1 = ordered_gridposition(pattern_number1)
    position2 = ordered_gridposition(pattern_number2)
    rowdist = position1[0] - position2[0]
    coldist = position1[1] - position2[1]
    pos_distance = np.sqrt((rowdist)**2+(coldist)**2)
    return(pos_distance)
    
def grid_distance_stats(grid):    #Calculate sum,avg,min of distances between adjacent patterns in a 3D grid
    distsum = 0 #initialize distance variable
    dist_vector = []
    for i in range(stim_number-1):   #For each pattern
        distsum += grid_distance(grid[:,:,i],grid[:,:,i+1]) 
        dist_vector = np.concatenate((dist_vector,grid_distance(grid[:,:,i],grid[:,:,i+1])),axis=0)
    min_dist = dist_vector.min()
    average_dist = dist_vector.mean()
    return(min_dist,average_dist)
     
def pos_distance_stats(grid_order):    #Calculate sum,avg,min of distances between adjacent patterns in a 3D  grid (0-indexed)
    distsum = 0 
    global dist_vector 
    dist_vector = []
    for i in range(stim_number-1):   #For each pattern
        distsum += position_distance(grid_order[i],grid_order[i+1]) 
        dist_vector.append(position_distance(grid_order[i],grid_order[i+1]))
    min_dist = np.min(dist_vector)
    average_dist = np.mean(dist_vector)
    return(min_dist,average_dist,dist_vector)

def reorder(grid,order):    #reorder patterns according to order vector
    global newgrid
    newgrid = np.zeros([len(grid[:,0,0]),len(grid[0,:,0]),len(grid[0,0,:])],dtype=int)
    counter = 0
    for i in order:
        newgrid[:,:,counter] = grid[:,:,i]
        counter += 1
    return(newgrid)
    
def matrixpos2gridpos(matrixpos):      #convert matrix position to ordered pattern position
    matrixpos -= 1 #Convert from 1-indexed to 0-indexed
    row_counter = math.floor(matrixpos/23)
    orderedpos = int((matrixpos-row_counter*11)/2)
    #orderedpos += 1 #Convert from 0-indexed to 1-indexed
    return(orderedpos)
    
def read_order(polygon_output):
    #Output order from polygon is 1-indexed
    #read_order("12x24grid_spaced_output.txt")
    order = list(np.loadtxt(polygon_output,dtype=int)[:,3])
    output_columns = list(np.loadtxt(polygon_output,dtype=int)[:,1])[0]
    grid_res = math.ceil(output_columns/12)
    global output_order
    output_order = []
    for i in range(288):
        order[i] -= 1 #Convert from 1-indexed to 0-indexed
        row_position = math.floor(order[i]/(grid_res*output_columns))
        column_position = (order[i]-(grid_res*output_columns*row_position))/grid_res
        gridpos = int(row_position*stim_columns + column_position)
        output_order.append(gridpos) 
    return(output_order)
    
    
# ---------------------------------------------------------------
"""                     Generate grid                         """
# ---------------------------------------------------------------

grid_block_size = grid_resolution + 1    #grid_block_size also includes the 1 blank space between each spot
row_number = grid_block_size*stim_rows - 1  #The -1 avoids the unnecessary last 0 in each row
column_number = grid_block_size*stim_columns - 1    #The -1 avoids the unnecessary bottom 0 in each column
stim_number = stim_rows * stim_columns
if all_stims == 0:
    ordered_grid = np.zeros([row_number,column_number,stim_number],dtype=int)
else:
    ordered_grid = np.zeros([row_number,column_number,stim_number+1],dtype=int)  

pattern_number = 0
#find top left pixel of each pattern spot:
for row in np.arange(0,len(ordered_grid[:,0,0]),grid_block_size):    #For each pattern row (skipping blank space rows)
    for column in np.arange(0,len(ordered_grid[0,:,0]),grid_block_size):  #For each pattern column (skipping blank space columns)
        #replace pixel values with 1s:
        for r in range(grid_resolution):
            for c in range(grid_resolution):
                ordered_grid[row+r,column+c,pattern_number] = 1  
        pattern_number += 1

    
#Reorder stims to maximise pattern separation. Iteratively create random vector such that euclidean distance >= min_distance
start = time.time()
random_order = [start_pattern]
ordered_list = list(np.arange(stim_number))
ordered_list.remove(start_pattern)
error_count = 0
counter = 0
for i in range(iterations):
    while len(ordered_list) > 0:
        next_position = rd.sample(ordered_list,1)[0]
        if min_distance < position_distance(random_order[counter],next_position) < max_distance: 
            random_order.append(next_position)
            ordered_list.remove(next_position)
            counter += 1
            error_count = 0 #reset error counter
        else:
            error_count += 1
            if error_count > len(ordered_grid):
                break   #if all numbers have been tried (avoid infinite loop)
    if len(ordered_list) == 0:
        print('Success!')
        break

end = time.time()
print(end - start)

if len(ordered_list)==0:
    maxseparated_order = random_order[:]
else:
    maxseparated_order = list(np.arange(stim_number))

        

# Remix patterns with smallest distance and largest distance to increase mindist
mindist = 0
for i in range(1000):
    mindist = max(mindist,pos_distance_stats(maxseparated_order)[0])
    smallest_distance_index = dist_vector.index(min(dist_vector))
    largest_distance_index = dist_vector.index(max(dist_vector))
    smallest_distance = maxseparated_order.pop(smallest_distance_index)
    maxseparated_order.insert(largest_distance_index,smallest_distance)
    if pos_distance_stats(maxseparated_order)[0] > mindist:
        print('Success!')
        break

print('mindist =',mindist)

#All orders are 0-indexed:
#maxseparated_order8 = [137, 18, 219, 134, 226, 125, 39, 146, 283, 60, 178, 31, 26, 168, 93, 154, 7, 123, 76, 5, 149, 278, 2, 185, 88, 29, 22, 186, 193, 109, 198, 68, 161, 243, 164, 79, 214, 64, 211, 55, 47, 210, 204, 113, 285, 33, 140, 196, 273, 277, 1, 225, 141, 238, 3, 220, 35, 218, 248, 57, 187, 73, 71, 41, 263, 235, 215, 121, 236, 252, 100, 270, 122, 256, 153, 15, 264, 224, 8, 61, 237, 56, 223, 274, 89, 267, 260, 119, 269, 130, 169, 67, 51, 180, 148, 208, 45, 160, 74, 30, 10, 151, 231, 52, 165, 233, 158, 247, 216, 266, 23, 242, 258, 156, 114, 244, 24, 4, 221, 103, 63, 202, 77, 234, 106, 14, 172, 107, 199, 286, 197, 150, 239, 138, 78, 131, 282, 162, 118, 49, 91, 254, 132, 246, 179, 44, 222, 117, 135, 203, 159, 19, 38, 48, 75, 272, 9, 188, 155, 58, 255, 86, 95, 99, 167, 25, 80, 127, 83, 85, 190, 253, 87, 177, 81, 230, 66, 128, 262, 104, 6, 142, 133, 175, 115, 209, 280, 152, 212, 21, 183, 90, 181, 227, 102, 0, 96, 276, 144, 69, 173, 143, 241, 126, 184, 189, 111, 250, 191, 62, 249, 170, 11, 182, 110, 82, 43, 171, 98, 166, 257, 46, 116, 240, 27, 124, 108, 53, 194, 105, 36, 174, 101, 279, 271, 65, 157, 206, 192, 32, 281, 72, 195, 200, 129, 16, 232, 17, 228, 136, 213, 112, 245, 97, 34, 37, 92, 40, 207, 94, 28, 145, 261, 50, 259, 251, 163, 13, 84, 201, 229, 12, 217, 139, 275, 147, 284, 205, 59, 268, 20, 176, 120, 287, 54, 70, 42, 265]
#maxseparated_order9 = [137, 213, 29, 240, 284, 10, 63, 144, 231, 158, 91, 234, 239, 13, 132, 4, 68, 119, 172, 78, 177, 259, 45, 139, 31, 82, 20, 260, 7, 72, 64, 190, 111, 236, 196, 79, 173, 247, 70, 102, 273, 159, 107, 229, 216, 3, 47, 141, 267, 61, 142, 147, 131, 270, 73, 8, 109, 189, 182, 263, 171, 251, 134, 200, 1, 98, 202, 96, 23, 246, 201, 254, 167, 262, 59, 252, 92, 207, 113, 49, 157, 232, 84, 258, 80, 135, 220, 120, 123, 106, 44, 186, 104, 99, 184, 41, 24, 281, 185, 127, 67, 222, 282, 206, 176, 194, 277, 165, 278, 180, 83, 6, 108, 279, 110, 38, 9, 145, 11, 168, 43, 255, 0, 218, 30, 166, 233, 28, 219, 52, 48, 105, 221, 95, 155, 223, 266, 101, 40, 274, 18, 193, 248, 33, 208, 275, 151, 51, 179, 118, 174, 228, 271, 138, 249, 241, 210, 21, 150, 203, 209, 264, 46, 242, 97, 66, 285, 204, 148, 116, 2, 133, 212, 86, 183, 39, 211, 130, 250, 197, 62, 199, 90, 276, 154, 170, 74, 56, 85, 26, 286, 112, 224, 128, 57, 124, 195, 117, 243, 175, 217, 122, 77, 230, 65, 205, 225, 287, 34, 245, 22, 93, 37, 163, 25, 42, 237, 32, 149, 15, 215, 115, 265, 156, 235, 192, 125, 178, 35, 257, 87, 129, 181, 88, 146, 269, 187, 121, 89, 272, 114, 60, 143, 16, 268, 161, 261, 12, 103, 256, 140, 50, 126, 27, 253, 55, 283, 188, 280, 191, 14, 54, 71, 17, 76, 226, 19, 94, 244, 214, 75, 160, 58, 153, 5, 69, 136, 53, 198, 227, 164, 169, 162, 81, 36, 238, 100, 152]
#maxseparated_order10 = [137, 167, 229, 61, 68, 140, 72, 189, 145, 210, 116, 175, 232, 18, 128, 10, 133, 205, 93, 287, 255, 23, 122, 285, 224, 67, 48, 29, 242, 69, 153, 262, 138, 6, 240, 209, 180, 71, 206, 114, 202, 99, 184, 81, 223, 152, 66, 199, 32, 151, 65, 166, 26, 190, 25, 148, 269, 38, 157, 143, 60, 146, 178, 160, 230, 129, 24, 188, 275, 22, 154, 98, 171, 243, 120, 134, 103, 241, 163, 270, 174, 181, 265, 50, 168, 198, 3, 258, 276, 40, 34, 277, 193, 149, 221, 97, 158, 279, 217, 212, 273, 200, 252, 73, 55, 177, 87, 56, 63, 108, 284, 194, 51, 187, 95, 183, 76, 235, 218, 111, 43, 176, 70, 253, 19, 231, 192, 89, 281, 135, 106, 215, 30, 130, 16, 274, 196, 2, 271, 227, 42, 144, 211, 268, 64, 150, 28, 219, 104, 256, 5, 84, 39, 127, 238, 53, 283, 86, 155, 172, 141, 182, 245, 170, 264, 208, 31, 14, 115, 96, 233, 109, 11, 207, 261, 36, 101, 83, 201, 88, 195, 80, 8, 126, 249, 100, 259, 62, 118, 282, 74, 15, 117, 280, 41, 165, 110, 267, 54, 213, 147, 225, 57, 220, 136, 156, 169, 58, 228, 17, 186, 52, 161, 244, 250, 191, 46, 254, 179, 92, 272, 266, 59, 185, 44, 107, 226, 1, 197, 123, 164, 78, 248, 13, 82, 21, 236, 27, 142, 125, 203, 246, 105, 286, 9, 132, 49, 214, 102, 121, 33, 131, 12, 234, 20, 119, 0, 216, 45, 85, 47, 75, 247, 4, 112, 204, 257, 37, 124, 94, 173, 79, 263, 7, 239, 90, 159, 222, 251, 113, 278, 162, 77, 35, 260, 139, 237, 91]
#maxseparated_order11 = [137, 22, 257, 51, 177, 279, 16, 268, 91, 267, 102, 191, 83, 200, 27, 46, 112, 5, 205, 269, 1, 36, 117, 247, 8, 161, 248, 74, 230, 111, 119, 162, 273, 37, 195, 127, 282, 228, 10, 245, 134, 229, 88, 31, 259, 135, 280, 97, 132, 263, 21, 94, 193, 81, 99, 179, 13, 175, 272, 124, 287, 242, 146, 25, 270, 49, 196, 141, 171, 35, 212, 9, 110, 42, 239, 120, 184, 214, 144, 234, 157, 249, 241, 176, 145, 54, 123, 17, 198, 116, 26, 125, 238, 57, 151, 68, 199, 275, 281, 190, 219, 78, 261, 158, 6, 189, 246, 168, 231, 19, 221, 84, 256, 4, 133, 164, 108, 122, 211, 71, 185, 69, 225, 33, 170, 2, 210, 28, 217, 235, 131, 271, 12, 152, 284, 140, 85, 255, 274, 147, 77, 159, 63, 59, 38, 264, 18, 243, 276, 48, 180, 138, 227, 169, 43, 61, 192, 92, 29, 101, 237, 107, 75, 70, 188, 11, 103, 250, 64, 118, 148, 202, 139, 285, 72, 142, 160, 87, 172, 115, 20, 50, 23, 62, 174, 156, 218, 106, 0, 286, 166, 209, 30, 260, 173, 96, 244, 126, 39, 222, 265, 104, 41, 181, 240, 73, 47, 15, 283, 82, 155, 252, 183, 216, 79, 197, 65, 251, 114, 204, 128, 226, 76, 194, 130, 44, 182, 201, 45, 206, 100, 258, 55, 215, 52, 24, 80, 167, 224, 253, 95, 121, 278, 7, 93, 203, 34, 186, 53, 223, 105, 232, 163, 254, 150, 277, 165, 233, 67, 207, 3, 236, 136, 178, 89, 153, 56, 109, 14, 149, 40, 213, 98, 187, 32, 143, 60, 90, 208, 113, 58, 129, 262, 266, 66, 220, 154, 86]
#maxseparated_order12 = [137, 22, 257, 51, 177, 16, 268, 91, 267, 102, 191, 83, 200, 27, 229, 46, 112, 5, 205, 1, 171, 214, 120, 269, 36, 117, 247, 8, 161, 248, 74, 230, 111, 273, 119, 279, 162, 37, 195, 31, 282, 127, 228, 10, 245, 134, 259, 135, 280, 97, 263, 132, 21, 193, 81, 281, 99, 184, 13, 175, 272, 124, 287, 145, 179, 88, 242, 146, 25, 270, 49, 196, 141, 35, 212, 9, 110, 190, 42, 239, 144, 234, 116, 249, 157, 94, 176, 54, 123, 231, 164, 17, 198, 26, 125, 238, 57, 151, 68, 275, 78, 261, 158, 6, 219, 108, 189, 61, 246, 29, 168, 19, 221, 84, 256, 4, 133, 241, 122, 211, 71, 185, 69, 225, 33, 170, 2, 210, 28, 217, 12, 235, 131, 271, 152, 284, 38, 140, 85, 255, 147, 274, 77, 159, 63, 276, 62, 172, 59, 264, 156, 18, 243, 48, 180, 138, 227, 87, 169, 43, 192, 92, 237, 101, 11, 107, 75, 188, 70, 250, 103, 199, 64, 118, 285, 148, 202, 72, 139, 50, 142, 20, 160, 23, 115, 218, 106, 0, 174, 286, 209, 30, 260, 173, 96, 244, 126, 39, 222, 104, 265, 153, 41, 181, 47, 240, 73, 166, 15, 216, 283, 82, 252, 155, 183, 178, 79, 197, 65, 251, 114, 204, 128, 226, 76, 194, 130, 182, 44, 201, 45, 206, 100, 258, 55, 215, 52, 224, 24, 80, 167, 90, 253, 95, 121, 278, 93, 7, 203, 34, 186, 53, 223, 105, 232, 163, 254, 150, 277, 89, 165, 233, 67, 207, 3, 236, 136, 56, 109, 14, 149, 40, 213, 98, 187, 32, 143, 60, 154, 208, 113, 58, 262, 129, 266, 66, 220, 86]
#maxseparated_order13 = [137, 22, 257, 51, 177, 16, 268, 91, 267, 102, 191, 83, 200, 27, 229, 46, 112, 5, 205, 1, 171, 214, 120, 24, 269, 36, 117, 247, 8, 161, 248, 74, 230, 111, 273, 119, 279, 162, 37, 195, 31, 282, 127, 228, 10, 245, 134, 259, 135, 280, 97, 263, 132, 21, 193, 81, 281, 99, 184, 13, 175, 272, 124, 287, 145, 179, 88, 242, 146, 25, 270, 49, 196, 141, 35, 212, 9, 110, 190, 42, 239, 144, 234, 116, 249, 157, 94, 176, 54, 123, 231, 164, 17, 198, 26, 125, 238, 57, 151, 68, 275, 78, 261, 158, 6, 219, 108, 189, 61, 246, 29, 168, 19, 221, 84, 256, 4, 133, 241, 122, 211, 71, 185, 69, 225, 33, 170, 2, 210, 28, 217, 12, 235, 131, 271, 152, 284, 38, 140, 85, 255, 147, 274, 77, 159, 63, 276, 62, 172, 59, 264, 156, 18, 243, 48, 180, 138, 227, 87, 169, 43, 192, 92, 237, 101, 11, 107, 75, 188, 70, 250, 103, 199, 64, 118, 285, 148, 202, 72, 139, 50, 142, 20, 160, 23, 115, 218, 106, 0, 174, 286, 209, 30, 260, 173, 96, 244, 126, 39, 222, 104, 265, 153, 41, 181, 47, 240, 73, 166, 15, 216, 283, 82, 252, 155, 183, 178, 79, 197, 65, 251, 114, 204, 128, 226, 76, 194, 130, 182, 44, 201, 45, 206, 100, 258, 55, 215, 52, 224, 80, 167, 90, 253, 95, 121, 278, 93, 7, 203, 34, 186, 53, 223, 105, 232, 163, 254, 150, 277, 89, 165, 233, 67, 207, 3, 236, 136, 56, 109, 14, 149, 40, 213, 98, 187, 32, 143, 60, 154, 208, 113, 58, 262, 129, 266, 66, 220, 86]
#maxseparated_order14 = [137, 22, 257, 51, 176, 16, 177, 268, 91, 267, 102, 191, 83, 200, 27, 229, 46, 112, 5, 205, 1, 171, 214, 120, 54, 180, 24, 232, 36, 117, 247, 8, 161, 248, 74, 230, 111, 273, 119, 279, 162, 37, 269, 31, 195, 282, 127, 228, 10, 245, 134, 259, 135, 280, 97, 263, 132, 21, 193, 81, 281, 99, 184, 13, 175, 272, 124, 287, 145, 179, 88, 242, 146, 25, 270, 49, 196, 141, 35, 212, 9, 110, 190, 42, 239, 144, 234, 116, 249, 157, 94, 123, 231, 164, 17, 198, 26, 125, 238, 57, 151, 68, 275, 78, 261, 158, 6, 219, 108, 189, 61, 246, 29, 168, 19, 221, 84, 256, 4, 133, 241, 122, 211, 71, 185, 69, 225, 33, 170, 2, 210, 28, 217, 105, 223, 12, 235, 131, 271, 152, 284, 38, 140, 85, 255, 147, 274, 77, 159, 63, 276, 62, 172, 59, 264, 156, 18, 243, 48, 138, 227, 87, 169, 43, 192, 92, 237, 101, 11, 107, 75, 188, 70, 250, 103, 199, 64, 118, 285, 148, 202, 72, 139, 50, 142, 20, 160, 23, 115, 218, 106, 0, 174, 286, 209, 30, 260, 173, 96, 244, 126, 39, 222, 104, 265, 153, 41, 181, 47, 240, 73, 166, 15, 216, 283, 82, 252, 155, 183, 178, 79, 197, 65, 251, 114, 204, 128, 226, 76, 194, 130, 182, 44, 201, 45, 206, 100, 258, 55, 215, 52, 224, 80, 167, 90, 253, 95, 121, 278, 93, 7, 203, 34, 186, 53, 163, 254, 150, 277, 89, 165, 233, 67, 207, 3, 236, 136, 56, 109, 14, 149, 40, 213, 98, 187, 32, 143, 60, 154, 208, 113, 58, 262, 129, 266, 66, 220, 86]
#maxseparated_order14p5 = [151, 22, 137, 257, 51, 176, 16, 177, 268, 91, 267, 102, 191, 83, 200, 27, 229, 46, 112, 5, 205, 1, 171, 214, 120, 54, 180, 24, 232, 36, 117, 247, 8, 161, 248, 74, 230, 111, 273, 119, 279, 162, 37, 113, 10, 269, 31, 195, 282, 127, 228, 57, 245, 134, 259, 135, 280, 97, 263, 132, 21, 193, 81, 281, 184, 13, 175, 272, 124, 287, 145, 179, 88, 242, 146, 25, 183, 155, 252, 143, 32, 270, 49, 141, 35, 212, 9, 110, 190, 42, 239, 144, 234, 116, 249, 157, 94, 123, 231, 164, 17, 198, 26, 125, 238, 68, 275, 78, 261, 158, 6, 219, 108, 189, 61, 246, 29, 168, 19, 221, 84, 256, 99, 4, 133, 241, 122, 211, 71, 185, 69, 225, 33, 170, 2, 210, 28, 217, 105, 223, 12, 235, 131, 271, 152, 284, 196, 38, 140, 85, 255, 147, 274, 77, 159, 63, 276, 62, 172, 59, 264, 156, 18, 243, 48, 138, 227, 87, 169, 43, 192, 92, 237, 101, 11, 107, 75, 188, 70, 250, 103, 199, 64, 118, 285, 148, 202, 72, 139, 50, 142, 20, 160, 23, 115, 218, 106, 0, 174, 286, 209, 30, 260, 173, 96, 244, 126, 39, 222, 104, 265, 153, 41, 181, 47, 240, 73, 166, 15, 216, 283, 82, 178, 79, 197, 65, 251, 114, 204, 128, 226, 76, 194, 130, 182, 44, 201, 45, 206, 100, 258, 55, 215, 52, 224, 80, 167, 90, 253, 95, 121, 278, 93, 7, 203, 34, 186, 53, 163, 254, 150, 277, 89, 165, 233, 67, 207, 3, 236, 136, 56, 109, 14, 149, 40, 213, 98, 187, 60, 154, 208, 58, 262, 129, 266, 66, 220, 86]
#maxseparated_order14p6 = [151, 22, 137, 257, 51, 176, 16, 177, 268, 91, 267, 102, 191, 83, 200, 27, 229, 46, 112, 5, 205, 1, 171, 214, 120, 54, 180, 159, 118, 199, 285, 77, 204, 89, 24, 114, 232, 32, 247, 117, 8, 161, 248, 74, 230, 111, 273, 119, 279, 113, 217, 150, 10, 165, 254, 37, 269, 162, 31, 195, 282, 127, 228, 57, 245, 134, 259, 135, 280, 97, 263, 132, 21, 193, 81, 281, 184, 13, 175, 272, 124, 287, 145, 179, 88, 242, 146, 25, 183, 155, 252, 143, 64, 270, 49, 141, 35, 212, 9, 110, 190, 42, 239, 144, 234, 116, 249, 157, 94, 123, 231, 164, 17, 198, 26, 125, 238, 68, 275, 78, 261, 158, 6, 219, 108, 189, 61, 246, 29, 168, 19, 221, 84, 256, 99, 4, 133, 241, 122, 211, 71, 185, 69, 225, 33, 170, 2, 210, 28, 223, 12, 105, 235, 131, 271, 152, 284, 196, 38, 140, 85, 255, 147, 274, 148, 276, 63, 36, 62, 172, 59, 264, 156, 18, 243, 48, 138, 227, 87, 169, 43, 192, 92, 237, 101, 11, 107, 75, 188, 70, 250, 103, 202, 72, 139, 50, 142, 20, 160, 23, 115, 218, 106, 0, 174, 286, 209, 30, 260, 173, 96, 244, 126, 39, 222, 104, 265, 153, 41, 181, 47, 240, 73, 166, 15, 216, 283, 82, 178, 79, 197, 65, 251, 128, 226, 76, 194, 130, 182, 44, 201, 45, 206, 100, 258, 55, 215, 52, 224, 80, 167, 90, 253, 95, 121, 278, 93, 7, 203, 34, 186, 53, 233, 67, 207, 3, 163, 277, 236, 136, 56, 109, 14, 149, 40, 213, 98, 187, 60, 154, 208, 58, 262, 129, 266, 66, 220, 86]
   
final_order = [137,22,257,51,177,16,268,91,267,102,191,83,200,27,229,46,112,5,205,1,171,269,36,117,247,8,161,248,74,230,111,273,119,279,162,37,195,282,127,228,10,245,134,259,135,280,97,263,132,21,193,81,281,99,184,13,175,272,124,287,145,179,88,242,146,25,270,49,196,31,141,35,212,9,110,190,42,239,120,214,144,234,116,249,157,94,176,54,123,231,164,17,198,26,125,238,57,151,68,275,78,261,158,6,219,108,189,61,246,29,168,19,221,84,256,4,133,241,122,211,71,185,69,225,33,170,2,210,28,217,12,235,131,271,152,284,38,140,85,255,147,274,77,159,63,276,62,172,59,264,156,18,243,48,180,138,227,87,169,43,192,92,237,101,11,107,75,188,70,250,103,199,64,118,285,148,202,72,139,50,142,20,160,23,115,218,106,0,209,286,30,260,173,96,244,126,39,222,104,265,153,41,181,47,240,73,174,15,283,82,252,155,183,166,216,178,79,197,65,251,114,204,128,226,76,194,130,182,44,201,45,206,100,258,55,215,52,224,24,80,167,90,253,95,121,278,93,7,203,34,186,53,223,105,232,163,254,150,277,89,165,233,67,207,3,236,136,56,109,14,149,40,213,98,187,32,143,60,154,208,113,58,262,129,266,66,220,86]
maxseparated_order = final_order

#Choose grid order that is most separated
if separation == 0:
        maxseparated_order = list(np.arange(stim_number))
        
reorder(ordered_grid, maxseparated_order)

if all_stims == 1:
    maxseparated_order.append(stim_number+1)
    for i in range(stim_number):
        newgrid[:,:,-1] += newgrid[:,:,i]

# ---------------------------------------------------------------
"""                     Generate grid                         """
# ---------------------------------------------------------------    

#Export matrix to a .txt file and append additional text to allow Polygon400 to generate a stimulus grid

gridstring = '''MightexVector1.0
# Mightex Vector file for 12 x 24 grid scan with space between stimuli.
# Minimum Distance = , Average Distance = 
# Type and BitDepth
Grid
1
# Columns and Rows\n''' +str(column_number) +'\n' + str(row_number)


for pattern in range(len(newgrid[0,0,:])):
    gridstring += ('\n#======== Pattern '+str(pattern+1)+', Row = '+str(ordered_gridposition(maxseparated_order[pattern]+1)[0])
                   +', Column = '+str(ordered_gridposition(maxseparated_order[pattern]+1)[1])+' ========\n')
    gridstring += 'Bin '
    patternstring = ''      # initialize pattern string
    for i in range(len(newgrid[:,0,0])):    #for every grid row
        rowstring = ' '.join(str(e)+'' for e in newgrid[i,:,pattern]) #add every colum value to the row
        if i>0:
            patternstring += '    '+ rowstring + ';\n'
        else:
            patternstring += rowstring + ';\n'
    gridstring += patternstring

counter = 0
filename = output_file
f = open(filename,'w')
f.write(gridstring)
f.close()