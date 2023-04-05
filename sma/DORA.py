#!/usr/bin/env python
# coding: utf-8

"""
@author: Jerry Wu
Last Update: 4.3.23
"""
# Note 1: This is the block with all the required imports

import scipy.stats as stats  # added to calculate z-score for Radius filtering
import mplcursors  # Jerry adds way to hover data points in matplotlib rather than plotly
import numpy as np
import pandas as pd
import os
import shutil
import plotly_express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib import ticker
import math
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from IPython import display
import random
import itertools
from matplotlib.widgets import Slider, Button
from IPython.core.display import HTML

# allows for large java HTML
matplotlib.rcParams['animation.embed_limit'] = 2**128 


# Note 1: DORA.find_center
# is a code that purely functions to run the centering algoritm on a csv of trajcetories.
# Inputs of this function are a LIST containing the variables:
# file_name, pixel_size, time_step, frame_start, frame_end, cmap, first_zero_end
# Ouptuts are a tuple center with the x and y coordinates of the center.
# Repeatedly run this code until you section of data to work with
'''
Template for the input parameters for DORA.find_center

# Define the below variables in the Variable Bank at the Top
# Relevant variables for this function to be  package into list:

# Initial Parameters is the Relevant Block

file_name = '00021.csv'  # 00086.csv or
pixel_size = 154  # in nanometers
time_step = 20  # miliseconds per frame in trajectory movie
frame_start = 0  # enter 0 to start from beginning of dataset
frame_end = 100  # enter -1 to end at the last value of the data set
cmap = "spring" # enter a color map string from this https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
exp_tag = "Glycerol_50_Laser_50" # a tag that caries the name of the experiment
first_zero_end = 'no'  # yes to cut off all values after first 0,0 = x,y
graph_centers = "no" #'yes' or 'no' to graphing the centers of the data.
save_plot = 'yes' 

# packaging all variables into list:
initial_parameters = [file_name, time_step, frame_start, frame_end, cmap, exp_tag, first_zero_end, graph_centers,save_plot]
center, data, ind_invalid_reading, data_back, my_rad_estimate = DORA.find_center(*initial_parameters)
'''


def find_center(*relevant_parameters):

    # unpackage our list of variables with their associated variables with variable names
    # unpack list into variables
    file_name, time_step, frame_start, frame_end, cmap, exp_tag, first_zero_end, graph_centers,save_plot = relevant_parameters

    # I will analyze the raw data from Ryan's code as pre_data and then covert that into two separate parts
    # 1) (data) the data formatted as arrays to be graphed [necessary numbers only]
    # 2) (data_back) the data formated as a Dataframe for record keeping [NaN placed where excluded values lie] (specifically they are excluded because of an invalid reading from Ryan's code)
    # read the csv file intended
    pre_data = pd.read_csv(file_name, header=None)
    # pre_data = pre_data.dropna()  # drop all NaN values?
    # create an array increasing in steps of
    pre_data['index'] = range(len(pre_data))
    pre_data.columns = ['X position',
                        'Y position', 'Intensity','index']  # label the columns
    pre_data = pre_data.iloc[:, [3, 0, 1 ,2]]  # reorganize the columns
    # create a boolean array of where 1s are when x position is 0 or invalid

    #section data from frame start to frame end
    pre_data = pre_data.iloc[frame_start:frame_end]

    # this is bc Ryan's code exports invalid readings as (0,0)
    ind_invalid_reading = pre_data['X position'] == 0

    if first_zero_end == 'yes':
        # run a boolean through the data set to find 0,0 points
        find_first_0 = (pre_data["X position"] == 0) & (
            pre_data["Y position"] == 0)
        # make an array can pre_x that holds all the x positions
        pre_x = pre_data["X position"].copy()
        # run the boolean through the pre_x variable, all 0,0 are true and stay.
        pre_x = pre_data[find_first_0]
        # the first 0,0 is the first item in this new pre_x
        my_first_0 = pre_x.index[0]
        frame_end = my_first_0  # set the index of the first 0,0 point as the end frame

    # SEPARATE data into front and back end (front==graphing ; back == tables)
    # if the index is not invalid (or valid) keep it and store in data
    data = pre_data[~ind_invalid_reading].copy()
    # section the pre data for all the invalid values
    data_back = pre_data[ind_invalid_reading].copy()

    # in data back develop a time colomn
    data_back['Time (ms)'] = data_back['index']*time_step

    # set all target x positions to NaN, if the reading was excluded due to invalid reading
    data_back['X position'] = np.nan
    # set all target y positions to NaN, if the reading was excluded due to invalid reading
    data_back['Y position'] = np.nan
    data_back['Excluded Type'] = 'Invalid Reading'

    ####################################### CENTERING ALGORITHM ###############################################

    # establish empty lists
    stand = []

    # find uniform guesses in range of max and min unaltered data values for y position
    # THE NUMBER OF UNIFORM GUESSES IS CURRENTLY HARD CODED AT 50 FOR X AND Y, CULMULITIVE 2,500
    guess_y = np.linspace(data.iloc[:, 2].max(
    ), data.iloc[:, 2].min(), 50)
    # put into list
    guess_y = guess_y.tolist()
    # find guesses for x position
    guess_x = np.linspace(data.iloc[:, 1].max(
    ), data.iloc[:, 1].min(), 50)
    guess_x = guess_x.tolist()

    # permute each x and y center guess together to create 10,000 unique center guesses
    center_guesses = list(itertools.product(guess_x, guess_y))
    # store center guesses in dataframe
    c = pd.DataFrame(center_guesses, columns=['X', 'Y'])
    # set up list to store average distances (radius) of circular trajectory path
    ave_distance = []
    # set up list to store standard deviation of distances to each point in the trajectory
    stand = []
    j = 0
    for j in range(len(c)):  # chnage to range(len(c))
        # find the distance between each point in a dataframe against guess[i]
        distance = np.power(
            ((data["X position"] - c['X'][j])**2 + (data["Y position"] - c['Y'][j])**2), 0.5)
        # store distances in a dataframe
        d = pd.DataFrame(distance, columns=['distance'])
        # find average of distances (this would be the radius)
        ave_d = d['distance'].mean(axis=0)
        # store all average distances from each guess[i] distance dataframes into list
        ave_distance.append(ave_d)
        # find standard deviation of center distance from each point in trajectory for each guess[i]
        std = d['distance'].std(axis=0)
        # store each standard deviation in a list
        stand.append(std)

        j += 1
    # put radius and std lists in a dataframe
    c['average_distance'] = ave_distance
    c['std'] = stand

    # this block finds the row with the lowest std, the corresponding radius and x,y coordinates for the center
    # want to return row with lowest std
    target_row = c['std'].idxmin()

    # x center guess with lowest std
    center_x = c.loc[target_row, 'X']
    # y center guess with lowest std
    center_y = c.loc[target_row, 'Y']
    # radius of trajectory
    dist = c.loc[target_row, 'average_distance']
    #Define dist as the estimated radius 
    my_rad_estimate = dist

    ###########GRAPHING BLOCK

    if graph_centers == "yes":

        # Our regularly scheduled 2D graphing program
        fig = plt.figure(figsize=(6, 6), dpi=100)
        # 121 # 1X1 grid plot 1, subplot(222) would be 2X2 grid plot 2, (223)--> 2X2 plot 3
        ax = fig.add_subplot(111)

        # Set up for color bar
        #HARD CODE:
        z_axis_label = "Frames" 

        #collect x and y values
        x =  data.iloc[:, 1] 
        y = data.iloc[:, 2]

        # A color bar associated with time needs two things c and cmap
        #these arguments go into ax.scatter as args

        # c (A scalar or sequence of n numbers to be mapped to colors using cmap and norm.)
        c = data["index"]

        #Make a ticks vector that spans the total number of frames
        # There is a bug because linspace doesn't understand what -1 is but the sequence does
        if frame_end == -1:   # negative 1
            last_frame = pre_data.iloc[:,0].iat[-1] #in the index column, give me the last valid value --> this is the max Frames
        else:
            last_frame = frame_end
        tix = np.linspace(frame_start,last_frame,8)
        tix_1 = np.round(tix,0)


        #scatter plot with a color vector
        p = ax.scatter(x, y, c=c, cmap = cmap,alpha=0.7)
        #add a vertical side bar that defines the color
        plt.colorbar(p, label=z_axis_label, shrink=.82, ticks = tix_1 )


        plt.axis('square')
        plt.xticks(rotation=45)

        # add a red dot to indicate center of trajectory
        ax.scatter(center_x, center_y, color='red')
        plt.text(x=center_x + 0.02, y=center_y + 0.02, s='algorithm centering')

        # add a red dot to indicate center of trajectory
        ax.scatter(center_x, center_y, color='red')
        plt.text(x=center_x + 0.02, y=center_y + 0.02, s='algorithm centering')

        # add a circle with center at our best guess and radius derived from our best guess
        circle = plt.Circle((center_x, center_y), dist, color='r', fill=False)
        ax.add_patch(circle)

        # Colorbar parameters below if we want one in the future

        # cbar = plt.colorbar(p, label= 'time' ,                asdfshrink= .82) #

        # #setting the ticks on the colorbar to span the length of the time column with 6 increments
        # cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1])

        # tix = np.linspace(0,len(data),6, dtype = int) # forces colorbar to show time in integers
        # tix_c = tix*20
        # cbar.set_ticklabels(tix_c)

        plt.axis('square')  # INTEGRAL to maintaining aspect ratio
        plt.xticks(rotation=45)
        ax.set_xlabel('X position (unaltered)', fontweight='bold', fontsize=14)
        ax.set_ylabel('Y position (unaltered)', fontweight='bold', fontsize=14)

        # plot title and font configurations

        # take the file name and separate from the extension
        # the first value in the tuple is the number
        # the second is .csv 
        # the number 00086.csv is the peak --> so this code takes the peak number
        pk = os.path.splitext(file_name)[0]

        graph_type = 'Algorithm_Center_Guess'

        # change title order!!! 
        list_of_strings = [graph_type, exp_tag]

        #in quotes is the the delimiter between the items in the string
        # by default it is a _ 
        my_title = "_".join(list_of_strings)

        plt.title(my_title, fontweight='bold', fontsize=16)

        plt.grid()

        if save_plot == "yes":
                # put title input and date time
                plt.savefig(my_title+"_centering.png")
                print("I have saved the picture with the name:")
                print(my_title+"_centering.png")
        plt.show()

    # #below is the key to maintaining the value of the center variable once the user satisfaction is achieved
    # global center

    center = (center_x, center_y)

    if graph_centers == "yes":
        print('The center is {0}'.format(center))

        print('If the center is satisfactory, change the find_center_coordinates parameter to no')
        print('If the center is unsatisfactory, adjust the frame_start and frame_end parameters and try again')


    return center, data, ind_invalid_reading, data_back, my_rad_estimate


# Note 2: DORA.downsample
# is a code that runs a down sampling method on the data
# it takes the following inputs:
# bin_size, processing, data, center, time_step, pixel_size, frame_start, frame_end
#      MOST are self explanatory, but for reference
#           "processing" is the type of processing to be run: 1)
#           "data" is the data from DORA.find_center


'''
# Relevant USER Input Parameters from the Variable Bank

bin_size = 20  # bin size for downsample/filter processing
processing = "none"  # enter downsample, moving average, or none

# DORA.down_sample(*downsample_parameters)
downsample_parameters = [bin_size, processing, data, center, time_step, pixel_size, frame_start, frame_end]
down_sampled_df = DORA.downsample(*downsample_parameters)

'''


def downsample(*downsample_parameters):
    # To hand our list input, we need to unpack it and redefine the values with the respective variable names
    bin_size, processing, data, center, time_step, pixel_size, frame_start, frame_end = downsample_parameters

    ### DATA NORMILIZATION AND UNIT ASSIGNMENT ###

    # # find the average of X and Y column respectively
    # ave = data.mean(axis=0) --> Jerry thinks this is antiquated code

    ## Step 1: center the data
    # substract averages from each column to find displacement, store into new columns
    data["X displacement (pixels)"] = data['X position'] - center[0]
    data["Y displacement (pixels)"] = data['Y position'] - center[1]
    ## Step 2: Convert pixels to nm
    data["X displacement (nm)"] = data['X displacement (pixels)']*pixel_size
    data["Y displacement (nm)"] = data['Y displacement (pixels)']*pixel_size
    ## Step 3; Create time step
    data["Time (ms)"] = data['index']*time_step

    # Recalculation of center using distance forumla -- Jerry
    # Radius Calculation from distance formula
    data['Radius (nm)'] = np.power(((data["X displacement (nm)"])
                                    ** 2 + (data["Y displacement (nm)"])**2), 0.5)

    # Z score calculation
    data['z-score Rad'] = stats.zscore(data["Radius (nm)"])

    # Angle Calculation

    # Radian to degree conversion factor
    r2d = 180/np.pi

    # Take Arc Tan function of x and y coord to get radius. Arctan 2 makes Quad 3 and 4 negative.
    data['Angle'] = -np.arctan2(data['Y displacement (nm)'],
                                data['X displacement (nm)'])*r2d

    # Make all negative Theta values positive equivalents
    data.loc[data.Angle < 0, ['Angle']] += 360


    ######################### PROCESSING BLOCK ##############################

    ### Moving Average ###

    # Simple Moving Average or "filter" dataframe:
    ma = pd.DataFrame(data.iloc[:, 0], columns=['index'])

    window = bin_size

    # for each column after the first (index) apply moving average filter
    for col in data.columns[1:]:
        ma[col] = data[col].rolling(window=window).mean()

    # Remove NaN's
    # In moving avg, all indices less than bin_size will be NaN 
    ma = ma.apply(pd.to_numeric, errors='coerce') 
    ma = ma.dropna()

    # Reset index
    ma = ma.reset_index(drop=True)

    #### Downsampling Dataframe ####

    # Create copy of data dataframe
    da = data.copy()

    # divide original index by sample size and round to nearest whole number to
    # achieve new index number underwhich the origial index is stored
    u = math.floor(frame_start/bin_size)
    v = math.floor(frame_end/bin_size)

    # define downsampling average function
    # This function taken from https://stackoverflow.com/questions/10847660/subsampling-averaging-over-a-numpy-array
    # allows us to downsample by averages over a set number
    # (change 'n' to the number of values you want to average over)

    def average_column(df, col_num, n):
        arr = df.iloc[:, col_num].values 
        end = n * int(len(arr)/n) 
        return np.mean(arr[:end].reshape(-1, n), 1)

    # Takes a column from our 'da' dataframe and runs the function over it
    # stores the new values in variables as an array (values in a row)


    # Iterate over columns of da except for the first column (index column)
    col_names = list(da.columns)
    averaged_cols = []
    for col_name in col_names:
        averaged_col = average_column(da, da.columns.get_loc(col_name), bin_size)
        averaged_cols.append(averaged_col)

    # Combine the averaged columns into a new dataframe
    dsa = pd.DataFrame(averaged_cols).T
    dsa.columns = da.columns

    # DETERMINE PROCESSING AND UNIT TYPE:
    # if more processing methods are to be added, an if statement must be
    # added with a key word to select that data frame
    # "df" becomes the variable used in the graphing block below
    if processing == "none":
        df = data
        return df
    if processing == "moving average":
        df = ma
        return df
    if processing == "downsample":
        df = dsa
        frame_start = math.floor(frame_start/bin_size)
        frame_end = math.floor(frame_end/bin_size)
        return df, frame_start, frame_end


# Note 3: DORA.graph(plot_type,*graph_parameters)
# This function graphs the following plots:
# Graphing options:
    # Trajectory Maps:
    # 2D: Colorful 2D visulization of the rotor from above
    # 3D: 2D plot but time is an axis

    # Grid plot
    # grid: a grid of little snippets of the data

    # Angular Analysis:

    #         By Jerry
    # radius_filter: Demarcate the excluded data points that will be eliminated from calculations
    # find_excluded_angle: Indicate excluded angles within angular_continuous by Jerry
    # angular_continuous_filtered: Angular Continuous recalculated with excluded points filtered. Excluded points skips indicated.
    # basal3: Graphs tailored for the basal graph analysis 3/14/2022
    # Angular Continuous with a downsampled curve as well. still has bugs with error labelling

    #         By Claire:  [NOT DONE]
    # angular: angle vs time, but it's not cummulative and resets at 360 to 0 (Claire)
    # angular_continuous: Claire's Calculation of a cummulative angle
    # find_excluded_angle_CR: Indicate erroneous angles within angular_continuous by Claire's calculations

    #          Animation   [DONE]
    # interactive: Interactive graph
    # animated: animated trajectory in notebook
    # HTML: Animated trajectory in a new window. May run better
'''
Template for DORA.graph(*INPUTS)

use the following for graph_parameters
#Graph Groupings:
# create a list of the acceptable groupings for the trajectory maps
trajectory_map = ["2D", "3D"]

# create a list of the acceptable groupings for the Angle Time grouping
AngleTime = ["radius_filter", "find_err_angle", "angular_continuous_filtered",
                "basal3", "angular", "angular_continuous", "find_err_angle_CR"]

# create a list of the acceptable groupings for the Animations Grouping

animations = ["interactive", "animated", "HTML"]

#Trajectory map parameters:
tajectory_map_parameters = [file_name, down_sampled_df, plot_type, display_center, expected_radius, x_axis_label, y_axis_label, z_axis_label, unit, 
pixel_min, pixel_max, axis_increment_nm, axis_increment_pixel, nm_min, nm_max, save_plot, frame_start, frame_end, time_step,cmap,exp_tag]

#Angle Versus Time (AVT or avt) parameters:
avt_parameters = [file_name, down_sampled_df, plot_type, display_center, ind_invalid_reading, rad_filter_type_upper,
                  rad_filter_type_lower, z_up, z_down, dist_high, dist_low, graph_style, bin_size, frame_start, frame_end,
                  display_center, exp_tag, x_axis_label, y_axis_label, z_axis_label, unit, pixel_min, pixel_max,
                  axis_increment_nm, axis_increment_pixel, nm_min, nm_max, save_plot, data_back, cmap, exp_tag] 

#Animated Parameters
animated_parameters = [file_name, down_sampled_df, plot_type, display_center, ind_invalid_reading, rad_filter_type_upper,
                  rad_filter_type_lower, z_up, z_down, dist_high, dist_low, graph_style, bin_size, frame_start, frame_end,
                  display_center, exp_tag, x_axis_label, y_axis_label, z_axis_label, unit, pixel_min, pixel_max,
                  axis_increment_nm, axis_increment_pixel, nm_min, nm_max, save_plot, data_back, cmap, exp_tag, frame_speed, tail_length] 

#Grid Parameters
grid_parameters = [file_name, down_sampled_df, plot_type, display_center, exp_tag, x_axis_label, y_axis_label, z_axis_label, unit, 
pixel_min, pixel_max, axis_increment_nm, axis_increment_pixel, nm_min, nm_max, save_plot, frame_start, frame_end, time_step,cmap,exp_tag, 
frames_per_plot, columns, fig_size_x, fig_size_y]


# #DORA.graph(plot_type,*relevant_parameters)

if plot_type == "2D" or plot_type == "3D":
    DORA.graph(plot_type,*tajectory_map_parameters)
if plot_type in animations:
    %matplotlib notebook
    DORA.graph(plot_type,*animated_parameters)
if plot_type == "grid":
    DORA.graph(plot_type, *grid_parameters)
if plot_type in AngleTime:
    DORA.graph(plot_type,*avt_parameters)
'''


def graph(plot_type, *graph_parameters):

    ###################### Which graph will be graphed?##########

    # graph groupings:

    # create a list of the acceptable groupings for the trajectory maps
    trajectory_map = ["2D", "3D"]

    # create a list of the acceptable groupings for the Angle Time grouping
    AngleTime = ["radius_filter", "find_excluded_angle", "angular_continuous_filtered",
                 "basal3", "angular", "angular_continuous", "find_excluded_angle_CR"]

    # create a list of the acceptable groupings for the Animations Grouping

    animations = ["interactive", "animated", "HTML","plotly_animated"]

    if plot_type in trajectory_map:

        #### Set up block#####

        # Accept my variables from graphing parameters
        [file_name, down_sampled_df, plot_type, display_center, expected_radius, x_axis_label, y_axis_label, z_axis_label, unit, 
        pixel_min, pixel_max, axis_increment_nm, axis_increment_pixel, nm_min, nm_max, save_plot, frame_start, frame_end, time_step,cmap,exp_tag] = graph_parameters

        # Claire's code accepts down_sampled_df as df
        df = down_sampled_df

        #####################Graphing data assignment block##############
        # Here the code determines the units of the graph, only for cartesian graphs
        if unit == "pixel":
            x = df["X displacement (pixels)"]
            y = df["Y displacement (pixels)"]
        if unit == "nm":
            x = df["X displacement (nm)"]
            y = df["Y displacement (nm)"]
        z = df["Time (ms)"]

        # graph either
        if plot_type == "2D":

            # Let the graphing beghin!!!
            fig = plt.figure(figsize=(6, 6), dpi=100)
            # this comand is here to take advantage of the "axes" plotting library
            ax = fig.add_subplot(111)

            # Set up for color bar
            #HARD CODE:
            z_axis_label = "Frames" 

            # A color bar associated with time needs two things c and cmap
            #these arguments go into ax.scatter as args

            # c (A scalar or sequence of n numbers to be mapped to colors using cmap and norm.)
            c = df["index"]

            #Make a ticks vector that spans the total number of frames
            # There is a bug because linspace doesn't understand what -1 is but the sequence does
            if frame_end == -1:   # negative 1
                last_frame = df["index"].iat[-1] #in the index column, give me the last valid value --> this is the max Frames
            else:
                last_frame = frame_end
            #tix = np.linspace(frame_start,last_frame,9)
            #tix_1 = np.round(tix,0)
            frame_step=int((last_frame-frame_start)/5)
    
            tix_1=np.arange(frame_start,last_frame,frame_step)

            #scatter plot with a color vector
            p = ax.scatter(x, y, c=c, cmap = cmap,alpha=0.7)
            #add a vertical side bar that defines the color
            plt.colorbar(p, label=z_axis_label, shrink=.82, ticks = tix_1 )


            plt.axis('square')
            plt.xticks(rotation=45)
            circle2 = plt.Circle((0, 0), expected_radius, color='m', fill=False)
            
            ax.add_patch(circle2)

            # display center
            if display_center == "yes":
                # in a centered graph, the center is actually(0,0)
                center1 = [0, 0]
                # plots center point as magenta X
                ax.scatter(0, 0, color='Magenta', marker="X", s=150)
                plt.text(x=center1[0] + 0.02,
                         y=center1[1] + 0.02, s='CENTER')

            # set graph limit conditions depending on unit specified
            if unit == "pixel":
                ax.set_xlim(pixel_min, pixel_max)
                ax.set_ylim(pixel_min, pixel_max)
                ax.yaxis.set_major_locator(
                    ticker.LinearLocator(axis_increment_pixel))
                ax.xaxis.set_major_locator(
                    ticker.LinearLocator(axis_increment_pixel))
            if unit == "nm":
                ax.set_xlim(nm_min, nm_max)
                ax.set_ylim(nm_min, nm_max)
                ax.yaxis.set_major_locator(
                    ticker.LinearLocator(axis_increment_nm))
                ax.xaxis.set_major_locator(
                    ticker.LinearLocator(axis_increment_nm))

            # Jerry Adds a hover cursor
            mplcursors.cursor(hover=True)
            mplcursors.cursor(highlight=True)

            # axis labels and font configurations
            ax.set_xlabel(x_axis_label, fontweight='bold', fontsize=14)
            ax.set_ylabel(y_axis_label, fontweight='bold', fontsize=14)

            # plot title and font configurations

            # take the file name and separate from the extension
            # the first value in the tuple is the number
            # the second is .csv 
            # the number 00086.csv is the peak --> so this code takes the peak number
            pk = os.path.splitext(file_name)[0]

            graph_type = '2D_Map'

            # change title order!!! 
            list_of_strings = [graph_type, exp_tag]

            #in quotes is the the delimiter between the items in the string
            # by default it is a _ 
            my_title = "_".join(list_of_strings)

            plt.title(my_title, fontweight='bold', fontsize=16)

            if save_plot == "yes":
                # put title input and date time
                plt.savefig(my_title+"_2D.png")
        if plot_type == "3D":
            # This block splices the segments between data points and assigns each segment to a color
            points = np.array([x, y, z]).transpose().reshape(-1, 1, 3)
            segs = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = Line3DCollection(segs, cmap=plt.get_cmap('cool'))
            lc.set_array(z)

            # This block plots the figure at a specified size, in 3D configuration, sets axis range, gathers the
            # colored segments from above, and assigns labels
            fig = plt.figure(figsize=(8, 8))
            ax = fig.gca(projection='3d')
            ax.set_zlim(min(z), max(z))

            # define graphing proportions according to unit used
            if unit == "pixel":
                ax.set_xlim(-1, 1)
                ax.set_ylim(-1, 1)
            if unit == "nm":
                ax.set_xlim(-150, 150)
                ax.set_ylim(-150, 150)

            ax.add_collection3d(lc, zs=z, zdir='z')

            # define title
            
            # take the file name and separate from the extension
            # the first value in the tuple is the number
            # the second is .csv 
            # the number 00086.csv is the peak --> so this code takes the peak number
            pk = os.path.splitext(file_name)[0]

            graph_type = '3D_Map'

            # change title order!!! 
            list_of_strings = [graph_type, exp_tag]

            #in quotes is the the delimiter between the items in the string
            # by default it is a _ 
            my_title = "_".join(list_of_strings)

            plt.title(my_title, fontweight='bold', fontsize=16)

            # set labels
            ax.set_xlabel(x_axis_label, fontweight='bold', fontsize=14)
            ax.set_ylabel(y_axis_label, fontweight='bold', fontsize=14)
            ax.set_zlabel(z_axis_label, fontweight='bold', fontsize=14)

            if save_plot == 'yes':
                plt.savefig(my_title+'_3D.png', dpi=300)

            plt.show()

    if plot_type in AngleTime:

        # Set up Block

        print("!entered angle time block!")
        # Accept my variables from graphing parameters
        (file_name, down_sampled_df, plot_type, display_center, ind_invalid_reading, rad_filter_type_upper,
         rad_filter_type_lower, z_up, z_down, dist_high, dist_low, graph_style, bin_size, frame_start, frame_end,
         display_center, title, x_axis_label, y_axis_label, z_axis_label, unit, pixel_min, pixel_max,
         axis_increment_nm, axis_increment_pixel, nm_min, nm_max, save_plot, data_back, cmap, exp_tag) = graph_parameters

        # import the 2 ways to analyze angular data: 1) Claire's way, 2) Jerry's Way [More Current]
        import AngleCalc

        # Gather inputs to Cacluate Angle under Jerry's Angle calculation paradigm
        inputs_avt_filter = [down_sampled_df, ind_invalid_reading, rad_filter_type_upper,
                             rad_filter_type_lower, z_up, z_down, dist_high, dist_low, bin_size]
        # Call and Run Jerry's Angle Calculation
        data, xy_goodbad, avt_good, avt_bad, data_fil_dsa, data_fil_down_bad, data_fil_up_bad = AngleCalc.avt_filter(
            *inputs_avt_filter)

        print("!Angle Calc completed!f")


        # Claire's code accepts down_sampled_df as df
        df = down_sampled_df

        #####################Graphing data assignment block##############
        # Here the code determines the units of the graph, only for cartesian graphs
        if unit == "pixel":
            x = df["X displacement (pixels)"]
            y = df["Y displacement (pixels)"]
        if unit == "nm":
            x = df["X displacement (nm)"]
            y = df["Y displacement (nm)"]
        z = df["Time (ms)"]
        
        # unpack list xy_goodbad into respective outputs
        x_good, y_good, x_bad, y_bad = xy_goodbad

        # Calculate The angle according to Claire's method
        df_theta = AngleCalc.avt_no_filter(
            down_sampled_df, frame_start, frame_end)

        
        t, theta, thetac = AngleCalc.avt_no_filter(
            down_sampled_df, frame_start, frame_end)

        # Let the graphing begin:

        if plot_type == "radius_filter":
            fig = plt.figure(figsize=(7, 6), dpi=100)
            # this comand is here to take advantage of the "axes" plotting library
            ax = fig.add_subplot(111)

            # Set up for color bar
            #HARD CODE:
            z_axis_label = "Frames" 

            # A color bar associated with time needs two things c and cmap
            #these arguments go into ax.scatter as args

            # c (A scalar or sequence of n numbers to be mapped to colors using cmap and norm.)
            c = df["index"]

            #Make a ticks vector that spans the total number of frames
            # There is a bug because linspace doesn't understand what -1 is but the sequence does
            if frame_end == -1:   # negative 1
                last_frame = df["index"].iat[-1] #in the index column, give me the last valid value --> this is the max Frames
            else:
                last_frame = frame_end
            #tix = np.linspace(frame_start,last_frame,9)
            #tix_1 = np.round(tix,0)
            frame_step=int((last_frame-frame_start)/5)
    
            tix_1=np.arange(frame_start,last_frame,frame_step)
            #scatter plot with a color vector
            p = ax.scatter(x, y, c=c, cmap = cmap,alpha=0.7)
            #add a vertical side bar that defines the color
            plt.colorbar(p, label=z_axis_label, shrink=.82, ticks = tix_1 )


            plt.axis('square')
            plt.xticks(rotation=45)

            # set graph limit conditions depending on unit specified
            if unit == "pixel":
                ax.set_xlim(pixel_min, pixel_max)
                ax.set_ylim(pixel_min, pixel_max)
                ax.yaxis.set_major_locator(
                    ticker.LinearLocator(axis_increment_pixel))
                ax.xaxis.set_major_locator(
                    ticker.LinearLocator(axis_increment_pixel))
            if unit == "nm":
                ax.set_xlim(nm_min, nm_max)
                ax.set_ylim(nm_min, nm_max)
                ax.yaxis.set_major_locator(
                    ticker.LinearLocator(axis_increment_nm))
                ax.xaxis.set_major_locator(
                    ticker.LinearLocator(axis_increment_nm))

            # display center
            if display_center == "yes":
                # in a centered graph, the center is actually(0,0)
                center1 = [0, 0]
                # plots center point as magenta X
                ax.scatter(0, 0, color='Magenta', marker="X", s=150)
                plt.text(x=center1[0] + 0.02,
                         y=center1[1] + 0.02, s='CENTER')

            # Jerry Adds a hover cursor
            mplcursors.cursor(hover=True)
            mplcursors.cursor(highlight=True)

            # axis labels and font configurations
            ax.set_xlabel(x_axis_label, fontweight='bold', fontsize=14)
            ax.set_ylabel(y_axis_label, fontweight='bold', fontsize=14)

            # plot title and font configurations
            # take the file name and separate from the extension
            # the first value in the tuple is the number
            # the second is .csv 
            # the number 00086.csv is the peak --> so this code takes the peak number
            pk = os.path.splitext(file_name)[0]

            graph_type = plot_type

            # change title order!!! 
            list_of_strings = [graph_type, exp_tag]

            #in quotes is the the delimiter between the items in the string
            # by default it is a _ 
            my_title = "_".join(list_of_strings)
            plt.title(my_title, fontweight='bold', fontsize=16)

            if save_plot == "yes":
                plt.savefig(my_title+"_2D.png")  # put title input and date time
        if plot_type == "find_excluded_angle":

            # data organization
            times = data["index"]
            conti_angle = data["Continuous Angle"]

            # setup fig and ax
            fig, ax = plt.subplots()

            # choose scatter plot or line plot
            if graph_style == 'scatter':
                ax.scatter(times, conti_angle, c='m', s=5)
            else:
                line, = plt.plot(times, conti_angle, 'm')

            # setting up for making vertical lines or indications of bad points (high and low points)
            ang_min = min(conti_angle)
            ang_max = max(conti_angle)

            # find and graph lower bad values
            bad_times_down = data_fil_down_bad["index"]
            ax.vlines([bad_times_down], ang_min, ang_max,
                      linestyles='dashed', colors='maroon')

            # find and graph upper bad values
            bad_times_up = data_fil_up_bad["index"]
            ax.vlines([bad_times_up], ang_min, ang_max,
                      linestyles='dashed', colors='tomato')

            # find and graph the invalid bad values
            invalid_times = data_back["index"]
            ax.vlines([invalid_times], ang_min, ang_max,
                      linestyles='dashed', colors='darkcyan')

            # Legend
            ax.legend(['Angle data', 'Below Lower Bound', 'Above Upper Bound',
                      'Invalid Readings'],bbox_to_anchor=(1.04,0), loc="lower left")

            # formatting
            plt.rcParams['figure.figsize'] = [13, 6]

            plt.xlabel('Frames')
            plt.ylabel('Angle Accumulation (degrees)')

            # Title configuration

            # take the file name and separate from the extension
            # the first value in the tuple is the number
            # the second is .csv 
            # the number 00086.csv is the peak --> so this code takes the peak number
            pk = os.path.splitext(file_name)[0]

            graph_type = 'Find_Excluded_Angle'

            # change title order!!! 
            list_of_strings = [graph_type, exp_tag]

            #in quotes is the the delimiter between the items in the string
            # by default it is a _ 
            my_title = "_".join(list_of_strings)
            my_title = graph_type + "\n" + my_title #add a line break between graph name and defining parameters
            plt.title(my_title)
            
            if frame_end == -1:   # negative 1
                last_frame = pre_data.iloc[:,0].iat[-1] #in the index column, give me the last valid value --> this is the max Frames
            else:
                last_frame = frame_end
                
            frame_step=int((last_frame-frame_start)/5)
    
            tix_1=np.arange(frame_start,last_frame,frame_step)
            plt.xlim(frame_start, last_frame+frame_step)
#             plt.ylim(-180,max(conti_angle)+180)
            ax.set_xticks(np.arange(frame_start,last_frame+frame_step, frame_step))
            ax.set_yticks(np.arange(-360, max(conti_angle)+180, 180))
            plt.xticks(rotation=-45)
            plt.grid()

            # hovering attempt 2
            # added by Jerry for Matplotlib compatible hovering
            mplcursors.cursor(hover=True)

            # Graph the newly calcuated Angular Continuous data, now filtered for good points only
        if plot_type == "angular_continuous_filtered":

            # data organization
            times = avt_good["Time (ms)"]
            conti_angle = avt_good["Continuous Angle"]

            # Graph a Scatter Plot otherwize the hover tool hovers to made up points
            fig, ax = plt.subplots()

            if graph_style == 'scatter':
                ax.scatter(times, conti_angle, c='m', s=5)
            else:
                line, = plt.plot(times, conti_angle, 'm')

            # setting up for making vertical lines or indications of bad points (high and low points)
            ang_min = min(conti_angle)
            ang_max = max(conti_angle)

            # find and graph lower bad values
            bad_times_down = data_fil_down_bad["Time (ms)"]
            ax.vlines([bad_times_down], ang_min, ang_max,
                      linestyles='dashed', colors='red')

            # find and graph upper bad values
            bad_times_up = data_fil_up_bad["Time (ms)"]
            ax.vlines([bad_times_up], ang_min, ang_max,
                      linestyles='dashed', colors='black')

            # find and graph the invalid bad values
            invalid_times = data_back["Time (ms)"]
            ax.vlines([invalid_times], ang_min, ang_max,
                      linestyles='dashed', colors='darkcyan')

            # Legend
            ax.legend(['Angle data', 'Below Lower Bound', 'Above Upper Bound',
                      'Invalid Readings'],bbox_to_anchor=(1.04,0), loc='lower right')

            # formatting
            plt.rcParams['figure.figsize'] = [20, 6]

            plt.xlabel('Time (ms)')
            plt.ylabel('Angle Accumulation (degrees)')

            # Title configuration

            graph_type = 'Accumulation of Angle (degrees) as a function of Time (ms)'

            #in quotes is the the delimiter between the items in the string
            # by default it is a 
            my_title = graph_type + "\n" + exp_tag #add a line break between graph name and defining parameters
            plt.title(my_title)

            plt.xlim(0, max(times)+1000)
            plt.ylim(ang_min-180, ang_max+180)
            ax.set_xticks(np.arange(0, max(times)+1000, 1000))
            ax.set_yticks(np.arange(ang_min-180, ang_max+180, 180))
            plt.xticks(rotation=-45)
            plt.grid()

            # hovering attempt 2
            # added by Jerry for Matplotlib compatible hovering
            mplcursors.cursor(hover=True)

        # Angular plot
        
        if plot_type == "angular":
            import plotly.express as px

            # Title configuration

            # take the file name and separate from the extension
            # the first value in the tuple is the number
            # the second is .csv 
            # the number 00086.csv is the peak --> so this code takes the peak number
            pk = os.path.splitext(file_name)[0]

            graph_type = "Angle vs Time"

            # change title order!!! 
            list_of_strings = [graph_type, exp_tag, pk]

            #in quotes is the the delimiter between the items in the string
            # by default it is a _ 
            my_title = " ".join(list_of_strings)

            
            fig = px.scatter(x = t, y = theta, title = my_title, width=1000, height=500, )
           
            fig.update_traces(hovertemplate='Time (ms): %{x} <br>Theta: %{y}')

            
        
            fig.update_layout(
            xaxis_title="Time (ms)",
            yaxis_title="Theta (degrees)"
            )
            

            fig.update_layout(
            yaxis = dict(
            tickmode = 'linear',
            tick0 = 360,
            dtick = 45))
            fig.show()
            
            if save_plot == "yes":
                fig.write_image(my_title +"_angular.png") 

           
            # Angular continuous plot
        if plot_type == "angular_continuous":
            
            import plotly.express as px

            # Title configuration

            # take the file name and separate from the extension
            # the first value in the tuple is the number
            # the second is .csv 
            # the number 00086.csv is the peak --> so this code takes the peak number
            pk = os.path.splitext(file_name)[0]

            graph_type = "Continuous Angular Rotation"

            # change title order!!! 
            list_of_strings = [graph_type, exp_tag, pk]

            #in quotes is the the delimiter between the items in the string
            # by default it is a _ 
            my_title = " ".join(list_of_strings)
            
            fig = px.line(x=t, y=thetac, title = my_title)
            fig.update_traces(hovertemplate='Time (ms): %{x} <br>Theta: %{y}')
            
            fig.update_layout(
            xaxis_title="Time (ms)",
            yaxis_title="Theta (degrees)")
            
            fig.update_layout(
            xaxis = dict(
            tickmode = 'linear',
            tick0 = 0,
            dtick = 1000))
            
            fig.update_layout(
            yaxis = dict(
            tickmode = 'linear',
            tick0 = -360,
            dtick = 180))
                        
            fig.show()
            if save_plot == "yes":
                fig.write_image(title+"_angular_continuous.png") 
    
    if plot_type in animations:
        
        # Accept my variables from graphing parameters
        (file_name, down_sampled_df, plot_type, display_center, ind_invalid_reading, rad_filter_type_upper,
                  rad_filter_type_lower, z_up, z_down, dist_high, dist_low, graph_style, bin_size, frame_start, frame_end,
                  display_center, exp_tag, x_axis_label, y_axis_label, z_axis_label, unit, pixel_min, pixel_max,
                  axis_increment_nm, axis_increment_pixel, nm_min, nm_max, save_plot, data_back, cmap, exp_tag, frame_speed, tail_length) = graph_parameters


        # Claire's code accepts down_sampled_df as df
        df = down_sampled_df

        #####################Graphing data assignment block##############
        # Here the code determines the units of the graph, only for cartesian graphs
        if unit == "pixel":
            x = df["X displacement (pixels)"]
            y = df["Y displacement (pixels)"]
        if unit == "nm":
            x = df["X displacement (nm)"]
            y = df["Y displacement (nm)"]
        z = df["Time (ms)"]

        if plot_type == "animated" or plot_type == 'HTML':
#             # allows for animation to animate
#             %matplotlib notebook
            #reassigning to a dataframe, setting up axis with empty list, aspect ratio
            coord = pd.DataFrame(x)
            coord.columns = ['x']
            coord['y'] = y
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sc = ax.scatter([], [])
            plt.axis('square')
            
            
            #Title configuration

            # take the file name and separate from the extension
            # the first value in the tuple is the number
            # the second is .csv 
            # the number 00086.csv is the peak --> so this code takes the peak number
            pk = os.path.splitext(file_name)[0]

            graph_type = 'Animated_2D_Plot'

            # change title order!!! 
            list_of_strings = [graph_type, exp_tag, pk]

            #in quotes is the the delimiter between the items in the string
            # by default it is a _ 
            my_title = "_".join(list_of_strings)
           

            # plot title and font configurations
            plt.title(my_title , fontweight = 'bold', fontsize = 16)
            
            # chosing units
            if unit == "pixel":
                    ax.set_xlim(pixel_min, pixel_max) 
                    ax.set_ylim(pixel_min, pixel_max)
                    ax.yaxis.set_major_locator(ticker.LinearLocator(axis_increment_pixel))# change to 5 for increments of .5
                    ax.xaxis.set_major_locator(ticker.LinearLocator(axis_increment_pixel))
                    ax.grid()
            if unit == "nm":
                    ax.set_xlim(nm_min, nm_max) 
                    ax.set_ylim(nm_min, nm_max)
                    ax.yaxis.set_major_locator(ticker.LinearLocator(axis_increment_nm))
                    ax.xaxis.set_major_locator(ticker.LinearLocator(axis_increment_nm))
                    ax.grid()
            ax.set_xlabel(x_axis_label, fontweight = 'bold', fontsize = 12)
            ax.set_ylabel(y_axis_label, fontweight = 'bold', fontsize = 12)

            # plot title and font configurations
            plt.title(my_title , fontweight = 'bold', fontsize = 16)





            # animation function feeds a window of dataframe values into the graphing function at a time,
            # iterates over user specified range in dataframe with user specified tail length
            # color of animation is also specified here
            def animate(count):
                sc.set_offsets(np.c_[coord.x.values[count-tail_length:count],coord.y.values[count-tail_length:count]])
                cmap = plt.cm.Greens
                norm = plt.Normalize(vmin=0, vmax=tail_length)
                z = np.array(range(tail_length))
                c = cmap(norm(z))
                sc.set_color(c)
                #button_ax = plt.axes([.78, .87, .1, .07]) # creates an outline for a potential button
               
                return sc

            ani = FuncAnimation(fig, animate, interval= frame_speed, frames = len(coord)) #frames=len(df)
            #ani.toggle(ax=button_ax)# potential button toggle for a potential button ;)
            if save_plot == 'yes':
                
                ani.save(my_title+'_animation_gif.gif', writer='pillow', fps=10, dpi=100)





            plt.tight_layout()
            plt.show()
            # With out the added if statement below, the animated plot will not animate 
            #(due to being a nested function)
        if plot_type == 'animated': 

                

            return ani
        # 
        if plot_type == 'HTML':
            plt.close('all')
            if save_plot == 'yes':
                with open(my_title+"_animation_html.html", "w") as f:
                    print(ani.to_jshtml(), file=f)
            return HTML(ani.to_jshtml())
            
#         if plot_type == "animated" or plot_type == 'HTML':
#             # allows for animation to animate
# #             %matplotlib notebook
#             #reassigning to a dataframe, setting up axis with empty list, aspect ratio
#             coord = pd.DataFrame(x)
#             coord.columns = ['x']
#             coord['y'] = y
#             fig = plt.figure()
#             ax = fig.add_subplot(111)
#             sc = ax.scatter([], [])
#             plt.axis('square')
            
#             # chosing units
#             if unit == "pixel":
#                     ax.set_xlim(pixel_min, pixel_max) 
#                     ax.set_ylim(pixel_min, pixel_max)
#                     ax.yaxis.set_major_locator(ticker.LinearLocator(axis_increment_pixel))# change to 5 for increments of .5
#                     ax.xaxis.set_major_locator(ticker.LinearLocator(axis_increment_pixel))
#                     ax.grid()
#             if unit == "nm":
#                     ax.set_xlim(nm_min, nm_max) 
#                     ax.set_ylim(nm_min, nm_max)
#                     ax.yaxis.set_major_locator(ticker.LinearLocator(axis_increment_nm))
#                     ax.xaxis.set_major_locator(ticker.LinearLocator(axis_increment_nm))
#                     ax.grid()
#             ax.set_xlabel(x_axis_label, fontweight = 'bold', fontsize = 12)
#             ax.set_ylabel(y_axis_label, fontweight = 'bold', fontsize = 12)
            
            
#             # Title configuration

#             # take the file name and separate from the extension
#             # the first value in the tuple is the number
#             # the second is .csv 
#             # the number 00086.csv is the peak --> so this code takes the peak number
#             pk = os.path.splitext(file_name)[0]

#             graph_type = 'Animated_2D_Plot'

#             # change title order!!! 
#             list_of_strings = [graph_type, exp_tag, pk]

#             #in quotes is the the delimiter between the items in the string
#             # by default it is a _ 
#             my_title = "_".join(list_of_strings)
           

#             # plot title and font configurations
#             plt.title(my_title , fontweight = 'bold', fontsize = 16)





#             # animation function feeds a window of dataframe values into the graphing function at a time,
#             # iterates over user specified range in dataframe with user specified tail length
#             # color of animation is also specified here
#             def animate(count):
#                 sc.set_offsets(np.c_[coord.x.values[count-tail_length:count],coord.y.values[count-tail_length:count]])
#                 cmap = plt.cm.Greens
#                 norm = plt.Normalize(vmin=0, vmax=tail_length)
#                 z = np.array(range(tail_length))
#                 c = cmap(norm(z))
#                 sc.set_color(c)
#                 #button_ax = plt.axes([.78, .87, .1, .07]) # creates an outline for a potential button
               
#                 return sc

#             ani = FuncAnimation(fig, animate, interval= frame_speed, frames = len(coord)) #frames=len(df)
#             #ani.toggle(ax=button_ax)# potential button toggle for a potential button ;)
#             if save_plot == 'yes':
                
#                 ani.save(my_title+'_animation_gif.gif', writer='pillow', fps=10, dpi=100)

#             plt.tight_layout()
#             plt.show()
#             # With out the added if statement below, the animated plot will not animate 
#             #(due to being a nested function)
#         if plot_type == 'animated': 

                

#             return ani
#         # 
#         if plot_type == 'HTML':
#             plt.close('all')
#             if save_plot == 'yes':
#                 with open(my_title+"_animation_html.html", "w") as f:
#                     print(ani.to_jshtml(), file=f)
#             return HTML(ani.to_jshtml())
        
        if plot_type == "interactive":

            #configure plot settings (currently called trace 1, may add more traces in the future)
            trace1=go.Scatter3d(x=x,
                                y=y,
                                z=z,
                                mode = "lines",
                                name = 'Original',
                                marker=dict(
                                    size=4,
                                    color='#e9ebf0',
                                    opacity=0.7,
                                    showscale=False,
                                    colorbar=dict(
                                        title='Time (ms)')),
                                line=dict(
                                    color='#e9ebf0',
                                    width=2))
            #assign traces
            fig = go.Figure(data=[trace1])
            
            # Title configuration

            # take the file name and separate from the extension
            # the first value in the tuple is the number
            # the second is .csv 
            # the number 00086.csv is the peak --> so this code takes the peak number
            pk = os.path.splitext(file_name)[0]

            graph_type = 'Interactive Plot'

            # change title order!!! 
            list_of_strings = [graph_type, exp_tag, pk]

            #in quotes is the the delimiter between the items in the string
            # by default it is a _ 
            my_title = "_".join(list_of_strings)
           

            # plot title and font configurations
            plt.title(my_title , fontweight = 'bold', fontsize = 16)

            
            
            #assign title
            
            fig.update_layout(title= my_title)
            #assign axis labels
            fig.update_layout(scene = dict(
                        xaxis_title= x_axis_label,
                        yaxis_title= y_axis_label,
                        zaxis_title= z_axis_label)) 

            #Here we can tweak the background color, grid color, and color of the origin for all axes/plane
            fig.update_layout(scene = dict(
                        xaxis = dict(
                             backgroundcolor="black",
                             gridcolor="gray",
                             showbackground=True,
                             zerolinecolor="white",),
                        yaxis = dict(
                            backgroundcolor="black",
                            gridcolor="gray",
                            showbackground=True,
                            zerolinecolor="white"),
                        zaxis = dict(
                            backgroundcolor="black",
                            gridcolor="gray",
                            showbackground=True,
                            zerolinecolor="white"),),
                      )

            #size and aspect ratio of the graph and the default camera zoom and angle 
            fig.update_layout(
            width=800,
            height=700,
            autosize=False,
            scene=dict(
            camera=dict(
                up=dict(
                    x=0,
                    y=0,
                    z=1
                ),
                eye=dict(
                    x=1,
                    y=2,
                    z=2,
                )
            ),
            aspectratio = dict( x=1, y=1, z=4 ),
            aspectmode = 'manual'
            ),
            )
            
            

            fig.show()

        
    if plot_type == "grid":

        #Receiving data Block
        [file_name, down_sampled_df, plot_type, display_center, exp_tag, x_axis_label, y_axis_label, z_axis_label, unit, 
        pixel_min, pixel_max, axis_increment_nm, axis_increment_pixel, nm_min, nm_max, save_plot, frame_start, frame_end, time_step,cmap,exp_tag, 
        frames_per_plot, columns, fig_size_x, fig_size_y] = graph_parameters

        #Changing names to accept DF
        df = down_sampled_df

        ######## SETUP FOR GRID PLOT

        #####################Graphing data assignment block##############
        # Here the code determines the units of the graph, only for cartesian graphs
        if unit == "pixel":
            x = df["X displacement (pixels)"]
            y = df["Y displacement (pixels)"]
        if unit == "nm":
            x = df["X displacement (nm)"]
            y = df["Y displacement (nm)"]
        z = df["Time (ms)"]
        
        #determine number of plots from amount of frames desired in each plot
        j = int(math.ceil(len(df)/frames_per_plot))
        
        if plot_type == 'grid':
            print(*['number of plots:',j])

        # Function for 2D plot parameters (called when user asks for multiple plots)
        # the grid_plot graph type runs best when this function is defined here and is called under plot_type == grid_plot if statement
        def do_plot(ax):
            #regular graphing parameters for 2D graph (color of scatter, size, shape, tick marks, etc.)
            colors = cm.spring(np.linspace(0, 1, len(z)))
            p=ax.scatter(x, y, c=colors)
            #fig = plt.figure(figsize=(6,6), dpi=100)
            tix = np.linspace(0,len(z),6)
            #tix_c = tix*time_step
            #cbar2.set_ticklabels(tix_c)
            plt.axis('square')
            plt.xticks(rotation=45)
            if unit == "pixel":
                ax.set_xlim(pixel_min, pixel_max) 
                ax.set_ylim(pixel_min, pixel_max)
                ax.yaxis.set_major_locator(ticker.LinearLocator(axis_increment_pixel))# change to 5 for increments of .5
                ax.xaxis.set_major_locator(ticker.LinearLocator(axis_increment_pixel))
                ax.grid()
            if unit == "nm":
                ax.set_xlim(nm_min, nm_max) 
                ax.set_ylim(nm_min, nm_max)
                ax.yaxis.set_major_locator(ticker.LinearLocator(axis_increment_nm))
                ax.xaxis.set_major_locator(ticker.LinearLocator(axis_increment_nm))
                ax.grid()
        
        i = 0
        dfs = np.array_split(df,j) # splits large dataframe into "j" equal dataframes
        #print (dfs[i]) #<--- command to print each dataframe # dataframe 0 is the first dataframe
        
        #this portion specifies subplot dimentions (N plots in 3 columns and amount of appropriate rows)
        cols = columns
        rows = int(math.ceil(j / cols)) #determining rows based on the number of graphs and columns

        gs = gridspec.GridSpec(rows, cols, wspace = .25, hspace = .25)# disallows overlap
        fig = plt.figure(figsize = (fig_size_x,fig_size_y))
        
        while i < j:
            
            if unit == "pixel":
                x = dfs[i]["X displacement (pixels)"]
                y = dfs[i]["Y displacement (pixels)"]
            if unit == "nm":
                x = dfs[i]["X displacement (nm)"]
                y = dfs[i]["Y displacement (nm)"]
            # x = dfs[i].iloc[:,x_unit] 
            # y = dfs[i].iloc[:,y_unit]
            z = dfs[i].iloc[:,7]
            ax = fig.add_subplot(gs[i])
            ax.set_title(dfs[i].index[0])
            do_plot(ax)
            i+= 1

        # plot title and font configurations

        # take the file name and separate from the extension
        # the first value in the tuple is the number
        # the second is .csv 
        # the number 00086.csv is the peak --> so this code takes the peak number
        pk = os.path.splitext(file_name)[0]

        graph_type = '2D_Map'

        # change title order!!! 
        list_of_strings = [graph_type, exp_tag]

        #in quotes is the the delimiter between the items in the string
        # by default it is a _ 
        my_title = "_".join(list_of_strings)

        if save_plot == "yes":
            
            framestr = '{}'.format(frames_per_plot)
            plt.savefig(my_title+'_'+framestr+'_frames_per_plot'+'_gridplot.png', dpi=300)

        plt.show()

# Note 4: DORA.table
# the following code creates the relevant data tables from the inputs detailed below in the template section
# the outputs are 4 data tables and, if you so choose, a saved .csv of a data table

'''
Template for DataTable inputs

table_parameters = [down_sampled_df, ind_invalid_reading, rad_filter_type_upper,
                    rad_filter_type_lower, z_up, z_down, dist_high, dist_low, bin_size, data_back,save_table,file_name]

data, avt_good, avt_bad, data_final_final = DORA.table(*table_parameters)

                    
'''


def table(*table_parameters):
    # Set up Block

    # Accept my variables from graphing parameters
    (down_sampled_df, ind_invalid_reading, rad_filter_type_upper,
     rad_filter_type_lower, z_up, z_down, dist_high, dist_low, bin_size, data_back, save_table, file_name) = table_parameters

    # import the 2 ways to analyze angular data: 1) Claire's way, 2) Jerry's Way [More Current]
    import AngleCalc

    # Gather inputs to Cacluate Angle under Jerry's Angle calculation paradigm
    inputs_avt_filter = [down_sampled_df, ind_invalid_reading, rad_filter_type_upper,
                         rad_filter_type_lower, z_up, z_down, dist_high, dist_low, bin_size]
    # Call and Run Jerry's Angle Calculation
    data, xy_goodbad, avt_good, avt_bad, data_fil_dsa, data_fil_down_bad, data_fil_up_bad = AngleCalc.avt_filter(
        *inputs_avt_filter)

    print("Checkpoint 1 PASSED: imported function parameters")

    ################################### [Final Data Table Assembly ] ######################################

    # Organzize Data Table with Final Filtered Data [Re insert sus points from lower and upper bound filtering]
    # slap all the bad data on the end of the good data
    data_final = pd.concat([avt_good, avt_bad])
    # sort by index so that values go back to where they are supposed to be:
    data_final = data_final.sort_values(by=["index"])

    # re insert sus points [re insert sus points from invalid]
    # Organzize Data Table with Front End data (data) and back end data (data_back)
    # slap all the bad data on the end of the good data
    data_final_final = pd.concat([data_back, data_final])
    # sort by index so that values go back to where they are supposed to be:
    data_final_final = data_final_final.sort_values(by=["index"])

    del data_final_final['X position']
    del data_final_final['Y position']

    # Label of the Data with either Normal, Upper bound, Lower Bound, Invalid Reading

    # Initialize the data table to be populated
    data_final_final["Excluded Type"] = 'None'
    # store this into a dummy vector
    dummy_vec = data_final_final["Excluded Type"].copy()

    # Set all indices of AVT_good to normal
    # Select all indicies that are NORMAL --> avt_good indicies
    ind_Normal = avt_good["index"].copy()
    dummy_vec[ind_Normal] = 'Normal'

    # Set all indicies of data_back to 'Invalid Reading' and put them in the dummy variable
    ind_IR = data_back["index"].copy()
    dummy_vec[ind_IR] = 'Invalid Reading'

    # Find bad lower bounds and index for them and set value to "below bound"
    ind_bad_down = data_fil_down_bad["index"].copy()
    dummy_vec[ind_bad_down] = 'Below Bound'

    # fFind bad upper bounds and index for them and set value to "upper bound"
    ind_bad_up = data_fil_up_bad["index"].copy()
    dummy_vec[ind_bad_up] = 'Above Bound'

    data_final_final["Excluded Type"] = dummy_vec
    data_final_final = data_final_final[[
        "index", "Time (ms)", "Angle", "Delta Angle", "Continuous Angle", "Excluded Type"]]


    print("Checkpoint 2 passed: all data is processed")

    if save_table == "yes":
        my_title = file_name + "_Final_DataTable.csv"
        data_final_final.to_csv(
            my_title, index="false")
        print("I have saved the table for " + file_name + " as " + my_title)

    return data, avt_good, avt_bad, data_final_final

# extract the relevant data from dataframe and convert into lists



# Note #5: DORA.collect_variable
# this function takes the input of the final data table that you want to analyze and takes that 
# data and stores it into its own separate csv 


def collect_variable(DataTable, col, file_name, sample_conditions, name_saving_folder):
    # this function operates on the final data exported by DORA.table
    # it takes 1) the data table 2) the column you want to operate on
    # and it gives you that column as a csv
    # 3) is the filename that this data came from
    # 4) are the name of the experimental set up so that you know which data this csv comes from
    # 5) name of the folder you want to save your outputted column .csv's into

    # Intialize my lists to store values
    collection_pot = []

    # Extract the Target Variable or col varible from both
    var_col = list(DataTable[col].copy())

    # Take extracted data and save into a my collection pot list
    collection_pot.append(var_col)

    # convert the extracted data into a dataframe
    df_collection_pot = pd.DataFrame(collection_pot)

    # set the name of the file to be saved
    my_file_name = sample_conditions + col + "_from_" + file_name + ".csv"

    # create a folder to save CSV into called "name_saving_folder"

    # if the folder does not exist already
    if not os.path.exists(name_saving_folder):
        # make a folder in the current path
        os.mkdir(name_saving_folder)
        # and report this to the user
        print("Directory ", name_saving_folder,  " Created ")
    else:
        # else tell user it exists already
        print("Directory ", name_saving_folder,  " already exists")

    # Where are we trying to save the file?
    # ANSWER: Current directory/name_saving_folder/my_file_name.csv

    # get current path
    path_OG = os.getcwd()

    # Create path to the saving folder's path
    my_saving_path = os.path.join(path_OG, name_saving_folder, my_file_name)

    # save the dataframe
    df_collection_pot.to_csv(my_saving_path)

    print('I have SAVED the file, {FileName}, in the directory {FolderName}'.format(
        FileName="my_file_name", FolderName="name_saving_folder"))

# Note 6: DORA.filter_centers
# this function takes a FOLDER of .csv files that are ALREADY SECTIONED trajectories 
# and outputs a inner folder and .csv of all the .csv's that pass the filter

# '''

# Parameters section:
# #universal parameters
# #get the name of your folder as a string and put r in front 
# folder_name = r"C:\Users\jerry\Desktop\Research\Kosuri\Rotor_Data_Interpretation\Jerry_Time_to_shine\DORA_Visualization-main\DORA_Visualization_Testing_3.31.2022\trajectory_bucket\26_laser_bucket"
# pixel_size = 154  # in nanometers (or nm per pixel or nm/px )
# time_step = 20  # miliseconds per frame in trajectory movie
# frame_start = 0  # enter 0 to start from beginning of dataset
# frame_end = -1  # enter -1 to end at the last value of the data set
# cmap = "spring" # enter a color map string from this https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
# exp_tag = "Laser_26" # a tag that caries the name of the experiment
# first_zero_end = 'no'  # yes to cut off all values after first 0,0 = x,y
# graph_centers = "no" #'yes' or 'no' to graphing the centers of the data. 
# save_plot = 'yes' 


# outter_parameters = [folder_name, pixel_size, time_step, frame_start, frame_end, cmap, exp_tag, first_zero_end, my_rad, graph_centers, save_plot]
# passing_files, Center_Checkpoint = DORA.filter_centers(*outter_parameters)

# '''

def filter_centers(*relevant_parameters):

    # unpackage our list of variables with their associated variables with variable names
    # unpack list into variables
    folder_name, pixel_size, time_step, frame_start, frame_end, cmap, exp_tag, first_zero_end, my_rad, graph_centers, save_plot = relevant_parameters

    #Create a for loop to run over a folder of csv's of .csv trajectories that meet the following criterion:
    # 1) filtered so that the relevant angular range is defined 

    # to iterate over a folder you need all the files of the folder first

    #Change the folder directory to be the current folder's 
    os.chdir(folder_name)

    #Take all files in the current folder(the one we just switched to) and store it in a list through which we will iterate
    my_files = os.listdir(os.getcwd())
    
    #Initialize Buckets to collect radii and csv_filenames
    rad_bucket = []
    csv_bucket = []

    for file_name in my_files:
        if '.csv' in file_name:
            
            #run DORA.find_center
            initial_parameters = [file_name, time_step, frame_start, frame_end, cmap, exp_tag, first_zero_end, graph_centers, save_plot]
            my_rad_estimate = find_center(*initial_parameters)[4]
            #collect radii
            rad_bucket.append(my_rad_estimate)
            #collect file_names
            csv_bucket.append(file_name)
            print(f"I have stored the radius of {file_name}")
    #make a dataframeout from csv_bucket
    Center_Checkpoint = pd.DataFrame(csv_bucket, columns=["File Name"]) 
    # add a column of data frame as rad_bucket
    Center_Checkpoint["Radius (pixels)"] = rad_bucket
    # add a column where radius is converted from pixels inot nm 
    Center_Checkpoint["Radius (nm)"] = Center_Checkpoint["Radius (pixels)"] * pixel_size #pixel size is nm/px so [nm = px * nm/px ]

    # unpack radius min and max into their own variables
    
    min_rad = my_rad[0]
    max_rad = my_rad[1]

    #store this in a separate checkpoint holder
    rad_nm = Center_Checkpoint["Radius (nm)"].copy()
    pnp_max_min = (min_rad < rad_nm) & ( rad_nm < max_rad)
    Center_Checkpoint["Passes Claire's Centering Algorithm"] = pnp_max_min

    #create a new folder
    name_saving_folder = "CSVs_that_PASS"
    # if the folder does not exist already
    if not os.path.exists(name_saving_folder):
        # make a folder in the current path
        os.mkdir(name_saving_folder)
        # and report this to the user
        print("Directory ", name_saving_folder,  " Created ")
    else:
        # else tell user it exists already
        print("Directory ", name_saving_folder,  " already exists")

    # Where are we trying to save the file?
    # ANSWER: Current directory/name_saving_folder/my_file_name.csv

    # get current path
    path_OG = os.getcwd()

    # create saving path 
    my_saving_path = os.path.join(path_OG, name_saving_folder)

    #make a dataframeout from csv_bucket
    Center_Checkpoint = pd.DataFrame(csv_bucket, columns=["File Name"]) 
    # add a column of data frame as rad_bucket
    Center_Checkpoint["Radius (pixels)"] = rad_bucket
    # add a column where radius is converted from pixels inot nm 
    Center_Checkpoint["Radius (nm)"] = Center_Checkpoint["Radius (pixels)"] * pixel_size #pixel size is nm/px so [nm = px * nm/px ]

    # unpack radius min and max into their own variables
    min_rad = my_rad[0]
    max_rad = my_rad[1]

    #store this in a separate checkpoint holder
    rad_nm = Center_Checkpoint["Radius (nm)"].copy()
    pnp_max_min = (min_rad < rad_nm) & ( rad_nm < max_rad)

    Center_Checkpoint["Passes Claire's Centering Algorithm"] = pnp_max_min
    # my_diameters = rad_bucket * pixel_size # radius (px) * pixel sie (nm/px) = nm


    #Give me only the rows in the df that pass, and then only the file names
    passing_files = Center_Checkpoint[pnp_max_min]["File Name"]

    #Create a CSV of the Dataframe that shows which pk's passed the test
    my_title = "Which_Peaks_Passed.csv"
    my_saving_path_title = os.path.join(path_OG, name_saving_folder, my_title)
    Center_Checkpoint.to_csv(my_saving_path_title, index="false")
    print(f"I have saved which folders have passed as datatable {my_title} ")

    #Creat a saving folder for the paths that pass
    name_saving_folder = "CSVs_that_PASS"
    # get current path
    path_OG = os.getcwd()
    #Add the name of the naming folder to create full path
    my_saving_path = os.path.join(path_OG, name_saving_folder)
    #for the current file,  f, in the files that passed
    for f in passing_files:
        shutil.copy(f, my_saving_path) #shuttle it to the name saving 
        print(f"I have moved {f} to {name_saving_folder}")

    return  passing_files, Center_Checkpoint


# purpose: When analyzing bring spots (non circular, just the grounded oligo), find center by binning the data and finding the brightest spot. This new spot is the center.
def find_center_hist_max(x, y,bin_num):
            # find the center of a plot using a low resolution histogram. Find the max value of that histogram. The corresponding x and y indicies become those of the center make the center
            H, xedges , yedges = np.histogram2d(x,y, bins = bin_num)

            #find the x and y index of the maximum value histogram 
            indMax = np.unravel_index(H.argmax(),H.shape)

            #Set the value of the max x and y indexes to be the new OR(OverRidden) center
            center_OR = [xedges[indMax[0]],yedges[indMax[1]]]

            plt.hist2d(x,y,bins = bin_num, cmap = "viridis")
            plt.scatter(x = center_OR[0], y= center_OR[1], color = "magenta")
            plt.text(x= center_OR [0] + 0.02, y= center_OR[1] + 0.02, s='CENTER')

            # Add title and axis names
            plt.title('2D Histogram of Centered Data')
            plt.xlabel('X pixels')
            plt.ylabel('Y pixels')
            plt.show()
            
            #return the center found
            return center_OR

    


            