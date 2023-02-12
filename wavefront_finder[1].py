# -*- coding: utf-8 -*-
"""
Created on Sat Oct 1 11:17:55 2022

Basic code formed with the function of plotting wave position against time for a 
spherical 3D wave. Here, the instantaneous position of the wavefront is found using 
a defined threshold to convert an array of K+ concentrations (within a radius) to 
a binary '1' (standing for the wavefront has passed this point) or a 0 if the wavefront
is yet to reach this radius. Once a 1 has been achieved, only the next radius is 
considered etc. 

@author: Harry Lipscomb 
"""

import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pds
from scipy.optimize import curve_fit

threshold = 3*(10**8)
number_headings = 1 # This is equal to number of rows prior to t = 0 K concentrations
number_initial_columns = 1 # This is equal to number of columns prior to K concentrations
radius_spacing_per_column = 1 # The radius change from one column to the next in micrometers
diffusion_coefficient_used = 0.1
number_bacteria_used = 3000
shape_used = 'Hemisphere'
duplicate_same_radii = True # Duplicate radii True if you want the same radii to be plot per time step
time_step = 1
intrinsic_uncertainty = 0.5


plt.rcParams.update({'font.size': 14})

path = 'C:/Fourth_Year/Masters_Project/Final Data/Menger Sponge/'

file_name = 'flatMenger.csv' # Use forward slashes here.
file = pds.read_csv(path + file_name, skiprows=number_headings) # Converts csv to a Pandas data frame
data_1 = file.to_numpy() # Converts a pandas data from to a numpy array

# read in more data files by repeating above 3 lines if necessary

def binary_wave_function(array, threshold, initial_columns):
    """
    Takes an array of concentrations, a threshold concentration and the number 
    of columns prior to this concentration data, and produces an array of 1s and 
    0s. A 1 represents an above threshold concentration within the respective 
    radii and a 0 represents a below threshold concentration. The array is sliced
    to only return those times where the initial wave has not yet reached the 
    edgeo of the area of interest. The function could be  developed to consider 
    a return or secondary wave.

    """
    [total_rows, total_columns] = np.shape(array)
    row = 0
    column = initial_columns
    for row in np.arange(0, total_rows, 1):
        for column in np.arange(0, total_columns, 1):
            if array[row, column]>threshold:
                array[row, column] = 1
            else:
                array[row, column] = 0
            
    # while (column<total_columns and row<total_rows-1):
    #     if array[row, column]>threshold:
    #         array[row, column] = 1
    #         column += 1
    #     else:
    #         array[row, column:] = 0
    #         column = initial_columns
    # row += 1
    time_reach_edge = array[row, 0]-1
    data_initial_wave = array[0:row+1, :]

    return data_initial_wave, time_reach_edge

def wavefront_finder(array, radius_spacing):
    """
    This function takes in an array made up of only binary values. It finds the 
    outermost '1' (where there is only 1s inside of this radii) and considers this
    the front of the wave. The function returns an array of times elapsed and a 
    corresponding array of values with the wavefront radii. The function could be 
    developed to consider a return or secondary wave.

    """
    times = array[:,0]
    wavefront_positions = np.array([0])
    row = 0
    column = 1
    radius = radius_spacing
    while row<np.shape(array)[0]-1: # While row being considered is less than the total number of rows
        while column<np.shape(array)[1]-1:    #
            if array[row, column] == 1:
                column += 1 
                if column > radius:
                    radius = column
            elif array[row, column] == 0:
                column += 1
            else:
                print("There has been an error in determining wavefront location.")
        wavefront_positions = np.append(wavefront_positions, radius-1)
        row += 1
        column = 1

    return times, wavefront_positions

def Data_Smoother(some_Array):
    
    some_Data = np.unique(some_Array, return_counts = True)
    boss_Array = np.array([])
    
    for index, count in enumerate(some_Data[1][:-1]):
        counter = 0
        if count == 1:
            boss_Array = np.append(boss_Array, some_Data[0][index])
        else:
            while counter < count:
                value_to_add = some_Data[0][index] + (some_Data[0][index+1] - some_Data[0][index])*counter/count
                boss_Array = np.append(boss_Array, value_to_add)
                counter += 1 
    
    no_times_final_value_added = some_Data[1][-1]
    final_count_to_add = np.ones(no_times_final_value_added)*some_Array[-1]
        
    boss_Array = np.append(boss_Array, final_count_to_add)
    print(boss_Array)
    
    return boss_Array

def standard_deviation_finder(arrays, time_step):
    lengths = np.array([])
    for array in arrays:
        lengths = np.append(lengths, len(array))
        
    longest_array_length = int(np.max(lengths))
    data = np.zeros((len(arrays), longest_array_length))
    
    for index, array in enumerate(arrays):
        for time, value in enumerate(array):
            try:
                data[index, time] = value
            except:
                data[index, time] = 0
    means = np.array([])
    std_results = np.array([])
    for i, time in enumerate(data[0, :]):
        
        results = data[:, i]
        mean = np.array([])
        results_cleaned = np.array([])
        for value in results:
            if value > 0:
                results_cleaned = np.append(results_cleaned, int(value))
        mean = np.mean(results_cleaned)
        std_result = np.std(results_cleaned)

        means = np.append(means, mean)
        std_results = np.append(std_results, std_result)
        
        means[0] = 0
        std_results[0] = 0

    longest_array_times = np.arange(0, longest_array_length, time_step)

    return means, std_results, longest_array_times

def func_to_fit(x, alpha, d, m):

    return ( m * x**alpha + d)

def red_chi_square(residuals, errors):
    
    chi_square = np.sum((residuals)**2 / errors**2)

    print(chi_square)
    reduced = chi_square / (len(residuals) - 3)
    
    return reduced

data_1, time_until_edge_1 = binary_wave_function(data_1, threshold, number_initial_columns)

times_1, wavefront_radii_1 = wavefront_finder(data_1, radius_spacing_per_column)

# more repeats for data taken can be included by repeated above two lines for different data (replace data_1)

if (duplicate_same_radii == False):
    wavefront_radii_1 = Data_Smoother(wavefront_radii_1)

all_arrays = [wavefront_radii_1**2] # add more repeats here if required

mean, std, long_t_array = standard_deviation_finder(all_arrays, time_step)

plt.figure()
plt.title('Sphere')
plt.ylabel(r'$<Y^{2}> [\mu m^2]$')
plt.xlabel(r"$\tau[s]$")
plt.plot(times_1, wavefront_radii_1**2, '.', label='Repeat 1')

plt.legend(fontsize=10)
plt.tight_layout()
plt.show()

left_area_interest = 4
right_area_interest = 84

# COMBINE UNCERTAINTIES IN QUADRATURE:

# TO FIND INTRINSIC UNCERTAINTIES: 

intrinsic_uncertainties = np.array([])
for index, value in enumerate(long_t_array):
    unc = 2*mean[index]*intrinsic_uncertainty / np.sqrt(5)
    intrinsic_uncertainties = np.append(intrinsic_uncertainties, unc)

# std here is the variance
for index, value in enumerate(std):
    std[index] = np.sqrt(value**2 + intrinsic_uncertainties[index]**2)
  
    
#std[0] = 0
uncertainties = np.sqrt(std)

times_fitted = long_t_array[left_area_interest: right_area_interest]
means_fitted = mean[left_area_interest: right_area_interest]
errors_fitted = uncertainties[left_area_interest: right_area_interest]

popt, pcov = curve_fit(func_to_fit, times_fitted, means_fitted, sigma=errors_fitted, maxfev=20000)
fit_curve = func_to_fit(times_fitted, popt[0], popt[1], popt[2])

residuals = means_fitted - fit_curve
chi_square_red = red_chi_square(residuals, errors_fitted) 

plt.figure()
plt.errorbar(long_t_array, mean, yerr=uncertainties, fmt='.', elinewidth=1, markersize=4, label='Data', color='k')
plt.xlabel(r'$\tau ~ [s]$')
plt.ylabel(r'$ \langle Y^2 \rangle ~ [\mu m^2]$')
plt.title('2D Menger Sponge')
plt.tick_params('x', labelbottom=True)

plt.legend(loc='upper left', fontsize=10)
plt.tight_layout()
plt.show()


