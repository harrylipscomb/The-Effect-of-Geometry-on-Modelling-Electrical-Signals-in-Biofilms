"""
Created on Sat Oct 5 12:11:44 2022

Code used to generate undulating surfaces using the diamond-square algorithm (DSA). A box-counting method is then used to determine the fractal dimension of the surfaces generated.
Fractal dimension, surface roughness (a parameter used within the DSA), root-mean square of the bacteria heights and kurtosis are all plotted together.

@author: Harry Lipscomb 
"""

import random
import numpy as np
import matplotlib.pyplot as plt

def main():
    k = 5 # size of box parameter
    n = 2**k +1
    rs = np.arange(0, 9, 0.5)
    repeats = 1
    show=False
    
    def fixed( d, i, j, v, offsets ):
        # For fixed bdries, all cells are valid. Define n so as to allow the
        # usual lower bound inclusive, upper bound exclusive indexing.
        n = d.shape[0]
        res, k = 0, 0
        for p, q in offsets:
            pp, qq = i + p*v, j + q*v
            if 0 <= pp < n and 0 <= qq < n:
                res += d[pp, qq]
                k += 1.0
        return res/k
    
    def periodic( d, i, j, v, offsets ):
        # For periodic bdries, the last row/col mirrors the first row/col.
        # Hence the effective square size is (n-1)x(n-1). Redefine n accordingly!
        n = d.shape[0] - 1
    
        res = 0
        for p, q in offsets:
            res += d[(i + p*v)%n, (j + q*v)%n]
        return res/4.0
    
    def single_diamond_square_step( d, w, s, avg ):
        # w is the dist from one "new" cell to the next
        # v is the dist from a "new" cell to the nbs to average over
        
        n = d.shape[0]
        v = w//2
        
        # offsets:
        diamond = [ (-1,-1), (-1,1), (1,1), (1,-1) ]
        square = [ (-1,0), (0,-1), (1,0), (0,1) ]
    
        # (i,j) are always the coords of the "new" cell
    
        # Diamond Step
        for i in range( v, n, w ):
            for j in range( v, n, w ):
                d[i, j] = avg( d, i, j, v, diamond ) + random.uniform(-s,s)
    
        # Square Step, rows
        for i in range( v, n, w ):
            for j in range( 0, n, w ):
                d[i, j] = avg( d, i, j, v, square ) + random.uniform(-s,s)
    
        # Square Step, cols
        for i in range( 0, n, w ):
            for j in range( v, n, w ):
                d[i, j] = avg( d, i, j, v, square ) + random.uniform(-s,s)
    
    def make_terrain( n, r, bdry ):
        # Returns an n-by-n landscape using the Diamond-Square algorithm, using 
        # roughness delta ds (0..1). bdry is an averaging fct, including the
        # bdry conditions: fixed() or periodic(). n must be 1+2**k, k integer.
        d = np.zeros( n*n ).reshape( n, n )
    
        w, s = n-1, r
        while w > 1:
            single_diamond_square_step( d, w, s, bdry )
            w //= 2
    
        return d
    
    def delete_row_and_column(terrain):
    
        return terrain[0:-1, 0:-1]
    
    def box_lengths_of_interest(n):
    
        boxes = np.array([])
    
        for value in np.arange(2, n):
            box = 2**value
            boxes = np.append(boxes, box)
    
        return boxes
    
    def check_box(box_size, z, terrain_of_interest):
    
        if box_size != 1:    
            for value in terrain_of_interest:
                for zs in value:
                    if z>=zs>(z-box_size):
                        return 1
            return 0
        else:
            if z>=terrain_of_interest[0]>(z-box_size):
                return 1
            return 0
    
    def boxes_overlapping_surface(box_sizes, terrain):
    
        boxes_overlapping = np.array([])
        length_terrain = len(terrain[:,0])
        for box_size in box_sizes:
            boxes_contained_for_box_size = 0
            x = 0
            y = 0
            z = length_terrain
            while x < length_terrain:
                while y < length_terrain:
                    while z>=(-3*length_terrain):
                        terrain_of_interest = terrain[x:int(x+box_size), y:int(y+box_size)]
                        is_contained = check_box(box_size, z, terrain_of_interest)
                        boxes_contained_for_box_size += is_contained
                        z -= box_size/10
                    y += box_size
                    y = int(y)
                    z = 3*length_terrain
                x += box_size
                x = int(x)
                y = int(0)
            boxes_overlapping = np.append(boxes_overlapping, boxes_contained_for_box_size)
    
        return boxes_overlapping
    
    def plot_results(x, y, r, show):
    
        parameters = np.polyfit(x,y,1)
    
        if show:
            fit_line = parameters[0]*x + parameters[1]
            plt.figure()
            plt.title(r"$log(N)$ vs $log(r)$ for $r={}$".format(r))
            plt.ylabel('$log(N)$')
            plt.xlabel('$log(r)$')
            plt.scatter(x, y)
            plt.plot(x, fit_line, '--')
            plt.show()
    
        return parameters
    
    def rms(terrain):
        average_line = np.average(terrain)
        shapes = np.shape(terrain)
        rows = shapes[0]
        columns = shapes[1]
        rms = 0
        for row in np.arange(0, rows, 1):
            for column in np.arange(0, columns, 1):
                rms += np.sqrt((terrain[row, column] - average_line)**2)
        rms = rms/(rows*columns)
        print(rms)
    
        return rms
    
    def kurtosis(terrain, rms):
        average_line = np.average(terrain)
        shapes = np.shape(terrain)
        rows = shapes[0]
        columns = shapes[1]
        kurtosis_number = 0
        for row in np.arange(0, rows, 1):
            for column in np.arange(0, columns, 1):
                kurtosis_number += np.sqrt((terrain[row, column] - average_line)**4)
        kurtosis_number = kurtosis_number/(rows*columns)
        kurtosis_number = kurtosis_number/(rms**4)
    
        return kurtosis_number
    
    def sub_main(n, r, fixed, show=False):
    
        terrain = make_terrain(n, r, fixed)
        terrain = delete_row_and_column(terrain)
        box_sizes = box_lengths_of_interest(k)
        rms_value = rms(terrain)
        kurtosis_value = kurtosis(terrain, rms_value)
        boxes_overlapping = boxes_overlapping_surface(box_sizes, terrain)
        x = np.log(1/box_sizes)
        y = np.log(boxes_overlapping)
        parameters = plot_results(x, y, r, show)
    
        return parameters[0], rms_value, kurtosis_value
    
    r_values = np.array([])
    
    D_means = np.array([])
    D_Stds = np.array([])
    
    rms_means = np.array([])
    rms_Stds = np.array([])
    
    kurtosis_means = np.array([])
    kurtosis_Stds = np.array([])
    
    for r in rs:
        r_values = np.append(r_values, r)
        i = 0
        D_values = np.array([])
        rms_array = np.array([])
        kurtosis_array = np.array([])
        while i<=(repeats-1):
            D, rms_value, kurtosis_value = sub_main(n, r, fixed, show)
            D_values = np.append(D_values, D)
            i += 1
            rms_array = np.append(rms_array, rms_value)
            kurtosis_array = np.append(kurtosis_array, kurtosis_value)
    
        D_value = np.average(D_values)
        D_Std = np.std(D_values)
        D_means = np.append(D_means, D_value)
        D_Stds = np.append(D_Stds, D_Std)
        
        rms_value = np.average(rms_array)
        rms_Std = np.std(rms_array)
        rms_means = np.append(rms_means, rms_value)
        rms_Stds = np.append(rms_Stds, rms_Std)
        
        kurtosis_value = np.average(kurtosis_array)
        kurtosis_Std = np.std(kurtosis_array)
        kurtosis_means = np.append(kurtosis_means, kurtosis_value)
        kurtosis_Stds = np.append(kurtosis_Stds, kurtosis_Std)
    
        print('Completed r = {}'.format(r))
    
    parameters = np.polyfit(r_values, D_means, 1)
    #best_fit_line = parameters[0]*r_values + parameters[1]
    
    plt.figure(figsize=(8,6))
    plt.ylabel('$D$', size=18)
    plt.xlabel('$r$', size=18)
    plt.title("Fractal Dimension against roughness parameter")
    plt.errorbar(r_values, D_means, yerr=D_Stds, color='k', marker='o', markersize=4)
    #plt.plot(r_values, best_fit_line, '--')
    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.savefig('r vs D for DSA', dpi=900)
    plt.legend()
    plt.show()
    
    
    plt.figure(figsize=(8,6))
    plt.ylabel('$R_{rms}$', size=18)
    plt.xlabel('$r$', size=18)
    plt.title('Root Mean Square vs Roughness parameter $r$', size=14)
    plt.errorbar(r_values, rms_means, yerr=rms_Stds, color='k', marker='o', markersize=4)
    #plt.plot(r_values, best_fit_line, '--')
    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.savefig('RMS vs r for DSA', dpi=900)
    plt.show()
    
    
    plt.figure(figsize=(8,6))
    plt.ylabel('$R_{ku}$', size=18)
    plt.xlabel('$r$', size=18)
    plt.title('Kurtosis vs Roughness parameter $r$', size=14)
    plt.errorbar(r_values, kurtosis_means, yerr=kurtosis_Stds, color='k', marker='o', markersize=4)
    #plt.plot(r_values, best_fit_line, '--')
    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.savefig('kurtosis vs r for DSA', dpi=900)
    plt.legend()
    plt.show()
    
    for index, r in enumerate(r_values):
        print("r: {}, D: {}, Unc_D: {}, rms: {}, Unc_rms {}, kurtosis: {}, Unc_kurtosis {} ".format(r_values[index], D_means[index], D_Stds[index],
                                                                                                    rms_means[index], rms_Stds[index],
                                                                              kurtosis_means[index], kurtosis_Stds[index]))
