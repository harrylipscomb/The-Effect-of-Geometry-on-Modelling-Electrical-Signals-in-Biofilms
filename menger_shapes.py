# -*- coding: utf-8 -*-
"""
Generates points within a menger sponge

Created on Tue Nov 29 11:40:27 2022

@author: Emily
"""

import numpy as np
import matplotlib.pyplot as plt

def menger(matrix, size, dimension):
    quotient, remainder = divmod(size, 3)
    
    if dimension == 2:

        if remainder == 0:
            for x in np.arange(0, size, quotient):
                for y in np.arange(0, size, quotient):
                    view = matrix[x:x + quotient, y:y + quotient]

                    if (x // quotient) % 3 == 1 and (y // quotient) % 3 == 1:
                        view *= 0

                    menger(view, quotient, 2)
                    
    if dimension == 3:
        if remainder == 0:
            for x in np.arange(0, size, quotient):
                for y in np.arange(0, size, quotient):
                    for z in np.arange(0, size, quotient):
                        view = matrix[x:x + quotient, y:y + quotient, z:z + quotient]

                        if ((x // quotient) % 3 == 1 and (y // quotient) % 3 == 1):
                            view *= 0
                        # Commment out lower if statements to have holes in one
                        # direction through the cube otherwise a fill menger 
                        # cube will be made
                        if ((x // quotient) % 3 == 1 and (z // quotient) % 3 == 1):
                            view *= 0
                        if ((y // quotient) % 3 == 1 and (z // quotient) % 3 == 1):
                            view *= 0

                        menger(view, quotient, 3)
        
                
def getCoordinates(matrix, dimension):
    shape = np.shape(matrix)
    
    xPositions = []
    yPositions = []
    zPositions = []
    
    if dimension == 2:
        for x in range(shape[0]):
            for y in range(shape[1]):
                if matrix[x,y] == 1:
                    xPositions = np.append(xPositions, x)
                    yPositions = np.append(yPositions, y)
                    zPositions = np.append(zPositions, 0)
                    
    if dimension == 3:
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    if matrix[x,y,z] == 1:
                        xPositions = np.append(xPositions, x)
                        yPositions = np.append(yPositions, y)
                        zPositions = np.append(zPositions, z)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(xPositions, yPositions, zPositions)
    plt.show()
    
    # Save positions in collumns in csv file
    positions = np.transpose([xPositions*2, yPositions*2, zPositions*2])

    
    path = "../rough Shapes/mengerShapes/"
    filename = input("Please enter a file name: ")
    np.savetxt(path + filename + ".csv", positions, delimiter=',')
    
    return positions
                
                

if __name__ == "__main__":

    SIZE = 27
    
    dimension = 3 # 2 for flat (2d) 3 for (3d)
    
    if dimension == 2:

        matrix = np.ones((SIZE, SIZE))

        menger(matrix, SIZE, dimension)

        plt.matshow(matrix)
        plt.colorbar()
        plt.show()
    
        positions = getCoordinates(matrix, 2)
        
    elif dimension == 3:
        
        matrix = np.ones((SIZE, SIZE, SIZE))
        
        menger(matrix, SIZE, dimension)
        
        positions = getCoordinates(matrix, 3)
        