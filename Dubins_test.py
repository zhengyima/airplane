from DubinsAirplaneFunctions import * 
from PlottingTools import plot3
import numpy as np
import time
import sys

pi = np.pi
dubins_case = 0
verbose_flag = 0
plot_flag = 16
R_min = 200
Gamma_max = pi/2

def get_dubin_L(x1,y1,z1,theta1,x2,y2,z2,theta2):
    # start_node  = np.array( [0,   0,   -100,   0*pi/180,    15] )
    # end_node = np.array( [0, 200,   -125, 270*pi/180,    15] )

    start_node  = np.array( [x1,   y1,   z1,   theta1,    200] )
    end_node = np.array( [x2, y2,   z2, theta2,    200] )
    # t0 = time.clock()
    DubinsAirplaneSolution = DubinsAirplanePath( start_node, end_node, R_min, Gamma_max )

    return DubinsAirplaneSolution['L'],DubinsAirplaneSolution