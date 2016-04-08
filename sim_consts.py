# Simulation constants:
import numpy as np

D2R = np.pi/180.0
R2D = 180.0/np.pi


R_E = 6378e3    # Earth radius in meters
H_IONO = 1e5    # Height of bottom of ionosphere, in meters
H_MAGNETO = 1e6 # Height of top of ionosphere, in meters

# For input power scaling:
Z0 = 377.0      # Free-space impedance
A = 5e3         # Constant for lightning power calculation
B = 1e5         # Constant for lightning power calculation
H_E = 5000      # Lightning incident height, in meters? (Confirm this pls)

# Scattering code simulation params:
T_MAX = 30
NUM_STEPS =  30000 #30000
T_STEP = T_MAX/NUM_STEPS


# Number of steps in fine-grid interpolation:
DIV_LAT_NUM = 5
DIV_FREQ_NUM = 3

# EA array grid settings:
# (EA = "Equal Area" -- defines a set of equal-area latitude slices
#  to compute scattering at. All values in degrees latitude.)
EALimS = -45.0  # Southern limit (degrees)
EALimN = 45.0   # Northern limit (degrees)
EAIncr = 1    # Step size (degrees)

# DL0: 
#DL0 = 6e-4    #(in L-shells)

# Size of bounding boxes for crossing detection:
# (Width around a field line, in L-shells, in which to consider a crossing)
L_MARGIN = 0.1 #(in L-shells)
