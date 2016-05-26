# Simulation constants:
import numpy as np

# Fundamental constants:
Q_EL = 1.602e-19                # electron charge 
M_EL = 9.1e-31                  # electron mass
E_EL = 5.105396765648739e5      # electron... uh... 
MU0  = np.pi*4e-7     
EPS0 = 8.854e-12
C    = 2.997956376932163e8      # Speed of light
B0   = 3.12e-5                  # Magnetic field strength (Teslas)

# Degrees, radians
D2R = np.pi/180.0
R2D = 180.0/np.pi

R_E = 6378e3    # Earth radius in meters
H_IONO = 1e5    # Height of bottom of ionosphere, in meters
H_MAGNETO = 1e6 # Height of top of ionosphere, in meters

# For input power scaling:
Z0 = 377.0      # Free-space impedance
A = 5e3         # Constant for lightning power calculation
B = 1e5         # Constant for lightning power calculation
H_E = 5000.0    # Lightning incident height, in meters? (Confirm this pls)

# Scattering code simulation params:
T_MAX = 20 #5.0        #30.0
NUM_STEPS =  2000 #30000
T_STEP = (1.0*T_MAX)/NUM_STEPS

# Degrees latitude around flash latitude that we'll use rays from:
LAT_SPREAD = 2

# Number of steps in fine-grid interpolation:
#DIV_LAT_NUM = 2 #6   # 2 does endpoints only
#DIV_FREQ_NUM = 11 #6

# Spacing between final time-frequency grids:
DT = T_STEP     # No extra interpolation here
#DF = 50         # Hz. REWORK THIS ONCE YOU ADD IN MULTIPLE FREQUENCIES AUSTIN
F_STEP = 40     # Hz. Separation in frequency. (i.e., do 1-hz interpolation)
LAT_STEP = 0.1    # Degrees interpolation of launch rays

# EA array grid settings:
# (EA = "Equal Area" -- defines a set of equal-area latitude slices
#  to compute scattering at. All values in degrees latitude.)
EALimS = -50.0  # Southern limit (degrees)
EALimN = 50.0   # Northern limit (degrees)
EAIncr = 1.0    # Step size (degrees)

# Size of bounding boxes for crossing detection:
# (Width around a field line, in L-shells, in which to consider a crossing)
L_MARGIN = 0.1 #6e-4  #0.1 #(in L-shells)


# Energy bins for resonant scattering calculation:
E_MIN = 1e1     # Minimum energy (ev)
E_MAX = 1e8     # Maximium energy (ev)
NUM_E = 512     # Number of bins in energy vector

E_BANDWIDTH = 0.3   # width plus/minus energy resonance to calculate for. Jacob hard-coded to 0.3.
#SQUARE = True   # 

# ---- Generated parameters
E_EXP_BOT = np.log10(E_MIN)
E_EXP_TOP = np.log10(E_MAX)
DE_EXP = ( (E_EXP_TOP - E_EXP_BOT)/NUM_E)

# Generate energy and velocity arrays
E_tot_arr = pow(10,E_EXP_BOT + DE_EXP*np.arange(0,NUM_E))
v_tot_arr = C*np.sqrt(1 - pow(E_EL/(E_EL + E_tot_arr),2))






