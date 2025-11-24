import numpy as np

# ------------------------------------------------------------
# USER PARAMETERS
# ------------------------------------------------------------
gridDelta = 0.85         # nm
xExtent   = 85           # nm
yExtent   = 50           # nm
substrateHeight = 50     # nm
holeRadius      = 20     # nm

# Angle Setup
angle_deg = 45.0
angle_rad = np.radians(angle_deg)

# Output file
filename = "Efield.csv"

# ------------------------------------------------------------
# BUILD GRID
# ------------------------------------------------------------
nx = int(xExtent / gridDelta)
ny = int(yExtent / gridDelta)

# Cell centers
xs = (np.arange(nx) + 0.5) * gridDelta
ys = (np.arange(ny) + 0.5) * gridDelta

# ------------------------------------------------------------
# FIELD CALCULATION
# ------------------------------------------------------------
def compute_E_vector(x, y):
    """
    Calculates E-field magnitude and splits it into vector components
    based on the incidence angle.
    """
    
    # 1. CALCULATE SHIFT (GEOMETRIC SLANT)
    # We assume the beam enters at the substrate surface (y = substrateHeight).
    # At 45 degrees, for every nm we go down (negative y direction), 
    # the beam shifts sideways.
    # If angle is +45 deg (coming from top-left), x increases as y decreases.
    
    dist_from_surface = substrateHeight - y
    
    # Calculate lateral shift based on depth and angle
    # tan(theta) = opp/adj = dx/dy  -> dx = dy * tan(theta)
    lateral_shift = dist_from_surface * np.tan(angle_rad)
    
    # The "effective" center of the beam shifts as we go deeper
    original_center = xExtent / 2.0
    beam_center_at_y = original_center - lateral_shift 

    # 2. MAGNITUDE CALCULATION (Similar to original, but shifted)
    E0 = 1.0
    
    # Check if we are inside the slanted beam
    inside_trench = abs(x - beam_center_at_y) < holeRadius

    # Vertical coordinate relative to substrate
    depth = y - substrateHeight

    # Decay logic
    if depth < 0:
        decay = np.exp(depth / 10.0) 
    else:
        decay = 1.0

    boost = 3.0 if inside_trench else 1.0
    
    E_magnitude = E0 * boost * decay

    # 3. VECTOR DECOMPOSITION
    # We need to split the magnitude into Ex and Ey.
    # Assuming the field points "down" into the substrate:
    # Ey component is negative (pointing down).
    # Ex component is positive (pointing right).
    
    Ex = E_magnitude * np.sin(angle_rad)
    Ey = -E_magnitude * np.cos(angle_rad) 
    Ez = 0.0

    return Ex, Ey, Ez

# ------------------------------------------------------------
# GENERATE FIELD FILE
# ------------------------------------------------------------
with open(filename, "w") as f:
    # Header (Optional - many solvers prefer raw data, 
    # but standard CSVs often have headers. Comment out if not needed.)
    # f.write("Ex,Ey,Ez\n") 
    
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            Ex, Ey, Ez = compute_E_vector(x, y)
            
            # Write vector components instead of just scalar E
            f.write(f"{Ex:.6f}, {Ey:.6f}, {Ez:.6f}\n")

print(f"Generated {filename} with {nx*ny} entries (Vector Format: Ex, Ey, Ez).")
print(f"Angle: {angle_deg} degrees.")