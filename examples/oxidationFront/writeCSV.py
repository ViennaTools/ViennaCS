import numpy as np
import sys

# ------------------------------------------------------------
# USER PARAMETERS
# ------------------------------------------------------------
gridDelta = 1.0         # nm (Matched to C++ default)
xExtent   = 150.0       # nm (Full length as requested)
yExtent   = 150.0       # nm (For 3D)
dimensions = 3          # 2 or 3

# Parse dimensions from command line
if len(sys.argv) > 1:
    try:
        dimensions = int(sys.argv[1])
    except ValueError:
        print(f"Invalid dimension argument: {sys.argv[1]}")
else:
    print("No dimension argument provided. Defaulting to 3D.")

# Output file
filename = "Efield.csv"

# ------------------------------------------------------------
# BUILD GRID
# ------------------------------------------------------------
nx = int(xExtent / gridDelta)
ny = int(yExtent / gridDelta)

# X: Lateral coordinates centered at 0 (-10 to 10)
xs = (np.arange(nx) + 0.5) * gridDelta - (xExtent / 2.0)
ys = (np.arange(ny) + 0.5) * gridDelta - (yExtent / 2.0)

# ------------------------------------------------------------
# FIELD CALCULATION
# ------------------------------------------------------------
def compute_E(x, y=0):
    """
    Calculates 1D E-field varying along X.
    Only Ey component is non-zero.
    """
    # Random Ey field
    Ey = np.random.uniform(0.0, 5.0)

    return Ey

# ------------------------------------------------------------
# GENERATE FIELD FILE
# ------------------------------------------------------------
print(f"Generating {filename} for {dimensions}D...")
if dimensions == 3:
    print(f"Grid: {nx} x {ny} points (2D plane X-Y)")
else:
    print(f"Grid: {nx} points (1D along X)")
print(f"X Range: [{xs[0]:.2f}, {xs[-1]:.2f}]")

with open(filename, "w") as f:
    if dimensions == 3:
        for y in ys:
            for x in xs:
                Ey = compute_E(x, y)
                f.write(f"{x:.6f}, {y:.6f}, {Ey:.6f}\n")
    else:
        for x in xs:
            Ey = compute_E(x)
            f.write(f"{x:.6f}, {Ey:.6f}\n")

print("Done.")