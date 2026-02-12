import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def plot_efield(filename):
    """
    Plots Electric Field data from a CSV file.
    Supports 2 column (X, E) and 3 column (X, Y, E) formats.
    """
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found.")
        return

    print(f"Reading {filename}...")
    try:
        data = np.genfromtxt(filename, delimiter=',')
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    if data.ndim == 1:
        data = data.reshape(1, -1)

    num_cols = data.shape[1]
    print(f"Data shape: {data.shape}")

    if num_cols == 2:
        # 1D Plot (X vs E)
        x = data[:, 0]
        e = data[:, 1]
        
        plt.figure(figsize=(10, 6))
        plt.plot(x, e, '-', linewidth=2)
        plt.title(f"Electric Field Profile (1D)\nSource: {filename}")
        plt.xlabel("Position X")
        plt.ylabel("E-Field Magnitude")
        plt.grid(True, alpha=0.3)
        
        out_name = os.path.splitext(filename)[0] + ".png"
        plt.savefig(out_name, dpi=300)
        print(f"Plot saved to {out_name}")
        # plt.show()

    elif num_cols >= 3:
        # 2D Plot (X, Y, E)
        x = data[:, 0]
        y = data[:, 1]
        e = data[:, 2]
        
        plt.figure(figsize=(10, 8))
        
        # Use tricontourf which works for both structured and unstructured grids
        contour = plt.tricontourf(x, y, e, levels=50, cmap='viridis')
        cbar = plt.colorbar(contour)
        cbar.set_label('E-Field Magnitude')
        
        plt.title(f"Electric Field Distribution (2D)\nSource: {filename}")
        plt.xlabel("Position X")
        plt.ylabel("Position Y")
        plt.axis('equal') # Ensure aspect ratio is correct for spatial data
        
        out_name = os.path.splitext(filename)[0] + ".png"
        plt.savefig(out_name, dpi=300)
        print(f"Plot saved to {out_name}")
        # plt.show()

    else:
        print(f"Error: Unexpected number of columns ({num_cols}). Expected 2 or 3.")

if __name__ == "__main__":
    # Default to Efield.csv if no argument provided
    target_file = "Efield.csv"
    if len(sys.argv) > 1:
        target_file = sys.argv[1]
    
    plot_efield(target_file)