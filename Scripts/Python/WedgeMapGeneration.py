import numpy as np
import mrcfile
import matplotlib.pyplot as plt

plt.clf()  # Clear current figure
plt.close()  # Close current figure if needed

# Load the slopeGrid array from your .npy file
wedgeSlopeGrid = np.load('Scripts/Python/WedgeSlopeGrid.npy')
wedgeSlopeGrid = wedgeSlopeGrid.reshape((56, 57, 83))

print(np.min(wedgeSlopeGrid), np.max(wedgeSlopeGrid))

# Check its shape
print(f"Slope grid shape: {wedgeSlopeGrid.shape}")

# Export a single Z-slice as an image using matplotlib
z_index = 1  # choose a slice to visualize
slice_data = wedgeSlopeGrid[:, :, z_index]

plt.imshow(slice_data, cmap='seismic', interpolation='nearest')
plt.colorbar(label='Slope')
plt.title(f"Slope Grid Slice at Z={z_index}")
plt.xlabel('Y index')
plt.ylabel('X index')
plt.show()

# Transpose from (X, Y, Z) to (Z, Y, X)
slope_grid_transposed = np.transpose(wedgeSlopeGrid, (2, 1, 0))

# Save as MRC
with mrcfile.new('Output/SlopeMaps/WedgeSlopeGrid.mrc', overwrite=True) as mrc:
    mrc.set_data(slope_grid_transposed.astype(np.float32))
    # Optional: set voxel size or origin if you know them
    mrc.voxel_size = (2, 2, 2)
    mrc.header.origin = (-43, -27, -53)  # in Å
    print("Saved slopeGrid.mrc successfully.")
