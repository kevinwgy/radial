import os;
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------
# Plot: Interpolation
# ---------------------------------------------
x = np.linspace(0, 1, 200);
y = np.linspace(0, 1, 200);
xy = np.array([[i,j] for j in y for i in x]);
np.savetxt('tmp_input_points.txt', xy);

subprocess.run(["../../radial", "2", "1", "training_data.txt", "tmp_input_points.txt", "40", "MultiQuadric", "tmp_rbf_result.txt"]);

z = np.genfromtxt("tmp_rbf_result.txt");
X, Y = np.meshgrid(x,y);
Z = z.reshape(200,200,3);
Z = Z[:,:,2];
plt.figure(figsize=(12,10));
plt.pcolormesh(X,Y,Z,vmin=-1.0,vmax=1.0,shading='auto');
plt.colorbar();
plt.savefig('rbf_interp.png')
plt.draw();

# clean-up
subprocess.run(["rm", "tmp_input_points.txt"]);
subprocess.run(["rm", "tmp_rbf_result.txt"]);
print("Deleted files tmp_input_points.txt and tmp_rbf_result.txt.");


# ---------------------------------------------
# Plot: Extrapolation
# ---------------------------------------------
x = np.linspace(-0.5, 1.5, 200);
y = np.linspace(-0.5, 1.5, 200);
xy = np.array([[i,j] for j in y for i in x]);
np.savetxt('tmp_input_points.txt', xy);

subprocess.run(["../../radial", "2", "1", "training_data.txt", "tmp_input_points.txt", "40", "MultiQuadric", "tmp_rbf_result.txt"]);

z = np.genfromtxt("tmp_rbf_result.txt");
X, Y = np.meshgrid(x,y);
Z = z.reshape(200,200,3);
Z = Z[:,:,2];
plt.figure(figsize=(12,10));
plt.pcolormesh(X,Y,Z,vmin=-1.0,vmax=1.0,shading='auto');
plt.colorbar();
plt.savefig('rbf_extrap.png')
plt.draw();

# clean-up
subprocess.run(["rm", "tmp_input_points.txt"]);
subprocess.run(["rm", "tmp_rbf_result.txt"]);
print("Deleted files tmp_input_points.txt and tmp_rbf_result.txt.");

plt.show();
