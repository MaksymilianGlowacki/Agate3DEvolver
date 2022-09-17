import numpy as np
import interpolation as intp
from skimage.measure import find_contours
from trimesh import Trimesh
from trimesh.repair import broken_faces
import os
import matplotlib.pyplot as plt
from tqdm import tqdm

img_dir = "AgateDataToContour"
images = []
with os.scandir(img_dir) as entries:
    for entry in entries:
        if entry.name[-4:] == ".png" and len(entry.name) == 6:
            images.append(entry.path)

images.sort()
im1 = plt.imread(images[0])[:, :, 3]
imsize = im1.shape


meshes = []
F = []
cont1r = np.roll(find_contours(plt.imread(images[0])[:, :, 3])[0], 1, axis=1)
cont1 = np.concatenate((cont1r, np.ones((len(cont1r), 1)) * 40), axis=1)
magic_point = np.mean(cont1, axis=0)
magic_point[2] -= 10
V = np.concatenate(([magic_point], cont1))
offset = 0
number_of_layers = len(images) - 1
for i in tqdm(range(number_of_layers)):
    cont2r = np.roll(find_contours(plt.imread(images[i + 1])[:, :, 3])[0], 1, axis=1)
    cont2 = np.concatenate((cont2r, np.ones((len(cont2r), 1)) * 20 * (i + 2) + 20), axis=1)

    first_layer = (i == 0)
    last_layer = (i == (number_of_layers - 1))

    v, f = intp.new_interpolate_between_two(cont1, cont2, 10, first_layer=first_layer, last_layer=last_layer)
    F.append(np.array(f) + offset)
    offset = len(V)
    V = np.concatenate((V, v))

    cont1 = v

final_mesh = Trimesh(vertices=V, faces=np.concatenate(F), validate=True, process=True)
final_mesh.fill_holes()
broken_faces(final_mesh, 'red')
final_mesh.export("AgateContour.ply")
