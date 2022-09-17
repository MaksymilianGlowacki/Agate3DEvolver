import matplotlib.pyplot as plt
from matplotlib import colors, cm
import numpy as np
import cv2 as cv

outer = np.load("OuterMask.npy")

i = 19
while i < 360:
    img = outer[i]
    background = cv.imread(f"AgateDataToContour/{(i + 1) // 20:02d}.png", cv.IMREAD_UNCHANGED)
    img2 = np.where(img > 0, np.power(np.sin(img), 2), 0)
    cmap: colors.Colormap = cm.get_cmap('ocean')
    img2 = cmap(img2)
    img = img.reshape((210, 300, 1))
    img = np.concatenate((img, img, img, img), axis=2)
    img2[:, :, 3] = 255
    img2[:, :, :-1] *= 255

    result = np.where(img > 33, img2, background)
    cv.imwrite(f"Pictures/{(i + 1) // 20:02d}.png", result)

    i += 20
