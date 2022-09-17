import matplotlib.pyplot as plt
import numpy as np
from skimage import measure

channel = np.load("InnerMask.npy")
raw = np.load("InnerRaw.npy")

i = 19
while i < 360:
    img = channel[i]
    img2 = raw[i]
    background = plt.imread(f"AgateDataToContour/{(i + 1) // 20:02d}.png")

    contours_ch = []
    isovalues = np.linspace(10, 40, 20)
    for val in isovalues:
        a = measure.find_contours(img, val)
        contours_ch.extend(a.copy())

    contours = []
    for val in isovalues:
        a = measure.find_contours(img2, val)
        contours.extend(a.copy())

    plt.imshow(background)
    for contour in contours:
        plt.plot(contour[:, 1], contour[:, 0], color="red", linewidth=0.5)
    for contour in contours_ch:
       plt.plot(contour[:, 1], contour[:, 0], color="blue", linewidth=0.5)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"Pictures4/{(i + 1) // 20:02d}.png", bbox_inches=0, dpi=300)
    plt.close()
    i += 20
