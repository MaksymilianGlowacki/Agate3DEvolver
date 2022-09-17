import numpy as np
import cv2 as cv
import os
from skimage.measure import find_contours


img_dir = "AgateDataToContour"
images_list = []
with os.scandir(img_dir) as entries:
    for entry in entries:
        if entry.name[-4:] == ".png":
            images_list.append(entry.path)
images_list.sort()

images = []
for entry in images_list:
    im = cv.imread(entry, cv.IMREAD_UNCHANGED)
    images.append(im.copy())
images = np.asarray(images)
matrix = np.load("Outer.npy")

num = 9
current_image = images[num]
current_image_copy = images[num]

def on_trackbar(val):
    val *= np.max((int(np.amax(matrix[num * 20 + 20])), 1)) / 100
    bcontours = find_contours(matrix[num * 20 + 20], val)
    global current_image_copy
    if len(bcontours) == 1:
        bcontours = np.roll(bcontours, 1, axis=2).astype(np.int32)
        current_image_copy = current_image.copy()
        cv.drawContours(current_image_copy, bcontours, 0, (0, 0, 255), 1)


cv.namedWindow("3D AgateEvolver", cv.WINDOW_GUI_NORMAL)
cv.createTrackbar("Surface value", "3D AgateEvolver", 0, 100, on_trackbar)
while True:
    cv.imshow("3D AgateEvolver", current_image_copy)
    button = cv.waitKey(30)
    if button == 27:
        break
    elif button == 46:
        num += 1
        num %= len(images)
        current_image = images[num]
        current_image_copy = images[num]
    elif button == 44:
        num -= 1
        num %= len(images)
        current_image = images[num]
        current_image_copy = images[num]
