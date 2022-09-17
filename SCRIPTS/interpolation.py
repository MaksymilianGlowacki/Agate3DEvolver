import numpy as np
from trimesh.proximity import signed_distance
from tqdm import tqdm


def distance(point1, point2):
    return np.sqrt(np.sum(np.power(point1 - point2, 2)))


def elongate(array):
    extras = []
    for index in range(array.shape[0] - 1):
        extras.append(np.mean(array[index:index + 2], axis=0))
    extras.append(np.mean(np.array([array[0], array[-1]]), axis=0))
    return np.concatenate((array, np.array(extras)), axis=1).reshape((-1, 3))


def new_interpolate_between_two(left_contour, right_contour, number_of_combinations=15, first_layer=False, last_layer=False):
    n = len(right_contour) // len(left_contour)
    m = (len(right_contour) % len(left_contour))
    if m > 0:
        combs = np.random.randint(len(left_contour), size=(number_of_combinations, m - 1))
        combs = np.concatenate((combs, np.zeros((number_of_combinations, 1), dtype=int)), axis=1)
    else:
        combs = np.random.randint(len(left_contour), size=(number_of_combinations, m))

    length_table = []
    for combination in combs:
        connections = np.ones(len(left_contour), dtype=int) * n
        for pos in combination:
            connections[pos] += 1

        lengths = []
        for conformation in range(len(right_contour)):
            index = 0
            total_length = 0
            rolled = np.roll(right_contour, conformation, axis=0)
            for count, point in enumerate(left_contour):
                total_length += distance(point, rolled[index - 1])
                for j in range(connections[count]):
                    total_length += distance(point, rolled[index + j])
                index += connections[count]
            lengths.append(total_length)
        length_table.append(lengths)

    optimal = np.argwhere(length_table == np.amin(length_table))[0]
    connections = np.ones(len(left_contour), dtype=int) * n
    for pos in combs[optimal[0]]:
        connections[pos] += 1
    right_contour = np.roll(right_contour, optimal[1], axis=0)

    vertices = np.concatenate((left_contour, right_contour))
    faces = []
    index = len(left_contour)

    faces.append([0, len(left_contour) - 1, len(vertices) - 1])
    faces.append([0, len(vertices) - 1, index])
    for j in range(connections[0] - 1):
        faces.append([0, index + j, index + j + 1])
    index += connections[0]

    for counter in range(1, len(left_contour)):
        faces.append([counter, counter - 1, index - 1])
        for j in range(connections[counter]):
            faces.append([counter, index + j - 1, index + j])
        index += connections[counter]

    faces = np.array(faces)

    if first_layer:
        magic_point = np.mean(left_contour, axis=0)
        magic_point[2] -= 20
        faces += 1
        vertices = np.concatenate(([magic_point], vertices))
        faces = np.concatenate((faces, [[0, len(left_contour), 1]]))
        for vertex in range(len(left_contour) - 1):
            faces = np.concatenate((faces, [[0, vertex + 1, vertex + 2]]))

    if last_layer:
        magic_point = np.mean(right_contour, axis=0)
        magic_point[2] += 10
        vertices = np.concatenate((vertices, [magic_point]))
        faces = np.concatenate((faces, [[len(vertices) - 1, len(vertices) - 2, len(left_contour) + first_layer]]))
        for vertex in range(len(right_contour) - 1):
            faces = np.concatenate((faces, [[len(vertices) - 1, len(vertices) - 3 - vertex, len(vertices) - 2 - vertex]]))

    return vertices[len(left_contour) + first_layer:], faces


def make_it_smooth(mesh, imsize=None):
    max_x = int(np.amax(mesh.vertices[:, 0]))
    max_y = int(np.amax(mesh.vertices[:, 1]))
    max_z = int(np.amax(mesh.vertices[:, 2]))

    if imsize is not None:
        matrix = np.zeros((max_z + 1, imsize[0], imsize[1]))
    else:
        matrix = np.zeros((max_z + 1, max_y, max_x))

    for z in tqdm(range(matrix.shape[0])):
        for y in range(matrix.shape[1]):
            for x in range(matrix.shape[2]):
                matrix[z, y, x] = signed_distance(mesh, [(x, y, z)])[0]
    np.save("AgateMatrix", matrix)
    return matrix
