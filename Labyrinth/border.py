import numpy as np
import random as rd

def generate_rectangle(x_len, y_len):
    center_x = x_len / 2
    center_y = y_len / 2
    ox = center_x + (0.5 - rd.random()) * x_len / 2
    oy = center_y + (0.5 - rd.random()) * y_len / 2

    angle = rd.random() * np.pi

    x_l, y_b = rd.random(), rd.random()
    x_r, y_t = rd.random(), rd.random()

    bounds = np.array([r([x_l, y_t], angle, [ox, oy]), r([x_r, y_t], angle, [ox, oy]),
                       r([x_r, y_b], angle, [ox, oy]), r([x_l, y_b], angle, [ox, oy])]), \
             np.array([r([x_l, y_b], angle, [ox, oy]), r([x_l, y_t], angle, [ox, oy]),
                       r([x_r, y_t], angle, [ox, oy]), r([x_r, y_b], angle, [ox, oy])])

    return bounds

def r(point, angle, center):
    return [center[0] + (point[0] * np.cos(angle) - point[1] * np.sin(angle)),
            center[1] + (point[0] * np.sin(angle) + point[1] * np.cos(angle))]
