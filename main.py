from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import random as rand
import numpy as np


class TubeBorder:
    radius = 0.1
    x_center = 0
    y_center = 0

    def __init__(self, radius=1, x_center=0, y_center=0):
        self.radius = radius
        self.x_center = x_center
        self.y_center = y_center

    def check_border_strike(self, x, y):
        return (x - self.x_center) ** 2 + (y - self.y_center) ** 2 > self.radius * 2

    def draw_tube(self, ax):
        phi = np.linspace(0, 6.28, 10000)
        x_of_tube = self.x_center + self.radius * np.cos(phi)
        y_of_tube = self.y_center + self.radius * np.sin(phi)
        ax.plot(x_of_tube, y_of_tube, marker='o')

    # def tube_strike_handler(self, ):


class Point2D(object):
    x0 = 0
    y0 = 0
    # z0 = 0
    Vx0 = 0
    Vy0 = 0
    point_radius = 0.1
    coord = (x0, y0)

    def __init__(self, x0, y0, v_x0, v_y0, r):
        self.x0 = x0
        self.y0 = y0
        self.coord = (x0, y0)
        self.Vx0 = v_x0
        self.Vy0 = v_y0
        self.velocity_vector = (v_x0, v_y0)
        self.point_radius = r

    def __init__(self, x0, y0, v_x0, v_y0):
        self.x0 = x0
        self.y0 = y0
        self.coord = (x0, y0)
        self.Vx0 = v_x0
        self.Vy0 = v_y0
        self.velocity_vector = (v_x0, v_y0)

    # def set_data(self, x_new, y_new):
    #    self.x0 = x_new
    #    self.y0 = y_new
    #    self.coord = (self.x0, self.y0)


def print_hi(name):
    print(f'{name}')


def __main__():
    number_of_points = int(input('Enter number of points '))
    points_tuple_x = [rand.random() for i in range(number_of_points)]
    points_tuple_y = [rand.random() for i in range(number_of_points)]
    points_tuple_vx0 = [rand.random() for i in range(number_of_points)]
    points_tuple_vy0 = [rand.random() for i in range(number_of_points)]
    # radius's_values = [rand.random() for i in range(number_of_points)]

    print(points_tuple_x)
    print(points_tuple_y)
    print(points_tuple_vx0)
    print(points_tuple_vy0)

    points = [Point2D(x0, y0, Vx0, Vy0) for x0, y0, Vx0, Vy0 in
              zip(points_tuple_x, points_tuple_y, points_tuple_vx0, points_tuple_vy0)]

    tube = TubeBorder(3, 1, 1)

    fig = plt.figure()
    # matplotlib.use('TkAgg')
    ax = fig.add_subplot(1, 1, 1)
    # ax.plot(points_tuple_x, points_tuple_y, marker='o') -- Draws points as a polyline
    drawn_points = []
    for point in points:
        p, = ax.plot(point.x0, point.y0, marker='o')  # Draws as separated points
        drawn_points.append(p)

    print(drawn_points)

    tube.draw_tube(ax)

    t = 0
    dt = 0.01

    def calculate_new_movement_equation(x, y, vx, vy, t_in):
        dx = vx
        dy = vy
        dvx = 0
        dvy = 0
        return [dx, dy, dvx, dvy]

    x_t0 = [point.x0 for point in points]
    y_t0 = [point.y0 for point in points]
    vx_t0 = [point.Vx0 for point in points]
    vy_t0 = [point.Vy0 for point in points]

    x_now = np.array(x_t0)
    y_now = np.array(y_t0)
    vx_now = np.array(vx_t0)
    vy_now = np.array(vy_t0)

    def calculate_new_points(i):
        global t, dt, x_now, y_now, vx_now, vy_now, drawn_points
        t += dt
        [d_x, d_y, d_vx, d_vy] = calculate_new_movement_equation(x_now, y_now, vx_now,
                                                                 vy_now, t)

        x_new = x_now + dt * d_x
        y_new = y_now + dt * d_y
        vx_new = vx_now + dt * d_vx
        vy_new = vy_now + dt * d_vy

        is_point_struck = tube.check_border_strike(x_new, y_new)
        print(is_point_struck)

        for i in range(is_point_struck):
            if is_point_struck[i]:
                vx_new[i] = 0
                vy_new[i] = 0

        x_now = x_new
        y_now = y_new
        vx_now = vx_new
        vy_now = vy_new

        for inner_point, xi, yi in zip(drawn_points, x_new, y_new):
            inner_point.set_data(xi, yi)

        return drawn_points

    anim = FuncAnimation(fig, calculate_new_points, interval=dt * 1000, blit=True)

    plt.show()


if __name__ == '__main__':
    __main__()
    print_hi('Finished')