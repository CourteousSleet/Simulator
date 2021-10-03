import numpy


class TubeBorder:
    radius = 0.1
    x_center = 0
    y_center = 0

    def __init__(self, radius=0.1, x_center=0, y_center=0):
        self.radius = radius
        self.x_center = x_center
        self.y_center = y_center

    def check_border_strike(self, x, y):
        return (x - self.x_center) ** 2 + (y - self.y_center) ** 2 > self.radius * 2

    def draw_tube_border(self, ax):
        phi = numpy.linspace(0, 6.28, 10000)
        x_of_tube = self.x_center + self.radius * numpy.cos(phi)
        y_of_tube = self.y_center + self.radius * numpy.sin(phi)
        ax.plot(x_of_tube, y_of_tube, 'black')


class Point2D:
    x0 = 0
    y0 = 0
    Vx0 = 0
    Vy0 = 0
    point_radius = 0.1
    coord = (x0, y0)
    plot_point = []

    def __init__(self, x0, y0, v_x0, v_y0):
        self.x0 = x0
        self.y0 = y0
        self.coord = (x0, y0)
        self.Vx0 = v_x0
        self.Vy0 = v_y0
        self.velocity_vector = (v_x0, v_y0)

    def draw_point(self, ax):
        # noinspection PyAttributeOutsideInit
        self.plot_point, = ax.plot(self.x0, self.y0, marker='o')

    def redraw_point(self, x, y):
        self.plot_point.set_data(x, y)


def implement_strike_to_points(x_new, y_new):
    strikes = numpy.zeros([len(x_new), len(x_new)], bool)
    for i in numpy.arange(len(x_new) - 1):
        for j in numpy.arange(i + 1, len(x_new)):
            strikes[i, j] = (x_new[i] - x_new[j]) ** 2 + (y_new[i] - y_new[j]) ** 2 < (2 * 0.1) ** 2

    return strikes
