from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import random

from system_classes import *


def __main__():
    number_of_points = int(input('Enter number of points: '))
    points_tuple_x = [random.uniform(-1, 3) for _ in range(number_of_points)]
    points_tuple_y = [random.uniform(-1, 3) for _ in range(number_of_points)]
    points_tuple_vx0 = [random.uniform(-1, 3) for _ in range(number_of_points)]
    points_tuple_vy0 = [random.uniform(-1, 3) for _ in range(number_of_points)]

    # print(points_tuple_x)
    # print(points_tuple_y)
    # print(points_tuple_vx0)
    # print(points_tuple_vy0)

    points = [Point2D(x0, y0, Vx0, Vy0) for x0, y0, Vx0, Vy0 in
              zip(points_tuple_x, points_tuple_y, points_tuple_vx0, points_tuple_vy0)]

    tube = TubeBorder(3, 1, 1)

    fig = plt.figure(num='Moving Points, Laboratory Work I, Ian Ilyasov')
    ax = fig.add_subplot(1, 1, 1)
    ax.axis('equal')
    # ax.plot(points_tuple_x, points_tuple_y, marker='o') -- Draws points as a polyline
    # drawn_points = []
    for point in points:
        point.draw_point(ax)

    tube.draw_tube_border(ax)

    global t, x_now, y_now, vx_now, vy_now
    t = 0
    dt = 0.01

    def calculate_new_movement_equation(vx, vy):
        dx = vx
        dy = vy
        dvx = 0
        dvy = 0  # -9.81
        return [dx, dy, dvx, dvy]

    x_t0 = [point.x0 for point in points]
    y_t0 = [point.y0 for point in points]
    vx_t0 = [point.Vx0 for point in points]
    vy_t0 = [point.Vy0 for point in points]

    x_now = numpy.array(x_t0)
    y_now = numpy.array(y_t0)
    vx_now = numpy.array(vx_t0)
    vy_now = numpy.array(vy_t0)

    def calculate_new_points(i):
        global t, x_now, y_now, vx_now, vy_now
        t += dt
        [d_x, d_y, d_vx, d_vy] = calculate_new_movement_equation(vx_now,
                                                                 vy_now)

        x_new = x_now + dt * d_x
        y_new = y_now + dt * d_y
        vx_new = vx_now + dt * d_vx
        vy_new = vy_now + dt * d_vy

        is_point_struck = tube.check_border_strike(x_new, y_new)

        for i in range(len(is_point_struck)):
            if is_point_struck[i]:
                n_x = (tube.x_center - x_now[i]) / (
                        ((tube.x_center - x_now[i]) ** 2 + (tube.y_center - y_now[i]) ** 2) ** 0.5)
                n_y = (tube.y_center - y_now[i]) / (
                        ((tube.x_center - x_now[i]) ** 2 + (tube.y_center - y_now[i]) ** 2) ** 0.5)
                vx_new[i] = -vx_now[i] * (n_x ** 2 - n_y ** 2) - 2 * vy_now[i] * n_x * n_y
                vy_new[i] = vy_now[i] * (n_x ** 2 - n_y ** 2) - 2 * vx_now[i] * n_x * n_y
                # print(vx_now[i], vy_now[i], vx_new[i], vy_new[i])
                # print(vx_now[i] * n_x + vy_now[i] * n_y, vx_new[i] * n_x + vy_new[i] * n_y)
                # print(vx_now[i] * n_y - vy_now[i] * n_x, vx_new[i] * n_y - vy_new[i] * n_x)
                x_new[i] = x_now[i] + dt * vx_new[i]
                y_new[i] = y_now[i] + dt * vy_new[i]

        is_point_struck_with_point = implement_strike_to_points(x_new, y_new)
        k = 0

        for i in numpy.arange(len(x_new) - 1):
            for j in numpy.arange(i + 1, len(x_new)):

                if is_point_struck_with_point[i, j]:
                    vx1, vy1 = vx_now[i], vy_now[i]
                    vx2, vy2 = vx_now[j], vy_now[j]

                    n_x = (x_now[j] - x_now[i]) / (((x_now[j] - x_now[i]) ** 2 + (y_now[j] - y_now[i]) ** 2) ** 0.5)
                    n_y = (y_now[j] - y_now[i]) / (((x_now[j] - x_now[i]) ** 2 + (y_now[j] - y_now[i]) ** 2) ** 0.5)

                    denominator = 2 * (n_x ** 2 + n_y ** 2)

                    vx_new1 = (n_x ** 2 * vx1 + n_x ** 2 * vx2 - n_x * n_y * vy1 + n_x * n_y * vy2 + 2 * n_y ** 2 * vx1 -
                              n_x * k * (n_x * vx1 - n_x * vx2 + n_y * vy1 - n_y * vy2)) / denominator

                    vx_new2 = (n_x ** 2 * vx1 + n_x ** 2 * vx2 + n_x * n_y * vy1 - n_x * n_y * vy2 + 2 * n_y ** 2 * vx2 +
                              n_x * k * (n_x * vx1 - n_x * vx2 + n_y * vy1 - n_y * vy2)) / denominator

                    vy_new1 = (n_y ** 2 * vy1 + n_y ** 2 * vy2 - n_x * vx1 * n_y + n_x * n_y * vx2 + 2 * n_x ** 2 * vy1 -
                              n_y * k * (n_x * vx1 - n_x * vx2 + n_y * vy1 - n_y * vy2)) / denominator

                    vy_new2 = (n_y ** 2 * vy1 + n_y ** 2 * vy2 + n_x * vx1 * n_y - n_x * n_y * vx2 + 2 * n_x ** 2 * vy2 +
                              n_y * k * (n_x * vx1 - n_x * vx2 + n_y * vy1 - n_y * vy2)) / denominator

                    x_new[i] = x_now[i] + dt * vx_new1
                    y_new[i] = y_now[i] + dt * vy_new1
                    x_new[j] = x_now[j] + dt * vx_new2
                    y_new[j] = y_now[j] + dt * vy_new2

        x_now = x_new
        y_now = y_new
        vx_now = vx_new
        vy_now = vy_new

        for inner_point, xi, yi in zip(points, x_new, y_new):
            inner_point.redraw_point(xi, yi)

        return [point_.plot_point for point_ in points]

    anim = FuncAnimation(fig, calculate_new_points, interval=dt * 1000, blit=True)

    plt.show()


if __name__ == '__main__':
    __main__()
