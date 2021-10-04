from PyQt5.QtWidgets import *
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt
import math
import numpy
import sympy
import FormOfTheTorpedo
import oceanwidget

global mass, radius, a, b, bigger_radius, k, drag_force, side_drag_force, wheel_force, rho, v_vpd, t, state_vector


class Aim:
    x = 0
    y = 0
    size = 10
    is_destroyed = False

    def __init__(self, x, y, s_):
        time = sympy.lambdify('t')
        self.x = sympy.lambdify(time, x)
        self.y = sympy.lambdify(time, y)
        self.size = s_
        self.drawn_aim = None
        phi = numpy.linspace(0, 6.28, 20)
        self.aim_x = self.size / 2 * numpy.sin(phi)
        self.aim_y = self.size / 2 * numpy.cos(phi)

    def draw(self, axes):
        self.drawn_aim = axes.plot(self.x + self.aim_x, self.y + self.aim_y)[0]


class Torpedo:
    x = 0
    y = 0
    phi = 0
    a = 1
    b = 1
    alpha_max = 0.8
    radius = 1
    sensor_r = -1
    sensor_phi = -1
    is_detected = False
    is_exploded = False
    sensor_radius_max = 50
    sensor_phi_max = 0.5

    def __init__(self, x, y, phi, a_, b_, r_):
        self.pitch = None
        self.x = x
        self.y = y
        self.phi = phi
        self.a = a_
        self.b = b_
        self.radius = r_
        self.is_detected = 0
        self.is_exploded = 0
        self.drawn_torpedo = None
        self.torpedo_x = numpy.array([b, b * 1.05, b, -a * 0.9, -a, -a * 0.9, b])
        self.torpedo_y = numpy.array([bigger_radius, 0, -bigger_radius, -bigger_radius, 0, bigger_radius, bigger_radius])

    def draw(self, axes):
        r_torpedo_x, r_torpedo_y = rot2d(self.torpedo_x, self.torpedo_y, self.phi)
        self.drawn_torpedo = axes.plot(self.x + r_torpedo_x, self.y + r_torpedo_y)[0]

    def redraw(self):
        if self.is_exploded:
            self.torpedo_x *= 1.002
            self.torpedo_y *= 1.004
            self.drawn_torpedo.set_data(self.x + self.torpedo_x, self.y + self.torpedo_y)
        else:
            r_torpedo_x, r_torpedo_y = rot2d(self.torpedo_x, self.torpedo_y, self.phi)
            self.drawn_torpedo.set_data(self.x + r_torpedo_x, self.y + r_torpedo_y)

    #
    #                                                                             **
    #                                                                            ****
    #               __________________________                                 /  **
    #              /                          \                               /))  --  beta
    #    -------------------------------------->---------------------------------------------
    #          (  /\__________________________/
    # alpha --  (/
    #           /
    def choose_course(self, aim_):
        r_vector_x = aim_.x - self.x
        r_vector_y = aim_.y - self.y
        math.atan2(r_vector_y, r_vector_x)

        beta = math.atan2(r_vector_y, r_vector_x) - self.phi
        alpha = self.pitch * beta
        if abs(alpha) > self.alpha_max:
            alpha = self.alpha_max
        elif abs(alpha) < -self.alpha_max:
            alpha = -self.alpha_max

        return alpha

    def explode(self, axes):
        self.is_exploded = True
        n = 6
        phi_1 = numpy.linspace(0, 6.28, n + 1)
        explosion_point_x_1 = self.a / 2 * numpy.sin(phi_1)
        explosion_point_y_1 = self.a / 2 * numpy.cos(phi_1)
        phi_2 = numpy.linspace(0, 6.28, 2 * n + 1)
        explosion_point_x_2 = self.a / 2 * numpy.sin(phi_2)
        explosion_point_y_2 = self.a / 2 * numpy.cos(phi_2)

        self.torpedo_x = numpy.hstack([explosion_point_x_1, explosion_point_x_2])
        self.torpedo_y = numpy.hstack([explosion_point_y_1, explosion_point_y_2])
        self.drawn_torpedo = \
        axes.plot(self.x + self.torpedo_x, self.y + self.torpedo_y, marker='o', color=[255, 248, 231])[0]

    def check_for_strike(self, aim_, axes):
        x_head = self.x + b * numpy.cos(self.phi)
        y_head = self.y + b * numpy.sin(self.phi)
        if (x_head - aim_.x) ** 2 + (y_head - aim_.y) ** 2 < (aim_.size / 2) ** 2:
            self.explode(axes)


def rot2d(x, y, phi):
    rot_x = x * numpy.cos(phi) - y * numpy.sin(phi)
    rot_y = x * numpy.sin(phi) + y * numpy.cos(phi)
    return rot_x, rot_y


def calculate_movement_equation(state_vector_, alpha):
    global mass, radius, a, b, bigger_radius, k, drag_force, side_drag_force, wheel_force, rho, v_vpd
    x, y, phi, v_x, v_y, omega = state_vector_
    dx = v_x
    dy = v_y
    d_phi = omega

    v_l = v_x * numpy.cos(phi) + v_y * numpy.sin(phi)
    v_b = -v_x * numpy.sin(phi) + v_y * numpy.cos(phi)

    s_l = 3.14 * bigger_radius ** 2
    s_b = (a + b) * 2 * bigger_radius
    s_r = bigger_radius ** 2 * k

    j = mass * (a ** 2 + b ** 2) / 6 * radius

    f_sl = rho * s_l * drag_force * v_l ** 2 * numpy.sign(v_l) / 2
    f_sb = rho * s_b * side_drag_force * v_b ** 2 * numpy.sign(v_b) / 2
    m_s = rho * s_b * wheel_force * omega ** 2 * (b + a) ** 2 * numpy.sign(omega) / 8

    f_dv = rho * 3.14 * bigger_radius ** 2 * v_vpd ** 2
    f_r = rho * s_r * v_vpd ** 2 * (numpy.sin(alpha)) ** 2 * numpy.sign(numpy.sin(alpha)) / 2

    d_v_x = (f_dv * numpy.cos(phi) + f_sb * numpy.sin(phi) - f_sl * numpy.cos(phi) + f_r * numpy.sin(
        phi - alpha)) / mass
    d_v_y = (f_dv * numpy.sin(phi) - f_sb * numpy.cos(phi) - f_sl * numpy.sin(phi) - f_r * numpy.cos(
        phi - alpha)) / mass
    d_omega = (a * f_r * numpy.cos(alpha) - m_s) / j

    return numpy.array([dx, dy, d_phi, d_v_x, d_v_y, d_omega])


class OceanWidget(QMainWindow, FormOfTheTorpedo.Ui_MainWindow):

    def __init__(self):
        QMainWindow.__init__(self)

        self.setupUi(self)

        self.setWindowTitle("Piece of shit")

        self.FireButton.clicked.connect(self.draw_the_thing)

        self.horizontalScrollBar.valueChanged.connect(self.value_change)
        self.addToolBar(NavigationToolbar(self.OceanWidget.canvas, self))

    def draw_the_thing(self):
        global mass, radius, a, b, bigger_radius, k, drag_force, side_drag_force, wheel_force, rho, v_vpd, t, state_vector
        # Alpha = self.AngleBar.value()/5000
        #     Mass parameters
        mass = float(self.m_edit.text())
        radius = float(self.radius_edit.text())

        #     Geometric parameters
        a = float(self.a_edit.text())
        b = float(self.b_edit.text())
        bigger_radius = float(self.bigger_radius_edit.text())
        k = float(self.k_edit.text())

        #     Forces parameters
        drag_force = float(self.drag_force_edit.text())
        side_drag_force = float(self.side_drag_force_edit.text())
        wheel_force = float(self.c_vr_edit.text())

        rho = float(self.rho_edit.text())

        v_vpd = float(self.v_vpd_edit.text())

        x0 = float(self.x0_edit.text())
        y0 = float(self.y0_edit.text())
        phi0 = float(self.phi0_edit.text())
        v_x0 = float(self.v_x0_edit.text())
        v_y0 = float(self.v_y0_edit.text())
        omega_0 = float(self.omega0_edit.text())

        self.OceanWidget.canvas.axes.clear()
        self.OceanWidget.canvas.axes.grid(True)
        self.OceanWidget.canvas.axes.axis('scaled')
        self.OceanWidget.canvas.axes.set(xlim=[0, (a + b) * 20], ylim=[-(a + b) * 10, (a + b) * 10])
        self.OceanWidget.canvas.axes.set_title('Water Surface')
        torpedo_x = numpy.array([b, b * 1.05, b, -a * 0.9, -a, -a * 0.9, b])
        torpedo_y = numpy.array([bigger_radius, 0, -bigger_radius, -bigger_radius, 0, bigger_radius, bigger_radius])
        # (np.random.uniform(0, (a + b) * 20),
        # np.random.uniform(-(a + b) * 10, (a + b) * 10), marker='o')[0]

        our_torpedo = Torpedo(x0, y0, phi0, a, b, bigger_radius)
        our_torpedo.draw(self.OceanWidget.canvas.axes)
        drawn_torpedo = our_torpedo.drawn_torpedo
        our_aim = Aim(50, 20, 10)
        our_aim.draw(self.OceanWidget.canvas.axes)
        drawn_aim = our_aim.drawn_aim
        # aim = self.OceanWidget.canvas.axes.plot(50, 20, marker='o')[0]

        self.OceanWidget.canvas.show()

        initial_vector = numpy.array([x0, y0, phi0, v_x0, v_y0, omega_0])
        state_vector = initial_vector

        t0 = 0
        dt = 0.01
        t = t0

        def animate(i):
            global t, state_vector
            alpha = self.AngleBar.value() / 180 * 3.14 / 6
            print(state_vector)
            t = t + dt

            d_vector = calculate_movement_equation(state_vector, alpha)
            state_vector = state_vector + dt * d_vector
            print(state_vector)
            r_torpedo_x, r_torpedo_y = rot2d(torpedo_x, torpedo_y, state_vector[2])
            drawn_torpedo.set_data(state_vector[0] + r_torpedo_x, state_vector[1] + r_torpedo_y)
            # aim.set_data(np.random.uniform(0, (a + b) * 20), np.random.uniform(-(a + b) * 10, (a + b) * 10))
            return [drawn_torpedo, drawn_aim]

        fig = self.OceanWidget.canvas.figure
        print(type(fig))
        animation = FuncAnimation(fig, animate, interval=10, blit=True)
        self.OceanWidget.canvas.draw()
        # plt.show()
