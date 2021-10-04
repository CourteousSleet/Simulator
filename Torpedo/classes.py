import math

from PyQt5.QtWidgets import *
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import numpy
import matplotlib.pyplot as plt
import oceanwidget
import OnlyPoint


class WidgetOfLife(QMainWindow, OnlyPoint.Ui_MainWindow):

    def __init__(self):
        QMainWindow.__init__(self)

        self.setupUi(self)

        self.setWindowTitle("Пошла!")

        self.ButtonOfTheStart.clicked.connect(self.start)

        # self.horizontalScrollBar.valueChanged.connect(self.valuechange)
        self.addToolBar(NavigationToolbar(self.WidgetOfLife.canvas, self))

    def start(self):
        print(1)
        r = self.RadiusBar.value() / 100
        print(1.5)
        self.WidgetOfLife.canvas.axes.clear()
        print(1.7)
        t0 = 0
        x0 = r * numpy.cos(t0)
        y0 = r * numpy.sin(t0)
        print(2)
        point = self.WidgetOfLife.canvas.axes.plot(x0, y0, marker='o')[0]
        print(3)
        self.WidgetOfLife.canvas.axes.grid(True)
        self.WidgetOfLife.canvas.axes.axis('scaled')
        self.WidgetOfLife.canvas.axes.set(xlim=[-1.1, 1.1], ylim=[-1.1, 1.1])
        self.WidgetOfLife.canvas.axes.set_title('Окружность')
        print(4)
        self.WidgetOfLife.canvas.draw()
        print(5)
        global t
        t = t0
        print(5.1)
        print(self.Dt_Edit.text().strip())
        print(type(self.Dt_Edit.text()))
        dt = float(self.Dt_Edit.text())
        print(6)

        def animate(i):
            print(7)
            global t, r
            t += dt
            r = self.RadiusBar.value() / 100

            x = r * numpy.cos(t)
            y = r * numpy.sin(t)

            point.set_data(x, y)

            return [point]

        anim = FuncAnimation(self.WidgetOfLife.canvas.figure, animate,
                             interval=10, blit=True)
        self.WidgetOfLife.canvas.draw()


class Aim:
    x = 0
    y = 0
    size = 10
    is_destroyed = 0

    def __init__(self, x, y, s_):
        self.x = x
        self.y = y
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
    is_detected = 0
    is_exploded = 0
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
        self.torpedo_y = numpy.array([R, 0, -R, -R, 0, R, R])

    def draw(self, axes):
        r_torpedo_x, r_torpedo_y = rot2d(self.torpedo_x, self.torpedo_y, self.phi)
        self.drawn_torpedo = axes.plot(self.x + r_torpedo_x, self.y + r_torpedo_y)[0]

    def redraw(self, axes):
        r_torpedo_x, r_torpedo_y = rot2d(self.torpedo_x, self.torpedo_y, self.phi)
        self.drawn_torpedo = axes.plot(self.x + r_torpedo_x, self.y + r_torpedo_y)

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

    def check_for_strike(self, aim_):
        if ()

def rot2d(x, y, phi):
    rot_x = x * numpy.cos(phi) - y * numpy.sin(phi)
    rot_y = x * numpy.sin(phi) + y * numpy.cos(phi)
    return rot_x, rot_y


def calculate_movement_equation(state_vector, alpha):
    global m, r, a, b, R, k, Cl, Cb, Cvr, rho, Vvpd
    x, y, phi, v_x, v_y, omega = state_vector
    dx = v_x
    dy = v_y
    d_phi = omega

    v_l = v_x * numpy.cos(phi) + v_y * numpy.sin(phi)
    v_b = -v_x * numpy.sin(phi) + v_y * numpy.cos(phi)

    s_l = 3.14 * R ** 2
    s_b = (a + b) * 2 * R
    s_r = R ** 2 * k

    j = m * (a ** 2 + b ** 2) / 6 * r

    f_sl = rho * s_l * Cl * v_l ** 2 * numpy.sign(v_l) / 2
    f_sb = rho * s_b * Cb * v_b ** 2 * numpy.sign(v_b) / 2
    m_s = rho * s_b * Cvr * omega ** 2 * (b + a) ** 2 * numpy.sign(omega) / 8

    f_dv = rho * 3.14 * R ** 2 * Vvpd ** 2
    f_r = rho * s_r * Vvpd ** 2 * (numpy.sin(alpha)) ** 2 * numpy.sign(numpy.sin(alpha)) / 2

    d_v_x = (f_dv * numpy.cos(phi) + f_sb * numpy.sin(phi) - f_sl * numpy.cos(phi) + f_r * numpy.sin(phi - alpha)) / m
    d_v_y = (f_dv * numpy.sin(phi) - f_sb * numpy.cos(phi) - f_sl * numpy.sin(phi) - f_r * numpy.cos(phi - alpha)) / m
    d_omega = (a * f_r * numpy.cos(alpha) - m_s) / j

    return numpy.array([dx, dy, d_phi, d_v_x, d_v_y, d_omega])


global t


class OceanWidget(QMainWindow, FormOfTheTorpeda.Ui_MainWindow):

    def __init__(self):
        QMainWindow.__init__(self)

        self.setupUi(self)

        self.setWindowTitle("Творение")

        self.FireButton.clicked.connect(self.HereAreWeGo)

        # self.horizontalScrollBar.valueChanged.connect(self.valuechange)
        self.addToolBar(NavigationToolbar(self.OceanWidget.canvas, self))

    def HereAreWeGo(self):
        global m, r, a, b, R, k, Cl, Cb, Cvr, rho, Vvpd, t
        # Alpha = self.AngleBar.value()/5000
        print(0)
        #     Параметры массы
        m = float(self.m_Edit.text())
        r = float(self.r_Edit.text())

        #     Геометрические параметры
        a = float(self.a_Edit.text())
        b = float(self.b_Edit.text())
        R = float(self.R_Edit.text())
        k = float(self.k_Edit.text())

        Cl = float(self.Cl_Edit.text())
        Cb = float(self.Cb_Edit.text())
        Cvr = float(self.Cvr_Edit.text())

        rho = float(self.rho_Edit.text())

        Vvpd = float(self.Vvpd_Edit.text())

        x0 = float(self.x0_Edit.text())
        y0 = float(self.y0_Edit.text())
        phi0 = float(self.phi0_Edit.text())
        Vx0 = float(self.Vx0_Edit.text())
        Vy0 = float(self.Vy0_Edit.text())
        Omega0 = float(self.Omega0_Edit.text())

        print(1)

        self.OceanWidget.canvas.axes.clear()
        self.OceanWidget.canvas.axes.grid(True)
        self.OceanWidget.canvas.axes.axis('scaled')
        self.OceanWidget.canvas.axes.set(xlim=[0, (a + b) * 20], ylim=[-(a + b) * 10, (a + b) * 10])
        self.OceanWidget.canvas.axes.set_title('Водная гладь')
        # TorpedaX = numpy.array([b, b*1.05, b, -a*0.9, -a, -a*0.9, b])
        # TorpedaY = numpy.array([R, 0, -R, -R, 0, R, R])
        print(1.1)
        # (numpy.random.uniform(0, (a + b) * 20),
        # numpy.random.uniform(-(a + b) * 10, (a + b) * 10), marker='o')[0]

        print(1.3)

        Our_Torpeda = Torpeda(x0, y0, phi0, a, b, R)
        Our_Torpeda.Draw(self.OceanWidget.canvas.axes)
        Drawed_Torpeda = Our_Torpeda.drawn_torpedo
        Our_Aim = Aim(50, 20, 10)
        Our_Aim.Draw(self.OceanWidget.canvas.axes)
        Drawed_Aim = Our_Aim.DrawedAim
        # aim = self.OceanWidget.canvas.axes.plot(50, 20, marker='o')[0]
        print(2)

        self.OceanWidget.canvas.show()
        print(3)

        global t, StateVector

        InitialVector = numpy.array([x0, y0, phi0, Vx0, Vy0, Omega0])
        StateVector = InitialVector

        t0 = 0
        dt = 0.01

        t = t0

        print(4)

        def anima(i):
            global t, StateVector
            Alpha = self.AngleBar.value() / 180 * 3.14 / 6
            print(StateVector)
            t = t + dt

            dVector = calculate_movement_equation(StateVector, t, Alpha)
            StateVector = StateVector + dt * dVector
            print(StateVector)
            RTorpedaX, RTorpedaY = rot2d(TorpedaX, TorpedaY, StateVector[2])
            Drawed_Torpeda.set_data(StateVector[0] + RTorpedaX, StateVector[1] + RTorpedaY)
            # aim.set_data(numpy.random.uniform(0, (a + b) * 20), numpy.random.uniform(-(a + b) * 10, (a + b) * 10))
            return [Drawed_Torpeda, Drawed_Aim]

        fig = self.OceanWidget.canvas.figure
        print(type(fig))
        anim = FuncAnimation(fig, anima, interval=10, blit=True)
        self.OceanWidget.canvas.draw()
        # plt.show()
