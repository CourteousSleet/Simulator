from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import FormOfTheTorpeda
import oceanwidget
import math
import numpy as np
import random as rd
from matplotlib.animation import FuncAnimation
import sympy as sp

class Aim():
    x_now = 0
    y_now = 0
    S = 10
    #DrawedAim
    IsDestroyed = 0

    def __init__(self, str_x, str_y, S):
        st = sp.Symbol('t')
        self.x = sp.lambdify(st, str_x)
        self.y = sp.lambdify(st, str_y)
        self.x_now = self.x(0)
        self.y_now = self.y(0)
        self.S = S
        phi = np.linspace(0, 6.28, 20)
        self.AimX = self.S/2*np.sin(phi)
        self.AimY = self.S/2*np.cos(phi)

    def Draw(self, axes):
        self.DrawedAim = axes.plot(self.x_now+self.AimX, self.y_now+self.AimY)[0]

    def Replace(self, t):
        self.x_now = self.x(t)
        self.y_now = self.y(t)

    def Redraw(self):
        self.DrawedAim.set_data(self.x_now + self.AimX, self.y_now + self.AimY)


class Torpeda():
    x = 0
    y = 0
    phi = 0
    a = 1
    b = 1
    R = 1
    K=1
    Alpha_max=0.8
    Sensor_r = -1
    Sensor_phi = -1
    IsDetecting = 0
    IsExploded = 0
    Sensor_r_max = 50
    Sensor_phi_max = 0.5

    def __init__(self, x, y, phi, a, b, R):
        self.x = x
        self.y = y
        self.phi = phi
        self.a = a
        self.b = b
        self.R = R
        self.IsDetecting = 0
        self.IsExploded = 0
        self.TorpedaX = np.array([b, b * 1.05, b, -a * 0.9, -a, -a * 0.9, b])
        self.TorpedaY = np.array([R, 0, -R, -R, 0, R, R])

    def Draw(self, axes):
        RTorpedaX, RTorpedaY = Rot2D(self.TorpedaX, self.TorpedaY, self.phi)
        self.DrawedTorpeda = axes.plot(self.x + RTorpedaX, self.y + RTorpedaY)[0]

    def Redraw(self):
        if self.IsExploded:
            self.TorpedaX *= 1 + (3*self.a - np.max(self.TorpedaX))/(a*100)
            self.TorpedaY *= 1 + (3*self.a - np.max(self.TorpedaX))/(a*100)
            self.DrawedTorpeda.set_data(self.x + self.TorpedaX, self.y + self.TorpedaY)
        else:
            RTorpedaX, RTorpedaY = Rot2D(self.TorpedaX, self.TorpedaY, self.phi)
            self.DrawedTorpeda.set_data(self.x + RTorpedaX, self.y + RTorpedaY)

    def CourseChose(self, Aim):
        R_vectorX = Aim.x_now - self.x
        R_vectorY = Aim.y_now - self.y
        beta = math.atan2(R_vectorY, R_vectorX) - self.phi
        if self.K*beta > self.Alpha_max:
            Alpha = self.Alpha_max
        elif self.K*beta < -self.Alpha_max:
            Alpha = - self.Alpha_max
        else:
            Alpha = self.K*beta
        return Alpha

    def ExplodingCondition(self, Aim):
        x_nose = self.x + b*np.cos(self.phi)
        y_nose = self.y + b*np.sin(self.phi)
        if (x_nose - Aim.x_now)**2 + (y_nose - Aim.y_now)**2 < (Aim.S/2)**2:
            return 1
        else:
            return 0

    def Explode(self, axes):
        self.IsExploded = 1
        n = 6
        Phi1 = np.linspace(0, 6.28, n+1)
        VzrivX1 = self.a/4*np.sin(Phi1)
        VzrivY1 = self.a/4*np.cos(Phi1)
        Phi2 = np.linspace(0, 6.28, 2*n + 1)
        VzrivX2 = self.a/2 * np.sin(Phi2)
        VzrivY2 = self.a/2 * np.cos(Phi2)
        self.TorpedaX = np.hstack([VzrivX1, VzrivX2])
        self.TorpedaY = np.hstack([VzrivY1, VzrivY2])
        self.DrawedTorpeda = axes.plot(self.x + self.TorpedaX, self.y + self.TorpedaY, marker='o', color=[1, 0.3, 0])[0]



def Rot2D(X,Y,phi):
    RotX = X*np.cos(phi) - Y*np.sin(phi)
    RotY = X*np.sin(phi) + Y*np.cos(phi)
    return RotX, RotY


def MoveEquation(StateVector,t,Alpha):
    global m, r, a, b, R, k, Cl, Cb, Cvr, rho, Vvpd, waterVx, waterVy
    x, y, phi, Vx, Vy, Omega = StateVector

    dx = Vx
    dy = Vy
    dphi = Omega

    sumVx = Vx - waterVx(x, y)
    sumVy = Vy - waterVy(x, y)
    Vl =  sumVx * np.cos(phi) + sumVy * np.sin(phi)
    Vb = -sumVx * np.sin(phi) + sumVy * np.cos(phi)

    Sl = 3.14 * R ** 2
    Sb = (a + b) * 2 * R
    Sr = R ** 2 * k

    J = m * (a ** 2 + b ** 2) / 6 * r

    Fsl = rho * Sl * Cl * Vl ** 2 * np.sign(Vl) / 2
    Fsb = rho * Sb * Cb * Vb ** 2 * np.sign(Vb) / 2
    Ms = rho * Sb * Cvr * Omega ** 2 * (b + a) ** 2 * np.sign(Omega) / 8

    Fdv = rho * 3.14 * R ** 2 * Vvpd ** 2
    Fr = rho * Sr * Vvpd ** 2 * (np.sin(Alpha)) ** 2 * np.sign(np.sin(Alpha)) / 2

    dVx = (Fdv * np.cos(phi) + Fsb * np.sin(phi) - Fsl * np.cos(phi) + Fr * np.sin(phi - Alpha)) / m
    dVy = (Fdv * np.sin(phi) - Fsb * np.cos(phi) - Fsl * np.sin(phi) - Fr * np.cos(phi - Alpha)) / m
    dOmega = (a * Fr * np.cos(Alpha) - Ms) / J

    return np.array([dx, dy, dphi, dVx, dVy, dOmega])
global t

class OceanWidget(QMainWindow, FormOfTheTorpeda.Ui_MainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.setWindowTitle("Творение")
        self.FireButton.clicked.connect(self.HereAreWeGo)
        self.addToolBar(NavigationToolbar(self.OceanWidget.canvas, self))

    def HereAreWeGo(self):
        global m, r, a, b, R, k, Cl, Cb, Cvr, rho, Vvpd, waterVx, waterVy, t

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
        x = sp.Symbol('x')
        y = sp.Symbol('y')
        waterVx = sp.lambdify([x, y], self.waterVx_Edit.text())
        waterVy = sp.lambdify([x, y], self.waterVy_Edit.text())

        xvalues = np.arange(0, 110, 10)
        yvalues = np.arange(-50, 55, 10)

        X, Y = np.meshgrid(xvalues, yvalues)
        u = waterVx(X, Y)
        v = waterVy(X, Y)

        x0 = float(self.x0_Edit.text())
        y0 = float(self.y0_Edit.text())
        phi0 = float(self.phi0_Edit.text())
        Vx0 = float(self.Vx0_Edit.text())
        Vy0 = float(self.Vy0_Edit.text())
        Omega0 = float(self.Omega0_Edit.text())

        self.OceanWidget.canvas.axes.clear()
        self.OceanWidget.canvas.axes.grid(True)
        self.OceanWidget.canvas.axes.axis('scaled')
        self.OceanWidget.canvas.axes.set(xlim=[0, (a+b)*20], ylim=[-(a+b)*10, (a+b)*10])
        self.OceanWidget.canvas.axes.set_title('Водная гладь')
        self.OceanWidget.canvas.axes.quiver(X, Y, u, v, scale=6, units='xy')


        Our_Torpeda = Torpeda(x0, y0, phi0, a, b, R)
        Our_Torpeda.Draw(self.OceanWidget.canvas.axes)
        Drawed_Torpeda = Our_Torpeda.DrawedTorpeda
        Our_Aim = Aim('50+20*cos(t)', '20+30*sin(t)', 10)
        Our_Aim.Draw(self.OceanWidget.canvas.axes)
        Drawed_Aim = Our_Aim.DrawedAim

        self.OceanWidget.canvas.show()

        global t, StateVector

        InitialVector = np.array([x0, y0, phi0, Vx0, Vy0, Omega0])
        StateVector = InitialVector

        t0 = 0
        dt = 0.01

        t=t0

        def anima(i):
            global t, StateVector
            t = t+dt
            if not Our_Torpeda.IsExploded:
                Alpha = Our_Torpeda.CourseChose(Our_Aim)

                dVector = MoveEquation(StateVector, t, Alpha)
                StateVector = StateVector + dt*dVector

                Our_Torpeda.x = StateVector[0]
                Our_Torpeda.y = StateVector[1]
                Our_Torpeda.phi = StateVector[2]

            Our_Aim.Replace(t)

            Our_Torpeda.Redraw()
            Our_Aim.Redraw()

            if Our_Torpeda.ExplodingCondition(Our_Aim) and not Our_Torpeda.IsExploded:
                Our_Torpeda.Explode(self.OceanWidget.canvas.axes)

            Drawed_Torpeda = Our_Torpeda.DrawedTorpeda
            return [Drawed_Torpeda,  Drawed_Aim]

        fig = self.OceanWidget.canvas.figure
        anim = FuncAnimation(fig, anima, interval=10, blit=True)
        self.OceanWidget.canvas.draw()



app = QApplication([])
window = OceanWidget()
window.show()
app.exec_()


