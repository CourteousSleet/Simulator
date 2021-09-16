import numpy as np


class TubeBorder():
    R = 1
    xC = 0
    yC = 0

    def __init__(self, R=1, xC=0, yC=0):
        self.R = R
        self.xC = xC
        self.yC = yC

    def BorderStrikeChecker(self, Xnew, Ynew):
        return (Xnew - self.xC) ** 2 + (Ynew - self.yC) ** 2 > self.R ** 2

    def DrawTubeBorder(self, ax):
        phi = np.linspace(0, 6.28, 100)
        XTube = self.xC + self.R * np.cos(phi)
        YTube = self.yC + self.R * np.sin(phi)
        ax.plot(XTube, YTube, 'black')


class Point2D():
    x0 = 0
    y0 = 0
    Vx0 = 0
    Vy0 = 0
    R = 0.1

    def __init__(self, x, y, Vx, Vy):
        self.x0 = x
        self.y0 = y
        self.Vx0 = Vx
        self.Vy0 = Vy

    def DrawPoint(self, ax):
        self.PlotPoint, = ax.plot(self.x0, self.y0, marker='o')

    def ReDrawPoint(self, x, y):
        self.PlotPoint.set_data(x, y)


def PointPointStrike(Xnew, Ynew):
    Strikes = np.zeros([len(Xnew), len(Xnew)], bool)
    for i in np.arange(len(Xnew) - 1):
        for j in np.arange(i + 1, len(Xnew)):
            Strikes[i, j] = (Xnew[i] - Xnew[j]) ** 2 + (Ynew[i] - Ynew[j]) ** 2 < (2 * 0.1) ** 2

    return Strikes
