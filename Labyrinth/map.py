import matplotlib.pyplot as plt
import numpy as np
import math


# Храню график сразу в карте. Дает преимущства:
#
# -   Разгружает main
# -   Упрощает работу в других классах, имеющих этот как поле
# -   Позволяет в main работать одновременно с несколькими картами (например, многокомнатная квартира)
#
# Недостатки:
#
# -  Любое изменение из другого класса изменит и карту
# - Нет выигрыша по памяти

class Map():
    def __init__(self, x0=0, xlim=4.5, y0=0, ylim=4.5, heads=None, tails=None):
        self.x0 = x0
        self.xlim = xlim
        self.y0 = y0
        self.ylim = ylim
        self.plot = plt.figure()
        self.ax = self.plot.add_subplot(111)
        self.ax.set_xlim(self.x0, self.xlim)
        self.ax.set_ylim(self.y0, self.ylim)
        self.boundaries = [Boundary(h, t) for h, t in zip(heads, tails)]
        self.heads = np.array(heads)
        self.tails = np.array(tails)

    def check_collision(self):
        return 0

    def draw_map(self):
        self.ax.cla()
        for bound in self.boundaries:
            self.ax.plot(bound.x, bound.y, color='black')
        self.ax.plot([self.x0, self.x0, self.xlim, self.xlim, self.x0],
                     [self.y0, self.ylim, self.ylim, self.y0, self.y0], color='black')
        return self.plot

    def add_bounds(self, h, t):
        result_h, result_t = np.zeros([len(self.heads) + len(h), 2]), np.zeros([len(self.tails) + len(t), 2])
        for i in range(len(self.heads)):
            result_h[i] = self.heads[i]
            result_t[i] = self.tails[i]
        for i in range(len(h)):
            self.boundaries.append(Boundary(h[i], t[i]))
            result_h[len(self.heads) + i] = h[i]
            result_t[len(self.tails) + i] = t[i]
        self.heads = result_h
        self.tails = result_t

    def divide_bounds(self, r):
        for h, t in zip(self.heads, self.tails):
            l = int(math.sqrt((h[0] - t[0]) ** 2 + (h[1] - t[1]) ** 2) * 2 / r)
            if l > 0:
                h_, t_ = np.zeros([l + 1, 2]), np.zeros([l + 1, 2])
                h_[0] = h
                for i in range(l):
                    dot = (h * (i + 1) + t * (l - i)) / (l + 1)
                    t_[i] = dot
                    h_[i + 1] = dot
                t_[-1] = t
            self.add_bounds(h_, t_)


class Boundary():
    def __init__(self, h, t):
        if h.any() != None and t.any() != None:
            self.h = h
            self.t = t
            self.x = [h[0], t[0]]
            self.y = [h[-1], t[-1]]
