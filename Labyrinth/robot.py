import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.animation import FuncAnimation


class RectangleGrid():
    def __init__(self, x, y, s, available=True):
        self.x, self.y, self.s, self.available = x, y, s, available
        self.rect = plt.Rectangle((x, y), s, s, linewidth=0.5, edgecolor='blue', fill=not available)
        self.att = None

    def check(self, heads, tails):
        for h, t in zip(heads, tails):

            if self.x <= h[0] and round(self.x + self.s, 3) > h[0]:
                if self.y <= h[1] and round(self.y + self.s, 3) > h[1]:
                    self.available = False
                    self.att = -1
            if self.x <= t[0] and round(self.x + self.s, 3) > t[0]:
                if self.y <= t[1] and round(self.y + self.s, 3) > t[1]:
                    self.available = False
                    self.att = -1

    def fill(self):
        self.rect = plt.Rectangle((self.x, self.y), self.s, self.s, linewidth=0.5, edgecolor='blue',
                                  fill=not self.available)

    def fill_att(self, color):
        self.rect = plt.Rectangle((self.x, self.y), self.s, self.s, linewidth=0.5, facecolor=color,
                                  edgecolor='blue', fill=True)


class Robot():
    def __init__(self, x=0, y=0, v=1, phi=0.5, omega=1, r=0.1, map_=None):
        self.x, self.y, self.v = x + r, y + r, v
        self.phi, self.omega, self.r = phi, omega, r
        self.map = map_
        self.grid = None  # Матрица квадратов карты
        self.pos = [x, y]
        self.drawed_robot = None

    def draw(self):
        self.drawed_robot = self.map.ax.plot(self.x, self.y, marker='o', color='green', markersize=10)[0]

    def pathfinder(self):
        return 0

    def map_grid(self):
        xs = np.arange(self.map.x0, self.map.xlim, 2 * self.r)
        ys = np.arange(self.map.y0, self.map.ylim, 2 * self.r)
        grid_matrix = [[None for _ in xs] for _ in ys]
        rects = []
        self.map.divide_bounds(self.r)
        for i in range(len(xs)):
            for j in range(len(ys)):
                r = RectangleGrid(round(xs[i], 3), round(ys[j], 3), round(2 * self.r, 3))
                rects.append(r)
                r.check(self.map.heads, self.map.tails)
                r.fill()
                self.map.ax.add_patch(r.rect)
                grid_matrix[i][j] = r
        self.grid = grid_matrix
        self.pos = [int(self.x * len(grid_matrix[1:]) / (self.map.xlim - self.map.x0)),
                    int(self.y * len(grid_matrix) / (self.map.ylim - self.map.y0))]
        print(self.pos)
        return self.map.plot, grid_matrix

    def move_up(self):
        if self.pos[1] < len(self.grid):
            self.y = round(self.y + 2 * self.r, 3)
            self.pos = [self.pos[0], self.pos[1] + 1]
            self.drawed_robot.set_data(self.x, self.y)

    def move_down(self):
        if self.pos[1] > 0:
            self.y = round(self.y - 2 * self.r, 3)
            self.pos = [self.pos[0], self.pos[1] - 1]
            self.drawed_robot.set_data(self.x, self.y)

    def move_left(self):
        if self.pos[0] > 0:
            self.x = self.x - 2 * self.r
            self.pos = [self.pos[0] - 1, self.pos[1]]
            self.drawed_robot.set_data(self.x, self.y)

    def move_right(self):
        if self.pos[0] < len(self.grid[0]):
            self.x = self.x + 2 * self.r
            self.pos = [self.pos[0] + 1, self.pos[1]]
            self.drawed_robot.set_data(self.x, self.y)

    def get_attainability(self, aim):
        # aim = [x, y]
        aim_pos = [int(aim[0] * len(self.grid[1]) / (self.map.xlim - self.map.x0)),
                   int(aim[1] * len(self.grid) / (self.map.ylim - self.map.y0))]
        self.grid[self.pos[0]][self.pos[1]].att = 0
        n = 0
        M = np.zeros([len(self.grid[0]), len(self.grid)])
        M[self.pos[0]][self.pos[1]] = 0
        while self.grid[aim_pos[0]][aim_pos[1]].att == None:
            for i in range(len(self.grid[0])):
                for j in range(len(self.grid[0])):
                    if self.grid[i][j].att == n:
                        if i + 1 <= len(self.grid[0]) - 1:
                            if self.grid[i + 1][j].att == None:  # 1) Сосед должен быть; 2)  Не помечен
                                self.grid[i + 1][j].att = n + 1
                                M[i + 1, j] = n + 1
                        if i - 1 >= 0:
                            if self.grid[i - 1][j].att == None:  # 1) Сосед должен быть; 2)  Не помечен
                                self.grid[i - 1][j].att = n + 1
                                M[i - 1, j] = n + 1
                        if j + 1 <= len(self.grid) - 1:
                            if self.grid[i][j + 1].att == None:  # 1) Сосед должен быть; 2)  Не помечен
                                self.grid[i][j + 1].att = n + 1
                                M[i, j + 1] = n + 1
                        if j - 1 >= 0:
                            if self.grid[i][j - 1].att == None:  # 1) Сосед должен быть; 2)  Не помечен
                                self.grid[i][j - 1].att = n + 1
                                M[i, j - 1] = n + 1
            n += 1

        for i in range(len(self.grid[0])):
            for j in range(len(self.grid[0])):
                if self.grid[i][j].att != None:
                    if self.grid[i][j].att >= 0:
                        self.grid[i][j].fill_att([1, 1 - self.grid[i][j].att / n, 1 - self.grid[i][j].att / n])
                        self.map.ax.add_patch(self.grid[i][j].rect)
        self.grid[aim_pos[0]][aim_pos[1]].fill_att([1, 1, 0])
        self.map.ax.add_patch(self.grid[aim_pos[0]][aim_pos[1]].rect)

        Trace = [aim_pos]
        for k in range(n):
            is_done = 0
            if Trace[-1][0] + 1 <= len(self.grid[0]) - 1:
                if self.grid[Trace[-1][0] + 1][Trace[-1][1]].att == n - 1 - k:
                    Trace = Trace + [[Trace[-1][0] + 1, Trace[-1][1]]]

                    self.grid[Trace[-1][0]][Trace[-1][1]].fill_att([1 * (n - 1 - k) / n, 1, 0])
                    self.map.ax.add_patch(self.grid[Trace[-1][0]][Trace[-1][1]].rect)
                    is_done = 1
            if Trace[-1][0] - 1 >= 0 and is_done == 0:
                if self.grid[Trace[-1][0] - 1][Trace[-1][1]].att == n - 1 - k:
                    Trace = Trace + [[Trace[-1][0] - 1, Trace[-1][1]]]
                    self.grid[Trace[-1][0]][Trace[-1][1]].fill_att([1 * (n - 1 - k) / n, 1, 0])
                    self.map.ax.add_patch(self.grid[Trace[-1][0]][Trace[-1][1]].rect)
                    is_done = 1
            if Trace[-1][1] + 1 <= len(self.grid) - 1 and is_done == 0:
                if self.grid[Trace[-1][0]][Trace[-1][1] + 1].att == n - 1 - k:
                    Trace = Trace + [[Trace[-1][0], Trace[-1][1] + 1]]
                    self.grid[Trace[-1][0]][Trace[-1][1]].fill_att([1 * (n - 1 - k) / n, 1, 0])
                    self.map.ax.add_patch(self.grid[Trace[-1][0]][Trace[-1][1]].rect)
                    is_done = 1
            if Trace[-1][1] - 1 >= 0 and is_done == 0:
                if self.grid[Trace[-1][0]][Trace[-1][1] - 1].att == n - 1 - k:
                    Trace = Trace + [[Trace[-1][0], Trace[-1][1] - 1]]
                    self.grid[Trace[-1][0]][Trace[-1][1]].fill_att([1 * (n - 1 - k) / n, 1, 0])
                    self.map.ax.add_patch(self.grid[Trace[-1][0]][Trace[-1][1]].rect)
                    is_done = 1
        Trace = Trace[::-1]
        print(Trace)
        print(len(Trace))
        global Kadr
        Kadr = 0

        def animation(i):
            global Kadr
            if Kadr < len(Trace):
                self.x = Trace[Kadr][0] * 2 * self.r + self.map.x0 + self.r
                self.y = Trace[Kadr][1] * 2 * self.r + self.map.x0 + self.r
                self.pos = [Trace[Kadr][0], Trace[Kadr][1]]
                self.drawed_robot.set_data(self.x, self.y)
            Kadr += 1
            return [self.drawed_robot]

        a = FuncAnimation(self.map.plot, animation, interval=100, frames=range(len(Trace)), blit=True)
        plt.show()
