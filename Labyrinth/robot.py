import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.animation import FuncAnimation
from scipy.interpolate import splprep, splev

from point import Point


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
    def __init__(self, x=0, y=0, vx=0, vy=0, phi=-1.57, omega=0, r=0.1, m =1, k=1, k1=1, k2=1, r_wheel=0.1, V=0.5, map_=None):
        self.m = m
        self.k = k
        self.k1 = k1
        self.k2 = k2
        self.r_wheel = r_wheel
        self.x, self.y, self.vx, self.vy = x + r, y + r, vx, vy
        self.phi, self.omega, self.r = phi, omega, r
        self.map = map_
        self.grid = None   # Матрица квадратов карты
        self.pos = [x, y]
        self.drawed_robot = None
        self.drawed_front = None
        self.V = V

    def draw(self):
        self.drawed_robot = self.map.ax.plot(self.x, self.y, marker='o', color='green', markersize=10)[0]
        self.drawed_front = self.map.ax.plot(self.x + 0.07 * np.cos(self.phi),
                                             self.y + 0.07 * np.sin(self.phi),
                                             marker='o', color='red', markersize=5)[0]
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
        self.pos = [int(self.x*len(grid_matrix[1:])/(self.map.xlim-self.map.x0)),
                    int(self.y*len(grid_matrix)/(self.map.ylim-self.map.y0))]
        print(self.pos)
        return self.map.plot, grid_matrix

    def get_coordinates(self, vector):
        x_scale = self.map.xlim / len(self.grid)
        y_scale = self.map.ylim / len(self.grid[0])

        return [(vector[0] + 0.5) * x_scale, (vector[1] + 0.5) * y_scale]

    def move_equations(self, M1, M2):
        dx = self.vx
        dy = self.vy
        dphi = self.omega
        V_forward = dx * np.cos(self.phi) + dy * np.sin(self.phi)

        J = self.m * (self.r ** 2) / 2
        F_resist = self.k * V_forward
        M_resist = self.k * dphi

        F1 = M1 / self.r_wheel
        F2 = M2 / self.r_wheel
        dV_forward = ((F1 + F2 - F_resist)) / self.m
        dVx = dV_forward * np.cos(self.phi) - V_forward * self.omega * np.sin(self.phi) / self.m
        dVy = dV_forward * np.sin(self.phi) + V_forward * self.omega * np.cos(self.phi) / self.m
        dOmega = (F1 * self.r - F2 * self.r - M_resist) / J

        return np.array([dx, dy, dphi, dVx, dVy, dOmega])

    def choose_course(self, point):
        r_x = point.x - self.x
        r_y = point.y - self.y

        beta = math.atan2(r_y, r_x) - self.phi

        while abs(beta) > np.pi:
            beta += 2 * np.pi * (-np.sign(beta))

        r_module = np.sqrt(r_x**2 + r_y**2)
        subtraction_M = self.k1 * beta

        if abs(beta) < np.pi / 2:
            sum_M = self.k2 * (r_module) * np.cos(beta)
        else:
            sum_M = 0

        M1 = (sum_M + subtraction_M) / 2
        M2 = (sum_M - subtraction_M) / 2

        return M1, M2

    def get_attainability(self,aim):
        aim_pos = [int(aim[0] * len(self.grid[1]) / (self.map.xlim - self.map.x0)),
                  int(aim[1] * len(self.grid) / (self.map.ylim - self.map.y0))]
        self.grid[self.pos[0]][self.pos[1]].att = 0
        n=0
        M = np.zeros([len(self.grid[0]),len(self.grid)])
        M[self.pos[0]][self.pos[1]] = 0
        while self.grid[aim_pos[0]][aim_pos[1]].att == None:
            for i in range(len(self.grid[0])):
                for j in range(len(self.grid[0])):
                    if self.grid[i][j].att == n:
                        if i+1 <= len(self.grid[0])-1:
                            if self.grid[i+1][j].att == None:    # 1) Сосед должен быть; 2)  Не помечен
                                self.grid[i + 1][j].att=n+1
                                M[i+1,j]=n+1
                        if i-1 >= 0:
                            if self.grid[i-1][j].att == None:    # 1) Сосед должен быть; 2)  Не помечен
                                self.grid[i - 1][j].att=n+1
                                M[i - 1, j] = n + 1
                        if j+1 <= len(self.grid)-1:
                            if self.grid[i][j+1].att == None:    # 1) Сосед должен быть; 2)  Не помечен
                                self.grid[i][j+1].att=n+1
                                M[i, j+1] = n + 1
                        if j-1 >= 0:
                            if self.grid[i][j-1].att == None:    # 1) Сосед должен быть; 2)  Не помечен
                                self.grid[i][j-1].att=n+1
                                M[i, j-1] = n + 1
            n+=1

        for i in range(len(self.grid[0])):
            for j in range(len(self.grid[0])):
                if self.grid[i][j].att != None:
                    if self.grid[i][j].att >= 0:
                        self.grid[i][j].fill_att([1,1-self.grid[i][j].att/n,1-self.grid[i][j].att/n])
                        self.map.ax.add_patch(self.grid[i][j].rect)
        self.grid[aim_pos[0]][aim_pos[1]].fill_att([1,1,0])
        self.map.ax.add_patch(self.grid[aim_pos[0]][aim_pos[1]].rect)

        Trace = [aim_pos]
        for k in range(n):
            is_done = 0
            if Trace[-1][0] + 1 <= len(self.grid[0]) - 1:
                if self.grid[Trace[-1][0] + 1][Trace[-1][1]].att == n - 1 - k:
                    Trace = Trace + [[Trace[-1][0] + 1, Trace[-1][1]]]

                    self.grid[Trace[-1][0]][Trace[-1][1]].fill_att([1*(n - 1 - k)/n, 1, 0])
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

        corner_indices = np.zeros(len(Trace))

        for i in range(1, len(Trace)-1):
            xP, yP = Trace[i - 1][0] == Trace[i][0], Trace[i - 1][1] == Trace[i][1]
            xN, yN = Trace[i + 1][0] == Trace[i][0], Trace[i + 1][1] == Trace[i][1]

            is_corner = xP != xN and yP != yN
            if is_corner:
                corner_indices[i] = 1

        filtered = [Trace[i] for (i, n) in enumerate(corner_indices) if n == 0]
        vectors = [self.get_coordinates(elem) for elem in filtered]

        x, y = zip(*vectors)
        dt = 0.001
        T_max = len(Trace[:][0]) / self.V
        tck, u = splprep([x, y], u=None, s=0.0, per=0)
        u_new = np.linspace(u.min(), u.max(), round(T_max/dt))
        x_new, y_new = splev(u_new, tck, der=0)

        print(len(x_new))
        self.map.ax.plot(x_new, y_new, 'w-')

        point = Point(x_new, y_new)
        point.draw(self.map.ax)
        dt = 0.001
        global Frame
        Frame=100
        def animation(i):
            global Frame
            if Frame<len(x_new):
                M1, M2 = self.choose_course(point)
                dVector = self.move_equations(M1, M2)
                state_vector = [self.x, self.y, self.phi, self.vx, self.vy, self.omega]

                self.x, self.y, self.phi, \
                self.vx, self.vy, self.omega = state_vector + dt * dVector

                self.drawed_robot.set_data(self.x, self.y)
                self.drawed_front.set_data(self.x + 0.07 * np.cos(self.phi),
                                           self.y + 0.07 * np.sin(self.phi))

                point.replace(Frame)
                point.redraw()
            Frame+=1
            return [self.drawed_robot, self.drawed_front, point.drawed_point]
        a = FuncAnimation(self.map.plot, animation, interval=2, frames=range(len(x_new)), blit=True)
        plt.show()
