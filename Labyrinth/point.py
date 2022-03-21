import numpy as np
import sympy as sp

class Point():
    x_now = 0
    y_now = 0

    def __init__(self, x_path, y_path):
        self.x_path = x_path
        self.y_path = y_path
        self.x = self.x_path[0]
        self.y = self.y_path[0]

    def draw(self, axes):
        self.drawed_point = axes.plot(self.x, self.y, marker='o', color='blue', markersize=3)[0]

    def replace(self, s):
        self.x = self.x_path[s]
        self.y = self.y_path[s]

    def redraw(self):
        self.drawed_point.set_data(self.x, self.y)