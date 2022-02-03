import matplotlib.pyplot as plt
import random as rd
import numpy as np
import math
import copy
# import plotly
# import plotly.graph_obj
from matplotlib.animation import FuncAnimation
from map import *
from robot import *

# Не лучший способ задания преград
if __name__ == "__main__":
    x = Map(heads=np.array([[1.5, 2.7], [2.5, 2.2], [1.8, 2.2], [2.4, 4.3], [2, 1]]),
            tails=np.array([[2.4, 4.3], [2, 1], [1.5, 2.7], [2.5, 2.2], [1.8, 2.2]]),
            xlim=5, ylim=5, x0=0, y0=0)
    x.add_bounds(np.array([[3, 0.5], [3.5, 2], [4.5, 0.5], [3, 0.5]]),
                 np.array([[3.5, 2], [4.5, 0.5], [3, 0.5], [3.5, 2]]))
    x.draw_map()
    r = Robot(map_=x, x=0.0, y=1.2)
    r.draw()
    r.map_grid()
    # r.move_up()
    # r.move_left()
    # r.move_down()
    r.move_right()
    r.get_attainability([4.5, 1.5])
    plt.show()
