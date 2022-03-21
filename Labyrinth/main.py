import matplotlib.pyplot as plt
import random as rd
import numpy as np
import math
import copy
import border
# import plotly
# import plotly.graph_obj
from matplotlib.animation import FuncAnimation
from map import *
from robot import *


        # Не лучший способ задания преград
if __name__ == "__main__":
    x_lim=6
    y_lim=6
    x = Map(heads=np.array([]),
            tails=np.array([]),
            xlim=x_lim,ylim=y_lim,x0=0,y0=0)

    for i in range(0, 20):
        x.add_bounds(*border.generate_rectangle(x_lim, y_lim))

    x.draw_map()
    r = Robot(map_ = x,x=0.0, y=1.2, k=2, k1=15, k2=2)
    r.draw()
    r.map_grid()

    r.get_attainability([4.5, 3.5])
    plt.show()