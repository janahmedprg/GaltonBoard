import numpy as np
import math
import matplotlib.pyplot as plt
import random

nrows = 3
ncols = 5
Z = np.zeros((3,5)).reshape(nrows, ncols)
Z[1][1]=5
Z[2][1]=5
Z[2][2]=5
Z[2][3]=5
x = np.arange(ncols + 1)
y = np.arange(nrows + 1)

fig, ax = plt.subplots()
ax.pcolormesh(x, y, Z, shading='flat', vmin=Z.min(), vmax=Z.max())


def _annotate(ax, x, y, title):
    # this all gets repeated below:
    X, Y = np.meshgrid(x, y)
    ax.plot(X.flat, Y.flat, 'o', color='m')
    ax.set_xlim(-0.7, 5.2)
    ax.set_ylim(-0.7, 3.2)
    ax.set_title(title)

_annotate(ax, x, y, "shading='flat'")

plt.show()
