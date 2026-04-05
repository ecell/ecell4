"""Regression test for https://github.com/ecell/ecell4/issues/100

fig.gca(projection='3d') was removed in matplotlib 3.6+.
"""
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ecell4.plotting._matplotlib import __prepare_mplot3d_with_matplotlib


def test_prepare_mplot3d_returns_3d_axes():
    wrange = {'x': (0, 1), 'y': (0, 1), 'z': (0, 1)}
    fig, ax = __prepare_mplot3d_with_matplotlib(
        wrange, figsize=6, grid=True, wireframe=False, angle=None, noaxis=False)
    assert isinstance(ax, Axes3D)
    plt.close(fig)
