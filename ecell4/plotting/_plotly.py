import random
import types

from ..util.session import load_world

from ._core import eval_key
from .styles import plotly_color_scale

__all__ = [
    'plot_number_observer_with_plotly',
    'plot_world_with_plotly',
    'plot_trajectory_with_plotly',
    ]

def init_notebook_mode(connected=False):
    import plotly.offline
    plotly.offline.init_notebook_mode()
    try:
        import google.colab
    except ImportError:
        pass
    else:
        import plotly.io as pio
        pio.renderers.default = "colab"

def plot_number_observer(
        *args, x=None, y=None, step=False, layout=None, **kwargs):
    """
    Generate a plot from NumberObservers and show it on IPython notebook
    with plotly.
    Require plotly and numpy.

    Parameters
    ----------
    obs : NumberObserver (e.g. FixedIntervalNumberObserver)
    step : bool, optional
        Piece-wise constant curve. False for default.
    layout : dict, optional
        The custom properties for layout.
        See also https://plot.ly/python/reference/#layout

    Examples
    --------
    >>> plot_number_observer(obs1)
    >>> plot_number_observer(obs1, obs2, obs3)
    >>> plot_number_observer(obs1, lambda t: 30)

    """
    import numpy
    import plotly
    import plotly.graph_objs as go

    x_key, y_keys = x, y

    color_scale = plotly_color_scale()

    if y_keys is not None and isinstance(y_keys, str):
        y_keys = (y_keys, )

    init_notebook_mode()
    fig = go.Figure()

    data, xdata = None, None
    for obs in args:
        if isinstance(obs, types.FunctionType):
            if xdata is None:
                raise ValueError("A function must be given after an observer.")
            y = [obs(xi) for xi in xdata]
            label = obs.__name__
            showlegend = (label not in color_scale.get_config())
            trace = go.Scatter(
                    x=xdata, y=y, name=label, line_color=color_scale.get_color(label),
                    legendgroup=label, showlegend=showlegend)
            fig.add_trace(trace)
            continue

        data = numpy.array(obs.data()).T

        if x_key is not None:
            xdata, _ = eval_key(x_key, obs.targets(), data)
        else:
            xdata = data[0]

        if y_keys is not None:
            targets_ = []
            for serial in y_keys:
                data_, err_ = eval_key(serial, obs.targets(), data)
                targets_.append((serial, data_, err_))
        else:
            targets_ = [(sp.serial(), data[idx + 1], None) for idx, sp in enumerate(obs.targets())]

        for label, data_, _ in targets_:
            showlegend = (label not in color_scale.get_config())
            trace = go.Scatter(
                    x=xdata, y=data_, name=label, line_shape=('linear' if not step else 'hv'),
                    line_color=color_scale.get_color(label), legendgroup=label,
                    showlegend=showlegend)
            fig.add_trace(trace)

    if x_key is not None:
        xaxis_title = 'The Number of Molecules [{}]'.format(x_key)
    else:
        xaxis_title = 'Time'
    layout_ = dict(xaxis_title=xaxis_title, yaxis_title="The Number of Molecules")
    if layout is not None:
        layout_.update(layout)
    fig.update(dict(layout=layout_))

    plotly.offline.iplot(fig)

plot_number_observer_with_plotly = plot_number_observer

def stl2mesh3d(filename, **kwargs):
    import plotly.graph_objs as go
    import numpy

    init_notebook_mode()

    from stl import mesh  # numpy-stl
    stl_mesh = mesh.Mesh.from_file(filename)

    p, q, r = stl_mesh.vectors.shape
    vertices, ixr = numpy.unique(stl_mesh.vectors.reshape(p * q, r), return_inverse=True, axis=0)
    I = numpy.take(ixr, [3 * k for k in range(p)])
    J = numpy.take(ixr, [3 * k + 1 for k in range(p)])
    K = numpy.take(ixr, [3 * k + 2 for k in range(p)])
    x, y, z = vertices.T
    # colorscale= [[0, '#e5dee5'], [1, '#e5dee5']]

    mesh3D = go.Mesh3d(
        x=x, y=y, z=z, i=I, j=J, k=K,
        flatshading=True, intensity=z, showscale=False, **kwargs)
    return mesh3D

def plot_stl(filename, layout=None, **kwargs):
    import plotly
    import plotly.graph_objs as go

    mesh3D = stl2mesh3d(filename, **kwargs)
    layout_ = dict(scene_aspectmode='data', margin=dict(l=0, r=0, b=0, t=0))
    if layout is not None:
        layout_.update(layout)
    layout_ = go.Layout(**layout_)
    fig = go.Figure(data=[mesh3D], layout=layout_)
    plotly.offline.iplot(fig)


def plot_world(world, species_list=None, max_count=1000, marker=None, layout=None, stl=None):
    """
    Generate a plot from received instance of World and show it on IPython notebook.

    Parameters
    ----------
    world : World or str
        World to render. A HDF5 filename is also acceptable.
    species_list : array of string, default None
        If set, plot_world will not search the list of species.
    max_count : Integer, default 1000
        The maximum number of particles to show for each species.
        None means no limitation.
    marker : dict, optional
        The custom properties for markers.
        See also https://plot.ly/python/marker-style/
    layout : dict, optional
        The custom properties for layout.
        See also https://plot.ly/python/reference/#layout

    """
    if isinstance(world, str):
        world = load_world(world)

    if species_list is None:
        species_list = [sp.serial() for sp in world.list_species()]
        species_list.sort()

    from ecell4_base.core import Species

    positions = {}
    for serial in species_list:
        x, y, z = [], [], []
        particles = world.list_particles_exact(Species(serial))
        if max_count is not None and len(particles) > max_count:
            particles = random.sample(particles, max_count)
        for pid, p in particles:
            pos = p.position()
            x.append(pos[0])
            y.append(pos[1])
            z.append(pos[2])

        positions[serial] = (x, y, z)

    import plotly
    import plotly.graph_objs as go

    init_notebook_mode()

    marker_ = dict(size=6, line=dict(color='rgb(204, 204, 204)', width=1),
            opacity=0.9, symbol='circle')
    if marker is not None:
        marker_.update(marker)

    traces = []

    if stl is not None:
        traces.extend(stl2mesh3d(filename, opacity=0.3) for filename in stl)

    for serial, (x, y, z) in positions.items():
        trace = go.Scatter3d(
            x=x, y=y, z=z, mode='markers',
            marker=marker_, name=serial)
        traces.append(trace)

    layout_ = dict(scene_aspectmode='data', margin=dict(l=0, r=0, b=0, t=0))
    if layout is not None:
        layout_.update(layout)
    layout_ = go.Layout(**layout_)
    fig = go.Figure(data=traces, layout=layout_)
    plotly.offline.iplot(fig)

plot_world_with_plotly = plot_world

def plot_trajectory(
        obs, max_count=10, line=None, layout=None, stl=None, **kwargs):
    """
    Generate a plot from received instance of TrajectoryObserver and show it
    on IPython notebook.

    Parameters
    ----------
    obs : TrajectoryObserver
        TrajectoryObserver to render.
    max_count : Integer, default 10
        The maximum number of particles to show. If None, show all.
    line : dict, optional
        The custom properties for lines.
        See also https://plot.ly/python/line-style/
    layout : dict, optional
        The custom properties for layout.
        See also https://plot.ly/python/reference/#layout

    """
    import numpy

    import plotly
    import plotly.graph_objs as go
    init_notebook_mode()

    line_ = line or {}

    data = obs.data()
    if max_count is not None and len(data) > max_count:
        data = random.sample(data, max_count)

    traces = []

    if stl is not None:
        traces.extend(stl2mesh3d(filename, opacity=0.5) for filename in stl)

    for i, trajectory in enumerate(data):
        trajectory = numpy.array([tuple(pos) for pos in trajectory]).T
        trace = go.Scatter3d(
            x=trajectory[0], y=trajectory[1], z=trajectory[2],
            line=line_, mode='lines', name=str(i))
        traces.append(trace)

    layout_ = dict(scene_aspectmode='data', margin=dict(l=0, r=0, b=0, t=0))
    if layout is not None:
        layout_.update(layout)
    layout_ = go.Layout(**layout_)
    fig = go.Figure(data=traces, layout=layout_)
    plotly.offline.iplot(fig)

plot_trajectory_with_plotly = plot_trajectory
