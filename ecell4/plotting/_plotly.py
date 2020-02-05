import random
import types

from ..util.session import load_world

from .styles import plotly_color_scale

__all__ = [
    'plot_number_observer_with_plotly',
    'plot_world_with_plotly',
    ]


def plot_number_observer(*args, step=False, layout=None, **kwargs):
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

    color_scale = plotly_color_scale()

    plotly.offline.init_notebook_mode()
    fig = go.Figure()

    data = None
    xidx = 0
    for obs in args:
        if isinstance(obs, types.FunctionType):
            if data is None:
                raise ValueError("A function must be given after an observer.")
            y = [obs(xi) for xi in data[xidx]]
            label = obs.__name__
            showlegend = (label not in color_scale.get_config())
            trace = go.Scatter(x=data[xidx], y=y, name=label, line_color=color_scale.get_color(label), legendgroup=label, showlegend=showlegend)
            fig.add_trace(trace)
            continue

        data = numpy.array(obs.data()).T

        targets = [sp.serial() for sp in obs.targets()]
        targets = list(enumerate(targets))

        for idx, serial in targets:
            showlegend = (serial not in color_scale.get_config())
            trace = go.Scatter(x=data[xidx], y=data[idx + 1], name=serial, line_shape=('linear' if not step else 'hv'), line_color=color_scale.get_color(serial), legendgroup=serial, showlegend=showlegend)
            fig.add_trace(trace)

    layout_ = dict(xaxis_title='Time', yaxis_title='The Number of Molecules')
    if layout is not None:
        layout_.update(layout)
    fig.update(dict(layout=layout_))

    plotly.offline.iplot(fig)

plot_number_observer_with_plotly = plot_number_observer

def plot_world(world, species_list=None, max_count=1000, marker=None, layout=None):
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
        The custom properties for marker.
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

    plotly.offline.init_notebook_mode()

    marker_ = dict(size=6, line=dict(color='rgb(204, 204, 204)', width=1),
            opacity=0.9, symbol='circle')
    if marker is not None:
        marker_.update(marker)

    data = []
    for serial, (x, y, z) in positions.items():
        trace = go.Scatter3d(
            x=x, y=y, z=z, mode='markers',
            marker=marker_, name=serial)
        data.append(trace)

    layout_ = dict(margin=dict(l=0, r=0, b=0, t=0))
    if layout is not None:
        layout_.update(layout)
    layout_ = go.Layout(**layout_)
    fig = go.Figure(data=data, layout=layout_)
    plotly.offline.iplot(fig)

plot_world_with_plotly = plot_world
