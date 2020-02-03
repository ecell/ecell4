import random

from ..util.session import load_world

from .styles import plotly_color_scale

__all__ = [
    'plot_number_observer_with_plotly',
    'plot_world_with_plotly',
    ]


def plot_number_observer(*args, **kwargs):
    import numpy

    import plotly
    import plotly.graph_objs as go

    step = kwargs.get('step', False)
    color_scale = plotly_color_scale()

    plotly.offline.init_notebook_mode()
    fig = go.Figure()

    for obs in args:
        data = numpy.array(obs.data()).T

        targets = [sp.serial() for sp in obs.targets()]
        targets = list(enumerate(targets))

        for idx, serial in targets:
            showlegend = (serial not in color_scale.get_config())
            trace = go.Scatter(x=data[0], y=data[idx + 1], name=serial, line_shape=('linear' if not step else 'hv'), line_color=color_scale.get_color(serial), legendgroup=serial, showlegend=showlegend)
            fig.add_trace(trace)

    fig.update(dict(layout=dict(xaxis_title='Time', yaxis_title='The Number of Molecules')))
    plotly.offline.iplot(fig)

plot_number_observer_with_plotly = plot_number_observer

def plot_world(world, species_list=None, max_count=1000):
    """
    Plot a World on IPython Notebook
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

    marker = dict(size=6, line=dict(color='rgb(204, 204, 204)', width=1),
                  opacity=0.9, symbol='circle')

    data = []
    for serial, (x, y, z) in positions.items():
        trace = go.Scatter3d(
            x=x, y=y, z=z, mode='markers',
            marker=marker, name=serial)
        data.append(trace)

    layout = go.Layout(margin=dict(l=0, r=0, b=0, t=0))
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.iplot(fig)

plot_world_with_plotly = plot_world
