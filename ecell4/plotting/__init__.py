"""ecell4.plotting: Visualizer of particles based on D3.js, THREE.js
and Elegans.
"""

import os.path
import base64
import copy
import random
import types

from ..util.session import load_world

from .styles import default_color_scale, attractive_mpl_color_scale
from ._core import get_range_of_world, get_range_of_trajectories, display_anim

from . import _matplotlib
BACKEND = _matplotlib


def __on_ipython_notebook():
    try:
        import IPython.terminal.interactiveshell
        if isinstance(get_ipython(), IPython.terminal.interactiveshell.TerminalInteractiveShell):
            return False
    except ImportError:
        return False
    except NameError:
        return False
    return True

def plot_number_observer(*args, backend=None, **kwargs):
    """
    Generate a plot from NumberObservers and show it.
    See plot_number_observer_with_matplotlib and _with_nya for details.

    Parameters
    ----------
    obs : NumberObserver (e.g. FixedIntervalNumberObserver)
    backend : str, optional
        backend. Either one of 'matplotlib', 'plotly' or 'elegans' is supported.

    Examples
    --------
    >>> plot_number_observer(obs1)
    >>> plot_number_observer(obs1, backend='plotly')

    """
    if backend == 'matplotlib':
        _matplotlib.plot_number_observer(*args, **kwargs)
    elif backend == 'plotly':
        from . import _plotly
        _plotly.plot_number_observer(*args, **kwargs)
    elif backend == 'elegans':
        from . import _elegans
        _elegans.plot_number_observer(*args, **kwargs)
    else:
        BACKEND.plot_number_observer(*args, **kwargs)

def plot_world(*args, backend=None, **kwargs):
    """
    Generate a plot from received instance of World and show it.
    See also plot_world_with_elegans and plot_world_with_matplotlib.

    Parameters
    ----------
    world : World or str
        World or a HDF5 filename to render.
    backend : str, optional
        backend. Either one of 'matplotlib', 'plotly' or 'elegans' is supported.

    Examples
    --------
    >>> plot_world(w)
    >>> plot_world(w, backend='plotly')

    """
    if backend == 'matplotlib':
        _matplotlib.plot_world(*args, **kwargs)
    elif backend == 'plotly':
        from . import _plotly
        _plotly.plot_world(*args, **kwargs)
    elif backend == 'elegans':
        from . import _elegans
        _elegans.plot_world(*args, **kwargs)
    else:
        BACKEND.plot_world(*args, **kwargs)

def plot_movie(*args, backend=None, **kwargs):
    """
    Generate a movie from received instances of World and show them.
    See also plot_movie_with_elegans and plot_movie_with_matplotlib.

    Parameters
    ----------
    worlds : list of World
        Worlds to render.
    backend : str, optional
        backend. Either one of 'matplotlib' or 'elegans' is supported.

    """
    if backend == 'matplotlib':
        _matplotlib.plot_movie(*args, **kwargs)
    elif backend == 'elegans':
        from . import _elegans
        _elegans.plot_movie(*args, **kwargs)
    else:
        BACKEND.plot_movie(*args, **kwargs)

def plot_trajectory(*args, backend=None, **kwargs):
    """
    Generate a plot from received instance of TrajectoryObserver and show it
    See also plot_trajectory_with_plotly, plot_trajectory_with_elegans
    and plot_trajectory_with_matplotlib.

    Parameters
    ----------
    obs : TrajectoryObserver
        TrajectoryObserver to render.
    backend : str, optional
        backend. Either one of 'matplotlib', 'plotly' or 'elegans' is supported.

    Examples
    --------
    >>> plot_trajectory(obs)

    """
    if backend == 'matplotlib':
        _matplotlib.plot_trajectory(*args, **kwargs)
    elif backend == 'plotly':
        from . import _plotly
        _plotly.plot_trajectory(*args, **kwargs)
    elif backend == 'elegans':
        from . import _elegans
        _elegans.plot_trajectory(*args, **kwargs)
    else:
        BACKEND.plot_trajectory(*args, **kwargs)

plot_movie_of_trajectory = _matplotlib.plot_movie_of_trajectory_with_matplotlib  # default

def logo(x=1, y=None):
    if not isinstance(x, int):
        x = 1
    else:
        x = min(10, max(1, x))
    if y is None or not isinstance(y, int):
        y = 1
    else:
        y = min(10, max(1, y))

    from IPython.core.display import display, HTML, Javascript

    template = """<script type="text/javascript">
    var y = 0;
    var running = false, stop = true;
    var base64a = ["%s", "%s", "%s", "%s", "%s",
        "%s", "%s", "%s", "%s", "%s",
        "%s", "%s", "%s", "%s", "%s"];
    var maxcnt = base64a.length;
    var timer_id;

    function move() {
        if (running)
        {
            y = (y + 1) %% maxcnt;
            var logos = document.getElementsByName('ecelllogo');
            for (var i = 0; i < logos.length; i++) {
                logos[i].src = "data:image/png;base64," + base64a[y + 1];
            }
            if (stop && y == maxcnt - 1) {
                // clearInterval(id);
                running = false;
                stop = true;
            }
        }
    }

    function action() {
        if (!stop) {
            stop = true;
        }
        else if (!running) {
            running = true;
            stop = false;
            if (timer_id != undefined) {
                clearInterval(timer_id);
            }
            timer_id = setInterval('move();', 120);
        }
    }
    </script>
    %s
    """

    filenames = [
       os.path.join(os.path.abspath(os.path.dirname(__file__)),
                    'templates/ecelllogo/logo%02d.png' % (i + 1))
       for i in range(15)]
    base64s = [
        base64.b64encode(open(filename, 'rt').read())
        for filename in filenames]
    img_html = ('<img name="ecelllogo" style="position:relative;'
                + ' left:0px;" alt="ecelllogo"'
                + ' src="data:image/png;base64,%s"' % (base64s[0])
                + ' onClick="action();" />')
    h = HTML(template % tuple(base64s + [("<p>%s</p>" % (img_html * x)) * y]))
    display(h)

def display_pdb(entity, width=400, height=400):
    from IPython.display import display, IFrame
    import ecell4.datasource.pdb as pdb
    entity_id = pdb.PDBDataSource.parse_entity(entity)
    if entity is None:
        raise ValueError('An invalid entity [{}] was given.'.format(repr(entity)))
    display(IFrame("https://gjbekker.github.io/molmil/#molmil.loadPDB('{}');".format(entity_id), width, height))

def __scatter_world_with_attractive_mpl(
        world, ax, species_list, marker_size, max_count, scale, **kwargs):
    from ecell4_base.core import Species
    color_scale = attractive_mpl_color_scale({})

    scatters, plots = [], []
    for i, name in enumerate(species_list):
        xs, ys, zs = [], [], []
        particles = world.list_particles_exact(Species(name))
        if max_count is not None and len(particles) > max_count:
            particles = random.sample(particles, max_count)
        for pid, p in particles:
            pos = p.position() * scale
            xs.append(pos[0])
            ys.append(pos[1])
            zs.append(pos[2])
        c = color_scale.get_color(name)
        opts = dict(marker='o', s=(2 ** marker_size), edgecolors='white', alpha=0.7)
        opts.update(kwargs)
        scatters.append(
            ax.scatter(xs, ys, zs, facecolor=c, label=name, **opts))
        # plots.extend(ax.plot([], [], 'o', c=c, markeredgecolor='white', label=name))  #XXX: A dirty hack to show the legends with keeping the 3d transparency effect on scatter
    return scatters, plots

def plot_movie_with_attractive_mpl(
        worlds, marker_size=6, figsize=6, grid=True,
        wireframe=False, species_list=None, max_count=None, angle=None, noaxis=False,
        interval=0.16, repeat_delay=3000, stride=1, rotate=None,
        legend=True, whratio=1.33, scale=1, output=None, crf=10, bitrate='1M', **kwargs):
    """
    Generate a move from the received list of instances of World,
    and show it on IPython notebook. This function may require ffmpeg.

    Parameters
    ----------
    worlds : list or FixedIntervalHDF5Observer
        A list of Worlds to render.
    marker_size : float, default 3
        Marker size for all species. Size is passed to scatter function
        as argument, s=(2 ** marker_size).
    figsize : float, default 6
        Size of the plotting area. Given in inch.
    species_list : array of string, default None
        If set, plot_world will not search the list of species.
    max_count : Integer, default None
        The maximum number of particles to show for each species.
        None means no limitation.
    angle : tuple, default None
        A tuple of view angle which is given as (azim, elev, dist).
        If None, use default assumed to be (-60, 30, 10).
    interval : Integer, default 0.16
        Parameters for matplotlib.animation.ArtistAnimation.
    stride : Integer, default 1
        Stride per frame.
    rotate : tuple, default None
        A pair of rotation angles, elev and azim, for animation.
        None means no rotation, same as (0, 0).
    legend : bool, default True
    whratio : float, default 1.33
        A ratio between figure width and height.
        Customize this to keep a legend within the figure.
    scale : float, default 1
        A length-scaling factor
    crf : int, default 10
        The CRF value can be from 4-63. Lower values mean better quality.
    bitrate : str, default '1M'
        Target bitrate
    output : str, default None
        An output filename. '.webm' or '.mp4' is only accepted.
        If None, display a movie on IPython Notebook.

    """
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from ecell4_base.core import Species, FixedIntervalHDF5Observer

    # print("Start generating species_list ...")

    if isinstance(worlds, FixedIntervalHDF5Observer):
        obs = worlds
        worlds = []
        for i in range(0, obs.num_steps(), stride):
            filename = obs.filename(i)
            if os.path.isfile(filename):
                worlds.append(load_world(filename))
            elif len(worlds) >0:
                worlds.append(worlds[-1])
    else:
        worlds = worlds[:: stride]

    if species_list is None:
        species_list = []
        for world in worlds:
            species_list.extend(
                [p.species().serial() for pid, p in world.list_particles()])
            species_list = sorted(
                set(species_list), key=species_list.index)  # XXX: pick unique ones

    # print("Start preparing mplot3d ...")

    fig, ax = __prepare_mplot3d_with_attractive_mpl(
        get_range_of_world(worlds[0], scale), figsize, grid, wireframe, angle,
        noaxis, whratio)

    from mpl_toolkits.mplot3d.art3d import juggle_axes

    def _update_plot(i, scatters, worlds, species_list):
        world = worlds[i]
        for i, name in enumerate(species_list):
            xs, ys, zs = [], [], []
            particles = world.list_particles_exact(Species(name))
            if max_count is not None and len(particles) > max_count:
                particles = random.sample(particles, max_count)
            for pid, p in particles:
                pos = p.position() * scale
                xs.append(pos[0])
                ys.append(pos[1])
                zs.append(pos[2])
            scatters[i]._offsets3d = juggle_axes(xs, ys, zs, 'z')

        if rotate is not None:
            ax.elev += rotate[0]
            ax.azim += rotate[1]

        fig.canvas.draw()

    # print("Start making animation ...")

    color_scale = attractive_mpl_color_scale({})
    scatters = []
    for i, name in enumerate(species_list):
        opts = dict(marker='o', s=(2 ** marker_size), edgecolors='white', alpha=0.7)
        opts.update(kwargs)
        scatters.append(
            ax.scatter(
                [], [], [], facecolor=color_scale.get_color(name), label=name, **opts))

    # if legend:
    #     ax.legend(loc='best', shadow=True)
    if legend is not None and legend is not False:
        legend_opts = {'loc': 'center left', 'bbox_to_anchor': (1.0, 0.5),
                       'shadow': False, 'frameon': False, 'fontsize': 'x-large',
                       'scatterpoints': 1}
        if isinstance(legend, dict):
            legend_opts.update(legend)
        ax.legend(**legend_opts)

    ani = animation.FuncAnimation(
        fig, _update_plot, fargs=(scatters, worlds, species_list),
        frames=len(worlds), interval=interval, blit=False)

    plt.close(ani._fig)

    # print("Start generating a movie ...")
    display_anim(ani, output, fps=1.0 / interval, crf=crf, bitrate=bitrate)

def plot_world_with_attractive_mpl(
        world, marker_size=6, figsize=6, grid=True,
        wireframe=False, species_list=None, max_count=1000, angle=None,
        legend=True, noaxis=False, whratio=1.33, scale=1.0, **kwargs):
    """
    Generate a plot from received instance of World and show it on IPython notebook.

    Parameters
    ----------
    world : World or str
        World to render. A HDF5 filename is also acceptable.
    marker_size : float, default 3
        Marker size for all species. Size is passed to scatter function
        as argument, s=(2 ** marker_size).
    figsize : float, default 6
        Size of the plotting area. Given in inch.
    species_list : array of string, default None
        If set, plot_world will not search the list of species.
    max_count : Integer, default 1000
        The maximum number of particles to show for each species.
        None means no limitation.
    angle : tuple, default None
        A tuple of view angle which is given as (azim, elev, dist).
        If None, use default assumed to be (-60, 30, 10).
    legend : bool, default True
    whratio : float, default 1.33
        A ratio between figure width and height.
        Customize this to keep a legend within the figure.
    scale : float, default 1
        A length-scaling factor

    """
    import matplotlib.pyplot as plt

    if species_list is None:
        species_list = [p.species().serial() for pid, p in world.list_particles()]
        species_list = sorted(
            set(species_list), key=species_list.index)  # XXX: pick unique ones

    fig, ax = __prepare_mplot3d_with_attractive_mpl(
        get_range_of_world(world, scale), figsize, grid, wireframe, angle,
        noaxis, whratio)
    scatters, plots = __scatter_world_with_attractive_mpl(
        world, ax, species_list, marker_size, max_count, scale, **kwargs)

    # if legend:
    #     ax.legend(handles=plots, labels=species_list, loc='best', shadow=True)
    if legend is not None and legend is not False:
        legend_opts = {'loc': 'center left', 'bbox_to_anchor': (1.0, 0.5),
                       'shadow': False, 'frameon': False, 'fontsize': 'x-large',
                       'scatterpoints': 1}
        if isinstance(legend, dict):
            legend_opts.update(legend)
        ax.legend(**legend_opts)
        # ax.legend(handles=plots, labels=species_list,  **legend_opts)

    plt.show()

def __prepare_mplot3d_with_attractive_mpl(
        wrange, figsize, grid, wireframe, angle, noaxis, whratio):
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.axis3d import Axis
    import matplotlib.pyplot as plt
    import itertools

    fig = plt.figure(figsize=(figsize * whratio, figsize))
    ax = plt.subplot(111, projection='3d')
    ax.set_facecolor('white')

    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        if wireframe:
            axis._axinfo['grid']['color'] = (0.9176470588235294, 0.9176470588235294, 0.9490196078431372, 1)
            axis._axinfo['grid']['linewidth'] = 0.6
        else:
            axis._axinfo['grid']['color'] = (1, 1, 1, 1)
            axis._axinfo['grid']['linewidth'] = 1.0

        for tick in axis.get_major_ticks():
            tick.label.set_fontsize(14)

    ax.set_xlim(*wrange['x'])
    ax.set_ylim(*wrange['y'])
    ax.set_zlim(*wrange['z'])
    ax.set_xlabel('X', fontsize=20, labelpad=12)
    ax.set_ylabel('Y', fontsize=20, labelpad=12)
    ax.set_zlabel('Z', fontsize=20, labelpad=12)

    for axis in (ax.w_xaxis, ax.w_yaxis, ax.w_zaxis):
        axis.line.set_color("white")
        axis.set_pane_color((0.9176470588235294, 0.9176470588235294, 0.9490196078431372, 0 if wireframe else 1))

    for line in itertools.chain(ax.get_xticklines(),
                                ax.get_yticklines(),
                                ax.get_zticklines()):
        line.set_visible(False)

    ax.grid(grid)

    if noaxis:
        ax.set_axis_off()

    if angle is not None:
        ax.azim, ax.elev, ax.dist = angle

    plt.subplots_adjust(left=0.0, right=1.0 / whratio, top=1.02, bottom=0.02)
    return (fig, ax)
