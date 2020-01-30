import os
import types
import random

from .styles import default_color_scale, matplotlib_color_scale

__all__ = [
    "plot_number_observer_with_matplotlib",
    "plot_world_with_matplotlib",
    "plot_trajectory_with_matplotlib",
    "plot_trajectory2d_with_matplotlib",
    "plot_movie_of_trajectory2d_with_matplotlib",
    "plot_movie_with_matplotlib",
    "plot_movie_of_trajectory_with_matplotlib",
    "plot_world2d_with_matplotlib",
    "plot_movie2d_with_matplotlib",
    ]

def plot_number_observer_with_matplotlib(*args, **kwargs):
    """
    Generate a plot from NumberObservers and show it on IPython notebook
    with matplotlib.

    Parameters
    ----------
    obs : NumberObserver (e.g. FixedIntervalNumberObserver)
    fmt : str, optional
    opt : dict, optional
        matplotlib plot options.

    Examples
    --------
    >>> plot_number_observer(obs1)
    >>> plot_number_observer(obs1, 'o')
    >>> plot_number_observer(obs1, obs2, obs3, {'linewidth': 2})
    >>> plot_number_observer(obs1, 'k-', obs2, 'k--')

    """
    import matplotlib.pylab as plt
    import numpy
    import collections

    special_keys = ("xlim", "ylim", "xlabel", "ylabel", "legend", "x", "y", "filename", "step")
    plot_opts = {key: value for key, value in kwargs.items()
                 if key not in special_keys}

    step = kwargs.get('step', None)
    if step is True:
        step = 'post'

    if 'axes.prop_cycle' in plt.rcParams.keys():
        color_cycle = [prop['color'] for prop in plt.rcParams['axes.prop_cycle']]
    else:
        color_cycle = plt.rcParams['axes.color_cycle']

    if "y" in kwargs.keys() and isinstance(kwargs["y"], str):
        kwargs["y"] = (kwargs["y"], )

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if len(args) > 1 and isinstance(args[1], str):
        if len(args) % 2 == 0:
            observers = [(args[i], args[i + 1]) for i in range(0, len(args), 2)]
        else:
            observers = [(args[i], args[i + 1]) for i in range(0, len(args) - 1, 2)]
            observers.append(args[-1], None)
    else:
        observers = [(obs, None) for obs in args]

    color_map = {}
    data, xidx = None, 0
    for obs, fmt in observers:
        if isinstance(obs, types.FunctionType):
            if data is None:
                raise ValueError("A function must be given after an observer.")
            y = [obs(xi) for xi in data[xidx]]
            opts = plot_opts.copy()
            label = obs.__name__
            opts["label"] = label
            if label not in color_map.keys():
                color_map[label] = color_cycle[len(color_map) % len(color_cycle)]
                opts["label"] = label
            opts["color"] = color_map[label]
            if fmt is None:
                ax.plot(data[xidx], y, **opts)
            else:
                ax.plot(data[xidx], y, fmt, **opts)
            continue

        data = numpy.array(obs.data()).T

        try:
            err = obs.error().T
        except AttributeError:
            err = None

        if "x" in kwargs.keys():
            targets = [sp.serial() for sp in obs.targets()]
            if kwargs["x"] not in targets:
                raise ValueError("[{0}] given as 'x' was not found.".fomrat(kwargs["x"]))
            xidx = targets.index(kwargs["x"]) + 1
        else:
            xidx = 0

        if "y" in kwargs.keys():
            targets = [sp.serial() for sp in obs.targets()]
            targets = [(targets.index(serial), serial)
                       for serial in kwargs["y"] if serial in targets]
        else:
            targets = [sp.serial() for sp in obs.targets()]
            targets = list(enumerate(targets))
            # targets.sort(key=lambda x: x[1])

        for idx, serial in targets:
            opts = plot_opts.copy()

            label = serial
            if len(label) > 0 and label[0] == '_':
                label = '$\_$' + label[1:]  # XXX: lazy escaping for a special character
            if label not in color_map.keys():
                color_map[label] = color_cycle[len(color_map) % len(color_cycle)]
                opts["label"] = label
            opts["color"] = color_map[label]

            if err is None:
                if step is None:
                    if fmt is None:
                        ax.plot(data[xidx], data[idx + 1], **opts)
                    else:
                        ax.plot(data[xidx], data[idx + 1], fmt, **opts)
                else:
                    if fmt is None:
                        ax.step(data[xidx], data[idx + 1], where=step, **opts)
                    else:
                        ax.step(data[xidx], data[idx + 1], fmt, where=step, **opts)
            else:
                if fmt is None:
                    ax.errorbar(data[xidx], data[idx + 1],
                        xerr=(None if xidx == 0 else err[xidx]), yerr=err[idx + 1],
                        **opts)
                else:
                    ax.errorbar(data[xidx], data[idx + 1],
                        xerr=(None if xidx == 0 else err[xidx]), yerr=err[idx + 1],
                        fmt=fmt, **opts)

    # if "legend" not in kwargs.keys() or kwargs["legend"]:
    #     ax.legend(*ax.get_legend_handles_labels(), loc="best", shadow=True)
    if "legend" not in kwargs.keys() or (kwargs["legend"] is not None and kwargs["legend"] is not False):
        legend_opts = {"loc": "best", "shadow": True}
        if "legend" in kwargs and isinstance(kwargs["legend"], dict):
            legend_opts.update(kwargs["legend"])
        ax.legend(*ax.get_legend_handles_labels(), **legend_opts)

    if "xlabel" in kwargs.keys():
        ax.set_xlabel(kwargs["xlabel"])
    elif "x" in kwargs.keys():
        ax.set_xlabel("The Number of Molecules [{0}]".format(kwargs["x"]))
    else:
        ax.set_xlabel("Time")
    if "ylabel" in kwargs.keys():
        ax.set_ylabel(kwargs["ylabel"])
    else:
        ax.set_ylabel("The Number of Molecules")
    if "xlim" in kwargs.keys():
        ax.set_xlim(kwargs["xlim"])
    if "ylim" in kwargs.keys():
        ax.set_ylim(kwargs["ylim"])

    if "filename" in kwargs.keys():
        plt.savefig(kwargs["filename"])
    else:
        plt.show()

def __prepare_mplot3d_with_matplotlib(
        wrange, figsize, grid, wireframe, angle, noaxis):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(figsize, figsize))
    ax = fig.gca(projection='3d')

    if wireframe:
        ax.w_xaxis.set_pane_color((0, 0, 0, 0))
        ax.w_yaxis.set_pane_color((0, 0, 0, 0))
        ax.w_zaxis.set_pane_color((0, 0, 0, 0))

    ax.grid(grid)
    ax.set_xlim(*wrange['x'])
    ax.set_ylim(*wrange['y'])
    ax.set_zlim(*wrange['z'])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    if noaxis:
        ax.set_axis_off()

    if angle is not None:
        ax.azim, ax.elev, ax.dist = angle

    return (fig, ax)

def __scatter_world_with_matplotlib(
        world, ax, species_list, marker_size, max_count, **kwargs):
    from ecell4_base.core import Species
    color_scale = matplotlib_color_scale()

    scatters, plots = [], []
    for i, name in enumerate(species_list):
        xs, ys, zs = [], [], []
        particles = world.list_particles_exact(Species(name))
        if max_count is not None and len(particles) > max_count:
            particles = random.sample(particles, max_count)
        for pid, p in particles:
            pos = p.position()
            xs.append(pos[0])
            ys.append(pos[1])
            zs.append(pos[2])
        c = color_scale.get_color(name)
        scatters.append(
            ax.scatter(
                xs, ys, zs,
                marker='o', s=(2 ** marker_size), lw=0, facecolor=c,
                label=name, **kwargs))
        plots.extend(ax.plot([], [], 'o', c=c, label=name))  #XXX: A dirty hack to show the legends with keeping the 3d transparency effect on scatter
    return scatters, plots

def __plot_trajectory_with_matplotlib(lines, ax, upto=None, **kwargs):
    color_scale = default_color_scale()
    plots = []
    for i, line in enumerate(lines):
        plots.append(
            ax.plot(line[0][: upto], line[1][: upto], line[2][: upto],
                label=i, color=color_scale.get_color(i), **kwargs)[0])
    return plots

def plot_world_with_matplotlib(
        world, marker_size=3, figsize=6, grid=True,
        wireframe=False, species_list=None, max_count=1000, angle=None,
        legend=True, noaxis=False, **kwargs):
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

    """
    import matplotlib.pyplot as plt

    if species_list is None:
        species_list = [p.species().serial() for pid, p in world.list_particles()]
        species_list = sorted(
            set(species_list), key=species_list.index)  # XXX: pick unique ones

    fig, ax = __prepare_mplot3d_with_matplotlib(
        __get_range_of_world(world), figsize, grid, wireframe, angle, noaxis)
    scatters, plots = __scatter_world_with_matplotlib(
        world, ax, species_list, marker_size, max_count, **kwargs)

    # if legend:
    #     ax.legend(handles=plots, labels=species_list, loc='best', shadow=True)
    if legend is not None and legend is not False:
        legend_opts = {"loc": "best", "shadow": True}
        if isinstance(legend, dict):
            legend_opts.update(legend)
        ax.legend(handles=plots, labels=species_list,  **legend_opts)

    plt.show()

def plot_trajectory_with_matplotlib(
        obs, max_count=10, figsize=6, legend=True, angle=None,
        wireframe=False, grid=True, noaxis=False, plot_range=None, **kwargs):
    """
    Generate a plot from received instance of TrajectoryObserver and show it
    on IPython notebook.

    Parameters
    ----------
    obs : TrajectoryObserver
        TrajectoryObserver to render.
    max_count : Integer, default 10
        The maximum number of particles to show. If None, show all.
    figsize : float, default 6
        Size of the plotting area. Given in inch.
    angle : tuple, default None
        A tuple of view angle which is given as (azim, elev, dist).
        If None, use default assumed to be (-60, 30, 10).
    legend : bool, default True
    plot_range : tuple, default None
        Range for plotting. A triplet of pairs suggesting (rangex, rangey, rangez).
        If None, the minimum volume containing all the trajectories is used.

    """
    import matplotlib.pyplot as plt

    data = obs.data()
    if max_count is not None and len(data) > max_count:
        data = random.sample(data, max_count)

    fig, ax = __prepare_mplot3d_with_matplotlib(
        __get_range_of_trajectories(data, plot_range),
        figsize, grid, wireframe, angle, noaxis)

    lines = []
    for i, y in enumerate(data):
        xarr, yarr, zarr = [], [], []
        for pos in y:
            xarr.append(pos[0])
            yarr.append(pos[1])
            zarr.append(pos[2])

        lines.append((xarr, yarr, zarr))

    __plot_trajectory_with_matplotlib(lines, ax, **kwargs)

    # if legend:
    #     ax.legend(loc='best', shadow=True)
    if legend is not None and legend is not False:
        legend_opts = {"loc": "best", "shadow": True}
        if isinstance(legend, dict):
            legend_opts.update(legend)
        ax.legend(**legend_opts)
    plt.show()

def __prepare_plot_with_matplotlib(
        wrange, figsize, grid, wireframe, noaxis):
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(figsize, figsize))
    ax = fig.gca()
    ax.set_aspect('equal')

    # if wireframe:
    #     ax.w_xaxis.set_pane_color((0, 0, 0, 0))
    #     ax.w_yaxis.set_pane_color((0, 0, 0, 0))
    #     ax.w_zaxis.set_pane_color((0, 0, 0, 0))

    ax.grid(grid)
    ax.set_xlim(*wrange['x'])
    ax.set_ylim(*wrange['y'])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    if noaxis:
        ax.set_axis_off()

    return (fig, ax)

def __plot_trajectory2d_with_matplotlib(lines, ax, upto=None, **kwargs):
    color_scale = default_color_scale()
    plots = []
    for i, line in enumerate(lines):
        plots.append(
            ax.plot(line[0][: upto], line[1][: upto],
                label=i, color=color_scale.get_color(i), **kwargs)[0])
    return plots

def plot_trajectory2d_with_matplotlib(
        obs, plane='xy', max_count=10, figsize=6, legend=True,
        wireframe=False, grid=True, noaxis=False, plot_range=None, **kwargs):
    """
    Make a 2D plot from received instance of TrajectoryObserver and show it
    on IPython notebook.

    Parameters
    ----------
    obs : TrajectoryObserver
        TrajectoryObserver to render.
    plane : str, default 'xy'
        'xy', 'yz', 'zx'.
    max_count : Integer, default 10
        The maximum number of particles to show. If None, show all.
    figsize : float, default 6
        Size of the plotting area. Given in inch.
    legend : bool, default True
    plot_range : tuple, default None
        Range for plotting. A triplet of pairs suggesting (rangex, rangey, rangez).
        If None, the minimum volume containing all the trajectories is used.

    """
    import matplotlib.pyplot as plt

    plane = plane.lower()
    if len(plane) != 2 or plane[0] not in ('x', 'y', 'z') or plane[1] not in ('x', 'y', 'z'):
        raise ValueError("invalid 'plane' argument [{}] was given.".format(repr(plane)))
    xidx = 0 if plane[0] == 'x' else (1 if plane[0] == 'y' else 2)
    yidx = 0 if plane[1] == 'x' else (1 if plane[1] == 'y' else 2)

    data = obs.data()
    if max_count is not None and len(data) > max_count:
        data = random.sample(data, max_count)

    wrange = __get_range_of_trajectories(data, plot_range)
    wrange = (wrange['x'], wrange['y'], wrange['z'])
    wrange = {'x': wrange[xidx], 'y': wrange[yidx]}
    fig, ax = __prepare_plot_with_matplotlib(
        wrange, figsize, grid, wireframe, noaxis)
    ax.set_xlabel(plane[0].upper())
    ax.set_ylabel(plane[1].upper())

    lines = []
    for i, y in enumerate(data):
        xarr, yarr, zarr = [], [], []
        for pos in y:
            xarr.append(pos[xidx])
            yarr.append(pos[yidx])

        lines.append((xarr, yarr))

    __plot_trajectory2d_with_matplotlib(lines, ax, **kwargs)

    # if legend:
    #     ax.legend(loc='best', shadow=True)
    if legend is not None and legend is not False:
        legend_opts = {"loc": "best", "shadow": True}
        if isinstance(legend, dict):
            legend_opts.update(legend)
        ax.legend(**legend_opts)
    plt.show()

def plot_movie_of_trajectory2d_with_matplotlib(
        obs, plane='xy', figsize=6, grid=True,
        wireframe=False, max_count=None, angle=None, noaxis=False,
        interval=0.16, repeat_delay=3000, stride=1, rotate=None,
        legend=True, output=None, crf=10, bitrate='1M', plot_range=None, **kwargs):
    """
    Generate a move from the received list of instances of World,
    and show it on IPython notebook. This function may require ffmpeg.

    Parameters
    ----------
    worlds : list or FixedIntervalHDF5Observer
        A list of Worlds to render.
    plane : str, default 'xy'
        'xy', 'yz', 'zx'.
    figsize : float, default 6
        Size of the plotting area. Given in inch.
    max_count : Integer, default None
        The maximum number of particles to show for each species.
        None means no limitation.
    interval : Integer, default 0.16
        Parameters for matplotlib.animation.ArtistAnimation.
    stride : Integer, default 1
        Stride per frame.
    legend : bool, default True
    output : str, default None
        An output filename. '.webm' or '.mp4' is only accepted.
        If None, display a movie on IPython Notebook.
    crf : int, default 10
        The CRF value can be from 4-63. Lower values mean better quality.
    bitrate : str, default '1M'
        Target bitrate
    plot_range : tuple, default None
        Range for plotting. A triplet of pairs suggesting (rangex, rangey, rangez).
        If None, the minimum volume containing all the trajectories is used.

    """
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from IPython.display import display, HTML
    from ecell4_base.core import Species, FixedIntervalHDF5Observer
    from .simulation import load_world
    import math

    # print("Taking all data ...")

    plane = plane.lower()
    if len(plane) != 2 or plane[0] not in ('x', 'y', 'z') or plane[1] not in ('x', 'y', 'z'):
        raise ValueError("invalid 'plane' argument [{}] was given.".format(repr(plane)))
    xidx = 0 if plane[0] == 'x' else (1 if plane[0] == 'y' else 2)
    yidx = 0 if plane[1] == 'x' else (1 if plane[1] == 'y' else 2)

    data = obs.data()
    if max_count is not None and len(data) > max_count:
        data = random.sample(data, max_count)

    lines = []
    num_frames = 0
    for i, y in enumerate(data):
        xarr, yarr, zarr = [], [], []
        for pos in y:
            xarr.append(pos[xidx])
            yarr.append(pos[yidx])

        lines.append((xarr, yarr))
        num_frames = max(num_frames, len(y))
    num_frames = int(math.ceil(float(num_frames) / stride))

    # print("Start preparing mplot3d ...")

    wrange = __get_range_of_trajectories(data, plot_range)
    wrange = (wrange['x'], wrange['y'], wrange['z'])
    wrange = {'x': wrange[xidx], 'y': wrange[yidx]}
    fig, ax = __prepare_plot_with_matplotlib(
        wrange, figsize, grid, wireframe, noaxis)
    ax.set_xlabel(plane[0].upper())
    ax.set_ylabel(plane[1].upper())

    def _update_plot(i, plots, lines):
        upto = i * stride
        for plot, line in zip(plots, lines):
            plot.set_xdata(line[0][: upto])
            plot.set_ydata(line[1][: upto])

        fig.canvas.draw()

    # print("Start making animation ...")

    plots = __plot_trajectory2d_with_matplotlib(lines, ax, 0, **kwargs)

    # if legend:
    #     ax.legend(loc='best', shadow=True)
    if legend is not None and legend is not False:
        legend_opts = {"loc": "best", "shadow": True}
        if isinstance(legend, dict):
            legend_opts.update(legend)
        ax.legend(**legend_opts)

    ani = animation.FuncAnimation(
        fig, _update_plot, fargs=(plots, lines),
        frames=num_frames, interval=interval, blit=False)

    plt.close(ani._fig)
    # print("Start generating a movie ...")
    display_anim(ani, output, fps=1.0 / interval, crf=crf, bitrate=bitrate)

def plot_movie_with_matplotlib(
        worlds, marker_size=3, figsize=6, grid=True,
        wireframe=False, species_list=None, max_count=None, angle=None, noaxis=False,
        interval=0.16, repeat_delay=3000, stride=1, rotate=None,
        legend=True, output=None, crf=10, bitrate='1M', **kwargs):
    """
    Generate a movie from the received list of instances of World,
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
    output : str, default None
        An output filename. '.webm' or '.mp4' is only accepted.
        If None, display a movie on IPython Notebook.
    crf : int, default 10
        The CRF value can be from 4-63. Lower values mean better quality.
    bitrate : str, default '1M'
        Target bitrate

    """
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from ecell4_base.core import Species, FixedIntervalHDF5Observer
    from .simulation import load_world

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

    fig, ax = __prepare_mplot3d_with_matplotlib(
        __get_range_of_world(worlds[0]), figsize, grid, wireframe, angle, noaxis)

    from mpl_toolkits.mplot3d.art3d import juggle_axes

    def _update_plot(i, scatters, worlds, species_list):
        world = worlds[i]
        for i, name in enumerate(species_list):
            xs, ys, zs = [], [], []
            particles = world.list_particles_exact(Species(name))
            if max_count is not None and len(particles) > max_count:
                particles = random.sample(particles, max_count)
            for pid, p in particles:
                pos = p.position()
                xs.append(pos[0])
                ys.append(pos[1])
                zs.append(pos[2])
            scatters[i]._offsets3d = juggle_axes(xs, ys, zs, 'z')

        if rotate is not None:
            ax.elev += rotate[0]
            ax.azim += rotate[1]

        fig.canvas.draw()

    # print("Start making animation ...")

    color_scale = matplotlib_color_scale()
    scatters = []
    for i, name in enumerate(species_list):
        scatters.append(
            ax.scatter([], [], [], marker='o', s=(2 ** marker_size),
                       lw=0, facecolor=color_scale.get_color(name), label=name))

    # if legend:
    #     ax.legend(loc='best', shadow=True)
    if legend is not None and legend is not False:
        legend_opts = {"loc": "best", "shadow": True}
        if isinstance(legend, dict):
            legend_opts.update(legend)
        ax.legend(**legend_opts)

    ani = animation.FuncAnimation(
        fig, _update_plot, fargs=(scatters, worlds, species_list),
        frames=len(worlds), interval=interval, blit=False)

    plt.close(ani._fig)
    # print("Start generating a movie ...")
    display_anim(ani, output, fps=1.0 / interval, crf=crf, bitrate=bitrate)

def plot_movie_of_trajectory_with_matplotlib(
        obs, figsize=6, grid=True,
        wireframe=False, max_count=None, angle=None, noaxis=False,
        interval=0.16, repeat_delay=3000, stride=1, rotate=None,
        legend=True, output=None, crf=10, bitrate='1M', plot_range=None, **kwargs):
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
    output : str, default None
        An output filename. '.webm' or '.mp4' is only accepted.
        If None, display a movie on IPython Notebook.
    crf : int, default 10
        The CRF value can be from 4-63. Lower values mean better quality.
    bitrate : str, default '1M'
        Target bitrate
    plot_range : tuple, default None
        Range for plotting. A triplet of pairs suggesting (rangex, rangey, rangez).
        If None, the minimum volume containing all the trajectories is used.

    """
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from ecell4_base.core import Species, FixedIntervalHDF5Observer
    from .simulation import load_world
    import math

    # print("Taking all data ...")

    data = obs.data()
    if max_count is not None and len(data) > max_count:
        data = random.sample(data, max_count)

    lines = []
    num_frames = 0
    for i, y in enumerate(data):
        xarr, yarr, zarr = [], [], []
        for pos in y:
            xarr.append(pos[0])
            yarr.append(pos[1])
            zarr.append(pos[2])

        lines.append((xarr, yarr, zarr))
        num_frames = max(num_frames, len(y))
    num_frames = int(math.ceil(float(num_frames) / stride))

    # print("Start preparing mplot3d ...")

    fig, ax = __prepare_mplot3d_with_matplotlib(
        __get_range_of_trajectories(data, plot_range),
        figsize, grid, wireframe, angle, noaxis)

    def _update_plot(i, plots, lines):
        upto = i * stride
        for plot, line in zip(plots, lines):
            plot.set_data(line[0][: upto], line[1][: upto])
            plot.set_3d_properties(line[2][: upto])

        if rotate is not None:
            ax.elev += rotate[0]
            ax.azim += rotate[1]

        fig.canvas.draw()

    # print("Start making animation ...")

    plots = __plot_trajectory_with_matplotlib(lines, ax, 0, **kwargs)

    # if legend:
    #     ax.legend(loc='best', shadow=True)
    if legend is not None and legend is not False:
        legend_opts = {"loc": "best", "shadow": True}
        if isinstance(legend, dict):
            legend_opts.update(legend)
        ax.legend(**legend_opts)

    ani = animation.FuncAnimation(
        fig, _update_plot, fargs=(plots, lines),
        frames=num_frames, interval=interval, blit=False)

    plt.close(ani._fig)
    # print("Start generating a movie ...")
    display_anim(ani, output, fps=1.0 / interval, crf=crf, bitrate=bitrate)

def plot_world2d_with_matplotlib(
        world, plane='xy', marker_size=3, figsize=6, grid=True,
        wireframe=False, species_list=None, max_count=1000, angle=None,
        legend=True, noaxis=False, scale=1.0, **kwargs):
    """
    Make a 2D plot from received instance of World and show it on IPython notebook.

    Parameters
    ----------
    world : World or str
        World to render. A HDF5 filename is also acceptable.
    plane : str, default 'xy'
        'xy', 'yz', 'zx'.
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
    scale : float, default 1
        A length-scaling factor

    """
    import matplotlib.pyplot as plt

    plane = plane.lower()
    if len(plane) != 2 or plane[0] not in ('x', 'y', 'z') or plane[1] not in ('x', 'y', 'z'):
        raise ValueError("invalid 'plane' argument [{}] was given.".format(repr(plane)))
    xidx = 0 if plane[0] == 'x' else (1 if plane[0] == 'y' else 2)
    yidx = 0 if plane[1] == 'x' else (1 if plane[1] == 'y' else 2)

    if species_list is None:
        species_list = [p.species().serial() for pid, p in world.list_particles()]
        species_list = sorted(
            set(species_list), key=species_list.index)  # XXX: pick unique ones

    wrange = __get_range_of_world(world, scale)
    wrange = (wrange['x'], wrange['y'], wrange['z'])
    wrange = {'x': wrange[xidx], 'y': wrange[yidx]}

    fig, ax = __prepare_plot_with_matplotlib(
        wrange, figsize, grid, wireframe, noaxis)
    scatters, plots = __scatter_world2d_with_matplotlib(
        world, (xidx, yidx), ax, species_list, marker_size, max_count, scale, **kwargs)
    ax.set_xlabel(plane[0].upper())
    ax.set_ylabel(plane[1].upper())

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

def __scatter_world2d_with_matplotlib(
        world, indices, ax, species_list, marker_size, max_count, scale, **kwargs):
    from ecell4_base.core import Species
    color_scale = matplotlib_color_scale()

    scatters, plots = [], []
    for i, name in enumerate(species_list):
        xs, ys = [], []
        particles = world.list_particles_exact(Species(name))
        if max_count is not None and len(particles) > max_count:
            particles = random.sample(particles, max_count)
        for pid, p in particles:
            pos = p.position() * scale
            xs.append(pos[indices[0]])
            ys.append(pos[indices[1]])
        c = color_scale.get_color(name)
        scatters.append(
            ax.scatter(
                xs, ys,
                marker='o', s=(2 ** marker_size), lw=0, facecolor=c,
                label=name, **kwargs))
    return scatters, plots

def plot_movie2d_with_matplotlib(
        worlds, plane='xy', marker_size=3, figsize=6, grid=True,
        wireframe=False, species_list=None, max_count=None, angle=None, noaxis=False,
        interval=0.16, repeat_delay=3000, stride=1, rotate=None,
        legend=True, scale=1, output=None, crf=10, bitrate='1M', **kwargs):
    """
    Generate a movie projected on the given plane from the received list
    of instances of World, and show it on IPython notebook.
    This function may require ffmpeg.

    Parameters
    ----------
    worlds : list or FixedIntervalHDF5Observer
        A list of Worlds to render.
    plane : str, default 'xy'
        'xy', 'yz', 'zx'.
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
    scale : float, default 1
        A length-scaling factor
    output : str, default None
        An output filename. '.webm' or '.mp4' is only accepted.
        If None, display a movie on IPython Notebook.
    crf : int, default 10
        The CRF value can be from 4-63. Lower values mean better quality.
    bitrate : str, default '1M'
        Target bitrate

    """
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from ecell4_base.core import Species, FixedIntervalHDF5Observer
    from .simulation import load_world

    plane = plane.lower()
    if len(plane) != 2 or plane[0] not in ('x', 'y', 'z') or plane[1] not in ('x', 'y', 'z'):
        raise ValueError("invalid 'plane' argument [{}] was given.".format(repr(plane)))
    xidx = 0 if plane[0] == 'x' else (1 if plane[0] == 'y' else 2)
    yidx = 0 if plane[1] == 'x' else (1 if plane[1] == 'y' else 2)

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

    wrange = __get_range_of_world(worlds[0], scale)
    wrange = (wrange['x'], wrange['y'], wrange['z'])
    wrange = {'x': wrange[xidx], 'y': wrange[yidx]}

    fig = plt.figure(figsize=(figsize, figsize))
    ax = fig.gca()

    color_scale = matplotlib_color_scale()

    def _update_plot(i, worlds, species_list):
        ax.cla()

        ax.set_aspect('equal')
        ax.grid(grid)
        ax.set_xlim(*wrange['x'])
        ax.set_ylim(*wrange['y'])
        ax.set_xlabel(plane[0].upper())
        ax.set_ylabel(plane[1].upper())

        if noaxis:
            ax.set_axis_off()

        _legend = False

        world = worlds[i]
        for i, name in enumerate(species_list):
            offsets = ([], [])
            particles = world.list_particles_exact(Species(name))
            if len(particles) == 0:
                continue
            _legend = True

            if max_count is not None and len(particles) > max_count:
                particles = random.sample(particles, max_count)
            for pid, p in particles:
                pos = p.position() * scale
                offsets[0].append(pos[xidx])
                offsets[1].append(pos[yidx])

            ax.scatter(
                offsets[0], offsets[1], marker='o', s=(2 ** marker_size),
                lw=0, facecolor=color_scale.get_color(name), label=name)

        if legend is not None and legend is not False and _legend:
            legend_opts = {"loc": "upper right", "shadow": True}
            if isinstance(legend, dict):
                legend_opts.update(legend)
            ax.legend(**legend_opts)

        fig.canvas.draw()

    ani = animation.FuncAnimation(
        fig, _update_plot, fargs=(worlds, species_list),
        frames=len(worlds), interval=interval, blit=False)

    plt.close(ani._fig)
    display_anim(ani, output, fps=1.0 / interval, crf=crf, bitrate=bitrate)
