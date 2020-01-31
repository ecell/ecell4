"""ecell4.plotting: Visualizer of particles based on D3.js, THREE.js
and Elegans.
"""

import os
import uuid
import json
import base64
import copy
import random
import types
from tempfile import NamedTemporaryFile

from .styles import default_color_scale, attractive_mpl_color_scale

from ._core import get_range_of_world, get_range_of_trajectories
from ._matplotlib import *

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

def plot_number_observer(*args, **kwargs):
    """
    Generate a plot from NumberObservers and show it.
    See plot_number_observer_with_matplotlib and _with_nya for details.

    Parameters
    ----------
    obs : NumberObserver (e.g. FixedIntervalNumberObserver)
    interactive : bool, default False
        Choose a visualizer. If False, show the plot with matplotlib.
        If True (only available on IPython Notebook), show it with nyaplot.

    Examples
    --------
    >>> plot_number_observer(obs1)
    >>> plot_number_observer(obs1, interactive=True)

    """
    interactive = kwargs.pop('interactive', False)
    if interactive:
        plot_number_observer_with_nya(*args, **kwargs)
    # elif __on_ipython_notebook():
    #     kwargs['to_png'] = True
    #     plot_number_observer_with_nya(*args, **kwargs)
    else:
        if kwargs.pop('to_png', None) is not None:
            #XXX: Remove an option available only on nyaplot for the consistency
            import warnings
            warnings.warn(
                "An option 'to_png' is not available with matplotlib. Just ignored.")
        plot_number_observer_with_matplotlib(*args, **kwargs)

def plot_world(*args, **kwargs):
    """
    Generate a plot from received instance of World and show it.
    See also plot_world_with_elegans and plot_world_with_matplotlib.

    Parameters
    ----------
    world : World or str
        World or a HDF5 filename to render.
    interactive : bool, default True
        Choose a visualizer. If False, show the plot with matplotlib.
        If True (only available on IPython Notebook), show it with elegans.

    Examples
    --------
    >>> plot_world(w)
    >>> plot_world(w, interactive=False)

    """
    interactive = kwargs.pop('interactive', True)
    if interactive:
        plot_world_with_elegans(*args, **kwargs)
    else:
        plot_world_with_matplotlib(*args, **kwargs)

def plot_movie(*args, **kwargs):
    """
    Generate a movie from received instances of World and show them.
    See also plot_movie_with_elegans and plot_movie_with_matplotlib.

    Parameters
    ----------
    worlds : list of World
        Worlds to render.
    interactive : bool, default True
        Choose a visualizer. If False, show the plot with matplotlib.
        If True (only available on IPython Notebook), show it with elegans.

    """
    interactive = kwargs.pop('interactive', False)
    if interactive:
        plot_movie_with_elegans(*args, **kwargs)
    else:
        plot_movie_with_matplotlib(*args, **kwargs)

def plot_trajectory(*args, **kwargs):
    """
    Generate a plot from received instance of TrajectoryObserver and show it
    See also plot_trajectory_with_elegans and plot_trajectory_with_matplotlib.

    Parameters
    ----------
    obs : TrajectoryObserver
        TrajectoryObserver to render.
    interactive : bool, default True
        Choose a visualizer. If False, show the plot with matplotlib.
        If True (only available on IPython Notebook), show it with elegans.

    Examples
    --------
    >>> plot_trajectory(obs)
    >>> plot_trajectory(obs, interactive=False)

    """
    interactive = kwargs.pop('interactive', True)
    if interactive:
        plot_trajectory_with_elegans(*args, **kwargs)
    else:
        plot_trajectory_with_matplotlib(*args, **kwargs)

def plot_number_observer_with_nya(obs, config=None, width=600, height=400, x=None, y=None, to_png=False):
    """
    Generate a plot from NumberObservers and show it on IPython notebook
    with nyaplot.

    Parameters
    ----------
    obs : NumberObserver (e.g. FixedIntervalNumberObserver)
    config : dict, optional
        A config data for coloring. The dictionary will be updated during this plot.
    width : int, optional
    height : int, optional
    x : str, optional
        A serial for x-axis. If None, x-axis corresponds time.
    y : str or list of str
        Serials for y axis.

    """
    config = config or {}

    from IPython.core.display import display, HTML
    import numpy

    config = {}
    color_scale = default_color_scale(config=config)

    data1, data2 = [], []
    data = numpy.array(obs.data())

    if x is None:
        xidx = 0
    else:
        tmp = [sp.serial() for sp in obs.targets()]
        if x not in tmp:
            raise ValueError("[{0}] given as 'x' was not found.".fomrat(x))
        xidx = tmp.index(x) + 1

    if y is None:
        targets = [sp.serial() for sp in obs.targets()]
        targets = list(enumerate(targets))
        # targets.sort(key=lambda x: x[1])
    else:
        if isinstance(y, str):
            y = (y, )
        targets = [sp.serial() for sp in obs.targets()]
        targets = [(targets.index(serial), serial)
                   for serial in y if serial in targets]

    for line in data:
        tmp = {"x": line[xidx]}
        for i, (idx, serial) in enumerate(targets):
            tmp["y{0}".format(i + 1)] = line[idx + 1]
        data1.append(tmp)
    for i, (idx, serial) in enumerate(targets):
        label = serial
        tmp = {"type": "line", "data": "data1",
               "options": {"x": "x", "y": "y{0}".format(i + 1),
                           "stroke_width": 2, "title": label,
                           "color": color_scale.get_color(label)}}
        data2.append(tmp)

    xmin, xmax = data.T[xidx].min(), data.T[xidx].max()
    yview = data.T.take([idx + 1 for idx, serial in targets], axis=0)
    ymin, ymax = yview.min(), yview.max()

    model = {
        "data": {"data1": data1},
        "panes": [{"type": 'rectangular',
                   "diagrams": data2,
                   "options": {"width": width, "height": height, "xrange": [xmin, xmax],
                               "yrange": [ymin, ymax], "legend": True, "zoom": True}}]}
    model_id = 'viz{0:s}'.format(str(uuid.uuid4()))
    display(HTML(generate_html(
        {'model': json.dumps(model), 'model_id': model_id, 'to_png': json.dumps(to_png)},
        'templates/nya.tmpl')))

def __parse_world(
        world, radius=None, species_list=None, max_count=None,
        predicator=None):
    """
    Private function to parse world. Return infomation about particles
    (name, coordinates and particle size) for each species.

    """
    from ecell4_base.core import Species

    if species_list is None:
        species_list = [
            p.species().serial() for pid, p in world.list_particles()]
        species_list = sorted(
            set(species_list), key=species_list.index)  # XXX: pick unique ones

    species = []
    for name in species_list:
        particles = [
            {'pos': p.position(), 'r': p.radius()}
            for pid, p in world.list_particles(Species(name))
            if predicator is None or predicator(pid, p)]
        # particles = [
        #     {'pos': p.position(), 'r': p.radius()}
        #     for pid, p in world.list_particles()
        #     if (p.species().serial() == name and
        #         (predicator is None or predicator(pid, p)))]

        if len(particles) == 0:
            continue

        if max_count is not None and len(particles) > max_count:
            particles = random.sample(particles, max_count)

        data = {
            'x': [p['pos'][0] for p in particles],
            'y': [p['pos'][1] for p in particles],
            'z': [p['pos'][2] for p in particles]
        }

        # assume that all particles belong to one species have the same radius
        r = max([p['r'] for p in particles]) if radius is None else radius
        r = r if r > 0 else min(world.edge_lengths()) * 0.005
        size = 30.0 / max(world.edge_lengths()) * r

        species.append({
            'name': name,
            'data': data,
            'size': size
        })

    return species

def plot_movie_with_elegans(
        worlds, radius=None, width=500, height=500, config=None, grid=False,
        species_list=None):
    """
    Generate a movie from received instances of World and show them
    on IPython notebook.

    Parameters
    ----------
    worlds : list of World
        Worlds to render.
    radius : float, default None
        If this value is set, all particles in the world will be rendered
        as if their radius are the same.
    width : float, default 500
        Width of the plotting area.
    height : float, default 500
        Height of the plotting area.
    config : dict, default {}
        Dict for configure default colors. Its values are colors unique
        to each speices. The dictionary will be updated during this plot.
        Colors included in config dict will never be used for other speices.
    species_list : array of string, default None
        If set, plot_movie will not search the list of species

    """
    config = config or {}

    from IPython.core.display import display, HTML
    from jinja2 import Template

    data = {}
    sizes = {}
    for i, world in enumerate(worlds):
        species = __parse_world(world, radius, species_list)
        for species_info in species:
            if data.get(species_info['name']) is None:
                data[species_info['name']] = []
            data[species_info['name']].append({
                'df': species_info['data'],
                't': i
            })
            sizes[species_info['name']] = species_info['size']

    options = {
        'player': True,
        'autorange': False,
        'space_mode': 'wireframe',
        'grid': grid,
        'range': get_range_of_world(worlds[0])
    }

    model_id = '"movie' + str(uuid.uuid4()) + '"'
    color_scale = default_color_scale(config=config)

    display(HTML(generate_html({
        'model_id': model_id,
        'names': json.dumps(list(data.keys())),
        'data': json.dumps(list(data.values())),
        'colors': json.dumps([color_scale.get_color(name)
                              for name in data.keys()]),
        'sizes': json.dumps([sizes[name] for name in data.keys()]),
        'options': json.dumps(options)
    }, 'templates/movie.tmpl')))

def plot_world_with_elegans(
        world, radius=None, width=350, height=350, config=None, grid=True,
        wireframe=False, species_list=None, debug=None, max_count=1000,
        camera_position=(-22, 23, 32), camera_rotation=(-0.6, 0.5, 0.6),
        return_id=False, predicator=None):
    """
    Generate a plot from received instance of World and show it on IPython notebook.
    This method returns the instance of dict that indicates color setting
    for each speices. You can use the dict as the parameter of plot_world,
    in order to use the same colors in another plot.

    Parameters
    ----------
    world : World or str
        World or a HDF5 filename to render.
    radius : float, default None
        If this value is set, all particles in the world will be rendered
        as if their radius are the same.
    width : float, default 350
        Width of the plotting area.
    height : float, default 350
        Height of the plotting area.
    config : dict, default {}
        Dict for configure default colors. Its values are colors unique
        to each speices. The dictionary will be updated during this plot.
        Colors included in config dict will never be used for other speices.
    species_list : array of string, default None
        If set, plot_world will not search the list of species.
    max_count : Integer, default 1000
        The maximum number of particles to show for each species.
    debug : array of dict, default []
        *** EXPERIMENTAL IMPRIMENTATION ***
        Example:
        >> [{'type': 'box', 'x': 10, 'y': 10, 'z': 10, 'options': {'width': 1, 'height': 1}}]
        type: 'box', 'plane', 'sphere', and 'cylinder'
        x, y, z: float
        options:
            box: width, height, depth
            plane: width, height
            sphere: radius
            cylinder: radius, height
    camera_position : tuple, default (-22, 23, 32)
    camera_rotaiton : tuple, default (-0.6, 0.5, 0.6)
        Initial position and rotation of camera.
    return_id : bool, default False
        If True, return a model id, which is required for `to_png` function.

    """
    config = config or {}

    from IPython.core.display import display, HTML
    from .simulation import load_world

    if isinstance(world, str):
        world = load_world(world)

    species = __parse_world(world, radius, species_list, max_count, predicator)
    color_scale = default_color_scale(config=config)
    plots = []

    for species_info in species:
        plots.append({
            'type': 'Particles',
            'data': species_info['data'],
            'options': {
                'name': species_info['name'],
                'color': color_scale.get_color(species_info['name']),
                'size': species_info['size']
            }
        })

    if debug is not None:
        data = {'type': [], 'x': [], 'y': [], 'z': [], 'options': []}
        for obj in debug:
            for k, v in obj.items():
                data[k].append(v)

        plots.append({
            'type': 'DebugObject',
            'data': data,
            'options': {}
        })

    model = {
        'plots': plots,
        'options': {
            'world_width': width,
            'world_height': height,
            'range': get_range_of_world(world),
            'autorange': False,
            'grid': grid,
            'save_image': True
            # 'save_image': False
        }
    }

    if wireframe:
        model['options']['space_mode'] = 'wireframe'

    model_id = '"viz' + str(uuid.uuid4()) + '"'
    display(HTML(generate_html(
        {'model': json.dumps(model), 'model_id': model_id,
        'px': camera_position[0], 'py': camera_position[1], 'pz': camera_position[2],
        'rx': camera_rotation[0], 'ry': camera_rotation[1], 'rz': camera_rotation[2]},
        'templates/particles.tmpl')))

    if return_id:
        return model_id

def to_png(plot_id):
    from IPython.display import display, HTML
    my_uuid = "\"png" + str(uuid.uuid4()) + "\""

    js = """
<script>
 function searchCell(uuid){
   var n = IPython.notebook.ncells();
   for(var i=0; i<n; i++){
     var cell = IPython.notebook.get_cell(i);
     if(typeof cell.output_area != "undefined"){
       var outputs = cell.output_area.outputs.filter(function(out){
console.log("Hi!");
         var html = out.data["text/html"];
         if(typeof html == "undefined")return false;
         if(html.includes(uuid))return true;
         return false;
       });
       if(outputs.length>0)return cell;
     }
   }
   return null;
 }

 var vis_id = %s;
 var my_uuid = %s;
 var vis_div = d3.select("#" + vis_id);
 var my_div =  d3.select("#" + my_uuid);

 var canvas = vis_div.select("canvas").node();
 var context = canvas.getContext("experimental-webgl", {preserveDrawingBuffer: true});
 var uri = canvas.toDataURL('image/png');

 my_div.append("img").attr("src", uri);

 window.setTimeout(function(){
 if(typeof window.IPython != "undefined"){
   try{
     var html = my_div.node().outerHTML;
     var cell = searchCell(my_uuid);
     if(cell == null)throw new Error("The cell whose id is " + my_uuid + " not found.");
     cell.output_area.outputs[0].data["text/html"] = html;
   }
   catch(e){
     console.warn("Maybe the front-end API of Jupyter has changed. message:" + e.message);
   }
 }
}, 0);
 
</script>
<div id=%s></div>
    """%(plot_id, my_uuid, my_uuid)
    display(HTML(js))


def plot_dense_array(
        arr, length=256, ranges=None, colors=("#a6cee3", "#fb9a99"), grid=False, camera_position=(-22, 23, 32), camera_rotation=(-0.6, 0.5, 0.6)):
    """
    Volume renderer

    Parameters
    ----------
    arr : list of numpy.array
        i.e. [array([[1,2,3], [2,3,4]]), array([[1,2,3]])]
    ranges : list of tuple
        ranges for x, y, and z axis
        i.e. [(-100, 100), (-100, 100), (-100, 100)]
    colors : list of string
        colors for species
    length : int
        length of the texture
        256 or 64
    camera_position : tuple, default (-22, 23, 32)
    camera_rotaiton : tuple, default (-0.6, 0.5, 0.6)
        Initial position and rotation of camera.

    """
    import numpy
    from PIL import Image
    from base64 import b64encode
    from tempfile import TemporaryFile
    from math import sqrt
    from IPython.core.display import display, HTML
    from functools import reduce

    # unfold 3d box into 2d grid
    def unfold(arr, dtype=None):
        dtype = arr.dtype if dtype is None else dtype
        i = sqrt(arr.shape[2])
        f_per_row, f_per_column = i, i
        # single channel (luminance)
        try:
            depth, height, width = arr.shape[:]
            arr = arr.reshape((depth*height, width))
            new_arr = numpy.empty((height*f_per_column, width*f_per_row), dtype=dtype)
        # multi channel (RGB)
        except ValueError:
            depth, height, width, channel = arr.shape
            arr = arr.reshape((depth*height, width, channel))
            new_arr = numpy.empty((height*f_per_column, width*f_per_row, channel), dtype=dtype)
        for h in range(0, int(f_per_column)):
            for w in range(0, int(f_per_row)):
                val = arr[(f_per_row*h+w)*height : (f_per_row*h+w+1)*height]
                new_arr[h*height : (h+1)*height, w*width : (w+1)*width] = val
        return new_arr

    def hist(arr, ranges, length, color):
        # create sample
        hist, bins = numpy.histogramdd(arr, bins=tuple([length]*3), range=tuple(ranges))
        # standardize value
        colors = [int(color[1:][i*2:(i+1)*2], 16) for i in range(0, 3)]
        len1d = reduce(lambda val, memo: memo*val, hist.shape, 1)
        arr = [((val/numpy.max(hist))*(hist.copy())).reshape(len1d) for val in colors]
        # add blue and green
        return numpy.array(arr, dtype=numpy.int8).transpose().reshape(tuple(list(hist.shape) + [3]))
    ranges = ranges if ranges is not None else [(numpy.min(a), numpy.max(a)) for a in numpy.array(arr).reshape((sum(map(len, arr)), 3)).transpose()]

    hist_arr = [hist(a, ranges, length, colors[i]) for i, a in enumerate(arr)]
    compressed = reduce(lambda p, n: p+n, hist_arr)

    img = Image.fromarray(unfold(compressed), "RGB")
    fp = TemporaryFile("r+b")
    img.save(fp, "PNG")
    fp.seek(0)
    encoded_url = "data:image/png;base64," + b64encode(fp.read())

    model = {
        'plots': [{
            'type': 'Volume',
            'data': encoded_url,
            'options': {
                'name': "",
                'width': length,
                'height': length,
                'depth': length,
                'f_per_row': sqrt(length),
                'f_per_column': sqrt(length)
            }
        }],
        'options': {
            'grid': grid,
            'save_image': True
        }
    }

    model_id = '"viz' + str(uuid.uuid4()) + '"'
    display(HTML(generate_html(
        {'model': json.dumps(model), 'model_id': model_id,
        'px': camera_position[0], 'py': camera_position[1], 'pz': camera_position[2],
        'rx': camera_rotation[0], 'ry': camera_rotation[1], 'rz': camera_rotation[2]},
        'templates/particles.tmpl')))

def generate_html(keywords, tmpl_path, package_name='ecell4.util'):
    """
    Generate static html file from JSON model and its own id.

    Parameters
    ----------
    model : dict
        JSON model from which ecell4.viz generates a plot.
    model_id : string
        Unique id for the plot.

    Returns
    -------
    html :
        A HTML object
    """
    from jinja2 import Template

    import pkgutil
    template = Template(pkgutil.get_data(package_name, tmpl_path).decode())
    # path = os.path.abspath(os.path.dirname(__file__)) + tmpl_path
    # template = Template(open(path).read())
    html = template.render(**keywords)
    return html


def plot_trajectory_with_elegans(
        obs, width=350, height=350, config=None, grid=True, wireframe=False,
        max_count=10, camera_position=(-22, 23, 32), camera_rotation=(-0.6, 0.5, 0.6),
        plot_range=None):
    """
    Generate a plot from received instance of TrajectoryObserver and show it
    on IPython notebook.

    Parameters
    ----------
    obs : TrajectoryObserver
        TrajectoryObserver to render.
    width : float, default 350
        Width of the plotting area.
    height : float, default 350
        Height of the plotting area.
    config : dict, default {}
        Dict for configure default colors. Its values are colors unique
        to each particle. The dictionary will be updated during this plot.
        Colors included in config dict will never be used for other particles.
    camera_position : tuple, default (-30, 31, 42)
    camera_rotaiton : tuple, default (-0.6, 0.5, 0.6)
        Initial position and rotation of camera.
    plot_range : tuple, default None
        Range for plotting. A triplet of pairs suggesting (rangex, rangey, rangez).
        If None, the minimum volume containing all the trajectories is used.

    """
    config = config or {}

    from IPython.core.display import display, HTML

    color_scale = default_color_scale(config=config)
    plots = []

    xmin, xmax, ymin, ymax, zmin, zmax = None, None, None, None, None, None

    data = obs.data()
    if max_count is not None and len(data) > max_count:
        data = random.sample(data, max_count)

    for i, y in enumerate(data):
        xarr, yarr, zarr = [], [], []
        for pos in y:
            xarr.append(pos[0])
            yarr.append(pos[1])
            zarr.append(pos[2])

        if xmin is None:
            if len(y) > 0:
                xmin, xmax = min(xarr), max(xarr)
                ymin, ymax = min(yarr), max(yarr)
                zmin, zmax = min(zarr), max(zarr)
        else:
            xmin, xmax = min([xmin] + xarr), max([xmax] + xarr)
            ymin, ymax = min([ymin] + yarr), max([ymax] + yarr)
            zmin, zmax = min([zmin] + zarr), max([zmax] + zarr)

        name = str(i + 1)
        c = color_scale.get_color(name)
        plots.append({
            'type': 'Line',
            'data': {'x': xarr, 'y': yarr, 'z': zarr},
            'options': {
                'name': name,
                'thickness': 2,  # XXX: 'thikness' doesn't work on Windows
                'colors': [c, c]}
        })

    if plot_range is None:
        if xmin is None:
            xmin, xmax, ymin, ymax, zmin, zmax = 0, 1, 0, 1, 0, 1

        max_length = max(xmax - xmin, ymax - ymin, zmax - zmin)
        rangex = [(xmin + xmax - max_length) * 0.5,
                  (xmin + xmax + max_length) * 0.5]
        rangey = [(ymin + ymax - max_length) * 0.5,
                  (ymin + ymax + max_length) * 0.5]
        rangez = [(zmin + zmax - max_length) * 0.5,
                  (zmin + zmax + max_length) * 0.5]
        wrange = {'x': rangex, 'y': rangey, 'z': rangez}
    else:
        wrange = get_range_of_trajectories(None, plot_range)

    model = {
        'plots': plots,
        'options': {
            'world_width': width,
            'world_height': height,
            'range': wrange,
            'autorange': False,
            'grid': grid,
            'save_image': True
        }
    }

    if wireframe:
        model['options']['space_mode'] = 'wireframe'

    model_id = '"viz' + str(uuid.uuid4()) + '"'
    display(HTML(generate_html(
        {'model': json.dumps(model), 'model_id': model_id,
        'px': camera_position[0], 'py': camera_position[1], 'pz': camera_position[2],
        'rx': camera_rotation[0], 'ry': camera_rotation[1], 'rz': camera_rotation[2]},
        'templates/particles.tmpl')))

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

def anim_to_html(anim, filename=None, fps=6, crf=10, bitrate='1M'):
    VIDEO_TAG = """<video controls>
     <source src="data:video/x-webm;base64,{0}" type="video/webm">
     Your browser does not support the video tag.
    </video>"""
    import base64

    if not hasattr(anim, '_encoded_video'):
        if filename is None:
            f = NamedTemporaryFile(suffix='.webm', delete=False)
            filename = f.name
            f.close()
            # anim.save(filename, fps=fps, extra_args=['-vcodec', 'libvpx'])
            anim.save(filename, fps=fps, codec='libvpx', extra_args=['-auto-alt-ref', '0', '-crf', str(crf), '-b:v', bitrate])
            # anim.save(filename, writer='mencoder', fps=fps, extra_args=['-lavcopts', 'vcodec=libvpx'])
            video = open(filename, "rb").read()
            os.remove(filename)
            # with NamedTemporaryFile(suffix='.webm') as f:
            #     anim.save(f.name, fps=fps, extra_args=['-vcodec', 'libvpx'])
            #     video = open(f.name, "rb").read()
        else:
            with open(filename, 'w') as f:
                anim.save(f.name, fps=fps, extra_args=['-vcodec', 'libvpx'])
                video = open(f.name, "rb").read()
        # anim._encoded_video = video.encode("base64")
        anim._encoded_video = base64.encodestring(video).decode('utf-8')
    return VIDEO_TAG.format(anim._encoded_video)

def display_anim(ani, output=None, fps=6, crf=10, bitrate='1M'):
    if output is None:
        from IPython.display import display, HTML
        display(HTML(anim_to_html(ani, output, fps=fps, crf=crf, bitrate=bitrate)))
    elif os.path.splitext(output.lower())[1] == '.webm':
        ani.save(output, fps=fps, codec='libvpx', extra_args=['-auto-alt-ref', '0', '-crf', str(crf), '-b:v', bitrate])
    elif os.path.splitext(output.lower())[1] == '.mp4':
        ani.save(output, fps=fps, codec='mpeg4', extra_args=['-crf', str(crf), '-b:v', bitrate])
    else:
        raise ValueError(
            "An output filename is only accepted with extension '.webm' or 'mp4'.")

plot_movie_of_trajectory = plot_movie_of_trajectory_with_matplotlib  # default

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
    from .simulation import load_world
    import os.path

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

def plot_world_with_plotly(world, species_list=None, max_count=1000):
    """
    Plot a World on IPython Notebook
    """
    if isinstance(world, str):
        from .simulation import load_world
        world = load_world(world)

    if species_list is None:
        species_list = [sp.serial() for sp in world.list_species()]
        species_list.sort()

    import random
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

def display_pdb(entity, width=400, height=400):
    from IPython.display import display, IFrame
    import ecell4.datasource.pdb as pdb
    entity_id = pdb.PDBDataSource.parse_entity(entity)
    if entity is None:
        raise ValueError('An invalid entity [{}] was given.'.format(repr(entity)))
    display(IFrame("https://gjbekker.github.io/molmil/#molmil.loadPDB('{}');".format(entity_id), width, height))
