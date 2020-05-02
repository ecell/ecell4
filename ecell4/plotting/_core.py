import os
import os.path
from tempfile import NamedTemporaryFile

def get_range_of_world(world, scale=1.0):
    edge_lengths = world.edge_lengths() * scale
    max_length = max(tuple(edge_lengths))

    rangex = [(edge_lengths[0] - max_length) * 0.5,
              (edge_lengths[0] + max_length) * 0.5]
    rangey = [(edge_lengths[1] - max_length) * 0.5,
              (edge_lengths[1] + max_length) * 0.5]
    rangez = [(edge_lengths[2] - max_length) * 0.5,
              (edge_lengths[2] + max_length) * 0.5]

    return {'x': rangex, 'y': rangey, 'z': rangez}

def get_range_of_trajectories(data, plot_range=None):
    from ecell4_base.core import Real3

    if plot_range is None:
        if len(data) == 0:
            xmin, xmax, ymin, ymax, zmin, zmax = 0, 1, 0, 1, 0, 1
        else:
            xmin, xmax, ymin, ymax, zmin, zmax = None, None, None, None, None, None

            for i, traj in enumerate(data):
                xarr, yarr, zarr = [], [], []
                for pos in traj:
                    xarr.append(pos[0])
                    yarr.append(pos[1])
                    zarr.append(pos[2])

                if xmin is None:
                    if len(traj) > 0:
                        xmin, xmax = min(xarr), max(xarr)
                        ymin, ymax = min(yarr), max(yarr)
                        zmin, zmax = min(zarr), max(zarr)
                else:
                    xmin, xmax = min([xmin] + xarr), max([xmax] + xarr)
                    ymin, ymax = min([ymin] + yarr), max([ymax] + yarr)
                    zmin, zmax = min([zmin] + zarr), max([zmax] + zarr)

        max_length = max(xmax - xmin, ymax - ymin, zmax - zmin)
        rangex = [(xmin + xmax - max_length) * 0.5,
                  (xmin + xmax + max_length) * 0.5]
        rangey = [(ymin + ymax - max_length) * 0.5,
                  (ymin + ymax + max_length) * 0.5]
        rangez = [(zmin + zmax - max_length) * 0.5,
                  (zmin + zmax + max_length) * 0.5]

        return {'x': rangex, 'y': rangey, 'z': rangez}
    elif isinstance(plot_range, dict):
        return plot_range
    elif isinstance(plot_range, (list, tuple)):
        if len(plot_range) != 3:
            raise ValueError(
                'The size of plot_range [{}] must be 3.'.format(len(plot_range)))
        elif (isinstance(plot_range[0], (list, tuple)) and
                isinstance(plot_range[1], (list, tuple)) and
                isinstance(plot_range[2], (list, tuple))):
            return {'x': plot_range[0], 'y': plot_range[1], 'z': plot_range[2]}
        else:
            return {'x': (0, plot_range[0]),
                    'y': (0, plot_range[1]),
                    'z': (0, plot_range[2])}
    elif isinstance(plot_range, Real3):
        return {'x': (0, plot_range[0]),
                'y': (0, plot_range[1]),
                'z': (0, plot_range[2])}
    else:
        raise ValueError(
            'plot_range must be list, tuple or dict. [{}] was given.'.format(
                repr(plot_range)))

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

from ..util.model_parser import Visitor, dispatch, RATELAW_RESERVED_FUNCTIONS, RATELAW_RESERVED_CONSTANTS
import ecell4_base

class EvalVisitor(Visitor):

    def __init__(self, variables):
        self.variables = variables

    def visit_const(self, obj):
        if str(obj) not in RATELAW_RESERVED_CONSTANTS:
            raise ValueError()
        return RATELAW_RESERVED_CONSTANTS[str(obj)]

    def visit_species(self, obj):
        pttrn = ecell4_base.core.Species(str(obj))
        return sum(ecell4_base.core.count_species_matches(pttrn, key) * value for key, value in self.variables.items())

    def visit_func(self, obj, *args):
        print(obj)
        assert obj._elems[0].name in RATELAW_RESERVED_FUNCTIONS
        return RATELAW_RESERVED_FUNCTIONS[obj._elems[0].name](*args)

    def visit_expression(self, obj, *args):
        return obj._execute(*args)

def eval_key(observable, targets, data, err=None):
    assert len(data) == len(targets) + 1
    assert err is None or len(err) == len(data)

    from ecell4.util.decorator_base import just_parse
    parsed = just_parse().eval(observable)

    import numpy
    variables = dict(zip(targets, data[1: ]))

    #XXX: The second return value is for standard deviations (err)
    return (dispatch(parsed, EvalVisitor(variables)), None)
