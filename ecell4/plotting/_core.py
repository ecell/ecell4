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

def eval_key(observable, targets, data, err=None):
    assert len(data) == len(targets) + 1
    assert err is None or len(err) == len(data)

    if observable in targets:
        idx = targets.index(observable) + 1
        data_ = data[idx]
        err_ = err[idx] if err is not None else None
    else:
        err_ = None
        import numpy
        locals_ = {
            'sqrt': numpy.sqrt, 'cbrt': numpy.cbrt, 'square': numpy.square, 'abs': numpy.fabs,
            'max': numpy.fmax, 'min': numpy.fmin, 'exp': numpy.exp, 'log': numpy.log,
            'pow': numpy.power, 'sin': numpy.sin, 'cos': numpy.cos, 'tan': numpy.tan,
            }
        try:
            data_ = eval(observable, {"__builtins__": None},
                    dict(locals_, **{serial: data[idx + 1] for idx, serial in enumerate(targets)}))
            # if err is not None:
            #     err_ = eval(observable, {"__builtins__": None},
            #             dict(locals_, **{serial: err[idx + 1] for idx, serial in enumerate(targets)}))
        except:
            raise ValueError("Unknown variable [{0}] was given.".format(observable))
    return (data_, err_)
