from __future__ import print_function

import os
import os.path
import logging
import tempfile
import pickle
import inspect
import textwrap
import re
import types
import itertools
import binascii
import multiprocessing
import copy
import csv

from . import sge
from . import slurm


def run_ensemble(target, jobs, repeat=1, method=None, **kwargs):
    """
    Evaluate the given function with each set of arguments, and return a list of results.

    Parameters
    ----------
    target : function
        A function to be evaluated. The function must accepts three arguments,
        which are a list of arguments given as `jobs`, a job and task id (int).
    jobs : list
        A list of arguments passed to the function.
    repeat : int, optional
        A number of tasks. Repeat the evaluation `n` times for each job.
        1 for default.
    method : str, optional
        The way for running multiple jobs.
        Choose one from 'serial', 'multiprocessing', 'sge', 'slurm', 'azure'.
        Default is None, which works as 'serial'.

    Returns
    -------
    results : list
        A list of results. Each element is a list containing `repeat` results.

    See Also
    --------
    ecell4.extra.ensemble.run_serial
    ecell4.extra.ensemble.run_sge
    ecell4.extra.ensemble.run_slurm
    ecell4.extra.ensemble.run_multiprocessing
    ecell4.extra.ensemble.run_azure

    """
    config = None
    config_filename = 'config.yaml'
    if os.path.isfile(config_filename):
        import yaml
        try:
            from yaml import CLoader as Loader
        except ImportError:
            from yaml import Loader
        with open(config_filename) as f:
            config_ = yaml.load(f.read(), Loader=Loader)
        if 'ensemble' in config_:
            config = config_['ensemble']

    if method is None:
        if config is None or 'method' not in config:
            method = 'serial'  # default
        else:
            method = config['method']
    method = method.lower()

    if config is not None and method in config:
        kwargs_ = config[method]
        kwargs = dict(kwargs_, **kwargs)

    if method == "serial":
        return run_serial(target, jobs, n=repeat, **kwargs)
    elif method == "sge":
        return run_sge(target, jobs, n=repeat, **kwargs)
    elif method == "slurm":
        return run_slurm(target, jobs, n=repeat, **kwargs)
    elif method == "multiprocessing":
        return run_multiprocessing(target, jobs, n=repeat, **kwargs)
    elif method == "azure":
        return run_azure(target, jobs, n=repeat, **kwargs)

    raise ValueError(
        'Argument "method" must be either "serial", "multiprocessing", "slurm", "sge" or "azure".')

def run_serial(target, jobs, n=1, **kwargs):
    """
    Evaluate the given function with each set of arguments, and return a list of results.
    This function does in series.

    Parameters
    ----------
    target : function
        A function to be evaluated. The function must accepts three arguments,
        which are a list of arguments given as `jobs`, a job and task id (int).
    jobs : list
        A list of arguments passed to the function.
    n : int, optional
        A number of tasks. Repeat the evaluation `n` times for each job.
        1 for default.

    Returns
    -------
    results : list
        A list of results. Each element is a list containing `n` results.

    Examples
    --------
    >>> jobs = ((1, 'spam'), (2, 'ham'), (3, 'eggs'))

    >>> target = lambda args, job_id, task_id: (args[1] * args[0])
    >>> run_serial(target, jobs)
    [['spam'], ['hamham'], ['eggseggseggs']]

    >>> target = lambda args, job_id, task_id: "{:d} {}".format(task_id, args[1] * args[0])
    >>> run_serial(target, jobs, n=2)
    [['1 spam', '2 spam'], ['1 hamham', '2 hamham'], ['1 eggseggseggs', '2 eggseggseggs']]

    >>> seeds = genseeds(3)
    >>> def target(arg, job_id, task_id):
    ...     from ecell4.extra.ensemble import getseed
    ...     return getseed(arg, task_id)
    >>> run_serial(target, (seeds, ), n=3)  # doctest: +SKIP
    [[127152315, 2028054913, 253611282]]

    See Also
    --------
    ecell4.extra.ensemble.run_serial
    ecell4.extra.ensemble.run_sge
    ecell4.extra.ensemble.run_slurm
    ecell4.extra.ensemble.run_multiprocessing
    ecell4.extra.ensemble.run_azure

    """
    return [[target(copy.deepcopy(job), i + 1, j + 1) for j in range(n)] for i, job in enumerate(jobs)]

def run_multiprocessing(target, jobs, n=1, nproc=None, **kwargs):
    """
    Evaluate the given function with each set of arguments, and return a list of results.
    This function does in parallel by using `multiprocessing`.

    Parameters
    ----------
    target : function
        A function to be evaluated. The function must accepts three arguments,
        which are a list of arguments given as `jobs`, a job and task id (int).
    jobs : list
        A list of arguments passed to the function.
        All the argument must be picklable.
    n : int, optional
        A number of tasks. Repeat the evaluation `n` times for each job.
        1 for default.
    nproc : int, optional
        A number of cores available once.
        If nothing is given, all available cores are used.

    Returns
    -------
    results : list
        A list of results. Each element is a list containing `n` results.

    Examples
    --------
    >>> jobs = ((1, 'spam'), (2, 'ham'), (3, 'eggs'))

    >>> target = lambda args, job_id, task_id: (args[1] * args[0])
    >>> run_multiprocessing(target, jobs, nproc=2)
    [['spam'], ['hamham'], ['eggseggseggs']]

    >>> target = lambda args, job_id, task_id: "{:d} {}".format(task_id, args[1] * args[0])
    >>> run_multiprocessing(target, jobs, n=2, nproc=2)
    [['1 spam', '2 spam'], ['1 hamham', '2 hamham'], ['1 eggseggseggs', '2 eggseggseggs']]

    See Also
    --------
    ecell4.extra.ensemble.run_serial
    ecell4.extra.ensemble.run_sge
    ecell4.extra.ensemble.run_slurm
    ecell4.extra.ensemble.run_multiprocessing
    ecell4.extra.ensemble.run_azure

    """
    def consumer(f, q_in, q_out):
        while True:
            val = q_in.get()
            if val is None:
                q_in.task_done()
                break
            i, x = val
            res = (i, f(*x))
            q_in.task_done()
            q_out.put(res)

    def mulpmap(f, X, nproc):
        nproc = nproc or multiprocessing.cpu_count()

        q_in = multiprocessing.JoinableQueue()
        q_out = multiprocessing.Queue()
        workers = [multiprocessing.Process(target=consumer, args=(f, q_in, q_out), daemon=True) for _ in range(nproc)]
        sent = [q_in.put((i, x)) for i, x in enumerate(X)]
        num_tasks = len(sent)
        [q_in.put(None) for _ in range(nproc)]  #XXX: poison pill
        [w.start() for w in workers]
        # [w.join() for w in workers]
        q_in.join()
        res = [q_out.get() for _ in range(num_tasks)]
        return [x for (_, x) in sorted(res)]

    res = mulpmap(
        target, ((job, i + 1, j + 1) for (i, job), j in itertools.product(enumerate(jobs), range(n))), nproc)
    return [res[i: i + n] for i in range(0, len(res), n)]

def run_sge(target, jobs, n=1, nproc=None, path='.', delete=True, wait=True, environ=None, modules=(), **kwargs):
    """
    Evaluate the given function with each set of arguments, and return a list of results.
    This function does in parallel on the Sun Grid Engine einvironment.

    Parameters
    ----------
    target : function
        A function to be evaluated. The function must accepts three arguments,
        which are a list of arguments given as `jobs`, a job and task id (int).
        This function can not be a lambda.
    jobs : list
        A list of arguments passed to the function.
        All the argument must be picklable.
    n : int, optional
        A number of tasks. Repeat the evaluation `n` times for each job.
        1 for default.
    nproc : int, optional
        A number of cores available once.
        If nothing is given, it runs with no limit.
    path : str, optional
        A path for temporary files to be saved. The path is created if not exists.
        The current directory is used as its default.
    delete : bool, optional
        Whether it removes temporary files after the successful execution.
        True for default.
    wait : bool, optional
        Whether it waits until all jobs are finished. If False, it just submits jobs.
        True for default.
    environ : dict, optional
        An environment variables used when running jobs.
        "PYTHONPATH" and "LD_LIBRARY_PATH" is inherited when no `environ` is given.
    modules : list, optional
        A list of module names imported before evaluating the given function.
        The modules are loaded as: `from [module] import *`.

    Returns
    -------
    results : list
        A list of results. Each element is a list containing `n` results.

    Examples
    --------
    >>> jobs = ((1, 'spam'), (2, 'ham'), (3, 'eggs'))

    >>> def target(args, job_id, task_id):
    ...     return (args[1] * args[0])
    ...
    >>> run_sge(target, jobs, nproc=2, path='.tmp')
    [['spam'], ['hamham'], ['eggseggseggs']]

    >>> def target(args, job_id, task_id):
    ...     return "{:d} {}".format(task_id, args[1] * args[0])
    ...
    >>> run_sge(target, jobs, n=2, nproc=2, path='.tmp')
    [['1 spam', '2 spam'], ['1 hamham', '2 hamham'], ['1 eggseggseggs', '2 eggseggseggs']]

    See Also
    --------
    ecell4.extra.ensemble.run_serial
    ecell4.extra.ensemble.run_sge
    ecell4.extra.ensemble.run_slurm
    ecell4.extra.ensemble.run_multiprocessing
    ecell4.extra.ensemble.run_azure

    """
    logging.basicConfig(level=logging.DEBUG)

    if isinstance(target, types.LambdaType) and target.__name__ == "<lambda>":
        raise RuntimeError("A lambda function is not accepted")

    src = textwrap.dedent(inspect.getsource(target)).replace(r'"', r'\"')
    if re.match('[\s\t]+', src.split('\n')[0]) is not None:
        raise RuntimeError(
            "Wrong indentation was found in the source translated")

    if not os.path.isdir(path):
        os.makedirs(path)  #XXX: MYOB

    if environ is None:
        environ = {}
        keys = ("LD_LIBRARY_PATH", "PYTHONPATH")
        for key in keys:
            if key in os.environ.keys():
                environ[key] = os.environ[key]
        if "PYTHONPATH" in environ.keys() and environ["PYTHONPATH"].strip() != "":
            environ["PYTHONPATH"] = "{}:{}".format(os.getcwd(), environ["PYTHONPATH"])
        else:
            environ["PYTHONPATH"] = os.getcwd()

    cmds = []
    pickleins = []
    pickleouts = []
    scripts = []
    for i, job in enumerate(jobs):
        (fd, picklein) = tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)
        with os.fdopen(fd, 'wb') as fout:
            pickle.dump(job, fout)
        pickleins.append(picklein)

        pickleouts.append([])
        for j in range(n):
            fd, pickleout = tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)
            os.close(fd)
            pickleouts[-1].append(pickleout)
        # pickleouts.append(
        #     [tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)[1]
        #      for j in range(n)])

        code = 'import sys\n'
        code += 'import os\n'
        code += 'import pickle\n'
        code += 'with open(\'{}\', \'rb\') as fin:\n'.format(picklein)
        code += '    job = pickle.load(fin)\n'
        code += 'pass\n'
        for m in modules:
            code += "from {} import *\n".format(m)
        code += src
        code += '\ntid = int(os.environ[\'SGE_TASK_ID\'])'
        code += '\nretval = {:s}(job, {:d}, tid)'.format(target.__name__, i + 1)
        code += '\nfilenames = {:s}'.format(str(pickleouts[-1]))
        code += '\npickle.dump(retval, open(filenames[tid - 1], \'wb\'))\n'

        (fd, script) = tempfile.mkstemp(suffix='.py', prefix='sge-', dir=path, text=True)
        with os.fdopen(fd, 'w') as fout:
            fout.write(code)
        scripts.append(script)

        cmd = '#!/bin/bash\n'
        for key, value in environ.items():
            cmd += 'export {:s}={:s}\n'.format(key, value)
        cmd += 'python3 {}'.format(script)  #XXX: Use the same executer, python
        # cmd += 'python3 -c "\n'
        # cmd += 'import sys\n'
        # cmd += 'import os\n'
        # cmd += 'import pickle\n'
        # cmd += 'with open(sys.argv[1], \'rb\') as fin:\n'
        # cmd += '    job = pickle.load(fin)\n'
        # cmd += 'pass\n'
        # for m in modules:
        #     cmd += "from {} import *\n".format(m)
        # cmd += src
        # cmd += '\ntid = int(os.environ[\'SGE_TASK_ID\'])'
        # cmd += '\nretval = {:s}(job, {:d}, tid)'.format(target.__name__, i + 1)
        # cmd += '\nfilenames = {:s}'.format(str(pickleouts[-1]))
        # cmd += '\npickle.dump(retval, open(filenames[tid - 1], \'wb\'))'
        # cmd += '" {:s}\n'.format(picklein)
        cmds.append(cmd)

    if isinstance(wait, bool):
        sync = 0 if not wait else 10
    elif isinstance(wait, int):
        sync = wait
    else:
        raise ValueError("'wait' must be either 'int' or 'bool'.")

    jobids = sge.run(cmds, n=n, path=path, delete=delete, sync=sync, max_running_tasks=nproc, **kwargs)

    if not (sync > 0):
        return None

    for jobid, name in jobids:
        outputs = sge.collect(jobid, name, n=n, path=path, delete=delete)
        for output in outputs:
            print(output, end='')

    retval = [[pickle.load(open(pickleout, 'rb')) for pickleout in tasks]
              for tasks in pickleouts]

    if delete:
        for tmpname in itertools.chain(pickleins, scripts, *pickleouts):
            os.remove(tmpname)

    return retval

def run_slurm(target, jobs, n=1, nproc=None, path='.', delete=True, wait=True, environ=None, modules=(), **kwargs):
    """
    Evaluate the given function with each set of arguments, and return a list of results.
    This function does in parallel with Slurm Workload Manager.

    Parameters
    ----------
    target : function
        A function to be evaluated. The function must accepts three arguments,
        which are a list of arguments given as `jobs`, a job and task id (int).
        This function can not be a lambda.
    jobs : list
        A list of arguments passed to the function.
        All the argument must be picklable.
    n : int, optional
        A number of tasks. Repeat the evaluation `n` times for each job.
        1 for default.
    nproc : int, optional
        A number of cores available once.
        If nothing is given, it runs with no limit.
    path : str, optional
        A path for temporary files to be saved. The path is created if not exists.
        The current directory is used as its default.
    delete : bool, optional
        Whether it removes temporary files after the successful execution.
        True for default.
    wait : bool, optional
        Whether it waits until all jobs are finished. If False, it just submits jobs.
        True for default.
    environ : dict, optional
        An environment variables used when running jobs.
        "PYTHONPATH" and "LD_LIBRARY_PATH" is inherited when no `environ` is given.
    modules : list, optional
        A list of module names imported before evaluating the given function.
        The modules are loaded as: `from [module] import *`.

    Returns
    -------
    results : list
        A list of results. Each element is a list containing `n` results.

    Examples
    --------
    >>> jobs = ((1, 'spam'), (2, 'ham'), (3, 'eggs'))

    >>> def target(args, job_id, task_id):
    ...     return (args[1] * args[0])
    ...
    >>> run_slurm(target, jobs, nproc=2, path='.tmp')
    [['spam'], ['hamham'], ['eggseggseggs']]

    >>> def target(args, job_id, task_id):
    ...     return "{:d} {}".format(task_id, args[1] * args[0])
    ...
    >>> run_slurm(target, jobs, n=2, nproc=2, path='.tmp')
    [['1 spam', '2 spam'], ['1 hamham', '2 hamham'], ['1 eggseggseggs', '2 eggseggseggs']]

    See Also
    --------
    ecell4.extra.ensemble.run_serial
    ecell4.extra.ensemble.run_sge
    ecell4.extra.ensemble.run_slurm
    ecell4.extra.ensemble.run_multiprocessing
    ecell4.extra.ensemble.run_azure

    """
    logging.basicConfig(level=logging.DEBUG)

    if isinstance(target, types.LambdaType) and target.__name__ == "<lambda>":
        raise RuntimeError("A lambda function is not accepted")

    src = textwrap.dedent(inspect.getsource(target)).replace(r'"', r'\"')
    if re.match('[\s\t]+', src.split('\n')[0]) is not None:
        raise RuntimeError(
            "Wrong indentation was found in the source translated")

    if not os.path.isdir(path):
        os.makedirs(path)  #XXX: MYOB

    if environ is None:
        environ = {}
        keys = ("LD_LIBRARY_PATH", "PYTHONPATH")
        for key in keys:
            if key in os.environ.keys():
                environ[key] = os.environ[key]
        if "PYTHONPATH" in environ.keys() and environ["PYTHONPATH"].strip() != "":
            environ["PYTHONPATH"] = "{}:{}".format(os.getcwd(), environ["PYTHONPATH"])
        else:
            environ["PYTHONPATH"] = os.getcwd()

    cmds = []
    pickleins = []
    pickleouts = []
    scripts = []
    for i, job in enumerate(jobs):
        (fd, picklein) = tempfile.mkstemp(suffix='.pickle', prefix='slurm-', dir=path)
        with os.fdopen(fd, 'wb') as fout:
            pickle.dump(job, fout)
        pickleins.append(picklein)

        pickleouts.append([])
        for j in range(n):
            fd, pickleout = tempfile.mkstemp(suffix='.pickle', prefix='slurm-', dir=path)
            os.close(fd)
            pickleouts[-1].append(pickleout)
        # pickleouts.append(
        #     [tempfile.mkstemp(suffix='.pickle', prefix='slurm-', dir=path)[1]
        #      for j in range(n)])

        code = 'import sys\n'
        code += 'import os\n'
        code += 'import pickle\n'
        code += 'with open(\'{}\', \'rb\') as fin:\n'.format(picklein)
        code += '    job = pickle.load(fin)\n'
        code += 'pass\n'
        for m in modules:
            code += "from {} import *\n".format(m)
        code += src
        code += '\ntid = int(os.environ[\'SLURM_ARRAY_TASK_ID\'])'
        # code += '\ntid = int(os.environ[\'SGE_TASK_ID\'])'
        code += '\nretval = {:s}(job, {:d}, tid)'.format(target.__name__, i + 1)
        code += '\nfilenames = {:s}'.format(str(pickleouts[-1]))
        code += '\npickle.dump(retval, open(filenames[tid - 1], \'wb\'))\n'

        (fd, script) = tempfile.mkstemp(suffix='.py', prefix='slurm-', dir=path, text=True)
        with os.fdopen(fd, 'w') as fout:
            fout.write(code)
        scripts.append(script)

        cmd = '#!/bin/bash\n'
        for key, value in environ.items():
            cmd += 'export {:s}={:s}\n'.format(key, value)
        cmd += 'python3 {}'.format(script)  #XXX: Use the same executer, python
        # cmd += 'python3 -c "\n'
        # cmd += 'import sys\n'
        # cmd += 'import os\n'
        # cmd += 'import pickle\n'
        # cmd += 'with open(sys.argv[1], \'rb\') as fin:\n'
        # cmd += '    job = pickle.load(fin)\n'
        # cmd += 'pass\n'
        # for m in modules:
        #     cmd += "from {} import *\n".format(m)
        # cmd += src
        # cmd += '\ntid = int(os.environ[\'SGE_TASK_ID\'])'
        # cmd += '\nretval = {:s}(job, {:d}, tid)'.format(target.__name__, i + 1)
        # cmd += '\nfilenames = {:s}'.format(str(pickleouts[-1]))
        # cmd += '\npickle.dump(retval, open(filenames[tid - 1], \'wb\'))'
        # cmd += '" {:s}\n'.format(picklein)
        cmds.append(cmd)

    if isinstance(wait, bool):
        sync = 0 if not wait else 10
    elif isinstance(wait, int):
        sync = wait
    else:
        raise ValueError("'wait' must be either 'int' or 'bool'.")

    #XXX: nproc only limits the maximum count for 'each' job, but not for the whole jobs.
    jobids = slurm.run(cmds, n=n, path=path, delete=delete, sync=sync, max_running_tasks=nproc, **kwargs)

    if not (sync > 0):
        return None

    for jobid, name in jobids:
        outputs = slurm.collect(jobid, name, n=n, path=path, delete=delete)
        for output in outputs:
            print(output, end='')

    retval = [[pickle.load(open(pickleout, 'rb')) for pickleout in tasks]
              for tasks in pickleouts]

    if delete:
        for tmpname in itertools.chain(pickleins, scripts, *pickleouts):
            os.remove(tmpname)

    return retval

def run_azure(target, jobs, n=1, nproc=None, path='.', delete=True, config=None, **kwargs):
    """
    Evaluate the given function with each set of arguments, and return a list of results.
    This function does in parallel with Microsoft Azure Batch.

    This function is the work in progress.
    The argument `nproc` doesn't work yet.
    See `ecell4.extra.azure_batch.run_azure` for details.

    See Also
    --------
    ecell4.extra.ensemble.run_serial
    ecell4.extra.ensemble.run_sge
    ecell4.extra.ensemble.run_slurm
    ecell4.extra.ensemble.run_multiprocessing
    ecell4.extra.ensemble.run_azure
    ecell4.extra.azure_batch.run_azure

    """
    import ecell4.extra.azure_batch as azure_batch
    return azure_batch.run_azure(target, jobs, n, path, delete, config)

def genseeds(n):
    """
    Return a random number generator seed for ensemble_simulations.
    A seed for a single run is given by ``getseed(rngseed, i)``.

    Parameters
    ----------
    n : int
        A size of the seed.

    Returns
    -------
    rndseed : bytes
        A random number seed for multiple runs.

    """
    return binascii.hexlify(os.urandom(4 * n))

def getseed(myseed, i):
    """
    Return a single seed from a long seed given by `genseeds`.

    Parameters
    ----------
    myseed : bytes
        A long seed given by `genseeds(n)`.
    i : int
        An index less than n.

    Returns
    -------
    rndseed : int
        A seed (less than (2 ** 31))

    """
    rndseed = int(myseed[(i - 1) * 8: i * 8], 16)
    rndseed = rndseed % (2 ** 31)  #XXX: trancate the first bit
    return rndseed


if __name__ == "__main__":
    # def myrun(job, job_id=0, task_id=0):
    #     import ecell4_base
    #     print("Hi, I'm in local!")
    #     print("My job id is {:d}, and my task id is {:d}.".format(job_id, task_id))
    #     print("My job is {:s}.".format(str(job)))
    #     return job['x'] + job['y']

    # jobs = [{'x': i, 'y': i ** 2} for i in range(1, 4)]
    # print(run_serial(myrun, jobs, n=2))
    # print(run_multiprocessing(myrun, jobs, n=2))
    # # print(run_sge(myrun, jobs, n=2, delete=False))
    # print(run_sge(myrun, jobs, n=2))

    from ecell4 import *
    from ecell4.extra import ensemble

    with reaction_rules():
        A + B == C | (0.01, 0.3)

    ensemble.ensemble_simulations(
        10.0, {'C': 60}, solver='gillespie', return_type='matplotlib',
        n=30, method='multiprocessing')
