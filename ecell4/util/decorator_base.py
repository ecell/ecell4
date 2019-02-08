import copy
import types
import warnings
import functools
import inspect
import re
import traceback

from . import parseobj


class Callback(object):
    """callback before the operations"""

    def __init__(self):
        pass

    def set(self):
        pass

    def get(self):
        return None

    def notify_unary_operations(self, obj):
        pass

    def notify_bitwise_operations(self, obj):
        pass

    def notify_comparisons(self, obj):
        pass

class JustParseCallback(Callback):

    def __init__(self):
        Callback.__init__(self)

        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_comparisons(self, obj):
        if isinstance(obj, parseobj.NeExp):
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)
        self.comparisons.append(obj)

class DispatcherCallback(Callback):

    def __init__(self, children):
        Callback.__init__(self)
        self.__children = children

    def set(self):
        for child in self.__children:
            child.set()

    def get(self):
        return [child.get() for child in self.__children]

    def notify_unary_operations(self, obj):
        for child in self.__children:
            child.notify_unary_operations(obj)

    def notify_bitwise_operations(self, obj):
        for child in self.__children:
            child.notify_bitwise_operations(obj)

    def notify_comparisons(self, obj):
        for child in self.__children:
            child.notify_comparisons(obj)

def callback_dispatcher(callback_classes):
    def wrapper():
        return DispatcherCallback([cls() for cls in callback_classes])
    return wrapper

class EchoCallback(Callback):

    def __init__(self):
        Callback.__init__(self)

    def set(self):
        pass

    def get(self):
        return None

    def notify_unary_operations(self, obj):
        stack = traceback.extract_stack()
        _, _, funcname, _ = stack[-2]
        self.say(funcname, obj)

    def notify_bitwise_operations(self, obj):
        stack = traceback.extract_stack()
        _, _, funcname, _ = stack[-2]
        self.say(funcname, obj)

    def notify_comparisons(self, obj):
        stack = traceback.extract_stack()
        _, _, funcname, _ = stack[-2]
        self.say(funcname, obj)

    def say(self, method, *args, **kwargs):
        print('{} [{}] method={}, args={}, kwargs={}'.format(self.__class__.__name__, id(self), method, args, kwargs))

class TransparentCallback(Callback):

    def __init__(self):
        Callback.__init__(self)

    def set(self):
        pass

    def get(self):
        return None

    def notify_unary_operations(self, obj):
        pass

    def notify_bitwise_operations(self, obj):
        pass

    def notify_comparisons(self, obj):
        pass

def keys_from_builtins(vardict):
    b = vardict['__builtins__']
    if isinstance(b, types.ModuleType):
        return dir(b)
    else:
        return dict(b).keys()

class ParseDecorator:

    def __init__(self, callback_class, func=None):
        self.__callback_class = callback_class

        self.__callback = None
        self.__newvars = {}

        if func is not None:
            self.__func = func
            functools.update_wrapper(self, func)
        else:
            self.__func = lambda *args, **kwargs: []

    def set_callback(self, callback=None):
        if callback is None:
            self.__callback = self.__callback_class()
        else:
            self.__callback = callback

    def wrapper(self, *args, **kwargs):
        try:
            vardict = copy.copy(self.__func.func_globals)
            func_code = self.__func.func_code
            name = self.__func.func_name
            defaults = self.__func.func_defaults
        except AttributeError:
            vardict = copy.copy(self.__func.__globals__)
            func_code = self.__func.__code__
            name = self.__func.__name__
            defaults = self.__func.__defaults__

        ignores = ("_", "__", "___", "_i", "_ii", "_iii",
            "_i1", "_i2", "_i3", "_dh", "_sh", "_oh")
        for ignore in ignores:
            if ignore in vardict.keys():
                del vardict[ignore]

        if "_eval" not in vardict.keys():
            vardict["_eval"] = self.__evaluate
        if "_callback" not in vardict.keys():
            vardict["_callback"] = self.__callback
        else:
            pass  #XXX: raise an exception?

        for k in func_code.co_names:
            if (not k in vardict.keys()
                and not k in keys_from_builtins(vardict)): # is this enough?
                vardict[k] = parseobj.AnyCallable(self.__callback, k)
        g = types.FunctionType(func_code, vardict, name=name, argdefs=defaults)
        with warnings.catch_warnings():
            # warnings.simplefilter("always")
            g(*args, **kwargs)
        return self.__callback.get()

    def __call__(self, *args, **kwargs):
        calling_frame = inspect.currentframe().f_back
        if ('_callback' in calling_frame.f_globals.keys()
            and isinstance(calling_frame.f_globals['_callback'], self.__callback_class)):
            self.set_callback(calling_frame.f_globals["_callback"])
        else:
            self.set_callback()
        retval = self.wrapper(*args, **kwargs)
        self.__callback = None
        return retval

    def __enter__(self):
        # print "ParseDecorator#__enter__"
        self.set_callback()
        calling_frame = inspect.currentframe().f_back
        vardict = copy.copy(calling_frame.f_globals)
        ignores = ("_", "__", "___", "_i", "_ii", "_iii",
            "_i1", "_i2", "_i3", "_dh", "_sh", "_oh")

        if "_eval" not in calling_frame.f_globals.keys():
            calling_frame.f_globals["_eval"] = self.__evaluate
            self.__newvars["_eval"] = None
        if "_callback" not in calling_frame.f_globals.keys():
            calling_frame.f_globals["_callback"] = self.__callback
            self.__newvars["_callback"] = None
        else:
            pass  #XXX: raise an exception?

        for k in calling_frame.f_code.co_names:
            if k in ('_eval', '_callback'):
                pass
            elif k in ignores:
                # print "WARNING: '%s' was overridden." % k
                calling_frame.f_globals[k] = parseobj.AnyCallable(self.__callback, k)
                self.__newvars[k] = vardict.get(k)
            elif (not k in vardict.keys()
                and not k in keys_from_builtins(vardict)):
                # print "WARNING: '%s' is undefined." % k
                calling_frame.f_globals[k] = parseobj.AnyCallable(self.__callback, k)
                self.__newvars[k] = None
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # print "ParseDecorator#__exit__", exc_type, exc_value, traceback
        if self.__callback is not None:
            if exc_type is None:
                # print self.__callback.get()
                self.__callback.set()

            self.__callback = None
            calling_frame = inspect.currentframe().f_back
            for k, v in self.__newvars.items():
                if v is None:
                    del calling_frame.f_globals[k]
                    # print "WARNING: '%s' was removed." % k
                else:
                    calling_frame.f_globals[k] = v
                    # print "WARNING: '%s' was recovered to be '%s'." % (k, v)

    def eval(self, expr, params=None):
        params = params or {}

        class AnyCallableLocals:

            def __init__(self, callback, locals):
                self.callback = callback
                self.locals = locals

            def __getitem__(self, key):
                if key in self.locals.keys():
                    return self.locals[key]
                return parseobj.AnyCallable(self.callback, key)

        l = locals()
        l.update(params)

        return eval(expr, globals(), AnyCallableLocals(self.__callback, l))

    def __evaluate(self, expr, params=None):
        return self.eval(expr, params)

just_parse = functools.partial(ParseDecorator, JustParseCallback)
