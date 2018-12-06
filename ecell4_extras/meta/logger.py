# coding: utf-8

import collections
from ecell4 import Species


class Logger:

    def __init__(self, coordinator, value):
        self.coordinator = coordinator

        self.__t = []
        self.__data = collections.OrderedDict()
        self.add(value)

    def add(self, value1, value2=None):
        if isinstance(value1, Species):
            if (value1, value2) in self.__data.keys():
                return
            self.__data[(value1, value2)] = [0] * len(self.__t)
        elif isinstance(value1, str):
            self.add(Species(value1), value2)
        elif isinstance(value1, collections.Iterable):
            for sp in value1:
                self.add(sp, value2)
        else:
            raise ValueError(
                'an invalid argument [{}] was given.'.format(repr(value1)))

    def log(self):
        self.__t.append(self.coordinator.t())
        for sp, w in self.__data.keys():
            if w is None:
                self.__data[(sp, w)].append(self.coordinator.get_value(sp))
            else:
                self.__data[(sp, w)].append(w.get_value_exact(sp))

    def savefig(
            self, filename="result.eps",
            legend=True, xlabel='Time', ylabel='The Number of Molecules'):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pylab as plt

        for (sp, w), data in self.__data.items():
            plt.plot(self.__t, data, '-', label=sp.serial())
        if legend:
            plt.legend(loc='best', shadow=True)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xlim(self.__t[0], self.__t[-1])
        plt.savefig(filename)
        plt.clf()

    def savetxt(self, filename="result.txt"):
        import numpy

        data = [self.__t]
        data.extend(self.__data.values())
        data = numpy.array(data).T
        header = 't {}'.format(' '.join(sp.serial() for sp, w in self.__data.keys()))
        numpy.savetxt(filename, data, header=header)
