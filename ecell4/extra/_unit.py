import itertools
import sys
import warnings
from logging import getLogger
_log = getLogger(__name__)

from ..util.parseobj import AnyCallable, ExpBase
from ..util.decorator_base import just_parse
from ..util import model_parser

import pint
from pint.quantity import _Quantity
from pint.unit import _Unit
# from pint.errors import UndefinedUnitError
from pint.errors import DimensionalityError


__all__ = [
    'getUnitRegistry', '_Quantity', '_Unit', 'wrap_quantity',
    'get_application_registry', 'check_model', 'DimensionalityMismatchError']

def wrapped_binary_operator(op1, op2):
    def wrapped(self, other):
        if isinstance(other, ExpBase):
            return op2(other, self)
        elif isinstance(other, AnyCallable):
            return op2(other._as_ParseObj(), self)
        return op1(self, other)
    return wrapped

def wrap_quantity(cls):
    cls.__add__ = wrapped_binary_operator(cls.__add__, ExpBase.__radd__)
    cls.__mul__ = wrapped_binary_operator(cls.__mul__, ExpBase.__rmul__)
    cls.__div__ = wrapped_binary_operator(cls.__div__, ExpBase.__rdiv__)
    cls.__truediv__ = wrapped_binary_operator(cls.__truediv__, ExpBase.__rtruediv__)
    cls.__pow__ = wrapped_binary_operator(cls.__pow__, ExpBase.__rpow__)
    return cls

def getUnitRegistry(length="meter", time="second", substance="item", volume=None, other=()):
    """Return a pint.UnitRegistry made compatible with ecell4.

    Parameters
    ----------
    length : str, optional
        A default unit for '[length]'. 'meter' is its default.
    time : str, optional
        A default unit for '[time]'. 'second' is its default.
    substance : str, optional
        A default unit for '[substance]' (the number of molecules). 'item' is its default.
    volume : str, optional
        A default unit for '[volume]'. Its default is None, thus '[length]**3'.
    other : tuple, optional
        A list of user-defined default units other than the above.

    Returns
    -------
    ureg : pint.UnitRegistry

    """
    ureg = pint.UnitRegistry()
    ureg.define('item = mole / avogadro_number')
    assert ureg.Quantity(1, 'item').check('[substance]')

    try:
        ureg.molar
    # except UndefinedUnitError:
    except AttributeError:
        # https://github.com/hgrecco/pint/blob/master/pint/default_en.txt#L75-L77
        ureg.define('[concentration] = [substance] / [volume]')
        ureg.define('molar = mol / (1e-3 * m ** 3) = M')

    base_units = [unit for unit in (length, time, substance, volume) if unit is not None]
    base_units.extend(other)
    _ = ureg.System.from_lines(
        ["@system local using international"] + base_units,
        ureg.get_base_units)
    ureg.default_system = 'local'

    wrap_quantity(ureg.Quantity)

    pint.set_application_registry(ureg)  # for pickling
    return ureg

def get_application_registry():
    """Just return `pint._APP_REGISTRY`."""
    return pint._APP_REGISTRY

def __check_dimensionality_of_quantity(q, dim, Quantity):
    import ecell4_base.core
    assert isinstance(q, (ecell4_base.core.Quantity_Real, ecell4_base.core.Quantity_Integer))
    if q.units == "":
        return (dim == "")
    q_ = Quantity(q.magnitude, q.units)
    return q_.check(dim)

def __is_structure_from_model(sp, m):
    for another in m.species_attributes():
        if (another.has_attribute("location")
                and another.get_attribute("location") == sp.serial()):
            return True
    return False

def __get_dimension_from_model(sp, m):
    from ecell4_base.core import get_dimension_from_model
    return get_dimension_from_model(sp, m)

def __get_dimensionality_of_reaction(obj, m):
    from ecell4_base.core import ReactionRule, Shape
    assert isinstance(obj, ReactionRule)
    if obj.has_descriptor():
        species_list = set(itertools.chain(obj.reactants(), obj.products()))
        dims = {sp.serial(): __get_dimension_from_model(sp, m) for sp in species_list}
        for serial in dims:
            if dims[serial] == Shape.UNDEF:
                raise RuntimeError(
                    f"The dimension of [{serial}] is undefined [{obj.as_string()}].")
        dim = max(int(dim) for dim in dims.values())
        #TODO: What if a structure in the list?
        return f"[substance]/[length]**{dim}/[time]"

    if len(obj.reactants()) > 2:
        dims = [__get_dimension_from_model(sp, m) for sp in obj.reactants()]
        if any(dim == Shape.UNDEF for dim in dims):
            raise RuntimeError(f"Dimension is undefined [{obj.as_string()}].")
        elif any(__is_structure_from_model(sp, m) for sp in obj.reactants()):
            raise RuntimeError(
                f"No structure is allowed in a reaction with more than 2 reactants [{obj.as_string()}].")
        elif len(set(dims)) != 1:
            raise RuntimeError(
                "All reactants must have the same dimension"
                " for a reaction with more than 2 reactants"
                f" [{obj.as_string()}].")
        return f"[length]**{int(tuple(set(dims))[0])}/[substance]/[time]"
    elif len(obj.reactants()) == 2:
        sp1, sp2 = obj.reactants()
        is_structure1, is_structure2 = __is_structure_from_model(sp1, m), __is_structure_from_model(sp2, m)
        assert not (is_structure1 and is_structure2)
        dim1, dim2 = __get_dimension_from_model(sp1, m), __get_dimension_from_model(sp2, m)
        if dim1 == Shape.UNDEF:
            raise RuntimeError(
                f"The dimension of [{sp1.serial()}] is undefined [{obj.as_string()}].")
        elif dim2 == Shape.UNDEF:
            raise RuntimeError(
                f"The dimension of [{sp2.serial()}] is undefined [{obj.as_string()}].")
        dim1, dim2 = int(dim1), int(dim2)
        if is_structure1:
            if dim1 > dim2:
                raise RuntimeError(f"Wrong dimensions [{dim1} > {dim2}] in [{obj.as_string()}].")
            return f"[length]**{max(1, dim2 - dim1)}/[time]"
        elif is_structure2:
            if dim1 < dim2:
                raise RuntimeError(f"Wrong dimensions [{dim1} < {dim2}] in [{obj.as_string()}].")
            return f"[length]**{max(1, dim1 - dim2)}/[time]"
        return f"[length]**{max(dim1, dim2)}/[substance]/[time]"
    elif len(obj.reactants()) == 1:
        sp1 = obj.reactants()[0]
        if __is_structure_from_model(sp1, m):
            dim = __get_dimension_from_model(sp1, m)
            if dim == Shape.UNDEF:
                raise RuntimeError(
                    f"The dimension of [{sp1.serial()}] is undefined [{obj.as_string()}].")
            return f"[substance]/[time]/[length]**{int(dim)}"
        return "1/[time]"
    elif len(obj.reactants()) == 0:
        assert len(obj.products()) != 0
        for sp in obj.products():
            if __is_structure_from_model(sp, m):
                raise RuntimeError(
                    f"No synthetic reaction of a structure [{sp.serial()}] is accepted [{obj.as_string()}].")
        dims = [__get_dimension_from_model(sp, m) for sp in obj.products()]
        if len(set(dims)) != 1:
            raise RuntimeError(
                f"Synthetic products must have the same dimension [{obj.as_string()}].")
        elif dims[0] == Shape.UNDEF:
            raise RuntimeError(
                f"The dimension is undefined [{obj.as_string()}].")
        dims = int(dims[0])
        return f"[substance]/[time]/[length]**{dims}"


class DimensionalityMismatchError(ValueError):
    pass

def check_model(m, Quantity=None):
    from ecell4_base.core import Model, Quantity_Real, Shape
    assert isinstance(m, Model)
    Quantity = Quantity if Quantity is not None else get_application_registry().Quantity

    attribute_dimensions = {'D': '[length]**2/[time]', 'radius': '[length]', 'dimension': ''}
    for sp in m.species_attributes():
        for key, dim in attribute_dimensions.items():
            if not sp.has_attribute(key):
                continue
            elif not __check_dimensionality_of_quantity(sp.get_attribute(key), dim, Quantity):
                raise DimensionalityMismatchError(
                    f"Attribute [{key}] in Species [{sp.serial()}] has wrong dimensionality."
                    f" [{dim}] is expected.")

    for rr in m.reaction_rules():
        dim = __get_dimensionality_of_reaction(rr, m)

        units = ''
        if not rr.has_descriptor():
            units = rr.get_k().units
        else:
            species_list = set(itertools.chain(rr.reactants(), rr.products()))
            dims = {sp.serial(): __get_dimension_from_model(sp, m) for sp in species_list}
            params = {}
            for serial in dims:
                assert dims[serial] != Shape.UNDEF  # Already checed in __get_dimensiononality_of_reaction
                #TODO: 1.0? The sloppy value can cause another error like ZeroDivisionError.
                params[serial] = Quantity(1.0, f"item / meter ** {int(dims[serial])}")

            desc = rr.get_descriptor()
            try:
                #XXX: Mathematical functions in the math module are not supported by pint
                params_ = dict(Quantity_Real=Quantity, Quantity_Integer=Quantity, **model_parser.RATELAW_RESERVED_FUNCTIONS, **model_parser.RATELAW_RESERVED_CONSTANTS)
                params_['pow'] = lambda base, exp: base ** exp  # numpy.power is not supported by pint
                params_.update(params)
                ret = just_parse().eval(desc.as_string(), params_)
            except DimensionalityError as e:
                #TODO: no units could be defined in the given law.
                raise DimensionalityMismatchError(f"Failed to evaluate [{desc.as_string()}] ({rr.as_string()}). {str(e)}")
            else:
                units = str(ret.units)

        if not __check_dimensionality_of_quantity(Quantity_Real(1.0, units), dim, Quantity):
            raise DimensionalityMismatchError(
                f"ReactionRule [{rr.as_string()}] has wrong dimensionality [{Quantity(1.0, units).dimensionality}]."
                f" [{dim}] is expected.")


if __name__ == '__main__':
    ureg = getUnitRegistry(length='micrometer', volume='liter')
    Q_ = ureg.Quantity
    q = Q_(1.0, "ml")
    print(q)
    print(q.to_base_units())
    q = Q_(1.0, "M")
    print(q)
    print(q.to_base_units())
