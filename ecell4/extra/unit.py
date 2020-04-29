from logging import getLogger
_log = getLogger(__name__)

try:
    import pint
except ImportError:
    HAS_PINT = False
    STRICT = False
    __all__ = ['use_strict']
    # _log.warn("No module named 'pint' required by '{}'".format(__name__))
else:
    HAS_PINT = True
    STRICT = True
    from ._unit import *
    __all__ = ['getUnitRegistry', 'get_application_registry', '_Quantity', '_Unit', 'wrap_quantity', 'check_model', 'DimensionalityMismatchError']

# def use_strict(strict=True):
#     global STRICT
#     STRICT = strict
