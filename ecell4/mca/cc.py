import numpy

from .network import generate_full_rank_matrix


def unscaled_control_coefficients(stoichiometry, elasticity):
    _, n = stoichiometry.shape

    # do Gaussian elimination,
    # and get reduced stoichiometry, kernel and link matrix
    link_matrix, kernel_matrix, independent_list = generate_full_rank_matrix(stoichiometry)
    reduced_matrix = numpy.take(stoichiometry, independent_list, 0)

    # constract Jacobian matrix from reduced, link matrix and elasticities,
    # M0 = N0 * epsilon * L
    epsilon_L = elasticity @ link_matrix
    jacobian = reduced_matrix @ epsilon_L

    # calculate unscaled concentration control coefficients
    # CS = -L * (M0)^(-1) * N0
    inv_jacobian = numpy.linalg.inv(jacobian)
    ccc = -link_matrix @ inv_jacobian
    ccc = ccc @ reduced_matrix

    # calculate unscaled flux control coefficients
    # CJ = I - epsilon * CS
    fcc = numpy.identity(n, dtype=numpy.float) + elasticity @ ccc

    return (ccc, fcc)

def invdiag(trace):
    '''
    return numpy.lib.twodim_base.diag(1.0 / trace)
    if there\'re zeros in the array, set zero for that
    trace: (array) one dimensional array
    return (matrix)
    '''
    inv_trace = numpy.zeros(len(trace), dtype=numpy.float)
    for i in range(len(trace)):
        if abs(trace[i]) > 0.0:
            inv_trace[i] = 1.0 / trace[i]
    return numpy.lib.twodim_base.diag(inv_trace)

def scale_control_coefficients(ccc, fcc, v, x):
    # calculate scaled concentration control coefficient
    # (scaled CS_ij) = E_j / S_i * (unscaled CS_ij)
    ccc = invdiag(x) @ ccc
    ccc = ccc @ numpy.lib.twodim_base.diag(v)
    # calculate scaled flux control coefficient
    # (scaled CJ_ij) = E_j / E_i * (unscaled CJ_ij)
    fcc = invdiag(v) @ fcc
    fcc = fcc @ numpy.lib.twodim_base.diag(v)
    return (ccc, fcc)

def scaled_control_coefficients(stoichiometry, elasticity, fluxes, x):
    ccc, fcc = unscaled_control_coefficients(stoichiometry, elasticity)
    ccc, fcc = scale_control_coefficients(ccc, fcc, fluxes, x)
    return (ccc, fcc)
