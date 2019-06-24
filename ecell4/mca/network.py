import numpy

# def get_incident_matrix(stoichiometry):
#     m, n = stoichiometry.shape
#     incident = numpy.zeros((m, n), dtype=int)
#     incident[stoichiometry > 0] = 1
#     return incident
# 
# def get_independent_groups(incident):
#     m, n = incident.shape
#     indices = [i for i in range(m) if any(incident[i] != 0)]
#     groups = []
#     while len(indices) > 0:
#         i = indices.pop(0)
#         grups.append([i])
#         row = incident[i].copy()
# 
#         c = 0
#         while c < len(indices):
#             j = indices[c]
#             if numpy.logical_and(row, incident[j]):
#                 row = numpy.logical_or(row, incident[j])
#                 groups[-1].append(indices.pop(c))
#             else:
#                 c += 1
#     return groups

def generate_full_rank_matrix(input_matrix):
    '''do Gaussian elimination and return the decomposed matrices
    input_matrix: (matrix)
    return (link_matrix, kernel_matrix, independent_list)
    '''
    m, n = input_matrix.shape
    reduced_matrix = input_matrix.copy()
    pivots = numpy.identity(m, dtype=numpy.float)

    dependent_list = list(range(m))
    independent_list = []
    skipped_list = []
    skipped_buffer = []

    for j in range(n):
        if len(dependent_list) == 0:
            break

        maxidx = dependent_list[0]
        maxelem = reduced_matrix[maxidx][j]
        for i in dependent_list:
            if abs(reduced_matrix[i][j]) > abs(maxelem):
                maxidx = i
                maxelem = reduced_matrix[i][j]

        if maxelem != 0:
            reduced_matrix[maxidx] /= maxelem
            pivots[maxidx] /= maxelem

            for i in range(m):
                if i != maxidx:
                    k = reduced_matrix[i][j]
                    reduced_matrix[i] -=  k * reduced_matrix[maxidx]
                    pivots[i] -=  k * pivots[maxidx]

            if len(skipped_buffer) > 0:
                skipped_list.extend(skipped_buffer)
                skipped_buffer = []

            dependent_list.remove(maxidx)
            independent_list.append(maxidx)
        else:
            skipped_buffer.append(j)

    assert len(dependent_list) + len(independent_list) == m

    link_matrix = numpy.identity(m, dtype=numpy.float)
    link_matrix[dependent_list] -= pivots[dependent_list]

    rank = len(independent_list)
    kernel_matrix = numpy.zeros((n, n - rank), dtype=numpy.float)
    parsed_rank = rank + len(skipped_list)
    reduced_matrix = numpy.take(reduced_matrix, list(range(parsed_rank, n)) + skipped_list, 1)

    cnt1 = 0
    cnt2 = 0
    for i in range(parsed_rank):
        if len(skipped_list) > cnt1 and skipped_list[cnt1] == i:
            kernel_matrix[i][n - parsed_rank + cnt1] = 1.0
            cnt1 += 1
        else:
            kernel_matrix[i][range(n - rank)] = -reduced_matrix[independent_list[cnt2]]
            cnt2 += 1

    for i in range(n - parsed_rank):
        kernel_matrix[i + parsed_rank][i] = 1.0

    independent_list = numpy.sort(independent_list)
    link_matrix = numpy.take(link_matrix, independent_list, 1)
    return (link_matrix, kernel_matrix, independent_list)
