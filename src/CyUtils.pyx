#cython: language_level=3
#cython : #cython: boundscheck=False, wraparound=False, nonecheck=False, profile=True

import numpy as np
#cimport numpy as np

# def low(np.ndarray[np.int_t, ndim=1] col):
#     '''col: 1-dimensional array.
    
#     Gets the index of the "lowest" element in col different from 0.
#     if col=0 then low = -1
#     '''
#     cdef int i, l = -1
#     cdef int N = col.shape[0]

#     # Use memoryview for faster access!
#     cdef long[:] column = col

#     for i in range(N):
#         if column[i] > 0:
#             l = i
#     return l


def low(col):
    '''col: 1-dimensional array.
    
    Gets the index of the "lowest" element in col different from 0.
    if col=0 then low = -1
    '''
    cdef int i, l = -1
    #cdef int N = col.shape[0]
    # try the old way
    cdef int N = len(col)

    # Use memoryview for faster access
    cdef long[:] column = col

    for i in range(N):
        if column[i] > 0:
            l = i
    return l