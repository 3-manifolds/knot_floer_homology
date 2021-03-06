# distutils: language = c++
# cython: language_level = 3

from cpython.ref cimport PyObject, Py_DECREF

cdef extern from "PyWrapper.h":
    cdef PyObject* PDCodeToMorse(const char *pd) except *
    cdef PyObject* PDCodeToHFK(const char *pd, int prime) except *

def _get_pd_string(pd):
    if hasattr(pd, 'PD_code'):
        pd = repr(pd.PD_code())
    if isinstance(pd, list):
        pd = repr(pd)
    return pd.encode('ascii')

def pd_to_morse(pd):
    """
    >>> pd = 'PD[(5,3,0,2),(1,5,2,4),(3,1,4,0)]'
    >>> morse = pd_to_morse(pd)
    >>> morse['girth']
    4
    """

    result = <object>PDCodeToMorse(_get_pd_string(pd))
    Py_DECREF(result)
    return result

def pd_to_hfk(pd_code, int prime = 2):
    """
    >>> pd = 'PD[(5,3,0,2),(1,5,2,4),(3,1,4,0)]'
    >>> HFK = pd_to_hfk(pd)
    >>> sorted(HFK['ranks'].keys())
    [(-1, -2), (0, -1), (1, 0)]
    >>> HFK['total_rank']
    3

    >>> pd = 'PD[(2,0,3,19),(6,2,7,1),(4,10,5,9),(15,8,16,9),(13,18,14,19),'
    >>> pd +=   '(7,14,8,15),(10,4,11,3),(17,12,18,13),(0,6,1,5),(11,16,12,17)]'
    >>> HFK = pd_to_hfk(pd)
    >>> HFK['total_rank']
    17
    >>> HFK['tau'], HFK['nu'], HFK['epsilon']
    (0, 0, 0)

    >>> pd_code = [(2,0,3,15), (0,6,1,5), (6,2,7,1), (3,11,4,10),
    ...            (11,5,12,4), (7,13,8,12), (13,9,14,8), (9,15,10,14)]
    >>> HFK = pd_to_hfk(pd_code)
    >>> sum(HFK['ranks'].values()) == HFK['total_rank']
    True
    >>> HFK['tau'], HFK['nu'], HFK['epsilon']
    (3, 3, 1)
    """
    result = <object>PDCodeToHFK(_get_pd_string(pd_code), prime)
    Py_DECREF(result)
    return result
