# distutils: language = c++
# cython: language_level = 3

from cpython.ref cimport PyObject, Py_DECREF
from libcpp cimport bool

cdef extern from "PyWrapper.h":
    cdef PyObject* PDCodeToMorse(const char *pd) except *
    cdef PyObject* PDCodeToHFK(const char *pd, int prime, bool complex) except *

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

def pd_to_hfk(pd_code, int prime=2, bool complex=False):
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

    If the parameter `complex` is set to True, then the simplified
    "UV = 0" knot Floer chain complex is returned. This complex is
    computed over the ring F[U,V]/(UV = 0), where F is the integers
    mod the chosen prime; this corresponds to only the horizontal and
    vertical arrows in the full knot Floer complex. The complex is
    specified by:

    * generators: a dictionary from the generator names to their
      (Alexander, Maslov) gradings.  The number of generators is
      equal to the total_rank.

    * differential: a dictionary whose value on (a, b) is an integer
      specifying the coefficient on the differential from generator a
      to generator b, where only nonzero differentials are
      recorded. (The coefficient on the differential is really an
      element of F[U,V]/(UV = 0), but the power of U or V can be
      recovered from the gradings on a and b so only the element of F
      is recorded.)

    For example, to compute the vertical differential, whose homology
    is HFhat(S^3), you can do:

    >>> data = pd_to_hfk(pd_code, complex=True)
    >>> gens, diff = data['generators'], data['differentials']
    >>> vert = {(i,j):diff[i, j] for i, j in diff
    ...                          if gens[i][1] == gens[j][1] + 1}
    >>> len(vert)
    2
    """
    result = <object>PDCodeToHFK(_get_pd_string(pd_code), prime, complex)
    Py_DECREF(result)
    return result
