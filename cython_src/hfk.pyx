# distutils: language = c++
from libc.stdlib cimport free
cdef extern void PDCodeToMorseAndHFK(char *pd, int prime, char** morse,
                                     char **hfk, char **error)

def pd_to_morse(pd):
    """
    >>> pd = 'PD(5,3,0,2),(1,5,2,4),(3,1,4,0)]'
    >>> morse = pd_to_morse(pd)
    >>> morse['girth']
    4
    """
    cdef char* morse
    cdef char* error
    cdef errorstring, morsestring
    if hasattr(pd, 'PD_code'):
        pd = 'PD' + repr(pd.PD_code())
    PDCodeToMorseAndHFK(pd.encode('ascii'), 2, &morse, NULL, &error)
    error_string = error.decode('ascii')
    free(error)
    morse_string = morse.decode('ascii')
    free(morse)
    if error_string:
        raise ValueError(error_string)
    else:
        return eval(morse_string)

def pd_to_hfk(pd, int prime=2):
    """
    >>> pd = 'PD(5,3,0,2),(1,5,2,4),(3,1,4,0)]'
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
    """

    cdef char* hfk
    cdef char* error
    cdef errorstring, morsestring
    if hasattr(pd, 'PD_code'):
        pd = 'PD' + repr(pd.PD_code())
    PDCodeToMorseAndHFK(pd.encode('ascii'), prime, NULL, &hfk, &error)
    error_string = error.decode('ascii')
    free(error)
    hfk_string = hfk.decode('ascii')
    free(hfk)
    if error_string:
        raise ValueError(error_string)
    else:
        result = {}
        exec('result.update(%s)'%hfk_string)
        return result
