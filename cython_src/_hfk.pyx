from libc.stdlib cimport free
cdef extern void PDCodeToMorseAndHFK(char *pd, int prime, char** morse, char **hfk, char **error)

def pd_to_morse(pd, int prime=2):
    cdef char* morse
    cdef char* error
    cdef errorstring, morsestring
    PDCodeToMorseAndHFK(pd.encode('ascii'), prime, &morse, NULL, &error)
    error_string = error.decode('ascii')
    free(error)
    morse_string = morse.decode('ascii')
    free(morse)
    if error_string:
        raise ValueError(error_string)
    else:
        return morse_string

def pd_to_hfk(pd, int prime=2):
    cdef char* hfk
    cdef char* error
    cdef errorstring, morsestring
    PDCodeToMorseAndHFK(pd.encode('ascii'), prime, NULL, &hfk, &error)
    error_string = error.decode('ascii')
    free(error)
    hfk_string = hfk.decode('ascii')
    free(hfk)
    if error_string:
        raise ValueError(error_string)
    else:
        return hfk_string
