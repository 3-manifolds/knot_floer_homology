from . import hfk
import os
import json

def regression_data():
    data_file = os.path.join(os.path.dirname(__file__), 'HFK_data.json')
    data = json.load(open(data_file))
    for d in data:
        ranks = d['ranks']
        for k in list(ranks.keys()):
            ranks[eval(k)] = ranks.pop(k)
    return data

def matches_saved_HFK_results():
    """
    >>> matches_saved_HFK_results()
    True
    """
    for datum in regression_data():
        name = datum.pop('name')
        pd = datum.pop('PD_code')
        result = hfk.pd_to_hfk(repr(pd))
        if datum != result:
            return False
        return True

if __name__ == '__main__':
    import doctest
    print('module tests:', doctest.testmod(hfk))
    print('regression tests:', doctest.testmod())
