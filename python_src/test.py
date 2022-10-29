from . import hfk
import os
import json
import sys

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
    failures = 0
    for datum in regression_data():
        name = datum.pop('name')
        pd = datum.pop('PD_code')
        result = hfk.pd_to_hfk(repr(pd))
        if datum != result:
            print('Regression failure: ' + name)
            failures += 1
                
    return failures == 0

if __name__ == '__main__':
    import doctest
    mod_tests = doctest.testmod(hfk)
    print('module tests:', mod_tests)
    reg_tests = doctest.testmod()
    print('regression tests:', reg_tests)
    if mod_tests.failed + reg_tests.failed > 0:
        sys.exit(1)
