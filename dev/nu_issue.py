import snappy
import collections
from spherogram.links.simplify import random_reverse_type_II

knots = snappy.HTLinkExteriors(cusps=1, alternating=False)
alt_knots = snappy.HTLinkExteriors(cusps=1, alternating=True)

def test(n):
    epsilons = collections.Counter()
    for i in range(n):
        K = knots.random().link()
        h = K.knot_floer_homology()
        m= K.mirror().knot_floer_homology()
        assert h['tau'] == -m['tau']
        if h['epsilon'] == -1:
            assert h['nu'] == h['tau'] + 1
            assert m['nu'] == m['tau']
        elif h['epsilon'] == 0:
            assert h['nu'] == h['tau'] and m['nu'] == m['tau']
        else:
            assert h['epsilon'] == 1
            assert m['nu'] == m['tau'] + 1
            assert h['nu'] == h['tau']
        epsilons.update([h['epsilon']])
        print(epsilons)

def test_fix(n):
    epsilons = collections.Counter()
    for i in range(n):
        M = alt_knots.random()
        K = M.link()
        L = K.copy()
        random_reverse_type_II(L, 'a', 'b')
        L._rebuild()
        h = K.knot_floer_homology()
        m = L.knot_floer_homology()
        assert h['tau'] == m['tau']
        assert h['nu'] == m['nu']
        assert h['epsilon'] == m['epsilon']
        epsilons.update([h['epsilon']])
        print(epsilons)
    
    
