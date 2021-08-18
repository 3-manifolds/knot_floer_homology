# Test passes for all primes < 2**15, fails for example with p=49157

import knot_floer_homology, spherogram
from sage.all import matrix, GF, prime_range

def test_vert(knot, p):
    data = knot_floer_homology.pd_to_hfk(K, prime=p, complex=True)
    gens, diff = data['generators'], data['differentials']
    vert = {(i,j):diff[i, j] for i, j in diff if gens[i][1] == gens[j][1] + 1}
    V = matrix(GF(p), len(gens), len(gens), vert, sparse=True)
    return (V*V).rank()

K = spherogram.Link('K15a5867')
for p in prime_range(0, 2**15):
    print(p, test_vert(K, p))

