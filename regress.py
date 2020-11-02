import spherogram, snappy, zs_hfk, knot_floer_homology

def old_way(knot_manifold):
    K = knot_manifold.link()
    ranks, props = zs_hfk.HFK(K)
    for new_key, old_key in [('L_space_knot', 'L_space'),
                             ('seifert_genus', 'genus'),
                             ('total_rank', 'rank')]:
        props[new_key] = props.pop(old_key)
    props['ranks'] = ranks
    props['modulus'] = 2
    return props

def new_way(knot_manifold):
    K = knot_manifold.link()
    return knot_floer_homology.pd_to_hfk(K)
