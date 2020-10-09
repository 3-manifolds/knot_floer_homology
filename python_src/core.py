import spherogram
from .hfk import pd_to_hfk

sample_output = """\
>>> import zs_hfk, spherogram
>>> L = spherogram.Link('K10n10')
>>> zs_hfk.HFK(L)
{'modulus': 2, 'ranks': {(-3, -2): 1, (-2, -2): 1, (-2, -1): 2, (-1, -1): 2, (-1, 0): 1, (0, 0): 3, (1, 1): 2, (1, 2): 1, (2, 2): 1, (2, 3): 2, (3, 4): 1}, 'total_rank': 17, 'seifert_genus': 3, 'fibered': True, 'L_space_knot': False, 'tau': 0, 'nu': 0, 'epsilon': 0}
"""

def HFK(knot, invariant=None):
    """
    >>> hfk = HFK('K8n3')
    >>> sum(hfk['ranks'].values()) == hfk['total_rank']
    True
    >>> hfk['tau'], hfk['nu'], hfk['epsilon']
    (3, 3, 1)
    """
    if not isinstance(knot, spherogram.Link):
        knot = spherogram.Link(knot)
    if len(knot.link_components) != 1 or knot.unlinked_unknot_components != 0:
        raise ValueError('Can only handle knots not links')
    hfk_data = pd_to_hfk('PD' + repr(knot.PD_code()))
    if isinstance(invariant, str):
        return hfk_data[invariant]
    else:
        return hfk_data
