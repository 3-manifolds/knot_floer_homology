import os
import subprocess
import tempfile
import re
import spherogram

executable_path = os.path.join(os.path.dirname(__file__), 'zs_hfk_binary')

sample_output = """\
Computations from file :  /tmp/tmptdcs_r95/knot.txt
with Coefficient = 2


Knot # 1
Ranks in Alexander, Maslov bigradings :
1   (-1,-2)
1   (0,-1)
1   (1,0)
Total rank : 3
Seifert genus : 1
Fibered : Yes
L-space knot : Yes
Tau : 1
Nu : 1
Epsilon : 1
"""

props_computed = [('rank', 'Total rank', int),
                  ('genus', 'Seifert genus', int),
                  ('fibered', 'Fibered', bool),
                  ('L_space', 'L-space knot', bool),
                  ('tau', 'Tau', int),
                  ('nu', 'Nu', int),
                  ('epsilon', 'Epsilon', int)]

def run_ZS_HFK_on_PD_data(PD_data_string, prime=2):
    """
    >>> trefoil = 'PD[(4,2,5,1),(2,6,3,5),(6,4,1,3)]'
    >>> run_ZS_HFK_on_PD_data(trefoil).strip().endswith('Epsilon : 1')
    True
    """
    tmpdir = tempfile.mkdtemp()
    prime = repr(prime)

    infile_name = os.path.join(tmpdir, 'knot.txt')
    infile = open(infile_name, 'w')
    infile.write(PD_data_string.strip() + '\n')
    infile.close()

    subprocess.check_call([executable_path, infile_name, prime],
                          stdout=subprocess.PIPE)
    data = open(infile_name + '.mod' + prime).read()
    return data
    
                           
def parse_raw_data(data):
    """
    >>> hom, props = parse_raw_data(sample_output)
    >>> sum(hom.values()) == props['rank']
    True
    """
    match = re.search("Alexander, Maslov bigradings :\n(.*)(Total rank.*)",
                      data, re.DOTALL)
    homology = dict()
    for rank, grading in re.findall(r'(\d)+\s+(\(.*?,.*?\))', match[1]):
        homology[eval(grading)] = int(rank)

    props = dict()
    for key, name, value_type in props_computed:
        props[key] = value_type(re.search(name + r'\s+:\s+(.*)', match[2])[1])
    return homology, props

def HFK(knot):
    """
    >>> hom, props = HFK('K8n3')
    >>> sum(hom.values()) == props['rank']
    True
    >>> props['tau'], props['nu'], props['epsilon']
    (3, 3, 1)
    """
    if not isinstance(knot, spherogram.Link):
        knot = spherogram.Link(knot)
    if len(knot.link_components) != 1 or knot.unlinked_unknot_components != 0:
        raise ValueError('Can only handle knots not links')
    data = run_ZS_HFK_on_PD_data('PD' + repr(knot.PD_code()))
    return parse_raw_data(data)
        
