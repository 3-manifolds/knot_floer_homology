Knot Floer Homology
===================

Python wrapper for Zoltán Szabó's `HFK Calculator`_.

Installing
----------

To install and test, do the following to install the current release
from `PyPI`_::

  python3 -m pip install knot_floer_homology
  python3 -m knot_floer_homology.test

You should see a result such as ``TestResults(failed=0, attempted=24)``.

Usage
-----

In Python, do::

  >>> import knot_floer_homology
  >>> PD = [(2,0,3,15),(0,6,1,5),(6,2,7,1),(3,10,4,11),(9,4,10,5),(7,12,8,13),(13,8,14,9),(11,14,12,15)]
  >>> knot_floer_homology.pd_to_hfk(PD)
  {'L_space_knot': False,
   'epsilon': 0,
   'fibered': True,
   'modulus': 2,
   'nu': 0,
   'ranks': {(-2, -2): 1, (-1, -1): 2, (0, 0): 3, (1, 1): 2, (2, 2): 1},
   'seifert_genus': 2,
   'tau': 0,
   'total_rank': 9}

It also accepts `Spherogram`_ knots as input::

  >>> import spherogram
  >>> L = spherogram.Link('K10n10')
  >>> ans = knot_floer_homology.pd_to_hfk(L)
  >>> ans['seifert_genus']
  3


License
-------

Copyright Zoltán Szabó, Marc Culler, Nathan M. Dunfield, and Matthias
Goerner, 2017-present.  This code is released under the `GNU General
Public License, version 2`_ or (at your option) any later version as
published by the Free Software Foundation.

.. _HFK Calculator: https://web.math.princeton.edu/~szabo/HFKcalc.html
.. _PyPI: https://pypi.org
.. _Spherogram: https://github.com/3-manifolds/Spherogram
.. _GNU General Public License, version 2: https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
