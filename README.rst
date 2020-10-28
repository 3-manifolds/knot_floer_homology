Knot Floer Homology
===================

Python wrapper for Zoltán Szabó's `HFK Calculator`_.  This is beta code and is under active development.  Interface details may change in the final release.

Installing
----------

To install and test, do the following::

  python3 -m pip install .
  python3 -m zs_hfk.test

You should see a result such as ``TestResults(failed=0, attempted=24)``.

Usage
-----

In Python, do::

  >>> import zs_hfk, spherogram
  >>> L = spherogram.Link('K10n10')
  >>> hfk = zs_hfk.HFK(L)
  >>> for key in hfk: print('%s: %s'%(key, hfk[key]))
  ...
  modulus: 2
  ranks: {(-3, -2): 1, (-2, -2): 1, (-2, -1): 2, (-1, -1): 2, (-1, 0):1,
  (0, 0): 3, (1, 1): 2, (1, 2): 1, (2, 2): 1, (2, 3): 2, (3, 4): 1}
  total_rank: 17
  seifert_genus: 3
  fibered: True
  L_space_knot: False
  tau: 0
  nu: 0
  epsilon: 0
  
License
-------

Copyright Zoltán Szabó, Marc Culler, Nathan M. Dunfield, and Matthias Goerner, 2017-present.  This code is released under the `GNU General Public License, version 2`_ or (at your option) any later version as published by the Free Software Foundation.

.. _HFK Calculator: https://web.math.princeton.edu/~szabo/HFKcalc.html
.. _GNU General Public License, version 2: https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
