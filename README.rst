ZS_HFK
======

Python wrapper for Zoltán Szabó's `HFK program`_.

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

.. _HFK Program: https://web.math.princeton.edu/~szabo/HFKcalc.html
