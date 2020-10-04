ZS_HFK
======

Python wrapper for Zoltán Szabó's `HFK program`_.

Installing
----------

For now, do the following::

  cd ComputeHFKv2
  bash compile.sh
  cd ../
  python3 -m pip install .
  python3 -m zs_hfk.test

You should see a result such as ``TestResults(failed=0,
attempted=7)``.

Usage
-----

In Python, do::

  >>> import zs_hfk, spherogram
  >>> L = spherogram.Link('K10n10')
  >>> hom, props = zs_hfk.HFK(L)
  >>> hom  # Grading is (Alexander, Maslov)
  {(-3, -2): 1, (-2, -2): 1, (-2, -1): 2,  (-1, -1): 2, (-1, 0): 1,
   (0, 0): 3, (1, 1): 2, (1, 2): 1, (2, 2): 1, (2, 3): 2, (3, 4): 1}
  >>> props
  {'rank': 17, 'genus': 3, 'fibered': True, 'L_space': True,
   'tau': 0, 'nu': 0, 'epsilon': 0}

.. _HFK Program: https://web.math.princeton.edu/~szabo/HFKcalc.html


