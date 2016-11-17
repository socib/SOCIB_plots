#!/usr/bin/python

if __name__ == '__main__':

  import sys
  import numpy as np
  from netCDF4 import Dataset
  t = 0
  s = 0
  for f in sys.argv[1:]:
    with Dataset(f, 'r') as d:
      p = d.variables['profile_index'][:]
      c = d.variables['conductivity'][:]
      g = ~np.isnan(c)
      m = np.int(np.floor(p[-1]))
      ok = np.zeros((m,), bool)
      for i in range(m):
          sel = g[p == i+1]
          ok[i] = np.sum(sel) > (0.3 * sel.size)
      n = np.sum(ok)
      s += m
      t += n
      print('{}: {} of {}'.format(f, n, m))
  print('total: {} of {}'.format(t, s))
