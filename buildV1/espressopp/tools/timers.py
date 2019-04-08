#  Copyright (C) 2012,2013,2016
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


import sys

def show(alltimers, precision=1):
  """
  Python functions to print timings from C++.
  """
  print 'timers (Run with overhead, Pair, FENE, Angle, Comm1, Comm2, Int1, Int2, Resort, Other, Overhead)'
  for rank in range(len(alltimers)):
      print 'RANK ' + str(rank) + ": " + str(alltimers[rank])
  
  num_fmt = '%' + str(precision + 4) + '.' + str(precision) + 'f'
  fmt1 = num_fmt + '\n'
  fmt2 = num_fmt + ' (' + num_fmt + ') \n'
  fmt3 = num_fmt + ' (' + num_fmt + '), MIN ' + num_fmt +', MAX ' + num_fmt + '\n'
  t_avg=[]
  t_min=[]
  t_max=[]
  nprocs = len(alltimers)
  for ntimer in xrange(11):
    t_avg.append(0.0)
    t_min.append(100000000)
    t_max.append(0.0)
    for k in xrange(nprocs):
      t_avg[ntimer] += alltimers[k][ntimer]
      t_min[ntimer] = min(t_min[ntimer], alltimers[k][ntimer])
      t_max[ntimer] = max(t_max[ntimer], alltimers[k][ntimer])
    t_avg[ntimer] /= nprocs

  sys.stdout.write('\n')
  sys.stdout.write('*  Times WITH overhead caused by duplicate meshes:\n')
  sys.stdout.write('Run (TOTAL)    AVG time (%) ' + fmt3 % (t_avg[0], 100, t_min[0], t_max[0]))
  run_wo_owh = t_avg[0] - t_avg[10]
  sys.stdout.write('Run w/o overh. AVG time (%) ' + fmt2 % (run_wo_owh, 100*run_wo_owh/t_avg[0]))
  sys.stdout.write('Overhead       AVG time (%) ' + fmt3 % (t_avg[10], 100*t_avg[10]/t_avg[0], t_min[10], t_max[10]))
  sys.stdout.write('*  Times WITHOUT overhead (comparable to standard espresso):\n')
  sys.stdout.write('Run    AVG time (%) ' + fmt2 % (run_wo_owh, 100*run_wo_owh/run_wo_owh))
  sys.stdout.write('Pair   AVG time (%) ' + fmt3 % (t_avg[1], 100*t_avg[1]/run_wo_owh, t_min[1], t_max[1]))
  sys.stdout.write('FENE   AVG time (%) ' + fmt3 % (t_avg[2], 100*t_avg[2]/run_wo_owh, t_min[2], t_max[2]))
  sys.stdout.write('Angle  AVG time (%) ' + fmt3 % (t_avg[3], 100*t_avg[3]/run_wo_owh, t_min[3], t_max[3]))
  sys.stdout.write('Comm1  AVG time (%) ' + fmt3 % (t_avg[4], 100*t_avg[4]/run_wo_owh, t_min[4], t_max[4]))
  sys.stdout.write('Comm2  AVG time (%) ' + fmt3 % (t_avg[5], 100*t_avg[5]/run_wo_owh, t_min[5], t_max[5]))
  sys.stdout.write('Int1   AVG time (%) ' + fmt3 % (t_avg[6], 100*t_avg[6]/run_wo_owh, t_min[6], t_max[6]))
  sys.stdout.write('Int2   AVG time (%) ' + fmt3 % (t_avg[7], 100*t_avg[7]/run_wo_owh, t_min[7], t_max[7]))
  sys.stdout.write('Resort AVG time (%) ' + fmt3 % (t_avg[8], 100*t_avg[8]/run_wo_owh, t_min[8], t_max[8]))
  sys.stdout.write('Other  AVG time (%) ' + fmt3 % (t_avg[9], 100*t_avg[9]/run_wo_owh, t_min[9], t_max[9]))
  sys.stdout.write('\n')
