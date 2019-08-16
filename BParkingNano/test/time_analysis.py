import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from pdb import set_trace
eps = 10**-7
from argparse import ArgumentParser
import numpy as np

parser = ArgumentParser()
parser.add_argument('inval', help='file path')
parser.add_argument('--tag', help='name to use')
args = parser.parse_args()

#Module Summary
infile = open(args.inval).read()
modules = infile.split('TimeReport ---------- Module Summary ---[Real sec]----')[1].split('T---Report end!')[0]
time_rep = infile.split('TimeReport ---------- Event  Summary')[1].split('\n\n')[0].split('\n')

module_times = []
for l in modules.split('\n'):
  line = l.strip()
  if not line: continue
  if line.endswith('Name'): continue
  info = line.split()
  module_times.append((info[-1], float(info[1])))

tot_time = sum(i for _, i in module_times)
from collections import defaultdict
groups = defaultdict(float)

# This part is ugly, but I have no better idea
for name, time in module_times:
  if 'gen' in name.lower():
    groups['GEN'] += time
  elif 'lhe' in name.lower():
    groups['GEN'] += time
  elif 'Table' in name:
    groups['Tables (not GEN)'] += time
  elif 'kee' in name.lower():
    groups['BToKee'] += time
  elif 'kmumu' in name.lower():
    groups['BToKmumu'] += time
  elif 'electron' in name.lower():
    groups['Electrons'] += time
  elif 'track' in name.lower():
    groups['Tracks'] += time
  elif 'muon' in name.lower():
    groups['Muons'] += time
  elif name.endswith('output'):
    groups['I/O'] += time
  else:
    groups['Other'] += time

import matplotlib.colors as colors
import matplotlib.cm as cmx
cm = plt.get_cmap('rainbow')
cNorm  = colors.Normalize(vmin=0, vmax=len(groups)-1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
cols = [scalarMap.to_rgba(i) for i in range(len(groups))]

def writer(pct):
  if pct < 5: return ''
  else: return '%.1f%%' % pct

plt.clf()
fig = plt.figure(
    figsize=(12, 6), 
)
plt.subplot(1, 2, 1)
vals = np.array(groups.values())/tot_time
wedges = plt.pie(vals, autopct = writer, colors = cols)
names = ['%s (%.1f%%)' % (n, float(p)*100/tot_time) for n, p in groups.iteritems()]
leg = plt.legend(
    wedges[0], names, loc = 5,
    bbox_to_anchor = (0.95, 0.5),
    mode="expand", borderaxespad=0., frameon=False
)
cpu = float(time_rep[1].split(' = ')[1])
wall = float(time_rep[2].split(' = ')[1])

title = 'Execution time: %.3f [s] (CPU / evt), %.3f [s] (Wall / evt)' % (cpu, wall)
plt.title(title)
fig.savefig('validation/timing_%s.png' % args.tag)
