import uproot
import pandas as pd
import numpy as np
import awkward
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('f_old', help='file path')
parser.add_argument('f_new', help='file path')
parser.add_argument('--legacy', action='store_true', help='compare against legacy version')
args = parser.parse_args()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from pdb import set_trace
eps = 10**-7

if not os.path.isdir('validation'):
  os.makedirs('validation')

legacy_mapping = { #mapping between George's Ntuples and Nano
  'nmuon' : 'nMuon',
  'muon_pt' : 'Muon_pt',
  'muon_eta' : 'Muon_eta',
  'muon_phi' : 'Muon_phi',
  'muon_charge' : 'Muon_charge',
  'muon_dxy' : 'Muon_dxy',
  'muon_edxy' : 'Muon_dxyErr',
  'muon_soft' : 'Muon_softId',
  #'muon_trgIndex' : '',
  
  'nelectron' : 'nElectron',
  'el_pt' : 'Electron_pt',
  'el_eta' : 'Electron_eta',
  'el_phi' : 'Electron_phi',
  'el_charge' : 'Electron_charge',
  'el_mva_unbiased' : 'Electron_unBiased',
  'el_islowpt' : 'Electron_isLowPt',
  # TODO: Add variables as they are validated and produced
}
# invert map as well
legacy_inverted = {v : k for k, v in legacy_mapping.iteritems()}

class NanoFrame(object):
  def __init__(self, infile, legacy = False):
    self.uf = uproot.open(infile)
    self.tt = self.uf['demo/mytree'] if legacy else self.uf['Events']
    self.legacy = legacy

  def __getitem__(self, key):
    return self.tt.array(legacy_inverted[key] if self.legacy else key)

  def keys(self):
    return legacy_mapping.keys() if self.legacy else self.tt.keys()

def byval_validation(v1, v2):
  if np.isinf(v2).any() or np.isnan(v2).any():
    v1 = v1[np.invert(np.isinf(v1) | np.isnan(v1).any())]
    v2 = v2[np.invert(np.isinf(v2) | np.isnan(v2).any())]

  try:
    if v1.dtype == 'bool' or np.issubdtype(v1.dtype, np.integer):
      return (v1 == v2).all()
    else:
      return ((np.abs(v1 - v2) / (abs(v1) + eps)) < 0.001).all()
  except ValueError:
    return False

def stat_validation(v1, v2, name = '', nbins = 20):
  if np.isinf(v2).any() or np.isnan(v2).any():
    print '\033[1;35m', name, '--> CONTAINS INFs/NANs!\033[0m'
    v1 = v1[np.invert(np.isinf(v1) | np.isnan(v1).any())]
    v2 = v2[np.invert(np.isinf(v2) | np.isnan(v2).any())]

  if v1.shape[0] == 0 and v2.shape[0] == 0:
    return True
  elif v1.shape[0] == 0 or v2.shape[0] == 0:
    return False

  M = max(v1.max(), v2.max())
  m = min(v1.min(), v2.min())
  m = m * 0.9 if m > 0 else m * 1.2
  M = M * 1.2 if M > 0 else M * 0.9
  if 'int' in str(v1.dtype):
    m = int(m) - 1
    M = int(M) + 1
    nbins = min(M - m, nbins*2)
  plt.clf()
  h1, _, _ = plt.hist(v1, range = (m,M), bins = nbins, label = 'old', histtype = 'step')
  h2, _, _ = plt.hist(v2, range = (m,M), bins = nbins, label = 'new', histtype = 'step')
  plt.legend(loc='best')
  plt.savefig('validation/%s.png' % name)
  plt.clf()
  return (h1 == h2).all()

old = NanoFrame(args.f_old, args.legacy)
new = NanoFrame(args.f_new)

#
# Size Checks
#
uf = new.uf
tt = new.tt

branches_and_size = {i.name : i.compressedbytes() for i in tt.allvalues()}
tot_branches = sum(branches_and_size.values())
n_entries = len(tt)
try:
  n_processed = int(uf['tag'].split('nevts:')[1])
except:
  n_processed = -1

from collections import defaultdict
groups = defaultdict(long)
for name, size in branches_and_size.iteritems():
    group = name.split('_')[0] if '_' in name else 'other'
    groups[group] += size

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

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
wedges = plt.pie(groups.values(), autopct = writer, colors = cols)
names = ['%s (%.1f%%)' % (n, float(p)*100/tot_branches) for n, p in groups.iteritems()]
leg = plt.legend(
    wedges[0], names, loc = 5,
    bbox_to_anchor = (0.95, 0.5),
    mode="expand", borderaxespad=0., frameon=False
)
title = 'Total size: %.3f kB / evt (%d events / %d processed)' % (tot_branches/(10.**3 * n_entries), n_entries, n_processed)
print title
plt.title(title)
fig.savefig('validation/size.png')

#
# Branch checks
#
old_k = set(old.keys())
new_k = set(new.keys())
intersection = old_k.intersection(new_k)

for branch in intersection:
  v_old = old[branch]
  v_new = new[branch]

  if hasattr(v_old, 'content'):
    v_old = v_old.content
    v_new = v_new.content

  stat_valid = stat_validation(v_old, v_new, branch)
  val_valid  = byval_validation(v_old, v_new)

  if val_valid and stat_valid:
    print '\033[1;32m', branch, '--> OK!\033[0m'
  elif stat_valid:
    print '\033[1;35m', branch, '--> FAILS BY VALUE CHECK ONLY!\033[0m'
  else:
    print '\033[1;31m', branch, '--> FAILS ALL CHECKS!\033[0m'
