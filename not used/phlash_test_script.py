#!/usr/bin/env python3.10

import phlash
import os.path
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import sys
import os

population=sys.argv[1]
folder_path=f"/scratch/perickso/private/ind_seq/popgen/phlash/{population}"
os.chdir(folder_path)

#from multiprocessing import Pool
#from multiprocessing import active_children
#from multiprocessing import cpu_count

#ncpus = int(os.environ['SLURM_CPUS_ON_NODE'])
#p = Pool(ncpus)

#with psmcfa files
#population="Kenya_2015_"

posterior_samples = phlash.psmc(["/scratch/perickso/private/ind_seq/popgen/phlash/Kenya_2015_/SRR11077345.autosome.psmcfa", "/scratch/perickso/private/ind_seq/popgen/phlash/Kenya_2015_/SRR11077344.autosome.psmcfa"])

times = np.array([dm.rescale(8.4e-9).eta.t[1:] for dm in posterior_samples])# choose a grid of points at which to evaluate the size history functions
T = np.geomspace(times.min(), times.max(), 1000)
Nes = np.array([dm.eta(T, Ne=True) for dm in posterior_samples])
median=np.median(Nes, axis=0)
lower_bound = np.percentile(Nes, 2.5, axis=0)
upper_bound = np.percentile(Nes, 97.5, axis=0)

df = pd.DataFrame()
df['time'] = T
df['median'] = median
df['lower2.5'] = lower_bound
df['upper97.5'] = upper_bound

df.to_csv(f"/scratch/perickso/private/ind_seq/popgen/phlash/{population}.autosome.phlash.csv")


#try with newer psmc files
