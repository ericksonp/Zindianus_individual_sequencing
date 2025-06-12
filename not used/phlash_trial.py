# conda activate phlash
# export PYTHONHOME=/usr/local/sw/anaconda/anaconda3/envs/phlash
# python3.10

import phlash
import os.path

onekg_base = "/scratch/perickso/private/ind_seq/phlash"

chr1 = phlash.contig("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gz", samples=["SRR11077296"], region="Scaffold_1:1-10000000")

chroms_1kg = []
for chrom in range(1, 5):
    chroms_1kg.append(
        phlash.contig("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gz", samples=["SRR11077296"], region=f"Scaffold_{chrom}:2100000-2200000")
    )
results = phlash.fit(chroms_1kg, mutation_rate=1.29e-8)
import matplotlib.pyplot as plt
import numpy as np

import pickle

with open('/scratch/perickso/private/ind_seq/phlash/data.pickle', 'wb') as f:
    pickle.dump(results, f)

import csv
with open("/scratch/perickso/private/ind_seq/phlash/data.csv", "w") as f:
    csv_writer = csv.writer(f)
    for mytuple in results:
        csv_writer.writerow(mytuple)

#with psmcfa files
posterior_samples = phlash.psmc(['/scratch/perickso/private/ind_seq/popgen/psmc/ZAP_3_C3.psmcfa', '/scratch/perickso/private/ind_seq/popgen/psmc/ZAP_3_C2.psmcfa'])
import pandas as pd


import matplotlib.pyplot as plt
import numpy as np
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

df.to_csv("/scratch/perickso/private/ind_seq/phlash/results_matrix_test.csv")


#try with newer psmc files

posterior_samples = phlash.psmc(['/scratch/perickso/private/ind_seq/popgen/phlash/Kenya_2015_/SRR11077344.psmcfa',
'/scratch/perickso/private/ind_seq/popgen/phlash/Kenya_2015_/SRR11077345.psmcfa', '/scratch/perickso/private/ind_seq/popgen/phlash/Kenya_2015_/SRR11077346.psmcfa' ])
