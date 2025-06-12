#!/usr/bin/env python3.10

import phlash
import os.path
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import sys
import os




def find_files(directory, pattern):
  """
  Lists files in a directory that match a given pattern.

  Args:
    directory: The directory to search in.
    pattern: The pattern to match (e.g., "*.txt", "image*").

  Returns:
    A list of file names that match the pattern.
  """
  matching_files = []
  return [filename for filename in os.listdir(directory) if pattern in filename]

if __name__ == '__main__':

    #get population from input
    population=sys.argv[1]
    #get path
    folder_path=f"/scratch/perickso/private/ind_seq/popgen/phlash/{population}"
    #change directory
    os.chdir(folder_path)

    #get list of file names
    files = find_files(folder_path, '.autosome.psmcfa')

    #run phlash
    posterior_samples = phlash.psmc(files)

    #summarize posteriors, use mutation rate from https://pmc.ncbi.nlm.nih.gov/articles/PMC3872194/
    times = np.array([dm.rescale(2.8e-9).eta.t[1:] for dm in posterior_samples])# choose a grid of points at which to evaluate the size history functions
    T = np.geomspace(times.min(), times.max(), 1000)
    Nes = np.array([dm.eta(T, Ne=True) for dm in posterior_samples])
    median=np.median(Nes, axis=0)
    lower_bound = np.percentile(Nes, 2.5, axis=0)
    upper_bound = np.percentile(Nes, 97.5, axis=0)

    #make a dataframe
    df = pd.DataFrame()
    df['time'] = T
    df['median'] = median
    df['lower2.5'] = lower_bound
    df['upper97.5'] = upper_bound

    #save file
    df.to_csv(f"/scratch/perickso/private/ind_seq/popgen/phlash/{population}.autosome.psmcfa.phlash.csv")

    #same posterior posterior_samples
    import pickle
    with open(f"/scratch/perickso/private/ind_seq/popgen/phlash/{population}.autosome.phlash.pscmfa.pickle", "wb") as file:
        pickle.dump(posterior_samples, file)
