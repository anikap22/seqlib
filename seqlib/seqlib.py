import numpy as np
import pandas as pd


class Seqlib:

    def __init__(self, ninds, nsites):
        self.ninds = ninds
        self.nsites = nsites
        self.seqs = self.simulate()
        # ...

    def mutate(self, base):
        diff = set("ACTG") - set(base)
        return np.random.choice(list(diff))

    # define simulate function to simulate sequences and mutations
    # define simulate function to simulate sequences and mutations
    def simulate(self):
        ninds = self.ninds
        nsites = self.nsites
        # create initial random sequence of length nsites
        oseq = np.random.choice(list("ACGT"), size=nsites)
        # define an array copying the sequences ninds times (each as a new row)
        arr = np.array([oseq for i in range(ninds)])
        # define whether any given site has a mutation using a binomial dist
        muts = np.random.binomial(1, 0.1, (ninds, nsites))
        # iterate through sites and insert newbase (randomly chosen other
        for col in range(nsites):
            newbase = self.mutate(arr[0, col])
            mask = muts[:, col].astype(bool)
            arr[:, col][mask] = newbase
        # make array from binomial distribution for sites to be missing
        missing = np.random.binomial(1, 0.1, (ninds, nsites))
        # replace value with N in sites to be missing
        arr[missing.astype(bool)] = "N"
        return arr

    # define function to remove columns with a missing value
    def filter_missing(self, maxfreq):
        arr = self.simulate()
        # calculate the fraction of rows with a missing base at that site
        freqmissing = np.sum(arr == "N", axis=0) / arr.shape[0]
        # return the array excluding the columns that had more missing than thr
        return arr[:, freqmissing <= maxfreq]

    # define function that takes array of sequences and its mutations through t
    def filter_maf(self, arr, minfreq):
        # arr = self.simulate()
        # calculates mutation freq (number of rows in a column that are not
        # equal to the first row over column length)
        freqs = np.sum(arr != arr[0], axis=0) / arr.shape[0]
        # make a copy of the frequency sequence
        maf = freqs.copy()
        # select sites where the freq is over .5 and replace with 1-that value
        maf[maf > 0.5] = 1 - maf[maf > 0.5]
        # return the array with all rows and columns where the maf value is
        return arr[:, maf > minfreq]

    # return maf as floats
    def maf(self):
        arr = self.simulate()
        # calculates mutation freq
        freqs = np.sum(arr != arr[0], axis=0) / arr.shape[0]
        # make a copy of the frequency sequence
        maf = freqs.copy()
        # select sites where the freq is over .5 and replace with 1-that value
        maf[maf > 0.5] = 1 - maf[maf > 0.5]
        return freqs

    # return filtered by missing and minfreq
    def filter(self, minfreq, maxfreq):
        maf = self.filter_maf(self.filter_missing(maxfreq), minfreq)
        return maf

    def calculate_statistics(self):
        arr = self.simulate()
        nd = np.var(arr == arr[0], axis=0).mean()
        mf = np.mean(np.sum(arr != arr[0], axis=0) / arr.shape[0])
        inv = np.any(arr != arr[0], axis=0).sum()
        var = arr.shape[1] - inv
        return pd.Series(
            {"mean nucleotide diversity": nd,
             "mean minor allele frequency": mf,
             "invariant sites": inv,
             "variable sites": var,
             })
