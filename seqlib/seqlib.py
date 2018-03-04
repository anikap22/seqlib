#!/usr/bin/env python

"""
seqlib library for class assignment
"""

import copy
import numpy as np
import pandas as pd


class Seqlib:

    def __init__(self, ninds, nsites):
        self.ninds = ninds
        self.nsites = nsites
        self.seqs = self._simulate()
        # ...

        ## store maf of the full seq array
        self.maf = self._get_maf()

    def _mutate(self, base):
        diff = set("ACTG") - set(base)
        return np.random.choice(list(diff))

    # define simulate function to simulate sequences and mutations
    # define simulate function to simulate sequences and mutations
    def _simulate(self):
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

    def _get_maf(self):
        "returns the maf of the full seqarray while not counting Ns"
        ## init an array to fill and iterate over columns
        maf = np.zeros(self.nsites)
        for col in range(self.nsites):
            ## select this column of bases
            thiscol = self.seqs[:, col]

            ## mask "N" bases and get new length
            nmask = thiscol != "N"
            no_n_len = np.sum(nmask)

            ## mask "N" bases and get the first base
            first_non_n_base = thiscol[nmask][0]

            ## calculate maf of "N" masked bases
            freq = np.sum(thiscol[nmask] != first_non_n_base) / no_n_len
            if freq > 0.5:
                maf[col] = 1 - freq
            else:
                maf[col] = freq
        return maf

    # define function to remove columns with a missing value
    def _filter_missing(self, maxmissing):
        # calculate the fraction of rows with a missing base at that site
        freqmissing = np.sum(self.seqs == "N", axis=0) / self.seqs.shape[0]
        return freqmissing > maxmissing

    # define function that takes array of sequences and its mutations through t
    def _filter_maf(self, minmaf):
        return self.maf < minmaf

    # return filtered by missing and minfreq
    def filter(self, minmaf, maxmissing):
        filter1 = self._filter_maf(minmaf)
        filter2 = self._filter_missing(maxmissing)
        fullfilter = filter1 + filter2
        return self.seqs[:, np.invert(fullfilter)]

    def filter_seqlib(self, minmaf, maxmissing):
        """
        Applies maf and missing filters to the array and returns a copy 
        of the seqlib object where the .seqs array has been filtered

        Parameters
        ----------
        minmaf: float
            The minimum minor allele frequency. Filter columns below this.
        maxmissing: float
            The maximum prop. missing data. Filter columns with prop Ns > this.
        """
        ## apply filters to get new array size
        newseqs = self.filter(minmaf, maxmissing)

        ## make a new copy of the seqlib object
        newself = copy.deepcopy(self)       
        newself.__init__(newseqs.shape[0], newseqs.shape[1]) 

        ## store the array (overwrite it)
        newself.seqs = newseqs

        ## call the _get_maf to match new array
        newself._get_maf()
        return newself

    def calculate_statistics(self):
        if self.seqs.size:
            nd = np.var(self.seqs == self.seqs[0], axis=0).mean()
            mf = np.mean(
                np.sum(self.seqs != self.seqs[0], axis=0) / self.seqs.shape[0])
            inv = np.all(self.seqs == self.seqs[0], axis=0).sum()
            var = self.seqs.shape[1] - inv
            return pd.Series(
                {"mean nucleotide diversity": nd,
                 "mean minor allele frequency": mf,
                 "invariant sites": inv,
                 "variable sites": var,
                })
        else:
            print("seqs array is empty")
