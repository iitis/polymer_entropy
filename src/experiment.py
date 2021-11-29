"""
A module holding Experiment class that represent entire dataset from a file along with basic operations on it.
"""
from collections import namedtuple
import os
from statistics import mean
from functools import lru_cache
from  scipy import stats
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from scipy.stats import entropy


class Experiment:
    """
    A class that represent an experiment - mucin behaviour when submerged to a salt solution.
    """
    ColumnMeaning = namedtuple('ColumnMeaning', 'description unit')
    columns = {
        0: ColumnMeaning('Time', 'ps'),
        1: ColumnMeaning('Total energy of the system', 'kJ/mol'),
        2: ColumnMeaning('Total bond energy', 'kJ/mol'),
        3: ColumnMeaning('Total angle energy', 'kJ/mol'),
        4: ColumnMeaning('Total dihedral energy', 'kJ/mol'),
        5: ColumnMeaning('Total planar energy', 'kJ/mol'),
        6: ColumnMeaning('Total Coulomb energy', 'kJ/mol'),
        7: ColumnMeaning('Total Van der Waals energy', 'kJ/mol'),
        8: ColumnMeaning('Solvent accessible surface of part A of albumin', 'A^2'),
        9: ColumnMeaning('Solvent accessible surface of part B of albumin', 'A^2'),
        10: ColumnMeaning('Solvent accessible surface of CS', 'A^2'),
        11: ColumnMeaning('Total number of hydrogen bonds', ''),
        12: ColumnMeaning('Total number of hydrogen bonds in part A of albumin', ''),
        13: ColumnMeaning('Total number of hydrogen bonds in part B of albumin', ''),
        14: ColumnMeaning('Total number of hydrogen bonds in CS', ''),
        15: ColumnMeaning('Total number of hydrogen bonds between A and B in albumin', ''),
        16: ColumnMeaning('Total number of hydrogen bonds between A in albumin and CS', ''),
        17: ColumnMeaning('Total number of hydrogen bonds between B in albumin and CS', ''),
        18: ColumnMeaning('RMSD CA', ''),
        19: ColumnMeaning('RMSDBackbone', ''),
        20: ColumnMeaning('RMSD HeavyAtoms', ''),
    }

    def __init__(self, filepath: str):
        self.dataframe = pd.read_csv(filepath, sep=';', skipfooter=2, engine='python')
        self.maxTau = 200

    def dropFirstObservations(self):
        """
        Drops first observations so that only points after stabilisation are further considered"
        """
        self.dataframe.drop(index=self.dataframe.index[0],axis=0,inplace=True)

    def plotColumns(self, ycol: int, xcol: int = 0, plotname: str = None):
        """
        Plots one column vs another (time by default).
        Stores to png file if plotname is provided.
        """
        myxlabel = f"{self.columns[xcol].description}[{self.columns[xcol].unit}]"
        myylabel = f"{self.columns[ycol].description}[{self.columns[ycol].unit}]"
        mytitle = f"{self.columns[ycol].description} vs {self.columns[xcol].description}"
        avg_value = self.dataframe.iloc[:, ycol].mean()
        if xcol == 0: #make time axis start at zero
            self.dataframe.plot(x=xcol, y=ycol, title=mytitle, xlabel=myxlabel, ylabel=myylabel, xlim=(0,100000), legend=False )
        else:
            self.dataframe.plot(x=xcol, y=ycol, title=mytitle, xlabel=myxlabel, ylabel=myylabel, legend=False)

        if plotname:
            plotFile = os.path.join('plots', f"{plotname}.png")
            plt.savefig(plotFile)
        else:
            plt.show()


    def plotHistogram2D(self, xcol: int, ycol: int, plotname: str = None):
        """
        Plots one column vs another 2D histogram
        """
        myxlabel = f"{self.columns[xcol].description}[{self.columns[xcol].unit}]"
        myylabel = f"{self.columns[ycol].description}[{self.columns[ycol].unit}]"
        mytitle = f"Histogram 2D"
        y = self.dataframe.iloc[:, ycol]
        x = self.dataframe.iloc[:, xcol]


        fig, ax = plt.subplots()
        plt.hist2d(x, y, bins=10, cmap=plt.cm.Reds)
        plt.colorbar(format = mtick.FormatStrFormatter('%.1e'))
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        if plotname:
            plotFile = os.path.join('plots', f"{plotname}.png")
            plt.savefig(plotFile)
        else:
            plt.show()

    def get_entropy(self, xcol: int, ycol: int):
        """
        computes entropy from histogram of xcol vs. ycol
        """
        y = self.dataframe.iloc[:, ycol]
        x = self.dataframe.iloc[:, xcol]

        h = np.histogram2d(x, y, bins=10)
        h_vec = np.concatenate(h[0])
        h_norm = h_vec/sum(h_vec)
        # use molar gas constant R = 8.314
        return 8.314*entropy(h_norm )
