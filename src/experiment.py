"""
A module holding Experiment class that represent entire dataset
 from a file along with basic operations on it.
"""
import os
from collections import namedtuple

import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import entropy


class Experiment:
    """
    A class that represent an experiment of entropy calculation
    from angles from the single file
    """

    ColumnMeaning = namedtuple("ColumnMeaning", "description unit")
    columns = {
        0: ColumnMeaning("Time", "ps"),
        1: ColumnMeaning("Total energy of the system", "kJ/mol"),
        2: ColumnMeaning("Total bond energy", "kJ/mol"),
        3: ColumnMeaning("Total angle energy", "kJ/mol"),
        4: ColumnMeaning("Total dihedral energy", "kJ/mol"),
        5: ColumnMeaning("Total planar energy", "kJ/mol"),
        6: ColumnMeaning("Total Coulomb energy", "kJ/mol"),
        7: ColumnMeaning("Total Van der Waals energy", "kJ/mol"),
        8: ColumnMeaning("Solvent accessible surface of part A of albumin", "A^2"),
        9: ColumnMeaning("Solvent accessible surface of part B of albumin", "A^2"),
        10: ColumnMeaning("Solvent accessible surface of CS", "A^2"),
        11: ColumnMeaning("Total number of hydrogen bonds", ""),
        12: ColumnMeaning("Total number of hydrogen bonds in part A of albumin", ""),
        13: ColumnMeaning("Total number of hydrogen bonds in part B of albumin", ""),
        14: ColumnMeaning("Total number of hydrogen bonds in CS", ""),
        15: ColumnMeaning(
            "Total number of hydrogen bonds between A and B in albumin", ""
        ),
        16: ColumnMeaning(
            "Total number of hydrogen bonds between A in albumin and CS", ""
        ),
        17: ColumnMeaning(
            "Total number of hydrogen bonds between B in albumin and CS", ""
        ),
        18: ColumnMeaning("RMSD CA", ""),
        19: ColumnMeaning("RMSDBackbone", ""),
        20: ColumnMeaning("RMSD HeavyAtoms", ""),
    }

    def __init__(self, filepath: str):
        self.dataframe = pd.read_csv(filepath, sep=";", skipfooter=2,
                                     engine="python")

    def dropFirstObservations(self):
        """
        Drops first observations so that only points after stabilisation are
        further considered
        """
        self.dataframe.drop(index=self.dataframe.index[0], axis=0,
                            inplace=True)

    def plotColumns(self, ycol: int, xcol: int = 0, plotname: str = None):
        """
        Plots one column vs another (time by default).
        Stores to png file if plot name is provided.
        """
        myxlab = f"{self.columns[xcol].description}[{self.columns[xcol].unit}]"
        myylab = f"{self.columns[ycol].description}[{self.columns[ycol].unit}]"
        mytitle = (f"{self.columns[ycol].description} vs {self.columns[xcol].description}")

        if xcol == 0:  # make time axis start at zero
            self.dataframe.plot(
                x=xcol,
                y=ycol,
                title=mytitle,
                xlabel=myxlab,
                ylabel=myylab,
                xlim=(0, 100000),
                legend=False,
            )
        else:
            self.dataframe.plot(
                x=xcol,
                y=ycol,
                title=mytitle,
                xlabel=myxlab,
                ylabel=myylab,
                legend=False,
            )

        if plotname:
            plotFile = os.path.join("plots", f"{plotname}.png")
            plt.savefig(plotFile)
        else:
            plt.show()

    def plotHistogram2D(self, xcol: int, ycol: int, plotname: str = None):
        """
        Plots one column vs another 2D histogram
        """
        myxlab = f"{self.columns[xcol].description}[{self.columns[xcol].unit}]"
        myylab = f"{self.columns[ycol].description}[{self.columns[ycol].unit}]"
        mytitle = f"Histogram 2D"
        y = self.dataframe.iloc[:, ycol]
        x = self.dataframe.iloc[:, xcol]

        fig, ax = plt.subplots()
        plt.hist2d(x, y, bins=10, cmap=plt.cm.Reds)
        plt.colorbar(format=mtick.FormatStrFormatter("%.1e"))
        plt.xlabel(myxlab)
        plt.ylabel(myylab)
        plt.title(mytitle)
        if plotname:
            plotFile = os.path.join("plots", f"{plotname}.png")
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
        h_norm = h_vec / sum(h_vec)
        # use molar gas constant R = 8.314
        return 8.314 * entropy(h_norm)


class SetOfExperiments:

    """
    this class represents a set of experiments of entropy calculation
     performed in a series of files for statistics
    """

    def __init__(self, partial_path: str, no_experiments: int = 12):
        self.partial_path = partial_path
        self.no_experiments = no_experiments

    def axis_descrtiption(self, myExperiment, xcol: int, ycol: int):
        """ this function will let us know which angle data
         we are dealing with """
        self.x_axis = f"{myExperiment.columns[xcol].description}"
        self.y_axis = f"{myExperiment.columns[ycol].description}"

    def hist_of_entropy(
        self, rest_of_path: str, xcol: int, ycol: int, plotname: str = None
    ):
        """ compute histogram of entropy over realisations """

        entropies = [0.0 for _ in range(self.no_experiments)]
        for i in range(self.no_experiments):
            path = self.partial_path
            path += str(i + 1)
            path += rest_of_path
            myExperiment = Experiment(path)
            self.axis_descrtiption(myExperiment, xcol, ycol)
            entropies[i] = myExperiment.get_entropy(xcol, ycol)

        if plotname:
            mytitle = f"{rest_of_path}"
            myxlabel = f"entropy{self.x_axis}{self.y_axis}"
            myylabel = "frequency"

            fig, ax = plt.subplots()
            plt.hist(entropies, bins=5)
            plt.title(mytitle[1:-4])
            plt.xlabel(myxlabel)
            plt.ylabel(myylabel)
            plt.show()
        return entropies

    def entropy_distribution_percentiles(
        self, chain_type: str, ions, xcol: int, ycol: int, plotname: str = None
    ):

        """  compute percentiles of the histogram of entropies """

        median_entropy = [0.0 for _ in ions]
        entropy_perc5 = [0.0 for _ in ions]
        entropy_perc95 = [0.0 for _ in ions]
        i = 0

        for ion in ions:

            rest_of_path = "_" + chain_type + "_" + ion + ".tab"
            entropies = np.array(self.hist_of_entropy(rest_of_path, xcol, ycol))
            median_entropy[i] = np.median(entropies)
            entropy_perc5[i] = np.percentile(entropies, 5)
            entropy_perc95[i] = np.percentile(entropies, 95)
            i += 1

        mytitle = f"{chain_type}"
        myylabel = f"entropy{self.x_axis}{self.y_axis}"
        myxlabel = "ions"

        plt.plot(ions, median_entropy, "o--", color="red", label="median")
        plt.plot(ions, entropy_perc5, ":", color="red",
                 label="5 and 95 percentile")
        plt.plot(ions, entropy_perc95, ":", color="red", label=None)
        plt.legend()
        plt.title(mytitle)
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        plt.show()

        return median_entropy
