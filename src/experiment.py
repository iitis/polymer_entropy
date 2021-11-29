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
        8: ColumnMeaning("ϕ₁₄ ", "deg"),
        9: ColumnMeaning("ψ₁₄ ", "deg"),
        10: ColumnMeaning("ϕ₁₃ ", "deg"),
        11: ColumnMeaning("ψ₁₃ ", "deg"),
    }

    columns_side = {
        0: ColumnMeaning("Time", "ps"),
        1: ColumnMeaning("Total energy of the system", "kJ/mol"),
        2: ColumnMeaning("Total bond energy", "kJ/mol"),
        3: ColumnMeaning("Total angle energy", "kJ/mol"),
        4: ColumnMeaning("Total dihedral energy", "kJ/mol"),
        5: ColumnMeaning("Total planar energy", "kJ/mol"),
        6: ColumnMeaning("Total Coulomb energy", "kJ/mol"),
        7: ColumnMeaning("Total Van der Waals energy", "kJ/mol"),
        8: ColumnMeaning("ϕ₁₄ ", "deg"),
        9: ColumnMeaning("ψ₁₄ ", "deg"),
        10: ColumnMeaning("ϕ₁₃ ", "deg"),
        11: ColumnMeaning("ψ₁₃ ", "deg"),
    }

    def __init__(self, filepath: str):
        self.dataframe = pd.read_csv(filepath, sep=";", skipfooter=2, engine="python")
        if "sidechain" in filepath:
            self.sidechain = True
        else:
            self.sidechain = False

    def correct_signs(self, x):
        """ if the sign of the angle jumps it is reversed"""
        η = 0.33
        x_cor = [el for el in x]

        for i in range(len(x_cor) - 1):
            if np.abs(x_cor[i] + x_cor[i + 1]) < η * np.abs(x_cor[i + 1]):

                x_cor[i + 1 : -1] = [-el for el in x_cor[i + 1 : -1]]

        if np.abs(x_cor[-1] + x_cor[-2]) < η * np.abs(x_cor[-1]):
            x_cor[-1] = -x_cor[-1]

        return x_cor

    def getColumnName(self, col: int):
        """ this is for analysis data """
        if self.sidechain:
            if col < 8:
                return self.columns[col].description, self.columns[col].unit
            else:
                offset = col % 4 + 8
                mer = col // 4 - 1
                desc = self.columns[offset].description
                desc += "mer " + str(mer) + " - mer " + str(mer + 1)
                unit = self.columns[offset].unit
                return desc, unit
        else:  # the rest will be changed
            if col < 8:
                return self.columns_side[col].description, self.columns_side[col].unit
            else:
                offset = col % 4 + 8
                mer = col // 4 - 1
                desc = self.columns_side[offset].description
                desc += "mer " + str(mer) + " - mer " + str(mer + 1)
                unit = self.columns_side[offset].unit
                return desc, unit

    def dropFirstObservations(self):
        """
        Drops first observations so that only points after stabilisation are
        further considered
        """
        self.dataframe.drop(index=self.dataframe.index[0], axis=0, inplace=True)

    def plotColumns(self, ycol: int, xcol: int = 0, plotname: str = None):
        """
        Plots one column vs another (time by default).
        Stores to png file if plot name is provided.
        """
        xname, xunit = self.getColumnName(xcol)
        yname, yunit = self.getColumnName(ycol)
        myxlab = f"{xname}[{xunit}]"
        myylab = f"{yname}[{yunit}]"
        mytitle = f"{xname} vs {yname}"

        y = self.dataframe.iloc[:, ycol]
        x = self.dataframe.iloc[:, xcol]

        y = self.correct_signs(y)

        plt.plot(x, y)
        plt.xlabel(myxlab)
        plt.ylabel(myylab)
        plt.title(mytitle)

        if plotname:
            plotFile = os.path.join("plots", f"{plotname}.png")
            plt.savefig(plotFile)
        else:
            plt.show()

    def plotHistogram2D(self, xcol: int, ycol: int, plotname: str = None):
        """
        Plots one column vs another 2D histogram
        """
        xname, xunit = self.getColumnName(xcol)
        yname, yunit = self.getColumnName(ycol)
        myxlab = f"{xname}[{xunit}]"
        myylab = f"{yname}[{yunit}]"
        mytitle = f"Histogram 2D"
        y = self.dataframe.iloc[:, ycol]
        y = self.correct_signs(y)

        x = self.dataframe.iloc[:, xcol]
        x = self.correct_signs(x)

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
        y = self.correct_signs(y)
        x = self.dataframe.iloc[:, xcol]
        x = self.correct_signs(x)

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
        xname, _ = myExperiment.getColumnName(xcol)
        yname, _ = myExperiment.getColumnName(ycol)
        self.x_axis = f"{xname}"
        self.y_axis = f"{yname}"

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
        self,
        chain_type: str,
        ion,
        xcol: int,
        ycol: int,
        plotname: str = None,
        no_mers: int = 24,
    ):

        """  compute percentiles of the histogram of entropies """

        first_mers = [mer for mer in range(no_mers - 1)]

        median_entropy = [0.0 for _ in range(no_mers - 1)]
        entropy_perc5 = [0.0 for _ in range(no_mers - 1)]
        entropy_perc95 = [0.0 for _ in range(no_mers - 1)]
        i = 0

        for mer in range(no_mers - 1):

            rest_of_path = "_" + chain_type + "_" + ion + ".tab"
            entropies = np.array(
                self.hist_of_entropy(rest_of_path, xcol + 4 * mer, ycol + 4 * mer)
            )
            median_entropy[i] = np.median(entropies)
            entropy_perc5[i] = np.percentile(entropies, 5)
            entropy_perc95[i] = np.percentile(entropies, 95)
            i += 1

        mytitle = f"file type {chain_type}, ion {ion}"
        myylabel = f"entropy{self.x_axis[0:4]}{self.y_axis[0:4]}"
        myxlabel = "first mer"

        plt.plot(first_mers, median_entropy, "o--", color="red", label="median")
        plt.plot(
            first_mers, entropy_perc5, ":", color="red", label="5 and 95 percentile"
        )
        plt.plot(first_mers, entropy_perc95, ":", color="red", label=None)
        plt.legend()
        plt.title(mytitle)
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        plt.show()

        return median_entropy
