"""
A module holding Experiment class that represent entire dataset
 from a file along with basic operations on it.
"""
from collections import namedtuple

import os
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
        8: ColumnMeaning("ϕ₁₄ ", "deg"),  # mer 1 - mer 2
        9: ColumnMeaning("ψ₁₄ ", "deg"),  # mer 1 - mer 2
        10: ColumnMeaning("ϕ₁₃ ", "deg"),  # mer 1 - mer 2
        11: ColumnMeaning("ψ₁₃ ", "deg"),  # mer 1 - mer 2
        # next 4  (12-15)  mer 2 - mer 3
        # up to mer 23 - mer 24
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
        8: ColumnMeaning("γ ", "deg"),  # mer 1 - mer 2
        # next mer 2 - mer 3
        # up to mer 23 - mer 24
        32: ColumnMeaning("ω ", "deg"),  # mer 1 - mer 2
        56: ColumnMeaning("δ", "deg"),  # mer 1 - mer 2
    }

    def __init__(self, filepath: str):
        self.dataframe = pd.read_csv(filepath, sep=";", skipfooter=2, engine="python")
        self.sidechain = "sidechain" in filepath


    def get_column_name(self, col: int):
        """ this is for analysis data """
        if not self.sidechain:
            if col < 8:
                return self.columns[col].description, self.columns[col].unit
            offset = col % 4 + 8
            mer = col // 4 - 1
            desc = self.columns[offset].description
            desc += f"mers {mer},{mer + 1}"
            unit = self.columns[offset].unit
            return desc, unit
        # this is the side chain
        if col < 8:
            desc = self.columns_side[col].description
            unit = self.columns_side[col].unit
            return desc, unit
        mer = (col - 8) % 24 + 1
        i = (col - 8) // 24
        i = i * 24 + 8
        desc = self.columns_side[i].description
        desc += f"mers {mer},{mer + 1}"
        unit = self.columns_side[i].unit
        return desc, unit

    def drop_first_observations(self):
        """
        Drops first observations so that only points after stabilisation are
        further considered
        """
        self.dataframe.drop(index=self.dataframe.index[0], axis=0, inplace=True)

    def plot_columns(self, ycol: int, plotname: str = None):
        """
        Plots one column vs another (time by default).
        Stores to pdf file if plot name is provided.
        """
        xcol = 0
        xname, xunit = self.get_column_name(xcol)
        yname, yunit = self.get_column_name(ycol)
        myxlab = f"{xname}  [{xunit}]"
        myylab = f"{yname}  [{yunit}]"
        if self.sidechain:
            mytitle = f"side chain {xname} vs {yname}"
        else:
            mytitle = f"{xname} vs {yname}"

        y = self.dataframe.iloc[:, ycol]
        x = self.dataframe.iloc[:, xcol]
        y = correct_signs(y)

        plt.plot(x, y)
        plt.xlabel(myxlab)
        plt.ylabel(myylab)
        plt.title(mytitle)

        if plotname:
            if self.sidechain:
                plotFile = f"{plotname}series_sidechain_{ycol}.pdf"
            else:
                plotFile = f"{plotname}series_{ycol}.pdf"
            plt.savefig(plotFile)
            plt.clf()
        else:
            plt.show()

    def plotHistogram2D(self, xcol: int, ycol: int, plotname: str = None):
        """
        Plots one column vs another 2D histogram
        """
        xname, xunit = self.get_column_name(xcol)
        yname, yunit = self.get_column_name(ycol)
        myxlab = f"{xname}   [{xunit}]"
        myylab = f"{yname}   [{yunit}]"
        if self.sidechain:
            mytitle = "side chain Histogram 2D"
        else:
            mytitle = "Histogram 2D"
        y = self.dataframe.iloc[:, ycol]
        y = correct_signs(y)

        x = self.dataframe.iloc[:, xcol]
        x = correct_signs(x)

        plt.subplots()
        plt.hist2d(x, y, bins=10, cmap=plt.cm.Reds)
        plt.colorbar(format=mtick.FormatStrFormatter("%.1e"))
        plt.xlabel(myxlab)
        plt.ylabel(myylab)
        plt.title(mytitle)
        if plotname:
            if self.sidechain:
                plotFile = f"{plotname}hist2D_side_{xcol}_{ycol}.pdf"
            else:
                plotFile = f"{plotname}hist2D_{xcol}_{ycol}.pdf"
            plt.savefig(plotFile)
            plt.clf()
        else:
            plt.show()

    def get_entropy(self, xcol: int, ycol: int):
        """
        computes entropy from histogram of xcol vs. ycol
        """
        self.drop_first_observations()
        y = self.dataframe.iloc[:, ycol]
        y = correct_signs(y)
        x = self.dataframe.iloc[:, xcol]
        x = correct_signs(x)

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
    magic_numbers = {'analysis': 4,
                     'sidechain' : 1}

    def __init__(self, partial_path: str, no_experiments: int = 12):
        self.partial_path = partial_path
        self.no_experiments = no_experiments

    def set_axis_description(self, myExperiment, xcol: int, ycol: int):
        """ this function will let us know which angle data
         we are dealing with """
        self.x_axis, _ = myExperiment.get_column_name(xcol)
        self.y_axis, _ = myExperiment.get_column_name(ycol)

    def hist_of_entropy(self, chain_type: str, ion: str, xcol: int, ycol: int, plotdir: str = None):
        """ compute histogram of entropy over realisations """

        entropies = []

        for i in range(self.no_experiments):
            path = f"{self.partial_path}{i+1}_{chain_type}_{ion}.tab"
            myExperiment = Experiment(path)
            self.set_axis_description(myExperiment, xcol, ycol) #FIXME sets multiple times
            entropies.append(myExperiment.get_entropy(xcol, ycol))

        if plotdir:
            xdesc = self.x_axis
            ydesc = self.y_axis
            plotFile = f"{plotdir}hist{xdesc[0:self.magic_numbers[chain_type]]}_{ydesc[0:self.magic_numbers[chain_type]]}.pdf"
            mytitle = f"{chain_type} {ion}"
            myxlabel = f"entropy {xdesc} vs. {ydesc}"
            myylabel = "frequency"

            fig, ax = plt.subplots()
            plt.hist(entropies, bins=5)
            plt.title(mytitle)
            plt.xlabel(myxlabel)
            plt.ylabel(myylabel)

            plt.savefig(plotFile)
            plt.clf()
        return entropies

    def entropy_distribution_percentiles(self, chain_type: str, ion, xcol: int, ycol: int, plotdir: str, no_mers: int = 23):
        """  compute percentiles of the histogram of entropies """
        first_mers = list(range(1, no_mers+1))

        median_entropy = []
        entropy_perc5 = []
        entropy_perc95 = []

        for mer in range(no_mers):
            entropies = np.array(self.hist_of_entropy(chain_type, ion, xcol + self.magic_numbers[chain_type] * mer, ycol + self.magic_numbers[chain_type] * mer))
            median_entropy.append(np.median(entropies))
            entropy_perc5.append(np.percentile(entropies, 5))
            entropy_perc95.append(np.percentile(entropies, 95))

        xdesc = self.x_axis[0:self.magic_numbers[chain_type]]
        ydesc = self.y_axis[0:self.magic_numbers[chain_type]]

        mytitle = f"{chain_type}, ion {ion}"
        myylabel = f"entropy  {xdesc} vs. {ydesc}"
        myxlabel = "number of first mer"

        plt.plot(first_mers, median_entropy, "o--", color="red", label="median")
        plt.plot(first_mers, entropy_perc5, ":", color="red", label="5 and 95 perc.")
        plt.plot(first_mers, entropy_perc95, ":", color="red", label=None)
        plt.legend()
        plt.title(mytitle)
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        plotFile = os.path.join(plotdir, f"entropy_{chain_type}_{ion}_{xdesc}_{ydesc}.pdf")
        plt.savefig(plotFile)
        plt.clf()

        return median_entropy

    def entropy_distribution_realisations(self, chain_type: str, ion, xcol: int, ycol: int, plotdir: str, no_mers: int = 23):
        """  compute percentiles of the histogram of entropies """
        no_struct = 12
        first_mers = list(range(1, no_mers+1))
        print("highest first mer", first_mers[-1])
        entropies = [[0.0 for i in range(no_struct)] for _ in range(no_mers)]

        for mer in range(no_mers):
            entropies[mer] = np.array(self.hist_of_entropy(chain_type, ion, xcol + self.magic_numbers[chain_type] * mer, ycol + self.magic_numbers[chain_type] * mer))

        xdesc = self.x_axis[0:self.magic_numbers[chain_type]]
        ydesc = self.y_axis[0:self.magic_numbers[chain_type]]

        mytitle = f"{chain_type}, ion {ion}"
        myylabel = f"entropy  {xdesc} vs. {ydesc}"
        myxlabel = "number of first mer"

        for j in range(no_struct):
            color = 5 / 6 * (1 - j / no_struct) * np.array((1, 1, 1))
            e = [entropies[i][j] for i in range(no_mers)]
            plt.plot(first_mers, e, label=f"{j+1}", color=color)

        plt.legend()
        plt.title(mytitle)
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        plotFile = os.path.join(plotdir, f"entropy_reals_{chain_type}_{ion}_{xdesc}_{ydesc}.pdf")
        plt.savefig(plotFile)
        plt.clf()

def correct_signs(series, thres: float = 0.5):
    """ if the sign of the angle is artificially reversed,
    we reverse it back"""
    series_corrected = np.array(series)

    for i in range(len(series_corrected) - 1):
        if np.abs(series_corrected[i] + series_corrected[i + 1]) < thres * np.abs(series_corrected[i + 1]):
            series_corrected[i + 1:] *= -1
    return series_corrected
