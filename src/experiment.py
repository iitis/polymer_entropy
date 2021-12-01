"""
A module holding Experiment class that represent entire dataset
 from a file along with basic operations on it.
"""
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
        self.dataframe = pd.read_csv(filepath, sep=";", skipfooter=2,
                                     engine="python")
        if "sidechain" in filepath:
            self.sidechain = True
        else:
            self.sidechain = False

    def correct_signs(self, x, thres: float = 0.5):
        """ if the sign of the angle is artificially reversed,
        we reverse it back"""
        x_cor = [el for el in x]

        for i in range(len(x_cor) - 1):
            if np.abs(x_cor[i] + x_cor[i + 1]) < thres * np.abs(x_cor[i + 1]):

                x_cor[i + 1: -1] = [-el for el in x_cor[i + 1: -1]]

        if np.abs(x_cor[-1] + x_cor[-2]) < thres * np.abs(x_cor[-1]):
            x_cor[-1] = -x_cor[-1]

        return x_cor

    def getColumnName(self, col: int):
        """ this is for analysis data """
        if not self.sidechain:
            if col < 8:
                return self.columns[col].description, self.columns[col].unit
            else:
                offset = col % 4 + 8
                mer = col // 4 - 1
                desc = self.columns[offset].description
                desc += f"mers {mer},{mer + 1}"
                unit = self.columns[offset].unit
                return desc, unit
        else:  # this is the mian chain
            if col < 8:
                desc = self.columns_side[col].description
                un = self.columns_side[col].unit
                return desc, un
            else:
                mer = (col - 8) % 24 + 1
                i = (col - 8) // 24
                i = i * 24 + 8
                desc = self.columns_side[i].description
                desc += f"mers {mer},{mer + 1}"
                unit = self.columns_side[i].unit
                return desc, unit

    def dropFirstObservations(self):
        """
        Drops first observations so that only points after stabilisation are
        further considered
        """
        self.dataframe.drop(index=self.dataframe.index[0], axis=0,
                            inplace=True)

    def plotColumns(self, ycol: int, plotname: str = None):

        xcol = 0
        """
        Plots one column vs another (time by default).
        Stores to pdf file if plot name is provided.
        """
        xname, xunit = self.getColumnName(xcol)
        yname, yunit = self.getColumnName(ycol)
        myxlab = f"{xname}  [{xunit}]"
        myylab = f"{yname}  [{yunit}]"
        if self.sidechain:
            mytitle = f"side chain {xname} vs {yname}"
        else:
            mytitle = f"{xname} vs {yname}"

        y = self.dataframe.iloc[:, ycol]
        x = self.dataframe.iloc[:, xcol]
        y = self.correct_signs(y)

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
        xname, xunit = self.getColumnName(xcol)
        yname, yunit = self.getColumnName(ycol)
        myxlab = f"{xname}   [{xunit}]"
        myylab = f"{yname}   [{yunit}]"
        if self.sidechain:
            mytitle = "side chain Histogram 2D"
        else:
            mytitle = "Histogram 2D"
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
        self.dropFirstObservations()
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
        self, rest_of_path: str, xcol: int, ycol: int, plotdir: str = None
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

        if plotdir:
            xdesc = self.x_axis
            ydesc = self.y_axis
            if myExperiment.sidechain:
                plotFile = f"{plotdir}hist{xdesc[0:1]}_{ydesc[0:1]}.pdf"
            else:
                plotFile = f"{plotdir}hist{xdesc[0:4]}_{ydesc[0:4]}.pdf"
            mytitle = f"{rest_of_path}"
            myxlabel = f"entropy {xdesc} vs. {ydesc}"
            myylabel = "frequency"

            fig, ax = plt.subplots()
            plt.hist(entropies, bins=5)
            plt.title(mytitle[1:-4])
            plt.xlabel(myxlabel)
            plt.ylabel(myylabel)

            plt.savefig(plotFile)
            plt.clf()
        return entropies

    def entropy_distribution_percentiles(
        self,
        chain_type: str,
        ion,
        xcol: int,
        ycol: int,
        plotdir: str,
        no_mers: int = 23,
    ):

        """  compute percentiles of the histogram of entropies """

        first_mers = [mer + 1 for mer in range(no_mers)]

        median_entropy = [0.0 for _ in range(no_mers)]
        entropy_perc5 = [0.0 for _ in range(no_mers)]
        entropy_perc95 = [0.0 for _ in range(no_mers)]

        for mer in range(no_mers):

            rest_of_path = "_" + chain_type + "_" + ion + ".tab"
            if chain_type == "analysis":
                entropies = np.array(
                    self.hist_of_entropy(rest_of_path,
                                         xcol + 4 * mer, ycol + 4 * mer
                                         )
                )
            elif chain_type == "sidechain":
                entropies = np.array(
                    self.hist_of_entropy(rest_of_path, xcol + mer, ycol + mer)
                )
            median_entropy[mer] = np.median(entropies)
            entropy_perc5[mer] = np.percentile(entropies, 5)
            entropy_perc95[mer] = np.percentile(entropies, 95)

        if chain_type == "sidechain":
            xdesc = self.x_axis[0:1]
            ydesc = self.y_axis[0:1]
        else:
            xdesc = self.x_axis[0:4]
            ydesc = self.y_axis[0:4]

        mytitle = f"{chain_type}, ion {ion}"
        myylabel = f"entropy  {xdesc} vs. {ydesc}"
        myxlabel = "number of first mer"

        plt.plot(
            first_mers, median_entropy, "o--", color="red", label="median"
        )
        plt.plot(
            first_mers, entropy_perc5, ":", color="red", label="5 and 95 perc."
        )
        plt.plot(first_mers, entropy_perc95, ":", color="red", label=None)
        plt.legend()
        plt.title(mytitle)
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        plotFile = f"{plotdir}entropy_{chain_type}_{ion}_{xdesc}_{ydesc}.pdf"
        plt.savefig(plotFile)
        plt.clf()

        return median_entropy

    def entropy_distribution_realisations(
        self,
        chain_type: str,
        ion,
        xcol: int,
        ycol: int,
        plotdir: str,
        no_mers: int = 23,
    ):

        """  compute percentiles of the histogram of entropies """

        no_struct = 12

        first_mers = [mer + 1 for mer in range(no_mers)]

        print("highest first mer", first_mers[-1])

        entropies = [[0.0 for i in range(no_struct)] for _ in range(no_mers)]
        # i = 0

        for mer in range(no_mers):

            rest_of_path = "_" + chain_type + "_" + ion + ".tab"

            if chain_type == "analysis":
                entropies[mer] = np.array(
                    self.hist_of_entropy(rest_of_path,
                                         xcol + 4 * mer, ycol + 4 * mer
                                         )
                )
            elif chain_type == "sidechain":
                entropies[mer] = np.array(
                    self.hist_of_entropy(rest_of_path, xcol + mer, ycol + mer)
                )

        if chain_type == "sidechain":
            xdesc = self.x_axis[0:1]
            ydesc = self.y_axis[0:1]
        else:
            xdesc = self.x_axis[0:4]
            ydesc = self.y_axis[0:4]

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
        plotFile = (
            f"{plotdir}entropy_reals_{chain_type}_{ion}_{xdesc}_{ydesc}.pdf"
        )
        plt.savefig(plotFile)
        plt.clf()
