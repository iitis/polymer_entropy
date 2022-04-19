"""
A module holding Experiment class that represent entire dataset
 from a file along with basic operations on it.
"""
import os
import glob
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import entropy

class ColumnMeaning:
    """Describes contents of a column"""
    def __init__(self, description: str, unit: str):
        self.description = description
        self.unit = unit

    def __str__(self):
        return f"{self.description} [{self.unit}]"

class Experiment:
    """
    A class that represent an experiment of entropy calculation
    from angles from the single file
    """

    plot_dpi = 350

    angles_config = {
        ('Albumin+HA','main chain'): ["ϕ₁₄","ψ₁₄","ϕ₁₃","ψ₁₃"],
        ('Albumin+HA','side chain'): ["γ","ω","δ"],
        ('Albumin+CS6','main chain'): ["ϕ₁₄","ψ₁₄","ϕ₁₃","ψ₁₃"],
    }

    initial_columns = [
        ColumnMeaning("Time", "ps"),                            #First column, index 0
        ColumnMeaning("Total energy of the system", "kJ/mol"),
        ColumnMeaning("Total bond energy", "kJ/mol"),
        ColumnMeaning("Total angle energy", "kJ/mol"),
        ColumnMeaning("Total dihedral energy", "kJ/mol"),
        ColumnMeaning("Total planar energy", "kJ/mol"),
        ColumnMeaning("Total Coulomb energy", "kJ/mol"),
        ColumnMeaning("Total Van der Waals energy", "kJ/mol"),  #Column index 7
    ]

    def __init__(self, filepath: str):
        self.complex = os.path.basename(filepath).split('_')[0]
        self.dataframe = pd.read_csv(filepath, sep=";", skipfooter=2, engine="python")
        self.chain = 'side chain' if 'sidechain' in filepath else 'main chain'
        self.angles = self.angles_config[(self.complex,self.chain)]
        self.no_mers = 23
        self.columns = self.initial_columns

        if self.complex == 'Albumin+HA':
            if self.chain == 'main chain':
                for mer in range(self.no_mers):
                    for angle in self.angles:
                        self.columns.append(ColumnMeaning(f"{angle} mers {mer+1}, {mer+2}", "deg"))
            elif self.chain == 'side chain':
                for angle in self.angles:
                    for mer in range(self.no_mers):
                        self.columns.append(ColumnMeaning(f"{angle} mers {mer+1}, {mer+2}", "deg"))

        if self.complex == 'Albumin+CS6':
            if self.chain == 'main chain':
                for mer in range(self.no_mers):
                    for angle in self.angles:
                        self.columns.append(ColumnMeaning(f"{angle} mers {mer+1}, {mer+2}", "deg"))
            elif self.chain == 'side chain':
                pass #not supported

    def get_colnum_by_meaning(self, meaning: str):
        """
        Returns column index when provided with column meaning
        """
        cols = [ x.description for x in self.columns]
        return cols.index(meaning)

    def drop_first_observations(self):
        """
        Drops first observations so that only points after stabilisation are
        further considered
        """
        self.dataframe.drop(index=self.dataframe.index[0], axis=0, inplace=True)

    def plot_columns(self, ycol: int, plotname: str = None):
        """
        Plots one column vs another (time by default).
        Stores to png file if plot name is provided.
        """
        xcol = 0
        xname = self.columns[xcol].description
        yname = self.columns[ycol].description
        mytitle = f"{self.chain} {xname} vs {yname}"

        y = self.dataframe.iloc[:, ycol]
        x = self.dataframe.iloc[:, xcol]
        y = correct_signs(y)

        plt.plot(x, y)
        plt.xlabel(self.columns[xcol])
        plt.ylabel(self.columns[ycol])
        plt.title(mytitle)

        if plotname:
            plot_filepath = f"{plotname}series_{self.chain.replace(' ','')}_{ycol}.png"
            plt.savefig(plot_filepath, dpi=self.plot_dpi)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_histogram_2d(self, xcolumn: str, ycolumn: str, plotname: str = None, numbins = 100):
        """
        Plots one column vs another 2D histogram
        """
        xcol = self.get_colnum_by_meaning(xcolumn)
        ycol = self.get_colnum_by_meaning(ycolumn)
        mytitle = f"{self.chain} Histogram 2D"
        y = self.dataframe.iloc[:, ycol]
        y = correct_signs(y)

        x = self.dataframe.iloc[:, xcol]
        x = correct_signs(x)

        plt.subplots()
        plt.hist2d(x, y, bins=numbins, range=[[-180,180],[-180,180]], cmap=plt.cm.Reds)
        plt.colorbar(format=mtick.FormatStrFormatter("%.1e"))
        plt.xlabel(self.columns[xcol])
        plt.ylabel(self.columns[ycol])
        plt.title(mytitle)
        if plotname:
            plot_filepath = f"{plotname}hist2D_{self.chain.replace(' ','')}_{xcol}_{ycol}_b{numbins}.png"
            plt.savefig(plot_filepath, dpi=self.plot_dpi)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_angle_histogram(self, angle_x, angle_y, plotname: str = None, numbins = 100):
        """
        Plots histogram fo angles for all mers in experiment (realisation)
        """
        x_data = np.array([])
        y_data = np.array([])

        x_columns_of_interest = [ f"{angle_x} mers {mer+1}, {mer+2}" for mer in range(self.no_mers) ]
        y_columns_of_interest = [ f"{angle_y} mers {mer+1}, {mer+2}" for mer in range(self.no_mers) ]

        x_cols = [ self.get_colnum_by_meaning(x) for x in x_columns_of_interest ]
        y_cols = [ self.get_colnum_by_meaning(x) for x in y_columns_of_interest ]

        for c in x_cols:
            x = self.dataframe.iloc[:, c]
            x = correct_signs(x)
            x_data = np.concatenate((x_data, x))

        for c in y_cols:
            y = self.dataframe.iloc[:, c]
            y = correct_signs(y)
            y_data = np.concatenate((y_data, y))

        plt.subplots()
        plt.hist2d(x, y, bins=numbins, range=[[-180,180],[-180,180]], cmap=plt.cm.Reds)
        plt.colorbar(format=mtick.FormatStrFormatter("%.1e"))
        plt.xlabel(angle_x)
        plt.ylabel(angle_y)
        plt.title("All mers angles")
        if plotname:
            plot_filepath = f"{plotname}hist2D_{self.chain.replace(' ','')}_{angle_x}_{angle_y}_b{numbins}.png"
            plt.savefig(plot_filepath, dpi=self.plot_dpi)
            plt.clf()
        else:
            plt.show()
        plt.close()

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

    def __init__(self, data_path: str, experiment_prefix: str, ion: str, chain: str):
        self.partial_path = os.path.join(data_path, f"{experiment_prefix}_")
        self.ion = ion
        assert chain in ['analysis', 'sidechain'], f"Incorrect chain type: {chain}"
        self.chain = chain
        experiment_file_names = glob.glob(f"{data_path}/{experiment_prefix}_*_{chain}_{ion}.tab")
        self.no_experiments = len(experiment_file_names)
        experiment_file_names = [ f"{data_path}/{experiment_prefix}_{number}_{chain}_{ion}.tab" for number in range(1,self.no_experiments+1) ]
        self.experiments = [Experiment(fp) for fp in experiment_file_names]
        self.plot_dpi = self.experiments[0].plot_dpi

    def hist_of_entropy(self, xcolumn: str, ycolumn: str, plotdir: str = None):
        """ compute histogram of entropy over realisations """
        xcol = self.experiments[0].get_colnum_by_meaning(xcolumn)
        ycol = self.experiments[0].get_colnum_by_meaning(ycolumn)
        entropies = [experiment.get_entropy(xcol, ycol) for experiment in self.experiments]

        if plotdir:
            plotFile = os.path.join(plotdir,f"hist{xcolumn}_{ycolumn}.png")
            mytitle = f"{self.chain} {self.ion}"
            mytitle = mytitle.replace("analysis", "main chain")
            myxlabel = f"entropy {xcolumn} vs. {ycolumn}"
            myylabel = "frequency"

            plt.subplots()
            plt.hist(entropies, bins=5)
            plt.title(mytitle)
            plt.xlabel(myxlabel)
            plt.ylabel(myylabel)

            plt.savefig(plotFile, dpi=self.plot_dpi)
            plt.close()
        return entropies

    def entropy_distribution_percentiles(self, angle1: str, angle2: str, plotdir: str):
        """  compute percentiles of the histogram of entropies """
        no_mers = self.experiments[0].no_mers
        first_mers = list(range(1, no_mers+1))

        median_entropy = []
        entropy_perc5 = []
        entropy_perc95 = []

        for mer in range(no_mers):
            entropies = np.array(self.hist_of_entropy(f"{angle1} mers {mer+1}, {mer+2}", f"{angle2} mers {mer+1}, {mer+2}"))
            median_entropy.append(np.median(entropies))
            entropy_perc5.append(np.percentile(entropies, 5))
            entropy_perc95.append(np.percentile(entropies, 95))

        mytitle = f"{self.chain}, ion {self.ion}"
        mytitle = mytitle.replace("analysis", "main chain")
        myylabel = f"entropy  {angle1} vs. {angle2}"
        myxlabel = "mers x, x+1"

        plt.plot(first_mers, median_entropy, "o--", color="red", label="median")
        plt.plot(first_mers, entropy_perc5, ":", color="red", label="5 and 95 perc.")
        plt.plot(first_mers, entropy_perc95, ":", color="red", label=None)
        plt.legend()
        plt.title(mytitle)
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        plot_filepath = os.path.join(plotdir, f"entropy_{self.chain}_{self.ion}_{angle1}_{angle2}.png")
        plt.savefig(plot_filepath, dpi=self.plot_dpi)
        plt.clf()

        return median_entropy

    def entropy_distribution_realisations(self, angle1: str, angle2: str, plotdir: str):
        """  compute percentiles of the histogram of entropies """
        no_mers = self.experiments[0].no_mers
        no_struct = self.no_experiments
        first_mers = list(range(1, no_mers+1))
        entropies = []

        for mer in range(no_mers):
            entropies.append(np.array(self.hist_of_entropy(f"{angle1} mers {mer+1}, {mer+2}", f"{angle2} mers {mer+1}, {mer+2}")))

        mytitle = f"{self.chain}, ion {self.ion}"
        mytitle = mytitle.replace("analysis", "main chain")
        myylabel = f"entropy  {angle1} vs. {angle2}"
        myxlabel = "mers x, x+1"

        for j in range(no_struct):
            color = 5 / 6 * (1 - j / no_struct) * np.array((1, 1, 1))
            e = [entropies[i][j] for i in range(no_mers)]
            plt.plot(first_mers, e, label=f"{j+1}", color=color)

        plt.legend()
        plt.title(mytitle)
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        plot_filepath = os.path.join(plotdir, f"entropy_reals_{self.chain}_{self.ion}_{angle1}_{angle2}.png")
        plt.savefig(plot_filepath, dpi=self.plot_dpi)
        plt.clf()

def correct_signs(series, thres: float = 0.5):
    """ if the sign of the angle is artificially reversed,
    we reverse it back"""
    series_corrected = np.array(series)

    for i in range(len(series_corrected) - 1):
        if np.abs(series_corrected[i] + series_corrected[i + 1]) < thres * np.abs(series_corrected[i + 1]):
            series_corrected[i + 1:] *= -1
    return series_corrected
