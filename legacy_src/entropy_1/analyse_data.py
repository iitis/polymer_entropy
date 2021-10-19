import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import math
from scipy.stats import entropy
import matplotlib.ticker as mtick
import matplotlib
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def read_txt(folder, file):
    """Reads .txt file with times and algles. Assume r = 101 time events.
    Returns vector of times from simulation (ts) and phi (psi) a r elements vector of vector of angles.
    Argumenst:
    - file: string
    -folder: string
    """
    r = 101
    t = ""
    phi = [[] for i in range(r)]
    psi = [[] for i in range(r)]
    ts = [[] for i in range(r)]
    k = -1
    for l in open(folder+file).readlines():
        if ('t' in l):
            k = k+1
            ts[k] = int(l.split("=")[1])*0.05
        else:
            x = l.split(",")
            if np.size(x) == 2:
                a = float(x[0])
                b = float(x[1])
                if (a != 0)*(b != 0):
                    phi[k].append(a)
                    psi[k].append(b)
    return np.array(phi), np.array(psi), np.array(ts)

def reshuffle_data(T, ks, data_folder):
    """
    Reads many realizations of data given a temperature T.

    Returns phis (psis) being in 3mode array with elements phis[i][j][k]
    where i numerates t form simulation (r+1 elements), j numerates realisations(1:10) and k numerates
    angles
    """
    r = 101
    phis = [[] for i in range(r)]
    psis = [[] for i in range(r)]
    tss = []
    for k in ks:
        phi, psi, ts = read_txt(data_folder, str(k)+"_"+str(T)+"K.txt")
        for i in range(101):
            phis[i].append(phi[i])
            psis[i].append(psi[i])
            tss = ts
    return phis, psis, tss


def plot_hist(Phii, Psii, name, hist_folder = ""):
    """
    Plot histogram for vector of angles Phii and Psii and save it in the folder.
    Returns float - entropy. Number of bins is set to 50.
    """
    x = np.histogram2d(Phii, Psii, bins=50)
    if hist_folder != "":
        fig, ax = plt.subplots()
        plt.hist2d(Phii, Psii, bins=50, normed = True, vmin = 0, vmax = 0.0016, cmap=plt.cm.Reds)
        plt.colorbar(format = mtick.FormatStrFormatter('%.1e'))
        plt.xlabel("$\phi$ - angle")
        plt.ylabel("$\psi$ - angle")
        plt.savefig(hist_folder+name+".pdf")
        plt.clf()
        with open(hist_folder+name+".txt", 'w') as f:
            f.write(str(x))
    y = np.concatenate(x[0])
    v = y/sum(y)
    return 8.314*entropy(v)

def save_entropies(ts, entropis, T, name, folder):
    """
     Save all entropies for various temperatures to the given folder"""
    for i in range(len(T)):
        save_entropy(ts, entropis[i], name+'_'+str(T[i])+'K', folder)


def save_entropy(ts, entropy, name, folder):
    """
    save entropy and simulation times in npz and txt formt to the given folder.
    """
    np.savez_compressed(folder+name+".npz", entropy=entropy, times = ts)
    with open(folder+name+".txt", 'w') as f:
        f.write(str(entropy))
        f.write(str(ts))



def compute_all(l_0, l, folder, hist_folder, data_folder, Ts, realizes):
    """
    compute and save all entropies and plot histrograms. Returns a vector of times included i.e. steps[l_0:l]
    end 3mode array of entropis all_ent, with elements all_ent[i][j][k] where i numerates
    i numerates temperature(0 for 300K, 1 for 303K ...), j numerates a time from  steps[l_0:l]
    (0 for l_0, 1 for l_0+1 ....) and k numerates a realisation (1:10).
    Save output in txt and npz formats
    Arguments:
        - l_0: start of the time period
        - l: (l-1) end of times period (l-1)
        - folder: folder where entropy is ploted and saved
        - hist folder: a folder where histograms are saved, if = "" no histograms are plotted (time is saved)
        - data_folder: a folder where data are,
        - Ts: a vector of temeratures (ints)
        - realizes: a vector of realizations (ints)
    """
    steps = []
    all_ent = []
    for T in Ts:
        phis, psis, steps = reshuffle_data(T, realizes, data_folder)
        entropies = []
        for i in range(l_0, l):
            print(i)
            ent = []
            j = 0
            for k in realizes:
                phi = phis[i][j]
                psi = psis[i][j]
                j = j + 1
                name = "wynik_"+str(k)+"_"+str(T)+"K_"+str(steps[i])+"ns"
                s = plot_hist(phi, psi, name, hist_folder)
                ent.append(s)
            entropies.append(ent)
        all_ent.append(entropies)
    name = 'Entropy_'+str(l_0)+'to'+str(l)
    save_entropies(steps[l_0:l], all_ent, Ts, name, folder)
    return steps[l_0:l], all_ent


def plot_entropies(Ts, name, folder):
    """
    Plot all entropies for various temperatures in T
    Save to given folder
    """
    nums = ["a)", "b)","c)","d)","e)","f)"]
    for i in range(len(Ts)):
        f = np.load(folder+name+"_"+str(Ts[i])+'K.npz')
        plot_entropy(f["times"], f["entropy"], Ts[i], name, nums[i], folder)



def plot_entropy(ts, entropy, T, name, f_num, folder):
    """
    Plot mean and mean +-sd of entropy (over realisations) vs simulation time in ts.
    Save to folder
    Arguments:
        - ts: a vector of times from simulation
        - entropy: a vector of vectors of entropy for a given temperature with elements:
            entropy[j][k] where j numerates a simulation time from ts and k numerates a realisation (1:10).
        - T: temperture (int)
        - name: a part of output file name
        - f_num: a number of figure for presentation (string)
        -folder: a folder where plots are saved
    """
    fig, ax = plt.subplots(figsize=(2.7, 2.2))
    plt.subplots_adjust(bottom = 0.16, top=0.92, left = 0.16, right = 0.92)
    plt.text(np.min(ts), 40.5, f_num)
    plt.text(np.mean(ts)-1, 41.1, 'T = '+str(T)+"K")
    σ = np.std(entropy, 1)
    μ = np.mean(entropy, 1)
    ax.plot(ts, μ, '-o', label = '$\mu $', linewidth = 1., markersize = 2.5, color = 'blue')
    ax.plot(ts, μ + σ, '--', linewidth = .75, label = '$\mu \pm \sigma$', color = 'red')
    ax.plot(ts, μ - σ, '--', linewidth = .75, color = 'red')
    ax.legend(ncol = 2)
    ax.set_ylim(36, 41)
    plt.ylabel('entropy', labelpad = 1.)
    plt.xlabel('simulation time [ns]', labelpad = -1)
    plt.savefig(folder+name+"_"+str(T)+"K.pdf")
    plt.clf()


# exemplary use of the code

Ts = [300, 303, 306, 309, 312, 315]
v = [i for i in range(1,11)]
data_folder = "./Entropy/"
output_folder = "./wyniki_entropy/"
# folder histograms are saved, if "" no histogram is ploted and time is saved
hist_folder = ""

# start and end (l-1) of the time period
l0 = 70
l = 100
s, r = compute_all(l0, l, output_folder, hist_folder, data_folder, Ts, v)
name = 'Entropy_'+str(l0)+'to'+str(l)
plot_entropies(Ts, name, output_folder)

l0 = 1
l = 31
s, r = compute_all(l0, l, output_folder, hist_folder, data_folder, Ts, v)
name = 'Entropy_'+str(l0)+'to'+str(l)
plot_entropies(Ts, name, output_folder)

# exemplary reading npz

file = output_folder+"Entropy_1to31_303K.npz"
f = np.load(file)
f["entropy"]
f["times"]
