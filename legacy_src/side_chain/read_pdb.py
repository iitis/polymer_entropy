import numpy as np
import os.path
import matplotlib.pyplot as plt
from itertools import groupby



def readxx(folder, file):
    ret = []
    for l in open(folder+file).readlines():
        nice_lines = l[19:-2]
        nl = nice_lines.split(" ")
        ret = np.concatenate((ret, nl), 0)
    return ret

def plotnobonds():
    data_folder = "./Wyniki_SideChain/"
    file = "sequence.txt"
    x = readxx(data_folder, file)
    z = np.unique(x)
    d = {y:list(x).count(y) for y in x}
    xx = [i+1 for i in range(np.size(z))]
    plt.plot(xx, [d[i] for i in z], "--d")
    plt.xticks(xx, z, rotation='vertical')
    plt.ylabel("number of bonds")
    plt.savefig("nobonds.pdf")
    plt.clf()

plotnobonds()


def count_angles(folder, file, polymer = "ALA"):
    """Reads .txt file with times and algles. Assume r = 101 time events.
    Returns

    Argumenst:
    - file: string
    - folder: string
    - polymer code: string
    """

    aminos = readxx(folder, "sequence.txt")
    r = 101
    t = ""
    numbers = [[] for i in range(r)]
    ts = [[] for i in range(r)]
    k = -1
    j = 0
    for l in open(folder+file).readlines():
        if ('t' in l):
            k = k+1
            ts[k] = int(l.split("=")[1])*0.05
            #print(j)
            j = 0
        else:
            x = l.split(",")
            if np.size(x) == 5:

                if aminos[j] == polymer:
                    #print(j)
                    numbers[k].append(sum([float(el) != 0 for el in x]))
                j = j+1
    print(j)
    return np.array(numbers), np.array(ts)




data_folder = "./Wyniki_SideChain/"
file = "1_side_300K.txt"
aminos = readxx(data_folder, "sequence.txt")
np.size(aminos)

n, t = count_angles(data_folder, file, "ALA")
"ALA"
n[0]


def read_pdb(folder, file):
    if not os.path.isfile('./seq/temp.txt'):
        print("no temp")
        nice_lines = [l for l in open(folder+file).readlines() \
            if ('ATOM' in l) and ('REMARK' not in l)]
        text = ''.join(nice_lines).replace("ATOM", "a ")
        with open('./seq/temp.txt', 'w') as f:
            f.write(text)
    data = np.genfromtxt('./seq/temp.txt', usecols=(1, 3, 5), dtype=object,
        converters={1: lambda x: int(x),
        3: lambda x: x.decode('utf-8'),
          5: lambda x: int(x)}, \
        names=('no', 'amin', 'series'))
    return data

def make_seq(polymer = 'PHE'):
    if not os.path.isfile("./seq/"+polymer+"seq.npz"):
        print("make seq")
        x = read_pdb("./Wyniki_SideChain/", "Albumin_303K_last.pdb")
        ins = [y[0] for y in list(filter(lambda y: y[1] == polymer, x))]
        np.savez_compressed("./seq/"+polymer+"seq.npz", index = ins)
        return ins
    else:
        d = np.load("./seq/"+polymer+"seq.npz", allow_pickle = True)
        return d["index"]



def no_aminos():
    amino = []
    no_aminos = []
    x = read_pdb("./Wyniki_SideChain/", "Albumin_303K_last.pdb")
    for amin in set(x["amin"]):
        amino.append(amin)
        no_aminos.append(np.size(make_seq(amin)))
    return amino, no_aminos






n[2]

if __name__ == "__main__":
    plotnobonds()

make_seq(polymer = 'PHE')


np.array([1,2,3]) > 0
