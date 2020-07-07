#!/usr/bin/env python
import numpy as np
import pysam
import os
import matplotlib.pyplot as plt
import sklearn.cluster as skc
from numba import njit
import datetime
import time



def get_chrlist(filename):
    # get the list of chromosomes
    # input：bam file; output：the list of chromosomes
    samfile = pysam.AlignmentFile(filename, "rb")
    List = samfile.references
    chrList = np.full(len(List), 0)
    for i in range(len(List)):
        chr = str(List[i]).strip('chr')
        if chr.isdigit():
            chrList[i] = int(chr)
    return chrList


def get_RC(filename, chrList, ReadCount):
    # get read count for each position
    # input：bam file, the list of chromosomes and rc array. output：rc array
    samfile = pysam.AlignmentFile(filename, "rb")
    for line in samfile:
        if line.reference_name:
            chr = line.reference_name.strip('chr')
            if chr.isdigit():
                num = np.argwhere(chrList == int(chr))[0][0]
                posList = line.positions
                ReadCount[num][posList] += 1
    return ReadCount


def get_mapq(filename, chrList, mapq):
    # get mapq for each position
    # input：bam file, the list of Chromosomes and rc array. output：rc array
    samfile = pysam.AlignmentFile(filename, "rb")
    for line in samfile:
        if line.reference_name:
            chr = line.reference_name.strip('chr')
            if chr.isdigit():
                num = np.argwhere(chrList == int(chr))[0][0]
                posList = line.positions
                mapq[num][posList] += line.mapq
    return mapq


def read_ref_file(filename, ref, num):
    # read reference file
    # input：fasta file, ref array and the chromosomes. output： ref array
    if os.path.exists(filename):
        print("Read reference file: " + str(filename))
        with open(filename, 'r') as f:
            line = f.readline()
            for line in f:
                linestr = line.strip()
                ref[num] += linestr
    else:
        print("Warning: can not open " + str(filename) + '\n')
    return ref


def ReadDepth(mapq,ReadCount, binNum, ref,binSize):
    # get read depth
    '''
       1. compute the mean of rc in each bin;
       2. count the number of 'N' in ref. If there is a 'N' in a bin，the rd is not counted;
       3. GC bias
    '''

    RD = np.full(binNum, 0.0)
    GC = np.full(binNum, 0)
    MQ = np.full(binNum, 0.0)
    pos = np.arange(1, binNum+1)
    for i in range(binNum):
        RD[i] = np.mean(ReadCount[i*binSize:(i+1)*binSize])
        MQ[i] = np.mean(mapq[i*binSize:(i+1)*binSize])
        cur_ref = ref[i*binSize:(i+1)*binSize]
        N_count = cur_ref.count('N') + cur_ref.count('n')
        if N_count == 0:
            gc_count = cur_ref.count('C') + cur_ref.count('c') + cur_ref.count('G') + cur_ref.count('g')
        else:
            RD[i] = -10000
            gc_count = 0
        GC[i] = int(round(gc_count / binSize, 3) * binSize)

    index = RD > 0
    RD = RD[index]
    GC = GC[index]
    MQ = MQ[index]
    pos = pos[index]
    RD = gc_correct(RD, GC)
    return pos, RD, MQ

def gc_correct(RD, GC):
    #  gc bias
    bincount = np.bincount(GC)
    global_rd_ave = np.mean(RD)
    for i in range(len(RD)):
        if bincount[GC[i]] < 2:
            continue
        mean = np.mean(RD[GC == GC[i]])
        RD[i] = global_rd_ave * RD[i] / mean
    return RD



def plot(pos, data):
    #plot data
    plt.scatter(pos, data, s=3, c="black")
    plt.xlabel("pos")
    plt.ylabel("rd")
    plt.show()


def plotRDMQ(RD, MQ):
    #plot RD and MQ
    x_value = range(1, len(RD) + 1)
    plt.scatter(x_value, MQ)
    plt.scatter(x_value, RD)
    plt.show()


def merge_bin(labels, maxbin, pos,binSize):
    # merge bin
    start = ends = step = 0
    td_range = []
    for i in range (len(labels)):
        if labels[i]==-1 and start == 0:
            start = pos[i]
            continue
        if labels[i] == -1 and start != 0:
            ends = pos[i]
            step = 0
            continue
        if labels[i] !=-1 and start != 0:
            if (i == len(labels)-1):
                td_range.append([start * binSize, ends * binSize])
            step = step + 1
            if step == maxbin:
                if (ends - start >= maxbin - 2):
                    td_range.append([start* binSize,ends* binSize])
                start = 0
                ends = 0
                step = 0


    return td_range


def get_exact_position(reference, binSize, bam, str1, td_range):
    # get exact position
    bp_exact_position = []
    not100Mperrange = np.empty(len(td_range), dtype=object)
    bf = pysam.AlignmentFile(bam, 'rb')
    maxbin = 1
    for i in range(len(td_range)):
        not100Mperrange[i] = []
        if td_range[i][0] - maxbin * binSize > td_range[i][1] + maxbin * binSize:
            continue

        for r in bf.fetch(reference.split('.')[0], td_range[i][0] - maxbin * binSize, td_range[i][1] + maxbin * binSize):
            if (r.cigarstring != str1 and r.cigarstring != None):
                not100Mperrange[i].append([r.reference_name, r.pos, r.cigarstring,r.isize,r.flag & 64,r.flag & 128])

    cigarresult = "bigrangeciagrresult"
    with open(cigarresult, 'w') as f1:
        for i in range(len(not100Mperrange)):
            start = 0
            end = 0
            
            f1.write("\nthis is " + str(i) + " big range:\n")
            f1.writelines(str(not100Mperrange[i]))
            for j in range(len(not100Mperrange[i])):
                pos = not100Mperrange[i][j][2].index('M')
                if pos == 4 or pos == 5:
                    start = not100Mperrange[i][j][1]
                else:
                    if pos == 1 or pos == 2:
                        end = not100Mperrange[i][j][1] + (int)(not100Mperrange[i][j][2][0: pos]) - 1

            if (start == 0 and end == 0):
                continue

            if (start == 0):
                start = td_range[i][0]

            if (end == 0):
                end = td_range[i][1]

            if (end - start > maxbin * binSize):
                bp_exact_position.append([start, end])

    return  bp_exact_position


def prox_L1(step_size: float, x: np.ndarray) -> np.ndarray:
    """
    L1 proximal operator
    """
    return np.fmax(x - step_size, 0) - np.fmax(- x - step_size, 0)



def prox_tv1d(step_size: float, w: np.ndarray) -> np.ndarray:
    """
    Computes the proximal operator of the 1-dimensional total variation operator.

    This solves a problem of the form

         argmin_x TV(x) + (1/(2 stepsize)) ||x - w||^2

    where TV(x) is the one-dimensional total variation

    Parameters
    ----------
    w: array
        vector of coefficients
    step_size: float
        step size (sometimes denoted gamma) in proximal objective function

    References
    ----------
    Condat, Laurent. "A direct algorithm for 1D total variation denoising."
    IEEE Signal Processing Letters (2013)
    """

    if w.dtype not in (np.float32, np.float64):
        raise ValueError('argument w must be array of floats')
    w = w.copy()
    output = np.empty_like(w)
    _prox_tv1d(step_size, w, output)
    return output



@njit
def _prox_tv1d(step_size, input, output):
    """low level function call, no checks are performed"""
    width = input.size + 1
    index_low = np.zeros(width, dtype=np.int32)
    slope_low = np.zeros(width, dtype=input.dtype)
    index_up  = np.zeros(width, dtype=np.int32)
    slope_up  = np.zeros(width, dtype=input.dtype)
    index     = np.zeros(width, dtype=np.int32)
    z         = np.zeros(width, dtype=input.dtype)
    y_low     = np.empty(width, dtype=input.dtype)
    y_up      = np.empty(width, dtype=input.dtype)
    s_low, c_low, s_up, c_up, c = 0, 0, 0, 0, 0
    y_low[0] = y_up[0] = 0
    y_low[1] = input[0] - step_size
    y_up[1] = input[0] + step_size
    incr = 1

    for i in range(2, width):
        y_low[i] = y_low[i-1] + input[(i - 1) * incr]
        y_up[i] = y_up[i-1] + input[(i - 1) * incr]

    y_low[width-1] += step_size
    y_up[width-1] -= step_size
    slope_low[0] = np.inf
    slope_up[0] = -np.inf
    z[0] = y_low[0]

    for i in range(1, width):
        c_low += 1
        c_up += 1
        index_low[c_low] = index_up[c_up] = i
        slope_low[c_low] = y_low[i]-y_low[i-1]
        while (c_low > s_low+1) and (slope_low[max(s_low, c_low-1)] <= slope_low[c_low]):
            c_low -= 1
            index_low[c_low] = i
            if c_low > s_low+1:
                slope_low[c_low] = (y_low[i]-y_low[index_low[c_low-1]]) / (i-index_low[c_low-1])
            else:
                slope_low[c_low] = (y_low[i]-z[c]) / (i-index[c])

        slope_up[c_up] = y_up[i]-y_up[i-1]
        while (c_up > s_up+1) and (slope_up[max(c_up-1, s_up)] >= slope_up[c_up]):
            c_up -= 1
            index_up[c_up] = i
            if c_up > s_up + 1:
                slope_up[c_up] = (y_up[i]-y_up[index_up[c_up-1]]) / (i-index_up[c_up-1])
            else:
                slope_up[c_up] = (y_up[i]-z[c]) / (i-index[c])

        while (c_low == s_low+1) and (c_up > s_up+1) and (slope_low[c_low] >= slope_up[s_up+1]):
            c += 1
            s_up += 1
            index[c] = index_up[s_up]
            z[c] = y_up[index[c]]
            index_low[s_low] = index[c]
            slope_low[c_low] = (y_low[i]-z[c]) / (i-index[c])
        while (c_up == s_up+1) and (c_low>s_low+1) and (slope_up[c_up]<=slope_low[s_low+1]):
            c += 1
            s_low += 1
            index[c] = index_low[s_low]
            z[c] = y_low[index[c]]
            index_up[s_up] = index[c]
            slope_up[c_up] = (y_up[i]-z[c]) / (i-index[c])

    for i in range(1, c_low - s_low + 1):
        index[c+i] = index_low[s_low+i]
        z[c+i] = y_low[index[c+i]]
    c = c + c_low-s_low
    j, i = 0, 1
    while i <= c:
        a = (z[i]-z[i-1]) / (index[i]-index[i-1])
        while j < index[i]:
            output[j * incr] = a
            output[j * incr] = a
            j += 1
        i += 1
    return


@njit
def prox_tv1d_cols(stepsize, a, n_rows, n_cols):
    """apply prox_tv1d along columns of the matri a
    """
    A = a.reshape((n_rows, n_cols))
    out = np.empty_like(A)
    for i in range(n_cols):
        _prox_tv1d(stepsize, A[:, i], out[:, i])
    return out.ravel()


@njit
def prox_tv1d_rows(stepsize, a, n_rows, n_cols):
    """apply prox_tv1d along rows of the matri a
    """
    A = a.reshape((n_rows, n_cols))
    out = np.empty_like(A)
    for i in range(n_rows):
        _prox_tv1d(stepsize, A[i, :], out[i, :])
    return out.ravel()


def main(params):

    start = time.time()
    binSize = int(params[0])
    bam = params[3]
    str1 = params[4]
    chrList = get_chrlist(bam)
    chrNum = len(chrList)
    refList = [[] for i in range(chrNum)]
    for i in range(chrNum):
        #reference = "chr21.fa"
        reference = params[2]
        refList = read_ref_file(reference, refList, i)

    chrLen = np.full(chrNum, 0)
    for i in range(chrNum):
        chrLen[i] = len(refList[i])
    print("Read bam file:", bam)
    ReadCount = np.full((chrNum, np.max(chrLen)), 0)
    ReadCount = get_RC(bam, chrList, ReadCount)
    mapq = np.full((chrNum, np.max(chrLen)), 0.0)
    mapq = get_mapq(bam, chrList, mapq)
    for i in range(chrNum):
        binNum = int(chrLen[i]/binSize)+1
        pos, RD, MQ = ReadDepth(mapq[0], ReadCount[0], binNum, refList[i],binSize)
        RDMQ = np.full((len(RD), 2), 0.0)
        RDMQ[:,0] = RD
        RDMQ[:,1] = MQ
        alpha = float(params[1])

        res = prox_tv1d(alpha, RDMQ[:,0])
        RDMQ[:,0] = res
        res = prox_tv1d(alpha, RDMQ[:,1])
        RDMQ[:,1] = res
        X = RDMQ
        eps = 0.7
        min_samples = 4
        maxbin = 3

        #algorithm = 'kd_tree'
        #db = skc.DBSCAN(eps, min_samples, algorithm).fit(X)
        db = skc.DBSCAN(eps, min_samples, algorithm='kd_tree').fit(X)        

        labels = db.labels_
        print('The label:', 'eps is',eps, 'min_samples is ',min_samples)
        n_noise_ = list(labels).count(-1)
        
        raito = len(labels[labels[:] == -1]) / len(labels)
        
        td_range = merge_bin(labels, maxbin, pos,binSize)
        #binfile = "binresult"
        binfile = params[5]
        with open (binfile,'w') as fbn:
            fbn.write("\neps is:")
            fbn.write(str(eps))
            fbn.write("\nmin_sample is ")
            fbn.write(str(min_samples))
            fbn.write('\nEstimated number of noise points:')
            fbn.write(str(n_noise_))
            print('Bin number is')
            fbn.write('\nBin number is:')
            for i in range(len(labels)):
                if labels[i]== -1:
                    print(pos[i])
                    fbn.write('\n')
                    fbn.write(str(pos[i]))
            fbn.write('\nBin number over.\n')
        print('Bin number over')
        
        bp_exact_position = get_exact_position(reference, binSize, bam, str1, td_range)

        #middle_result_file = "middle_result"
        middle_result_file = params[6]
        with open (middle_result_file,'w') as fbn:
            fbn.write('After merge, the possible region:\n')
            for line in td_range:
                fbn.write(str(line) + '\n')
            fbn.write("final result is:\n")
            for line in bp_exact_position:
                fbn.write(str(line) + '\n')

        #final_result_file = "final_result"
        final_result_file = params[7]
        np.savetxt(final_result_file, bp_exact_position,fmt='%d')


    end = time.time()
    print(" ** the run time of is: ", end-start, " **")



