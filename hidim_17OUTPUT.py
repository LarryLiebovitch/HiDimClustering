# DEF FU7CTIONS FOR HIDIM
# hidim_16.py
# starting to add some calculation pieces
# FOR THIS PROGRAM USE DATA5.TXT and the smaller DATA6.txt
# Doing the SORT and COUNT aa GEMINI
# NOW trying the MASKING, for real FIRST TIME
# ADDING THE LAST PIECE (of this set) THE BIGLOOP
#COMPUTATIONS similar to hdim_13 but a few info added to TITLE

import numpy as np
from scipy import stats
import math as math
import matplotlib.pyplot as plt
import time
from datetime import datetime
import matplotlib.pyplot as plt
import random


def filein(fname):
    numlines = 0
    xin = []
    f = open(fname, 'r')
    for line in f:
        # print (line, end='')
        xin.append(line)
        numlines = numlines+1
    f.close()
    return xin, numlines
    data, numlines = filein(fname)
    # print ('\ndata\n',data,'\nlines',numlines)
    x = ['0' for i in range(numlines)]
    for i in range(numlines):
        x[i] = data[i].replace('\n', '')
        x[i] = eval(x[i])
    return x, numlines

# reads any text table with n columns
# elements in a row separated by a TAB
# rows separated by a RETURN


def getxn(fname):
    data, numlines = filein(fname)
    dataline = ['0' for i in range(numlines)]
    for i in range(numlines):
        x = data[i]
        y = x.split('\t')
        y[-1] = y[-1].replace('\n', '')
        dataline[i] = y
    # print ('\n\nascii-input',dataline)
    xdata = dataline[:]
    for i in range(numlines):
        inline = len(dataline[i])
        for j in range(inline):
            if xdata[i][j] != '':
                xdata[i][j] = eval(xdata[i][j])
            else:
                xdata[i][j] = None
    # print ('dataline',dataline)
    # print ('xdata',xdata)
    return xdata, numlines

def outmake(x):
    outstring=''
    personso=len(x)
    timeslotso=len(x[0])
    for i in range (personso):
        for j in range (timeslotso):
            outstring=outstring + str(x[i][j]) + '\t'
        outstring=outstring+'\n'
    outstring=outstring[0:-1]
    return outstring

def fileout (filename,filedata):
    f2=open(filename,'w')
    f2.write(filedata)
    f2.close()


def binary_summary_report(arr):
    """
    Takes a 2D binary array, sorts it, and prints a summary table.
    """
    nrows, ncols = arr.shape

    # 1. Convert to Decimal using Dot Product (Fastest)
    # Using 'object' dtype handles up to 100+ columns if necessary,
    # but for < 64 bits, standard int64 is used automatically.
    powers = 1 << np.arange(ncols)[::-1]
    decimal_values = arr.dot(powers)

    # 2. Get unique values and their frequencies
    unique_vals, counts = np.unique(decimal_values, return_counts=True)

    # 3. Reconstruct Binary Rows for display (Vectorized)
    # This creates a bitmask to turn decimals back into [1, 0, 1, 0] format
    unique_binary_rows = ((unique_vals[:, None] & (
        1 << np.arange(ncols)[::-1])) > 0).astype(int)

    # 4. Print the Summary Table
    header = f"{'Unique Decimal':<15} | {
        'Binary Equivalent':<25} | {'Count (Frequency)':<15}"
    print(header)
    print("-" * len(header))

    for i in range(len(unique_vals)):
        # Convert to list for cleaner printing
        bin_row = unique_binary_rows[i].tolist()
        print(f"{unique_vals[i]:<15} | {str(bin_row):<25} | {counts[i]:<15}")

    # Return the dictionary as requested previously
    return dict(zip(unique_vals, counts))


def draw_flexible_box_grid(n, label_list, datain, now):
    total_boxes = 2**n

    # 1. Determine Grid Dimensions (Width >= Height)
    width = 2**((n + 1) // 2)
    height = 2**(n // 2)

    # 2. Handle Custom Labels
    # Convert to list and pad with empty strings if the user provided too few
    full_labels = list(label_list)
    if len(full_labels) < total_boxes:
        full_labels += [""] * (total_boxes - len(full_labels))

    # Take only the amount needed for the grid
    display_labels = full_labels[:total_boxes]

    # 3. Reshape into the grid matrix
    grid_labels = np.array(display_labels).reshape(height, width)

    fig, ax = plt.subplots(figsize=(width * 1.5, height * 1.5))

    for i in range(height):
        for j in range(width):
            # Coordinates for top-down numbering
            rect = plt.Rectangle((j, height - i - 1), 1, 1,
                                 fill=False, edgecolor='black', linewidth=2)
            ax.add_patch(rect)

            # Place the custom label
            ax.text(j + 0.5, height - i - 0.5, grid_labels[i, j],
                    ha='center', va='center',
                    fontsize=12 if n < 5 else 9)

    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_aspect('equal')
    ax.axis('off')
    # plt.title(f"Custom Grid: n={n} ({width}x{height})", pad=20)
    # newtitle = datain+'    '+f"total_boxes={total_boxes}"+ '    ' +now
    newtitle = programname+'    '+datain+'    '+f"boxes={total_boxes}"+\
        '    ' +titlenote+'    '+now
    plt.title(newtitle, pad=20)
    plt.show()


programname='hidim_14.py'
# inout data
# cols separated by TAB
# reos sepsarated by RETURN
# floatihng point numers [0,1]
inputname = input('\ninput file namem (e.g. data6.txt)')
titlenote=input('\ninput title note (e.g. noise=0.7)')
data, ndata = getxn(inputname)
print('\ndata= ', data)
adata = np.array(data)
print('\nadata.shape (row,col)', adata.shape)
[nrow, ncol] = adata.shape
bdata = np.zeros((nrow, ncol), dtype=int)

# #make bdata[i,j] = 0 (if adata <= 0) = 1 (if adata > 0)
for irow in range(nrow):
    for icol in range(ncol):
        if adata[irow, icol] > 0.:
            bdata[irow, icol] = 1
print('\nbdata\n', bdata)

# DO HAMMING DISTANCE COLS
target_0 = np.zeros((nrow), dtype=int)
target_1 = np.ones((nrow), dtype=int)
# for icol in range (ncol)
ham = np.zeros((ncol), dtype=int)
for icol in range(ncol):
    col_test = bdata[:, icol].copy()
    # print ('\ni, col_test= ',icol,col_test)
    ham_0 = np.sum(col_test != target_0)
    ham_1 = np.sum(col_test != target_1)
    ham[icol] = min(ham_0, ham_1)
    print('\ntarget_0, target_1, ham     ', ham_0, ham_1, ham[icol])

# first sort ham
# want the position in ham to be the min value
hamlow = np.zeros((ncol), dtype=int)
hamlow = np.argsort(ham)
cdata = bdata[:, hamlow]  # pretty fancy, thank you Gemini
print('\nbdata\n', bdata)
print('\ncdata\n', cdata)

# THIS IS **TEST** OF ONLY <=63 COLS!!
# 1. Create a "Powers of 2" vector: [2^15, 2^14, ..., 2^0]
# Using bitwise shift (1 << n) is slightly faster than 2**n
powers = 1 << np.arange(ncol)[::-1]

# 2. Convert each row to a decimal integer using Matrix Multiplication
# This is the heavy lifting performed in optimized C/Fortran
decimal_values = cdata.dot(powers)

# 3. Get the sorted indices and reorder the original array
sorted_indices = np.argsort(decimal_values)
sorted_c = cdata[sorted_indices]
print('\nsorted_c\n', sorted_c)

# 4. Get Unique Decimal Values and their Counts
# return_counts=True makes this very fast
# unique_decimals=np.zeros((32),dtype=int)
unique_decimals, counts = np.unique(decimal_values, return_counts=True)
# labels=unique_decimals.astype(str).tolist()

# 5. Result Display
print("Sorted Array:\n", sorted_c)
print("\nUnique Decimal Values:", unique_decimals)
print("Frequency Counts:     ", counts)

# --- Example Usage ---
# Creating a dummy array with 100 rows and 5 columns
# data = np.random.randint(0, 2, size=(100, 5))
results = binary_summary_report(cdata)

# OK that's ALL sorted, let's now do the BIGLOOP
# OK MESSED UP dont need to run SORT but do need COURNTS
now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
for bigloop in range(6):  # SET TO DEGUB correct is range(6)
    n = 2**bigloop
    # edata = sorted_c[:, 0:bigloop]
    edata = cdata[:, 0:bigloop]
    mrow, mcol = edata.shape
    mpowers = 1 << np.arange(mcol)[::-1]
    mdecimal_values = edata.dot(mpowers)
#ADD FIRST COL 0, 1, 2, ETC
#ADD LAST COL OF THE MDECIMAL VALUES
    arow=np.arange(mrow) #ADD A COL OF THIS AS FIRST COL
    bigedata = np.column_stack((arow,edata, mdecimal_values))
    msorted_indices = np.argsort(mdecimal_values)
    sorted_e = edata[msorted_indices]
    sorted_bige=bigedata[msorted_indices]
    munique_decimals, mcounts = np.unique(mdecimal_values, return_counts=True)
    # 5. Result Display
    print("\nSorted Array:", sorted_e)
    print("\nUnique Decimal Values:", munique_decimals)
    print("\nFrequency Counts:     ", mcounts)
    print("\nBIG Array:\n", bigedata)
    print("\nSORTED BIG Array:\n", sorted_bige)
    results = binary_summary_report(edata)
# my_labels = ["Start", "A1", "B2", "C3", "D4", "E5", "F6", "End"]
    my_labels = mcounts.astype(str).tolist()
    draw_flexible_box_grid(bigloop, my_labels, inputname, now)

# for i5 in range (ntotal):
#     xout[i5][0]=x[i5]
#     xout[i5][1]=y[i5]
# zo=outmake(xout)
zo=outmake(bigedata)
outtest=input('\nSAVE OUTPUT TO FILE? (y/n)')
if outtest == 'y':
    filenameout=input ('Please give me filename, e.g. iris2out.txt ')
    fileout (filenameout,zo)
    print ('\nfile saved')
else:
    print ('\nNO file saved')

