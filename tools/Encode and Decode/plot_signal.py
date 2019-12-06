#!/usr/bin/env python2
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import subprocess, os, platform
import argparse
from matplotlib.ticker import FormatStrFormatter



def generate_rawlora_graph(filepath, sf, sr, output_pdf, vmin):
    filepath = os.path.realpath(filepath)
    fs=sr
    start_offset = 0
    length = os.path.getsize(filepath)
    test_file = np.fromfile(filepath, dtype=np.complex64)[start_offset:start_offset+length]

    matplotlib.rcParams.update({
    'font.size': 4,
    })
    fig = plt.figure(figsize=(10, 2),dpi=800)
    axis = plt.gca()
    axis.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    plt.axvline(x=176513/1e6, linewidth=0.4)
    plt.axvline(x=180609/1e6, linewidth=0.4)
    plt.axvline(x=184705/1e6, linewidth=0.4)
    plt.axvline(x=188802/1e6, linewidth=0.4)
    plt.axvline(x=192897/1e6, linewidth=0.4)
    plt.axvline(x=196994/1e6, linewidth=0.4)
    plt.axvline(x=201089/1e6, linewidth=0.4)
    plt.axvline(x=205185/1e6, linewidth=0.4)
    plt.text(0.3, 125000, "g", fontsize=14)


    symbol_len = 2**sf * 8
    axis.set_ylabel("Frequency (kHz)")
    axis.set_xlabel("Sample")

    plt.specgram(test_file, NFFT=256, Fs=fs, noverlap=64, cmap='plasma', vmin=vmin )
    axis.set_xticklabels(np.array(axis.get_xticks() * fs, dtype=int))
    axis.set_ylim([-500000*0.4, 500000*0.4])
    plt.yticks([ -125000, -62500,  0,  62500, 125000])
    axis.set_yticklabels(np.array(axis.get_yticks() / 1000.0, dtype=float))
    plt.tight_layout()
    plt.legend()
    with PdfPages(output_pdf) as pdf:
        pdf.savefig()

def main ():
    parser = argparse.ArgumentParser()
    parser.add_argument("filepath", help="abs or rel path to file / recording")
    parser.add_argument("sf", help="spreading factor of the signal")
    parser.add_argument("-s", "--sampling_rate",  default=1e6,
                        help="sampling rate at which signal was recorded, default=1000000")
    parser.add_argument("-p", "--power_min", 
                        help="ignore out signals with power below threshold, try values -40 til -80 ")
    args = parser.parse_args()

    filepath= args.filepath
    sf = int(args.sf)
    sr = int(args.sampling_rate)
    vmin = None
    if args.power_min:
        vmin = int(args.power_min)
    
    output_pdf = 'rawframe.pdf'

    generate_rawlora_graph(filepath=filepath, sf=sf, sr=sr, output_pdf=output_pdf, vmin=vmin)
    if platform.system() == 'Darwin':       # macOS
            subprocess.call(('open', output_pdf))
    elif platform.system() == 'Windows':    # Windows
        os.startfile(output_pdf)
    else:                                   # linux variants
            subprocess.call(('xdg-open', output_pdf))
    print ("Done")

if __name__ == "__main__":
    main()


