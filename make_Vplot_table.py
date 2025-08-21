import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool
import argparse
import progressbar

parser = argparse.ArgumentParser(description='Takes a single sorted bed file and outputs the midpoint and size of all fragments')
parser.add_argument("-f", '--file', type=str, required=True,  help='path to bed file')
parser.add_argument("-o", '--out', type=str, required=True,  help='path to output')

args = vars(parser.parse_args())

def write_vplot(starts, ends):
    # this creates a table of all fragments in sample with 2 columns: midpoint and size
    with open(args["out"], "w") as out:
        with progressbar.ProgressBar(max_value=(len(starts)/2)) as bar:

            for frag in range(0, len(starts), 2):
                # get end and start of the 4 positions that make up the 2 reads of a fragment
                maximum = max([starts[frag], starts[frag + 1], ends[frag], ends[frag +1]])
                minimum = min([starts[frag], starts[frag + 1], ends[frag], ends[frag +1]])
                size = maximum - minimum
                mid = int(minimum + size/2)

                # add line
                out.write(str(mid) + "\t" + str(size) + "\n")

                bar.update(frag/2)


path_to_bed = args["file"]

# read the table and extract the important columns
df = pd.read_table(path_to_bed, names = ["chr", "start", "end", "name", "score", "strand"])
starts = df["start"]
ends = df["end"]

# get the name of the sample
name = os.path.basename(path_to_bed).split(".")[0]

# # multiprocessing
# pool = Pool(os.cpu_count())
# pool.map(write_vplot, starts, ends)

# calculate
print("\nNow processing " + name + "\n")
write_vplot(starts, ends)
