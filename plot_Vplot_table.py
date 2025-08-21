import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns;  sns.set_style({'axes.linewidth': 2, 'axes.edgecolor':'black'})
import os

parser = argparse.ArgumentParser(description='Takes a single Vplot tab file created by make_Vplot_table.py and plots a vplot form it. The table should consist of 2 columns: the fragment midpoint and the size of the fragment.')
parser.add_argument("-f", '--file', type=str, required=True, help='path to tab file')
parser.add_argument("-o", '--out', type=str, required=True, help='base outfile name. For example the name of the sample')
parser.add_argument("-x", '--ext', type=str, nargs = "*", required=False, default=[".svg", ".png"], help='list of file extensions to save the output as.')
parser.add_argument("-b", "--binsize", type=int, default=100, help="Size of the bins the genome is binned into for Vplot plotting.")
parser.add_argument("-s", "--size_binsize", type=int, default=1, help="Size of the bins the sizes are binned into for Vplot plotting.")

args = vars(parser.parse_args())

sns.set(font_scale=2)
sns.set_style("ticks")

sample = os.path.basename(args["file"]).replace(".tab","")

# set genome size
if "ts1" in sample:
    genome_size = 30890
elif "wtc" in sample:
    genome_size = 30070
else:
    genome_size = 34062

df = pd.read_csv(args["file"], sep="\t", names=["midpoint", "size"])

# create genomic bins of a size to plot the vplot into
pos_bin_size = args["binsize"]
pos_bins = list(range(0, genome_size, pos_bin_size))
pos_bins.append(genome_size)
pos_labels = list(range(pos_bin_size, genome_size, pos_bin_size))
pos_labels.append(genome_size)
df["bin"] = pd.cut(df["midpoint"], bins=pos_bins, labels=pos_labels) # add to the dataframe

if args["size_binsize"] != 1:
    # also create bins to bin the sizes into
    size_bin_size = args["size_binsize"]
    size_bins = list(range(0,501, size_bin_size))
    size_labels = list(range(size_bin_size,501, size_bin_size))
    df["size_bin"] = pd.cut(df["size"], bins=size_bins, labels=size_labels) # add to the dataframe
    # count the number of combinations of sizebin and genomicbin
    df_counts = df.groupby(['size_bins','bin']).size()

else:
    # count the number of combinations of size and genomicbin
    df_counts = df.groupby(['size','bin']).size()


df_counts = df_counts.apply(np.log1p) # log of the counts

# replace outliers that lie outside of the 95 percentile with the 95 percentile
value = np.percentile(df_counts, 95)
for i in range(len(df_counts)):
    if df_counts.iloc[i] > value:
        df_counts.iloc[i] = value

heatmap_data = df_counts.unstack(level=0).T # convert it into a image/heatmap
heatmap_data = heatmap_data.reindex(sorted(heatmap_data.columns), axis=1) # sort by the column names
heatmap_data = heatmap_data.iloc[::-1] # turn the index upside down, so smaller sized fragments are at the bottom
heatmap_data = heatmap_data.reindex(list(range(500,0,-1)), fill_value = 0)
heatmap_data.replace(0, np.nan, inplace=True) # replace 0s with NaN to get a white background

# plot heatmap
figure, ax = plt.subplots(figsize=(20,5))
ax = sns.heatmap(heatmap_data, cmap="coolwarm", cbar_kws={'label': 'log of fragment count'})
# set spine visible
for _, spine in ax.spines.items():
    spine.set_visible(True)


# text
plt.title(sample)
plt.xlabel("genome")
plt.ylabel("fragment size")

# set xticks
xticks = list(range(0, genome_size, 5000))
xticks.append(genome_size)
xticks_perc = [x/pos_bin_size for x in xticks]
plt.xticks(xticks_perc, xticks)

# set y ticks
yticks = list(range(0, 501, 50))
yticks_pos = list(range(0, 501, 50))
plt.yticks(ticks = yticks_pos[::-1], labels = yticks)


for ext in args["ext"]:
    plt.savefig(args["out"] + ext, dpi=300, bbox_inches = "tight")
plt.close()
