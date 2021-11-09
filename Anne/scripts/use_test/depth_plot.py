#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt
from bokeh.io import export_png
from bokeh.layouts import row 
from bokeh.plotting import figure
import seaborn as sns
# total arguments
n = len(sys.argv)
print("Total arguments passed:", n)
print(sys.argv[0])
print(sys.argv[1])
print(sys.argv[2])
df = pd.read_csv(sys.argv[1], sep="\t", header=None, names=["Chr", "locus", "depth"])
print(df.head())
print(len(df))
get_rows = df #.head(200)
# bokeh
#print('bokeh')
#p = figure(title = "test")
#p.circle('locus','depth',source=get_rows, fill_alpha=0.2, size=10)
#export_png(p, filename=sys.argv[2])

#hist, edges = np.histogram(get_rows['depth'], density=True, bins=10)
#p2 = figure()
#p2.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color="white")
#export_png(p2, filename='test.png')

# seaborn
print('seaborn')
plt.figure(figsize=(10, 10))
sns.scatterplot(data=get_rows, x="locus", y="depth")
plt.savefig('seaborn_scatter.png')
plt.clf()

plt.figure(figsize=(15, 15))
sns.histplot(data=get_rows, x="depth")
plt.yscale('log')
plt.savefig('seaborn_hist.png')
plt.clf()

# matplotlib
print('matplotlib')
plt.figure(figsize=(5, 5))
plt.scatter(get_rows['locus'], get_rows['depth'], alpha=0.5)
plt.savefig('matplot_scatter.png')
plt.clf()

plt.figure(figsize=(12, 12))
plt.hist(get_rows['depth'], bins=50)
plt.yscale('log')
plt.savefig('matplot_hist.png')
plt.clf()