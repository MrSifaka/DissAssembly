import pysam
import seaborn as sns
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="Plot mapping quality score from bam file")

parser.add_argument('--bam', required=True, help='path to input bam file', action='store')
parser.add_argument('--png', required=True, help='path to output png file', action='store')

args=parser.parse_args()

output = args.png
input = args.bam

bam = pysam.AlignmentFile(input, "rb")

quals = [read.mapping_quality for read in bam.fetch()]

sns.set_style("dark")
sns.set_context("paper")

plt.figure(figsize=(10,6))

ax = sns.kdeplot(quals, shade = True)
ax.set_ylabel('Density')
ax.set_xlabel('Mapping Quality Score')

plt.savefig(output)
