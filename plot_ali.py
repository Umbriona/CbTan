import pandas

from helper import parse_fix_matrix, view_alignment
from Bio import AlignIO, SeqIO

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import panel as pn
aln = AlignIO.read('data/alignment_sequence/Tom_sequences_split_22_11_29.ali','fasta')
p = view_alignment(aln, plot_width=90)
pn.pane.Bokeh(p)