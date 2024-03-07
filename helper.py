import os, io, random
import string
import numpy as np

from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO

import panel as pn
import panel.widgets as pnw
pn.extension()

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

import logomaker as lm
import matplotlib.pyplot as plt
from  matplotlib.ticker import FuncFormatter

def parse_fix_matrix(fileLocation):


    dict_mat = {"id": [], "vec":[]}
    matx = []
    f = open(fileLocation).readlines()
    firstLine = f.pop(0)
    for line in f:
        line = line.strip().split(' ')
        dict_mat["id"].append(line[0]) 
        line.pop(0)

        line = [elem for elem in line if elem.strip()]
        assert (int(firstLine) == (len(line)))
        dict_mat["vec"].append(np.array(line))
        matx.append(line)
    # print (firstLine)
    
    return dict_mat, np.array(matx)

def atom2res(plddt_atom_vec, idx_vec):
    plddt_atom_vec = [float(rec) for rec in plddt_atom_vec.split(";") if rec != ""]
    idx_vec = [float(rec) for rec in idx_vec.split(";") if rec != ""]
    plddt_res=[]
    prev_idx=0
    for plddt_atom, res_idx_atom in zip(plddt_atom_vec, idx_vec):
        if res_idx_atom > prev_idx:
            plddt_res.append(plddt_atom)
            prev_idx=res_idx_atom
        else:
            continue
    return plddt_res
            
#test_idx_vec = "1;1;1;2;2;2;2;3;3;4;5;5;5;5;5;5;5;5;5;5;5;5;5"
#test_plddt_atom_vec = "10;10;10;20;20;20;20;20;20;40;50;50;50;50;50;50;50;50;50;50;50;50;50"


def make_seq(length=40):    
    return ''.join([random.choice(['A','C','T','G']) for i in range(length)])

def mutate_seq(seq):
    """mutate a sequence randomly"""
    seq = list(seq)
    pos = np.random.randint(1,len(seq),6)    
    for i in pos:
        seq[i] = random.choice(['A','C','T','G'])
    return ''.join(seq)

def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i.upper() for s in list(seqs) for i in s]
    color_lookup =  {'A':[200,200,200],
             'R':[ 20, 90,255],
             'N':[  0,220,220],
             'D':[230,230, 10],
             'C':[230,230,  0],
             'E':[230,230, 10],
             'Q':[  0,220,220],
             'G':[235,235,235],
             'H':[130,130,210],
             'I':[ 15,130, 15],
             'L':[ 15,130, 15],
             'V':[ 15,130, 15],
             'K':[ 20, 90,255],
             'M':[230,230,  0],
             'F':[ 50, 50,170],
             'P':[220,150,130],
             'S':[250,150,  0],
             'T':[250,150,  0],
             'W':[180, 90,180],
             'Y':[ 50, 50,170],
             '-':[0,0,0]}
 #   ASP, GLU	Bright Red	 	[230,230, 10]	E60A0A
#    CYS, MET	Yellow	 	[230,230,  0]	E6E600
#    LYS, ARG	Blue	 	[ 20, 90,255]	145AFF
#     SER, THR	Orange	 	[250,150,  0]	FA9600
#     PHE, TYR	Mid Blue	 	[ 50, 50,170]	3232AA
#ASN, GLN	Cyan	 	[  0,220,220]	00DCDC
#GLY	Light Grey	 	[235,235,235]	EBEBEB
#LEU, VAL, ILE	Green	 	[ 15,130, 15]	0F820F
#ALA	Dark Grey	 	[200,200,200]	C8C8C8
#TRP	Purple	 	[180, 90,180]	B45AB4
#HIS	Pale Blue	 	[130,130,210]	8282D2
#PRO	Flesh	 	[220,150,130]	DC9682
#Others	Tan	 	[190,160,110]	BEA06E




    colors = [color_lookup[i] for i in text]
    return colors

def muscle_alignment(seqs):
    """Align 2 sequences with muscle"""
    filename = 'temp.faa'
    SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=filename, out=name+'.txt')
    stdout, stderr = cline()
    align = AlignIO.read(name+'.txt', 'fasta')
    return align

def view_alignment(aln, fontsize="9pt", plot_width=800):
    """Bokeh sequence alignment view"""

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]    
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)    
    N = len(seqs[0])
    S = len(seqs)    
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto')
    if N>1000:
        viewlen=1000
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"

    #entire sequence view (no text, with zoom)
    p = figure(title=None, width= plot_width, height=50,
               x_range=x_range, y_range=(0,S), tools=tools,
               min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False  

    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, width=plot_width, height=plot_height,
                x_range=view_range, y_range=ids, tools="xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)          
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font="monospace",text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = gridplot([[p],[p1]], toolbar_location='below')
    return p


df_col = {"D": list(np.array([230,230, 10])/255),
          "E": list(np.array([230,230, 10])/255),
          "C": list(np.array([230,230,  0])/255),
          "M": list(np.array([230,230,  0])/255),
          "K": list(np.array([ 20, 90,255])/255),
          "R": list(np.array([ 20, 90,255])/255),
          "S": list(np.array([250,150,  0])/255),
          "T": list(np.array([250,150,  0])/255),
          "F": list(np.array([ 50, 50,170])/255), 
          "Y": list(np.array([ 50, 50,170])/255),
          
          "N": list(np.array([  0,220,220])/255),
          "Q": list(np.array([  0,220,220])/255),
          "G": list(np.array([200,200,200])/255),  # [235,235,235]
          "L": list(np.array([ 15,130, 15])/255),
          "V":  list(np.array([15,130, 15])/255),
          "I":  list(np.array([15,130, 15])/255),
          "A": list(np.array([160,160,160])/255),
          "W": list(np.array([180, 90,180])/255),
          "H": list(np.array([130,130,210])/255),
          "P": list(np.array([220,150,130])/255),
          "X": list(np.array([190,160,110])/255)}
def plotRegions(df, title, start, end , y_max = 4):
   # plt.figure(figsize=[20,5])
    end = start + df[start:end].shape[0]
    lm.Logo(df[start:end], figsize=[(end-start)/2,5], color_scheme = df_col)
    plt.title(f'Profile of {title}')
    plt.ylabel('Information (Bits)')
    #plt.ylabel('Probability')
    plt.xlabel('Position')
    plt.gca().set_xticks([i for i in range(start,end+1, 5)])
    print(f"Y_lim: {plt.gca().get_ylim()}")
    plt.gca().set_ylim([0, y_max])
    plt.savefig(f"results/{title}_pattern.png")
    return 0


def label_maker_TanCb_active(x):
    if x[0] != "S" or x[2] != "H":
        return "Not Tan"
    elif x[1] in ["D", "E"]:
        return "Acid"
    elif x[1] in ["N", "Q"]:
        return "Amide"
    else:
        return "Not Tan"

def label_maker_TanCb_active_why_not(x):

    string = "Missing: "
    if x[0] != "S":
        string = string + "".join("S, ")
    if x[2] != "H":
        string = string + "".join("H, ")
    if x[1] not in ["D", "E", "N", "Q"]:
        string = string + "".join("[D/E/N/Q]")

    if string == "Missing: ":
        return "Is Tan"
    if string[-2:] == ", ":
        return string[:-2]
    return string