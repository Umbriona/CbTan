import os
import pandas as pd
import numpy as np
from time import time

#from helper import parse_fix_matrix, view_alignment
from Bio import AlignIO, SeqIO

from multiprocessing import Pool

import pymol2

import pickle

import argparse 


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", type=str, help= "Input directory containing structures")
parser.add_argument("--selection", type=str, default=None, help= "" )
parser.add_argument("-o", "--output", type=str, help="output file name, output file will be a pickle file")
parser.add_argument("-t", "--threads", type=int, default=1)

with pymol2.PyMOL() as pymol:

    def load_structures(args):
        print(f"Loading structures from directory: {args.input}")
        list_structures = os.listdir(args.input)
        names_list = []
        count = 0
        for rec in list_structures:
            if rec[-4:] != ".pdb":
                continue
            count += 1
            pymol.cmd.load(f'{args.input}{rec}', f'{rec}')
            names_list.append(rec)
        print(f"Loaded {count} structures")
        return names_list

    def create_selections(names_list, args):
        print(f"Making selections for the {len(names_list)} structures based on file: {args.selection}")
        selection_df = pd.read_csv(args.selection)
        sel_names_list = []
        for idx, rec in enumerate(names_list):
            sel_start = selection_df.loc[selection_df["ID"]==rec.split(".")[0]]["Needle_start"].values 
            sel_end = selection_df.loc[selection_df["ID"]==rec.split(".")[0]]["Needle_stop"].values 
            print(f"index: {idx} Start: {sel_start[0]} End: {sel_end[0]}")
            pymol.cmd.select(f'{rec}_sel', f'{rec} resi {sel_start[0]}-{sel_end[0]}')
            sel_names_list.append(f'{rec}_sel')
        return sel_names_list

    def create_experiment_list(names):
        list_pairs=[]
        dict_structure_dist={}
        for idx_outer, rec_outer in enumerate(names):
            dict_structure_dist[f"{rec_outer}"]={}
            for idx_inner, rec_inner in enumerate(names):
                list_pairs.append([f"{rec_outer}", f"{rec_inner}"])
        print(f"Created {idx_inner * idx_outer} piars of structures to align")
        return list_pairs, dict_structure_dist

    def rmsd(rec):
        foo = rec[0]
        bar = rec[1]
        aln = pymol.cmd.cealign(foo, bar)
        return (foo, bar, aln["RMSD"])

    def get_rmsds(list_pairs, dict_structure_dist, args):
        print("Starting alignment process")
        p = Pool(args.threads)
        list_rmsd = p.map(rmsd, list_pairs)
        for rec in list_rmsd:
            dict_structure_dist[rec[0]][rec[1]]=rec[2]
        print("Alignments done!!!")
        return dict_structure_dist



    def save_dict(dict_structure_dist, args):
        with open(args.output, "wb") as file_writer:
            pickle.dump(dict_structure_dist, file_writer)
        print(f"Saved distances in file {args.output}")

    def align(args):
        dict_structure_dist = {}   

        names = load_structures(args)
        if args.selection != None:
            names = create_selections(names, args)

        list_pairs, dict_structure_dist = create_experiment_list(names)

        dict_structure_dist = get_rmsds(list_pairs, dict_structure_dist, args)
        save_dict(dict_structure_dist, args)


    def main(args):

        align(args)

        return 0


    if __name__ == "__main__":
        args = parser.parse_args()
        main(args)