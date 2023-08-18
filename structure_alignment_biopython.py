import os
import pandas as pd
import numpy as np
from time import time
import gc

#from helper import parse_fix_matrix, view_alignment
from Bio import AlignIO, SeqIO
from Bio.PDB import PDBParser, Superimposer, StructureAlignment

from multiprocessing import Pool
from numba import jit

import pickle

import argparse 


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", type=str, help= "Input directory containing structures")
parser.add_argument("--selection", type=str, default=None, help= "" )
parser.add_argument("--selection_type", type=str, required=True)
parser.add_argument("--alignment", type=str, required=True)
parser.add_argument("-o", "--output", type=str, help="output file name, output file will be a pickle file")
parser.add_argument("-t", "--threads", type=int, default=1)
parser.add_argument("--threshold", type=int, default=70)
parser.add_argument("--norm", action="store_true")

def load_alignment(args):
    align = AlignIO.read(args.alignment, "fasta")
    aliName2idx={rec.name:i for i, rec in enumerate(align)} 
    return align, aliName2idx

def load_structures(args, list_structures):
    print(f"Loading structures from directory: {args.input}")
    parser = PDBParser()
    #list_structures = os.listdir(args.input)
    #list_structures = 
    print(list_structures)
    models_dict = {}
    count = 0
    for rec in list_structures:
        #print(f'{args.input}/{rec}.pdb')
        #if rec[-4:] != ".pdb":
        #    continue
        models_dict[f'{rec}'] = parser.get_structure(f'{rec}', f'{args.input}/{rec}.pdb')
        print(f"Loaded structure: {f'{rec}'}")
        count += 1

    print(f"Loaded {count} structures")
    return models_dict
#@jit
def create_selection_paird_atom_list(align, aliName2idx,  model_s, model_m, selection_list, args):
    #print(f"Making selections for the {len(names_list)} structures based on file: {args.selection}")
    Atom_list_stationary = []
    Atom_list_move = []
    si = aliName2idx[model_s.id] #.split(".")[0]]
    sj = aliName2idx[model_m.id] #.split(".")[0]]
    #try:
    struct_align = StructureAlignment(align, model_s, model_m, si=si, sj=sj)
    #except:
    #    print(f"si is {si}")
    #    print(f"sj is {sj}")
    #    print(model_s)
    #    print(model_m)
    #    print(align)
    #    print(align[si])
    #    print(align[sj])
    #    exit()
    d_0, d_1 = struct_align.get_maps()
    for idx, rec in enumerate(d_0.items()):
        if args.selection_type == "All":
            cond = rec[1] != None
        else:
            cond = rec[1] != None and idx+1 in selection_list
        if cond:
            #print(rec[0])
            #help(rec[0])
            #exit()
            for atom_s in rec[0].get_atoms():
                
                if atom_s.get_name() == "CA":
                    #Atom_list_stationary.append(atom_s)
                    flag_s = atom_s.get_bfactor() >= args.threshold #70
                    break
                    
            #Atom_arr_stationary = np.stack(Atom_list_stationary)
            for atom_m in rec[1].get_atoms():
                if atom_m.get_name() == "CA":
                    #Atom_list_move.append(atom_m)
                    flag_m = atom_m.get_bfactor() >= args.threshold #70
                    break
                    
            if flag_s and flag_m:
                Atom_list_stationary.append(atom_s)
                Atom_list_move.append(atom_m)
                
    return Atom_list_stationary, Atom_list_move
#create_selection_paird_atom_list_jit = jit(nopython=False)(create_selection_paird_atom_list)
#@jit
def get_rmsd(Atom_list_stationary, Atom_list_move):
    super_imposer = Superimposer()
    super_imposer.set_atoms(Atom_list_stationary, Atom_list_move)
    
    if args.norm:
        rmsd = super_imposer.rms / len(Atom_list_stationary) # Normalise RMSD with length of attoms in alignment
    else:
        rmsd = super_imposer.rms    
    #del super_imposer
    #gc.collect()
    return rmsd
#get_rmsd_jit = jit(nopython=False)(get_rmsd)

def create_experiment_list(models_dict, align, aliName2idx, args):

    list_of_experiments=[]
    dict_structure_dist={}
    for idx_outer, rec_outer in enumerate(models_dict.items()):
        
        #print(f"{rec_outer}")
        if args.selection_type == "range":
            selection_df = pd.read_csv(args.selection)
            sel_start = selection_df.loc[selection_df["ID"]==rec_outer[0].split(".")[0]]["Needle_start"].values 
            sel_end = selection_df.loc[selection_df["ID"]==rec_outer[0].split(".")[0]]["Needle_stop"].values 
            selection_list = list(range(int(sel_start), int(sel_end)))
            dict_structure_dist[f"{rec_outer[0]}"]={} 
        elif args.selection_type == "list":
            selection_df = pd.read_csv(args.selection)
            sel_S = selection_df.loc[selection_df["ID"]==rec_outer[0].split(".")[0]]["S"].values
            sel_D = selection_df.loc[selection_df["ID"]==rec_outer[0].split(".")[0]]["D"].values 
            sel_H = selection_df.loc[selection_df["ID"]==rec_outer[0].split(".")[0]]["H"].values 
            selection_list = [sel_S, sel_D, sel_H]
            
            length_active_outer = len(selection_df.loc[selection_df["ID"]==rec_outer[0].split(".")[0]]["Active_triad"].values[0])
            if length_active_outer != 3:
                print(f'skipping calculating alignment for id {rec_outer[0].split(".")[0]}. active triad is {selection_df.loc[selection_df["ID"]==rec_outer[0].split(".")[0]]["Active_triad"].values[0]}') 
                continue
            dict_structure_dist[f"{rec_outer[0]}"]={} 
        elif args.selection_type == "All":
            selection_list = "All"
            dict_structure_dist[f"{rec_outer[0]}"]={} 
        else:
            raise ValueError('--selection_type must be either range list or All')
            
        for idx_inner, rec_inner in enumerate(models_dict.items()):
            
            if  args.selection_type == "list":
                length_active_inner = len(selection_df.loc[selection_df["ID"]==rec_inner[0].split(".")[0]]["Active_triad"].values[0])
                if length_active_inner < 3:
                    print(f"skipping calculating alignment between id {rec_outer[0].split('.')[0]} and {rec_inner[0].split('.')[0]}") 
                    continue
            
            
            experiment={}
            experiment["name"] = f"{rec_outer[0]}_{rec_inner[0]}"
            experiment["model_s"] = rec_outer[1]
            experiment["model_m"] = rec_inner[1]
            experiment["selection_list"] = selection_list
            experiment["alignment"] = align
            experiment["aliName2idx"] = aliName2idx
            list_of_experiments.append(experiment)
    print(f"Created {idx_inner * idx_outer} piars of structures to align")
    return list_of_experiments, dict_structure_dist

def calc_rmsd_select(experiment):
    #print(experiment)
    Atom_list_stationary, Atom_list_move = create_selection_paird_atom_list(align = experiment["alignment"],
                                     aliName2idx=experiment["aliName2idx"],
                                     model_s=experiment["model_s"],
                                     model_m=experiment["model_m"],
                                     selection_list=experiment["selection_list"],
                                     args=args)
    rmsd_selection = get_rmsd(Atom_list_stationary, Atom_list_move)
    return (rmsd_selection, experiment["name"]) 


def save_dict(dict_structure_dist, args):
    with open(args.output, "wb") as file_writer:
        pickle.dump(dict_structure_dist, file_writer)
    print(f"Saved distances in file {args.output}")

def main(args):
    align, aliName2idx = load_alignment(args)
    structure_dict = load_structures(args,  list(aliName2idx.keys()))
    
    list_pairs, dict_structure_dist = create_experiment_list(structure_dict, align, aliName2idx, args)
    start = time()
    rmsds=[]
    #list_pairs = list_pairs[:665]
    n_experiments = len(list_pairs)
    #print(len(list_pairs))
    #print(list_pairs)
    with Pool(10) as p:
        for i in range(n_experiments//200+1):
            rmsds += p.map(calc_rmsd_select, list_pairs[i*200:(i+1)*200])
            print(f"Done with {len(rmsds)} alignments")
    #for i in range(n_experiments):
    #    rmsds += calc_rmsd_select(list_pairs[i])
    #    print(f"Done with {len(rmsds)} alignments")
    print(rmsds)
    print(len(rmsds))
    for rec in rmsds:
        rmsd = rec[0]
        name_outer = rec[1].split("_")[0]
        name_inner = rec[1].split("_")[1]
        dict_structure_dist[name_outer][name_inner]=rmsd
    save_dict(dict_structure_dist, args)
    print(time()-start)
    return 0


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)