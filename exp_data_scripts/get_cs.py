#!/usr/bin/env python3
import sys
import mdtraj as md
import os
import shutil
import math
import numpy as np
import pandas
import textwrap
from Bio.PDB import *
from Bio.SeqUtils import seq1
import argparse
import multiprocessing

def return_cs(trj): 
    return md.chemical_shifts_spartaplus(trj, rename_HN=True)


PDB = sys.argv[1]
TRAJ = sys.argv[2]
offset = int(sys.argv[3])
chunksize = int(sys.argv[4])
nproc = int(sys.argv[5])
ofile_name = './cs_assignment.csv'
# read trajectory files and topology (from pdb)
# this requires mdtraj
# http://mdtraj.org/1.9.3/
trj = md.load(TRAJ, top=PDB)

# remove the pdb for all tasks except 0
# calculate number of residues


nres=[]
for res in trj.topology.residues: nres.append(res.resSeq)

log = open("log", "w")
log.write("** SYSTEM INFO **\n")
log.write("Structure filename: %s\n" % PDB)
log.write("Trajectory filenames: %s\n" % str(TRAJ))
log.write("Number of atoms: %d\n" % trj.n_atoms)
log.write("Number of residues: %d\n" % len(set(nres)))
log.write("Number of frames: %d\n" % trj.n_frames)
del(trj)
# Calculate chemical shifts
# This requires SPARTA+ installed and pandas
# https://spin.niddk.nih.gov/bax/software/SPARTA+/
# https://pandas.pydata.org
log.write("- Calculating chemical shifts\n")
with multiprocessing.Pool(processes=nproc) as Pool:
    cs=Pool.map(return_cs,md.iterload(TRAJ, top=PDB, chunk=chunksize))

# print the panda DataFrame to csv
nchunks = len(cs)
final = pandas.DataFrame()

for num in range(len(cs)):
    if num!=0:
        cs[num].columns += chunksize*(num+1)
    sys.stdout.write(f"\r{((num+1)/nchunks)*100}\% Done relabling frames.")
    sys.stdout.flush()

final = pandas.concat(cs, axis=1, verify_integrity=True)
print(final.head())

check_offset=True
with open(ofile_name,'w') as file_out:
    ofinal = final.to_csv().splitlines()[1:]
    len_ofinal = len(ofinal)
    for num, line in enumerate(ofinal):
        if check_offset:
            oset = offset - int(line.split(',')[0])
            check_offset = False
        rem = len(line.split(',')[0])
        resnr = str(int(line.split(',')[0]) + oset)
        file_out.write(resnr+line[rem:]+'\n')
        sys.stdout.flush()
        sys.stdout.write(f"\r{((num+1)/len_ofinal)*100}\% Done Writing to file {ofile_name}.")
        

#closing log file
log.write("ALL DONE!\n")
log.close()
