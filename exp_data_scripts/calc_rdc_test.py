import sys
import mdtraj as md
import os
import shutil
import math
import numpy as np
import textwrap
from Bio.PDB import *
from Bio.SeqUtils import seq1
import argparse


# create parser
parser = argparse.ArgumentParser(prog='python calc_exp_data.py', \
         formatter_class=argparse.RawDescriptionHelpFormatter, \
         epilog=textwrap.dedent('''\
Required software/libraries:
- Python 3.x: https://www.python.org
- catdcd: https://www.ks.uiuc.edu/Development/MDTools/catdcd
- SPARTA+: https://spin.niddk.nih.gov/bax/software/SPARTA+
- PALES: https://spin.niddk.nih.gov/bax/software/PALES
- mdtraj: http://mdtraj.org
- numpy: https://numpy.org
- pandas: https://pandas.pydata.org
- Biopython: https://biopython.org
 '''))
# define arguments
parser.add_argument('top',     metavar='structure.xxx',  nargs=1,   type=str, help='topology/structure file [pdb,gro,mae,...]')
parser.add_argument('traj',    metavar='traj.xxx',       nargs='*', type=str, help='list of trajectory files [xtc,trr,dcd,...]')
parser.add_argument('--itask', metavar='I', default=[0], nargs=1,   type=int, help='task id (for splitting trajectory)')
parser.add_argument('--ntask', metavar='N', default=[1], nargs=1,   type=int, help='number of tasks (for splitting trajectory)')

parser.add_argument('--debug',     action='store_true', default=False, help='keep temporary directories')

# define which experimental data will be calculated
parser.add_argument('--rdc',       action='store_true', default=False, help='calculate RDCs')
# parse arguments
args = parser.parse_args()

# name of the topology file
TOP_=vars(args)['top'][0]
# get PREFIX and TYPE
PREFIX_ = (os.path.splitext(TOP_)[0]).split("/")[-1]
TYPE_   = (os.path.splitext(TOP_)[1]).split(".")[1]
# define pdb name
PDB_=PREFIX_+".pdb"
# list of trajectory files
TRAJ_=vars(args)['traj']
# tasks info
ITASK_=vars(args)['itask'][0]
NTASK_=vars(args)['ntask'][0]

# create a working directory labelled by ITASK_
wdir="task-"+str(ITASK_)
os.mkdir(wdir)

# convert to pdb
# This requires catdcd installed:
# https://www.ks.uiuc.edu/Development/MDTools/catdcd/
if(TYPE_!="pdb"):
 os.system("catdcd -o "+wdir+"/"+PDB_+" -otype pdb -stype "+TYPE_+" -s "+TOP_+" -"+TYPE_+" "+TOP_)
else:
 os.system("cp "+TOP_+" "+wdir+"/"+PDB_)

# read trajectory files and topology (from pdb)
# this requires mdtraj
# http://mdtraj.org/1.9.3/
trj = md.load(TRAJ_, top=wdir+"/"+PDB_)
# remove the pdb for all tasks except 0
if(ITASK_!=0): os.remove(wdir+"/"+PDB_)

# slice the trajectory based on ITASK_ and NTASK_
n_frames=int(math.floor(float(trj.n_frames)/float(NTASK_)))
# set initial and final frame
first_frame = ITASK_ * n_frames
last_frame  = first_frame + n_frames
# adjust last task
if(ITASK_==NTASK_-1): last_frame = trj.n_frames
# do the actual slicing
trj=trj.slice(range(first_frame,last_frame), copy=False)

# calculate number of residues
nres=[]
for res in trj.topology.residues: nres.append(res.resSeq)

# Print information about the system to log file
log = open(wdir+"/log", "w")
log.write("** SYSTEM INFO **\n")
log.write("Structure filename: %s\n" % PDB_)
log.write("Trajectory filenames: %s\n" % str(TRAJ_))
log.write("Number of atoms: %d\n" % trj.n_atoms)
log.write("Number of residues: %d\n" % len(set(nres)))
log.write("Number of frames: %d\n" % trj.n_frames)
log.write("Starting frame: %d\n" % first_frame)
log.write("Last frame: %d\n" % last_frame)
log.write("Task id: %d\n" % ITASK_)
log.write("Number of tasks: %d\n" % NTASK_)

# Define format for output
fmt0='%d,'; fmt1=''
for i in range(0, trj.n_frames-1): fmt1+='%.4lf,'
fmt1+='%.4lf'


if(vars(args)['rdc']):
    log.write("- Calculating RDC\n")
    for t in range(0, trj.n_frames):
        # create a temporary directory
        tmpdir = wdir + "/tmp-" + str(t)
        os.mkdir(tmpdir)

        # save pdb file
        ipdb = tmpdir + "/out.pdb"
        trj[t].save_pdb(ipdb)

        # clean it - this requires Bio.PDB
        structure = PDBParser().get_structure('PDB', ipdb)
        # get sequence info (residue name and number)
        # WARNING: assuming pdb with 1 model and 1 chain
        resname = []
        resnum = []
        for i in structure[0].get_chains():
            for j in structure.get_residues():
                resname.append(j.get_resname())
                resnum.append(j.get_id()[1])

        # get sequence (one letter code)
        seq = seq1("".join(resname))

        # number of residues
        nres = len(resname)

        # sanity check of sequence length
        if(nres != len(seq)):
            print("Check length of the protein failed!")
            exit()

        # save clean pdb
        io = PDBIO()
        io.set_structure(structure)
        opdb = tmpdir + "/out-clean.pdb"
        io.save(opdb)

        # create PALES input file
        ifile = tmpdir + "/PALES_input.dat"
        with open(ifile, "w") as f:
            f.write("DATA SEQUENCE ")
            for i in range(0, nres):
                f.write("%s" % seq[i])
                # add a space every 10 residues, but not for last one
                if((i + 1) % 10 == 0 and i != (nres - 1)):
                    f.write(" ")
            f.write("\n\n")
            f.write("VARS   RESID_I RESNAME_I ATOMNAME_I RESID_J RESNAME_J ATOMNAME_J D      DD    W\n")
            f.write("FORMAT %5d     %6s       %6s        %5d     %6s       %6s    %9.3f   %9.3f %.2f\n")
            f.write("\n")
            for i in range(0, nres):
                f.write("%d %3s H %d %3s N 0 1.00 1.00\n" %
                        (resnum[i], resname[i], resnum[i], resname[i]))

        # run PALES on clean pdb
        rdc_frame = []
        expected_resnum = resnum[1:nres-1]
        # cycle on residues, except first and last
        for ires in range(1, nres - 1):
            # determine window
            l = 7
            h = 7
            if(ires < 7):
                l = ires
            if(ires > nres - 8):
                h = nres - 1 - ires
            # window is the minimum between l and h
            w = min(l, h)

            # output file
            ofile = tmpdir + "/" + str(resnum[ires]) + ".dat"

            # run PALES inside apptainer
            pales_cmd = (
                "apptainer run "
                "-B /dartfs-hpc/rc/lab/R/RobustelliP/ "
                "/dartfs-hpc/rc/home/5/f004nf5/labhome/SINGULARITY_BUILDS/pales/pales.sif "
                "-inD {ifile} -pdb {opdb} -r1 {r1} -rN {rN} -outD {ofile}"
            ).format(ifile=ifile, opdb=opdb, r1=resnum[ires - w], rN=resnum[ires + w], ofile=ofile)

            os.system(pales_cmd)

            # parse the output file and add to list of RDCs
            #for lines in open(ofile, "r").readlines():
                #riga = lines.strip().split()
                #if(len(riga) == 12 and riga[0].isdigit()):
                    #if(float(riga[0]) == resnum[ires]):
                        #rdc_frame.append(float(riga[8]))


        rdc_frame = []
        for rn in expected_resnum:
            value = 0.0000  # default if missing
            datfile = tmpdir + "/" + str(rn) + ".dat"
            if os.path.exists(datfile):
               with open(datfile, "r") as fh:
                   for line in fh:
                       riga = line.strip().split()
                       if len(riga) == 12 and riga[0].isdigit():
                          if int(riga[0]) == rn:
                             value = float(riga[8])
                             break
            rdc_frame.append(value)


        # delete the temporary directory
        if(vars(args)['debug'] == False):
            shutil.rmtree(tmpdir)

        # create/add to global numpy array (n_data, n_frames)
        if(t == 0):
            rdc = np.array(rdc_frame)
        else:
            rdc = np.column_stack((rdc, np.array(rdc_frame)))

    # save RDCs to file
    label = np.array(resnum[1:nres - 1])
    np.savetxt(wdir + "/RDC.csv", np.column_stack((label, rdc)),
               fmt=fmt0 + fmt1, header="resSeq,frame")

    # closing log file
log.write("ALL DONE!\n")
log.close()
