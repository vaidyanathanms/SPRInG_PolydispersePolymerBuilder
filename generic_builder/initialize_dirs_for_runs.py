# To generate initial directories and fill with initial files for
# lignin systems. Will run LigninBuilder and generate initial
# pdb/top/psf files

# Version: Mar-14-2021
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
import fileinput
from supp_initdirs import *
#------------------------------------------------------------------

# Version Info and # of input args for parsing Lignin Builder
print("Generating GROMACS run-time inputs")
print("Version: Mar-03-2021")
if len(sys.argv) == 2:
    inp_fyle = str(sys.argv[1])
else:
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()
#------------------------------------------------------------------

# Input Data
inp_type  = 'cosolvents'       # melts, solvents, cosolvents
biomass   = 'MYB' # name of the biomass type
disperse  = 'mono' # mono/poly; only for melts
o_sol_typ = 'EOH'  # prefix for solvent file for solvents/cosolvents
wat_type  = 'tip3p' # prefix for water file for coslvents
solv_name = o_sol_typ # change if prefix is different from name in PDB
wat_name  = 'TIP3_' # diff from prefix
run_arr   = [1,2,3] # run number
nchains   = 1     # number of polymer chains
npoly_res = 22  # number of polymer residues
n_orgsolv = 1000 # number of organic solvents
nwater    = 4000 # number of water molecules (for cosolvents)
box_dim   = 15 # box size for solvent only. cosolvent=+3
#------------------------------------------------------------------

# Directory Paths
main_dir  = os.getcwd() # current dir
scr_dir   = '/lustre/or-hydra/cades-bsd/v0e' # scratch dir

if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
scr_dir  = scr_dir + '/lignin'
if not os.path.isdir(scr_dir):
    os.mkdir(scr_dir)
#------------------------------------------------------------------

# Required GMX/sh Files
lig_fyles = ['make_genpsf.py','genconf.py','findmissingterms.py']
#------------------------------------------------------------------

#Main Code
for run_id in range(len(run_arr)):
    
    print( "Initializing for run_" + str(run_arr[run_id]))
      
    # Make directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        os.mkdir(head_dir)

    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        os.mkdir(poly_dir)

    if inp_type == 'melts':
        poly_dir = poly_dir + '/' + disperse
        if not os.path.isdir(poly_dir):
            os.mkdir(poly_dir)
    
    outdir = poly_dir + '/run_' + str(run_arr[run_id])
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    findir = set_working_dir(outdir,inp_type,o_sol_typ)
    supp_findir = findir + '/init_files'
    if not os.path.isdir(supp_findir):
        os.mkdir(supp_findir)

    # Retrieve case number for moving into the new directory
    casenum = retrieve_case_num(inp_fyle)
    # Run genconf.py
    run_genconf(inp_fyle,lig_fyles,casenum)
    # Check casenum directory is made
    init_dir = main_dir + '/casenum_' + str(casenum)
    if not os.path.isdir(init_dir):
        raise RuntimeError(init_dir + " does not exist!")
    os.chdir(init_dir)
    # Run all tcl steps from init_dir
    run_all_steps(init_dir)
    # Copy all pdb/topology files
    cpy_pdb_top(init_dir,findir,biomass)
    cpy_supp_files(init_dir,supp_findir)

    # End details and return directory handle
    print("End making files for run_" + str(run_arr[run_id]))
    os.chdir(main_dir)
