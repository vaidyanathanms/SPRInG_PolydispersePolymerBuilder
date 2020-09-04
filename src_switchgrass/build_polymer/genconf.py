#------------------------------------------------------------------
# Ver: Sept-04-2020
# Author: Vaidyanathan Sethuraman
# To generate the initial configuration file for lignin topology
# Use NAMD to run the script
# An iterative procedure is used

# References for lignin structure: 
# (A) Samuel et al., Frontiers in Energy Research, 1 (14) 2014
# (B) Yan et al., Biomass and Bioenergy 34, 48-53, 2010 

# Switchgrass variety: Alamo switchgrass

# G:S = 0.78, H = 2 (A); H:G:S = 27:41:32 (B)
# pCA:FA = 6:32 (A); pCA:FA = 1 (B)

#------------------------------------------------------------------

# Import modules

import os
import sys
import numpy
import re
import shutil
import glob

# Function definitions

from make_genpsf import psfgen_headers
from make_genpsf import monomer_ratios
from make_genpsf import psfgen_postprocess
from make_genpsf import gencpy

# Input data

casenum    = 1  # case number
biomas_typ  = 'switchgrass_' # type of biomass; end with _
deg_poly   = 18 # degree of polymerization (final)
swit_opt   = 1 # references, A,B (add more and hard code if necesary)
input_top = 'lignin.top' # topology file input
input_pdb  = 'G-b04L-G.pdb' # file input - dimer

# Output file names (will be generated automatically)

tcl_fname  = biomas_typ + str(casenum) + '.tcl' # outfile for tcl
pdbpsf_fname = biomas_typ + str(casenum)  #prefix for pdb/psf files

# Define empty arrays and default values

struct_arr = []
gsrat   = 0 # G:S
pfrat   = 0 # PCA:FA
num_php = 0 # num of PHP (H) monomers
num_gua = 0 # num of GUA (G) monomers
num_syr = 0 # num of SYR (S) monomers
num_pca = 0 # num of PCA monomers
num_fa  = 0 # num of FA monomers


# Throw exceptions

if not os.path.exists(input_top):
    exit('Topology file not found \n')

if not os.path.exists(input_pdb):
    exit('Dimer file not found \n')

optflag = monomer_ratios(swit_opt)
if optflag == -1:
    exit('ERROR: Unknown option for parameters \n')
    

# Get directory info

srcdir = os.getcwd()
outdir = cwd + str('/casenum_') + str(casenum) # out directory

# Create main directory and copy required files

if not os.path.isdir(outdir):
    os.mkdir(outdir)
    
gencpy(srcdir,destdir,input_top)
gencpy(srcdir,destdir,input_pdb)

# Open file and write headers

fmain = open(outdir + '/' + tcl_fname,'w')
psfgen_headers(fmain,input_top,pdbpsf_name)
monomer_ratios(swit_opt)
psfgen_postprocess(fmain,input_pdb)
         


# Generate segments



# Generate patches


# Run NAMD


