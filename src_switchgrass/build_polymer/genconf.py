#------------------------------------------------------------------
# Ver: Sept-04-2020
# Author: Vaidyanathan Sethuraman
# To generate the initial configuration file for lignin topology
# Use NAMD to run the script
# An iterative procedure is used

# References for lignin structure: 
# (A) Yan et al., Biomass and Bioenergy 34, 48-53, 2010 
# (B) Samuel et al., Frontiers in Energy Research, 1 (14) 2014


# Switchgrass variety: Alamo switchgrass

# H:G:S = 26:42:32 (A); G:S = 0.75 - 0.78, H = 2 (B)
# pCA:FA = 1 (A); pCA:FA = 6:32 (B) 

#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math

# Function definitions
from make_genpsf import gencpy
from make_genpsf import psfgen_headers
from make_genpsf import psfgen_postprocess
from make_genpsf import residue_ratios
from make_genpsf import linker_ratios
from make_genpsf import init_logwrite
from make_genpsf import cumul_probdist
from make_genpsf import create_segments
from make_genpsf import make_segments

# Input data
casenum    = 1  # case number
biomas_typ  = 'switchgrass_' # type of biomass; end with _
deg_poly   = 16 # degree of polymerization (final)
swit_opt   = 'A' # references, A,B (add more and hard code if necesary)
input_top = 'lignin.top' # topology file input
input_pdb  = 'G-bO4L-G.pdb' # file input - dimer
seg_name = 'swli' #name of segment: switchgrass lignin
num_chains = 1 # number of chains
tol = 0.1 # relative tolerance
maxatt = 500 # maximum attempts to obtain avg configuration

# Output file names (will be generated automatically)
tcl_fname  = biomas_typ + str(casenum) + '.tcl' # outfile for tcl
pdbpsf_name = biomas_typ + str(casenum)  #prefix for pdb/psf files
alllist_fname = 'alllist_' + str(casenum) + '.tcl' #all seg/link list
log_fname = 'log_' + str(casenum) + '.txt' #log file

# Read defaults and throw exceptions
if not os.path.exists(input_top):
    exit('Topology file not found \n')

if not os.path.exists(input_pdb):
    exit('Dimer file not found \n')

resperc_dict = residue_ratios(swit_opt) # residue % dictionary mode
if not bool(resperc_dict):
    exit('ERROR: Unknown option for monomer parameters \n')

linkperc = linker_ratios(swit_opt) # linker % dictionary mode
if not bool(linkperc):
    exit('ERROR: Unknown option for linker parameters \n')

# Get directory info
srcdir = os.getcwd()
outdir = srcdir + str('/casenum_') + str(casenum) # out directory

# Create main directory and copy required files
if not os.path.isdir(outdir):
    os.mkdir(outdir)

gencpy(srcdir,outdir,input_top)
gencpy(srcdir,outdir,input_pdb)

# Open file and write headers
fmain = open(outdir + '/' + tcl_fname,'w')
psfgen_headers(fmain,input_top,pdbpsf_name)
flist = open(outdir + '/' + alllist_fname,'w')
flist.write('; Contains all segments and linkers for NAMD files.\n')
flog = open(outdir + '/' + log_fname,'w')
init_logwrite(flog,casenum,biomas_typ,deg_poly,swit_opt,input_top\
              ,input_pdb,seg_name,num_chains)
print('Begin analysis for ',biomas_typ,', case_num: ', casenum)
            
# Number of iterations needed per chain
niter = int(math.log(deg_poly)/math.log(2)) #closest power of 2
flog.write('Number of iterations per chain: %d\n ' %(niter))

# Create cumulative probability distribution of segments
flog.write('Making cumulative distribution..\n')
cumul_resarr = cumul_probdist(resperc_dict)
flog.write(str(cumul_resarr))
    
# Set 2D default list and generate segments/linkers
res_list = [[] for i in range(num_chains)]
link_list = [[] for i in range(num_chains-1)]

# Create segments/links and check for avg probability 
flog.write('Creating residue list..\n')
res_list = create_segments(flist,deg_poly,num_chains,seg_name,\
                           resperc_dict,cumul_resarr,tol,maxatt,flog)


# Write segments according to iteration number
for iterval in range(niter):
    iter_num  = iterval + 1
    nmonsthisiter = pow(2,iter_num)
#    write_segments(fmain,iter_num,nmonsthisiter,seg_name,\
#                   res_cumulprob_list)
              
# Generate patches


# Run NAMD
psfgen_postprocess(fmain,input_pdb)
