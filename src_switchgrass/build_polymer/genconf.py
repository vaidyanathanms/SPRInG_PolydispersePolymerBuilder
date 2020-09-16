#------------------------------------------------------------------
# Ver: Sept-04-2020
# Author: Vaidyanathan Sethuraman
# To generate the initial configuration file for lignin topology
# Use NAMD to run the script
# An iterative procedure is used

# References for lignin structure: 
# (A) Yan et al., Biomass and Bioenergy 34, 48-53, 2010 
# (B) Samuel et al., Frontiers in Energy Research, 1 (14) 2014

# 'None' is a keyword reserved - DONT USE IT for PDB/PSF filenames.
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
from make_genpsf import patch_ratios
from make_genpsf import init_logwrite
from make_genpsf import cumul_probdist
from make_genpsf import create_segments
from make_genpsf import create_patches
from make_genpsf import write_multi_segments
from make_genpsf import write_segments_onego
from make_genpsf import run_namd

# Input data
casenum    = 1  # case number
fl_constraint = 1 # flag for reading constraints (constraints.inp)
biomas_typ  = 'switchgrass_' # type of biomass; end with _
deg_poly   = 17 # degree of polymerization (final)
swit_opt   = 'A' # references, A,B (add more and hard code if necesary)
seg_name = 'swli' #name of segment: switchgrass lignin
num_chains = 1 # number of chains
graft_opt = [1,'PCA','GOG'] # graft option, res, patch
tol = 0.1 # relative tolerance for residue/patch generation
maxatt = 500 # maximum attempts to obtain avg configuration
itertype  = 'multi' # O/p style: single-> one go. multi-> multi iter
iterinc   = 4 # iteration increments (for multi itertype)

# Input file names
input_ctr = 'constraints.inp' # patch constraint file input
input_top = 'lignin.top' # topology file input
input_pdb = 'G-bO4L-G.pdb' # file input - dimer

# Output file names (will be generated automatically)
pdbpsf_name = biomas_typ + str(casenum)  #prefix for pdb/psf files
reslist_fname = 'reslist_' + str(casenum) + '.tcl' #all res list
patch_fname = 'patchlist_' + str(casenum)+ '.tcl' #all patch list
log_fname = 'log_' + str(casenum) + '.txt' #log file

# Open log file
flog = open(outdir + '/' + log_fname,'w')
init_logwrite(flog,casenum,biomas_typ,deg_poly,swit_opt,input_top\
              ,input_pdb,seg_name,num_chains)
print('Begin analysis for: ',biomas_typ,', case_num: ', casenum)

# Read defaults and throw exceptions
if not os.path.exists(input_top):
    exit('Topology file not found \n')

if not os.path.exists(input_pdb):
    exit('Dimer file not found \n')

if fl_constraint:
    if not os.path.exists(input_ctr):
        exit('Constraint file not found \n')

# residue % dictionary mode
resperc_dict = residue_ratios(swit_opt) 
if not bool(resperc_dict):
    exit('ERROR: Unknown option for monomer parameters \n')

if graft_opt[0] == 1:
    flog.write('Grafting incorporated while building the chains\n')
else:
    flog.write('Linear chains are built\n')

# patches %dict mode
patchperc_dict = patch_ratios(swit_opt,graft_opt,resperc_dict) 
if not bool(patchperc_dict):
    exit('ERROR: Unknown option for patch parameters \n')

# Get directory info
srcdir = os.getcwd()
outdir = srcdir + str('/casenum_') + str(casenum) # out directory

# Create main directory and copy required files
if not os.path.isdir(outdir):
    os.mkdir(outdir)

gencpy(srcdir,outdir,input_top)
gencpy(srcdir,outdir,input_pdb)
gencpy(srcdir,outdir,input_ctr)

# Open file and write headers
fresin = open(outdir + '/' + reslist_fname,'w')
fresin.write(';# Contains all segments for NAMD files.\n')
fpatchin = open(outdir + '/' + patch_fname,'w')
fpatchin.write(';# Contains all patches for NAMD files.\n')
            
# Create cumulative probability distribution of segments/patches
flog.write('Making cumulative distribution for segments..\n')
print('Making cumulative distribution for segments..')
cumul_resarr = cumul_probdist(resperc_dict,flog)
flog.write(str(cumul_resarr)+'\n')

flog.write('Making cumulative distribution for patches..\n')
print('Making cumulative distribution for patches..')
cumul_patcharr = cumul_probdist(patchperc_dict,flog)
flog.write(str(cumul_patcharr)+'\n')
    
# Set 2D default list and generate segments/patches
res_list = [[] for i in range(num_chains)]
patch_list = [[] for i in range(num_chains-1)]

# Create segments and check for avg probability 
flog.write('Creating residue list..\n')
res_list = create_segments(fresin,deg_poly,num_chains,seg_name,\
                           resperc_dict,cumul_resarr,tol,maxatt,flog)

# Create patches with constraints and check for avg probability 
flog.write('Creating patches list..\n')
patch_list = create_patches(fpatchin,deg_poly,num_chains,seg_name,\
                           patchperc_dict,cumul_patcharr,tol,\
                           maxatt,flog,fl_constraint,input_ctr,\
                           res_list)

flog.write('Writing data to files \n')
flog.write('Output style %s\n' %(itertype))
print('Writing data to files..')
print('Output style: ', itertype)

if itertype == 'single':
    # Open single file and write all details
    tcl_fname  = biomas_typ + str(casenum) + '_allchains.tcl' 
    fmain = open(outdir + '/' + tcl_fname,'w')
    psfgen_headers(fmain,input_top,pdbpsf_name)
    flog.write('Writing config for n-segments: %d\n' %(deg_poly))
    print('Writing config for n-segments: ', deg_poly)
    write_segments_onego(fmain,deg_poly,num_chains,seg_name,\
                         res_list,patch_list)
    psfgen_postprocess(fmain,input_pdb,itertype,0,'None')

    #Exit and close file
    fmain.write('exit')
    fmain.close()

elif itertype == 'multi':

    for chcnt in range(num_chains):
        chnum = chcnt + 1
        flog.write('Writing chain number: ', chnum)
        tcl_fname  = biomas_typ + 'case_' + str(casenum) + \
                     '_chnum_' + str(chnum) +'.tcl' 
        fmain = open(outdir + '/' + tcl_fname,'w')

        flog.write('Iteration increment counter %d\n' %(iterinc))

        # Write segments according to iteration number
        iter_num = 1
        nmonsthisiter = iterinc
    
        while nmonsthisiter <= deg_poly:
            flog.write('Writing config for n-segments: %d\n' %(nmonsthisiter))
            print('Writing config for n-segments: ', nmonsthisiter)
            write_multi_segments(fmain,iter_num,nmonsthisiter,chnum,\
                                 num_chains,seg_name,res_list,patch_list)
            psfgen_postprocess(fmain,input_pdb,itertype,iter_num,seg_name)
            run_namd(fmain,'namd2', 'mini.conf', 'mini.out')
            iter_num  = iter_num + 1
            nmonsthisiter = nmonsthisiter + iterinc

        # Write the rest in one go
        if deg_poly%iterinc != 0:
            flog.write('Writing config for n-segments: %d\n' %(deg_poly))
            print('Writing config for n-segments: ', deg_poly)
            iter_num = iter_num + 1
            write_multi_segments(fmain,iter_num,deg_poly,chcnt,\
                                 num_chains,seg_name,res_list,patch_list)
            psfgen_postprocess(fmain,input_pdb,itertype,iter_num,seg_name)
            run_namd(fmain, 'namd2', 'mini.conf', 'mini.out')

        # Exit and close nchth file
        fmain.write('exit')
        fmain.close() 

else:
    exit('ERROR: Unknown output write style option: ' + itertype)


flog.write('Completed psf generation for casenum: %d\n' %(casenum))
print('Completed psf generation for casenum: ', casenum)

# Close files
fresin.close()
fpatchin.close()
flog.close()

