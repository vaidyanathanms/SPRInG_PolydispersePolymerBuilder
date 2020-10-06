#------------------------------------------------------------------
# Ver: Oct-05-2020
# Author: Vaidyanathan Sethuraman
# Fun def: to create melt using PACKMOl with input pdbs
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math

# General copy script
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR: ', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)
#---------------------------------------------------------------------
# Set defaults
def def_vals():
    return -1, 0, 0, 0, 0, 0, 0, 0, 0, 0
#---------------------------------------------------------------------

# Check all flags 
def check_all_flags(casenum,fpdbflag,ftopflag,fresflag,fpatflag,\
                    fl_constraint,fpresctr,fppctr,opt,ffflag,fnamd):
    outflag = 1
    if casenum == -1:
        print('ERR: Case number not input'); outflag = -1
    elif fpdbflag == 0 or ftopflag == 0:
        print('ERR: PDB/Topology file not entered'); outflag = -1
    elif ffflag == 0:
        print('ERR: Force field type not set: A,B, None');outflag = -1
    elif fresflag == 0 and (opt == 'none' or opt=='None'):
        print('ERR: Residue file/option input not entered'); outflag = -1
    elif fpatflag == 0 and (opt == 'none' or opt=='None'):
        print('ERR: Patch file/option not entered'); outflag = -1
    elif fl_constraint == 1:
        if fpresctr == 0 or fppctr == 0:
            print('ERR: Constraint files not given: constraint flag ON')
            outflag = -1
    elif fnamd == 0:
        print('WARNING: No NAMD file found')

    return outflag
#---------------------------------------------------------------------

# Define headers for packmol files
def packmol_headers(fin,topname,outname):
    fin.write(';# headers and inputs \n')
    fin.write('package require psfgen \n')
    topinp = '../' + topname
    fin.write('%s\t %s\n' %('topology',topinp))
    fin.write('%s\t %s\n' %('set outputname', outname))
#---------------------------------------------------------------------             

# Run generic namd script
def run_packmol(fin,execfyle,inpfyle,outfyle):
    fin.write(';# Run NAMD\n')
    fin.write('%s  %s  > %s\n' %(execfyle,inpfyle,outfyle))        
    fin.write(';# exit \n')
    fin.write(';# -------------------------------------\n')
    fin.write('\n')
#---------------------------------------------------------------------
