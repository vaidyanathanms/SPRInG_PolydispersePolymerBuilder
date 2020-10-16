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
from make_genpsf import * # function definitions

# Read input file
if len(sys.argv) != 2:
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()
print('Input file name: ', sys.argv[1])

# Set defaults
graft_opt = []; swit_opt = 'None' 
input_namd = 'none'; input_prm = 'none'
casenum,mono_deg_poly,num_chains,fpdbflag,ftopflag,fresflag,fpatflag,\
    fl_constraint,disperflag,fpresctr,fppctr,ffflag,fnamdflag,\
    pmolflag,cleanslate,packtol = def_vals()

# Read from file: see definitions/defaults at the end of the script
with open(sys.argv[1]) as farg:
    for line in farg:
        line = line.rstrip('\n')
        words = line.split()
        # call all functions
        if words[0] == 'case_num': 
            casenum = int(words[1])
        elif words[0] == 'biomass_type': 
            biomas_typ = words[1]
        elif words[0] == 'ff_type':
            swit_opt = words[1]; ffflag = 1
        elif words[0] == 'disperse':
            disper_fyle = words[1]; disperflag = 1
        elif words[0] == 'num_resids':
            mono_deg_poly = int(words[1])
        elif words[0] == 'num_chains':
            num_chains = int(words[1])
        elif words[0] == 'seg_name':
            seg_name = words[1]
        elif words[0] == 'grafting':
            if len(words) < 4 or (len(words)-2)%2 != 0:
                exit('Unknown number of graft options: ' + line)
            else:
                graft_opt.append(int(words[1]))
                for wcnt in range(len(words)-2):
                    graft_opt.append(words[wcnt+2])
        elif words[0] == 'tol':
            tol = float(words[1])
        elif words[0] == 'nattempts':
            maxatt = int(words[1])
        elif words[0] == 'op_style':
            itertype  = words[1]
            if itertype == 'multi':
                iterinc = int(words[2]) if len(words) == 3 \
                          else exit('Args for multi not found: '+line)
        elif words[0] == 'constraint_flag':
            fl_constraint = int(words[1])
        elif words[0] == 'patch_res_constraint':
            fpresctr = 1
            input_pres = words[1] if fl_constraint == 1 \
                        else exit('ERR: Constraint flag not set')
        elif words[0] == 'patch_patch_constraint':
            fppctr = 1
            input_pp = words[1] if fl_constraint == 1 \
                       else exit('Constraint flag not set')
        elif words[0] == 'pdb_ipfile':
            if len(words) != 3:
                exit('Unknown number of args: ' + line)
            else:
                input_pdb = words[1]; def_res = words[2]
                fpdbflag = 1
        elif words[0] == 'top_ipfile':
            input_top = words[1]; ftopflag = 1
        elif words[0] == 'resid_inp':
            resinpfyle = words[1]; fresflag = 1
        elif words[0] == 'patch_inp':
            patinpfyle = words[1]; fpatflag = 1
        elif words[0] == 'namd_inp':
            if len(words) != 3:
                exit('Unknown number of arguments: ' + line)
            else:
                input_namd = words[1]; input_prm = words[2]
                fnamdflag = 1
        elif words[0] == 'clean_directories':
            cleanslate = 1 if words[1] == 'Y' \
                      else 0 if words [1] == 'N' \
                           else exit('Unknown args: '+line)
        elif words[0] == 'gen_packmol':
            pmolflag = 1
            if len(words) != 3 and len(words) != 9:
                exit('Unknown number of arguments: '+ line)
            input_packmol = words[1]; trans_list = []
            if words[1] != 'DEF':
                packtol = float(words[2])
            if len(words) > 3:
                for k in range(6):
                    trans_list.append(words[3+k])            
        else:
            exit('Unknown keyword ' + str(words[0]))


outflag = check_all_flags(casenum,fpdbflag,ftopflag,fresflag,fpatflag\
                          ,fl_constraint,fpresctr,fppctr,swit_opt,\
                          ffflag,fnamdflag,disperflag,mono_deg_poly,num_chains)
if outflag == -1:
    exit()

# Output file names (will be generated automatically)
reslist_fname = 'reslist_' + str(casenum) + '.tcl' #all res list
patch_fname = 'patchlist_' + str(casenum)+ '.tcl' #all patch list
log_fname = 'log_' + str(casenum) + '.txt' #log file

# Get directory info
srcdir = os.getcwd()
head_outdir = srcdir + str('/casenum_') + str(casenum) # main outdir

# Create main directories and copy required files
if not os.path.isdir(head_outdir):
    os.mkdir(head_outdir)
elif cleanslate:
    print('Removing old: ', head_outdir)
    srcdir
    shutil.rmtree(head_outdir)
    os.mkdir(head_outdir)

print('Begin analysis for: ',biomas_typ,', case_num: ', casenum)

# Make monomer array for all chains
if disperflag:
    print('Polydispersity file: ', disper_fyle)
    deg_poly_all,pdival = make_polydisp_resids(disper_fyle,num_chains)
    if pdival == 0:
        exit()
else:
    deg_poly_all = [mono_deg_poly]*num_chains
    pdival = 1.0
    print('Monodispersed case')
print('Tot ch/res/pat/pdi',num_chains,sum(deg_poly_all),\
      sum(deg_poly_all)-num_chains,pdival)

# Open log file
flog = open(head_outdir + '/' + log_fname,'w')
init_logwrite(flog,casenum,biomas_typ,deg_poly_all,swit_opt,input_top\
              ,input_pdb,seg_name,num_chains,maxatt,tol,itertype\
              ,fl_constraint,resinpfyle,patinpfyle,disperflag,pdival)
    
# Read defaults and throw exceptions
if not os.path.exists(input_top):
    exit('Topology file not found \n')

if not os.path.exists(input_pdb):
    exit('Dimer file not found \n')

if fl_constraint:
    if not os.path.exists(input_pres) and not os.path.exists(input_pp):
        exit('No constraint file not found \n')

if fnamdflag == 1:
    if not os.path.exists(input_namd) or not os.path.exists(input_prm):
        exit('No NAMD file found \n')

# residue % dictionary mode
resperc_dict = residue_ratios(swit_opt,resinpfyle) 
if not bool(resperc_dict):
    exit('ERROR: Unknown option for monomer parameters \n')

if graft_opt[0] == 1:
    flog.write('Grafting incorporated while building the chains\n')
else:
    flog.write('Linear chains are built\n')

# patches %dict mode
patchperc_dict = patch_ratios(graft_opt,resperc_dict,swit_opt,patinpfyle) 
if not bool(patchperc_dict):
    exit('ERROR: Unknown option for patch parameters \n')

# copy all files
gencpy(srcdir,head_outdir,input_top)
gencpy(srcdir,head_outdir,input_pdb)
if fl_constraint:
    gencpy(srcdir,head_outdir,input_pres)
    gencpy(srcdir,head_outdir,input_pp)

# Open file and write headers
fresin = open(head_outdir + '/' + reslist_fname,'w')
fresin.write(';# Contains all segments for NAMD files.\n')
fpatchin = open(head_outdir + '/' + patch_fname,'w')
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

# check pdb file defaults
pdbfyleflag = check_pdb_defaults(input_pdb,def_res,seg_name)
if pdbfyleflag == -1:
    exit()

# Create residues and check for avg probability 
print('Generating residues..')
flog.write('Creating residue list..\n')
res_list = create_segments(fresin,deg_poly_all,num_chains,seg_name,\
                           resperc_dict,cumul_resarr,tol,maxatt,\
                           flog,graft_opt,def_res)
if res_list == -1:
    exit()


# Create patches with constraints and check for avg probability 
if fl_constraint:
    # Read patch-patch constraints (in one go)
    print('Reading patch-patch constraints..')
    flog.write('Reading patch-patch constraints..\n')
    flog.write('All constraints \n')
    ppctr_list = read_patch_incomp(input_pp)
    flog.writelines('\t'.join(str(jval) for jval in ival) +\
                    '\n' for ival in ppctr_list)

print('Generating patches..')
flog.write('Creating patches list..\n')
patch_list = create_patches(fpatchin,deg_poly_all,num_chains,seg_name,\
                            patchperc_dict,cumul_patcharr,tol,\
                            maxatt,flog,fl_constraint,input_pres,\
                            res_list,ppctr_list,graft_opt)

if patch_list == -1:
    exit()


# Write to file
flog.write('Writing data to files \n')
flog.write('Output style %s\n' %(itertype))
print('Writing data to files..')
print('Output style: ', itertype)

if pmolflag:
    flog.write('Initiating packmol files for %d\n..' %(casenum))
    print('Initiating packmol files for ', casenum)
    if input_packmol == 'DEF':
        packname = 'packmol_'+biomas_typ+'_case_'+str(casenum)+'.inp'
    else:
        packname = input_packmol
    fpack = open(head_outdir + '/' + packname,'w')        
    initiate_packmol(fpack,biomas_typ,num_chains,packtol)

for chcnt in range(num_chains):
    chnum = chcnt + 1
    flog.write('****Writing chain number: %d***\n' %(chnum))
    print('Writing chain number: ', chnum)

    # Make chain level output directory
    tcldir = head_outdir + '/all_tclfiles_chnum_' + str(chnum)
    if not os.path.isdir(tcldir):
        os.mkdir(tcldir)

    #prefix for pdb/psf/tcl files
    pdbpsf_name = biomas_typ + '_case_' + str(casenum) + \
                  '_chnum_' + str(chnum) 
    tcl_fname  =  pdbpsf_name +'.tcl' 
    fmain = open(tcldir + '/' + tcl_fname,'w')

    # Copy NAMD files
    gencpy(srcdir,tcldir,input_prm) 
    fr = open(input_namd,'r')
    fw = open('mini.conf','w')
    fid = fr.read().replace("py_inpname",pdbpsf_name)
    fw.write(fid)
    fr.close(); fw.close()
    gencpy(srcdir,tcldir,'mini.conf')

    deg_poly_this_chain = deg_poly_all[chcnt]

    if itertype == 'single':
        
        psfgen_headers(fmain,input_top,pdbpsf_name)
        flog.write('Writing config for  %d chains\n' %(num_chains))
        write_multi_segments(fmain,-1,deg_poly_this_chain,num_chains,chnum,\
                             seg_name,res_list,patch_list,graft_opt,\
                             deg_poly_this_chain)
        psfgen_postprocess(fmain,input_pdb,itertype,0,'None')

    elif itertype == 'multi':
        flog.write('Iteration increment counter %d\n' %(iterinc))
        # Write segments according to iteration number
        iter_num = 1
        nmonsthisiter = iterinc
    
        while nmonsthisiter <= deg_poly_this_chain:
            if iter_num == 1:
                psfgen_headers(fmain,input_top,pdbpsf_name)
            flog.write('Writing config for n-segments: %d\n' %(nmonsthisiter))
            write_multi_segments(fmain,iter_num,nmonsthisiter,num_chains,\
                                 chnum,seg_name,res_list,patch_list,\
                                 graft_opt,deg_poly_this_chain)
            psfgen_postprocess(fmain,input_pdb,itertype,\
                               iter_num,seg_name)
            if fnamdflag == 1:
                out_namd = 'mini' + str(iter_num) + '.out'
                run_namd(fmain,'namd2','mini.conf',out_namd)
            iter_num  += 1
            nmonsthisiter = nmonsthisiter + iterinc

        # Write the rest in one go
        if deg_poly_this_chain%iterinc != 0:
            flog.write('Writing config for n-segments: %d\n' \
                       %(deg_poly_this_chain))
            iter_num += 1
            write_multi_segments(fmain,iter_num,deg_poly_this_chain,\
                                 num_chains,chnum,seg_name,res_list,\
                                 patch_list,graft_opt,deg_poly_this_chain)
            psfgen_postprocess(fmain,input_pdb,itertype,\
                               iter_num,seg_name)

            if fnamdflag == 1:
                out_namd = 'mini' + str(iter_num) + '.out'
                run_namd(fmain,'namd2','mini.conf',out_namd)

    else:
        exit('ERROR: Unknown output write style option: ' + itertype)
  
    if pmolflag:
        make_packmol(fpack,pdbpsf_name,1,trans_list)

#Exit and close file
fmain.write('exit')
fmain.close()

if pmolflag:
    fpack.write('\n')
    fpack.write('#---End of PACKMOL file----\n')
    fpack.close()
flog.write('Completed psf generation for casenum: %d\n' %(casenum))
print('Completed psf generation for casenum: ', casenum)


# Close files
fresin.close()
fpatchin.close()
flog.close()



#-------------------------Defaults------------------------------------
# Input data
#casenum    = 1  # case number
#fl_constraint = 1 # flag for reading constraints (patch-patch/patch-res)
#biomas_typ  = 'switchgrass' # type of biomass
#deg_poly   = 22 # degree of polymerization (final)
#swit_opt   = 'A' # references, A,B (add more and hard code if necesary)
#seg_name = 'swli' #name of segment: switchgrass lignin
#num_chains = 10 # number of chains
#graft_opt = [1,'PCA','GOG'] # graft option, res, patch
#tol = 0.1 # relative tolerance for residue/patch generation
#maxatt = 500 # maximum attempts to obtain avg configuration
#itertype  = 'multi' # O/p style: single-> one go. multi-> multi iter
#iterinc   = 4 # iteration increments (for multi itertype)

# Input file names
#input_pres = 'constraints.inp' # patch-res constraint file input
#input_top = 'lignin.top' # topology file input
#input_pdb = 'G-bO4L-G.pdb' # file input - dimer
#input_pp  = 'patch_incomp.inp' # patch-patch incompatibility
