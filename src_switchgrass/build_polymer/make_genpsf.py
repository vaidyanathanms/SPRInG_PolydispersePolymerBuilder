#---------------------------------------------------------------------------
# Supporting scripts for making a general psf/pdb file
# Switchgrass variety: Alamo switchgrass

# H:G:S = 26:42:32 (A); G:S = 0.75 - 0.78, H = 2 (B)
# pCA:FA = 1 (A); pCA:FA = 6:32 (B)

# Import modules

import os
import sys
import numpy
import re
import shutil
import glob
import random
import collections

# General copy script
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR: ', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)


# Define headers for psf files
def psfgen_headers(fin,topname,outname):
    fin.write(';# headers and inputs \n')
    fin.write('package require psfgen \n')
    fin.write('%s\t %s\n' %('topology',topname))
    fin.write('%s\t %s\n' %('set outputname', outname))
              
# Details for closing input files
def psfgen_postprocess(fin,basic_pdb):
    # outname is already there. no need again
    fin.write(';# Writing output \n')
    fin.write('regenerate angles dihedrals \n')
    fin.write('coordpdb %s\t%\n' %(basic_pdb,';# dimer pdb for reference'))
    fin.write('guesscoord ;# guesses rest of the coordinates \n')
    fin.write('writepdb $outputname.pdb \n')
    fin.write('writepsf $outputnmae.psf \n')
    fin.write(';# exit')


# Define monomer ratios from literature    
def residue_ratios(opt):
# add monomer details
    frac_mons = collections.OrderedDict()
    
    if opt == 'A' or opt == 'a':

        print('here')
        # H:G:S = 26:42:32 (B); pCA:FA = 1

        frac_mons['PHP'] = 26/140 # % PHP (H) monomers
        frac_mons['GUA'] = 42/140 # % GUA (G) monomers
        frac_mons['SYR'] = 32/140 # % SYR (S) monomers
        frac_mons['PCA'] = 20/140 # % PCA monomers
        frac_mons['FEA'] = 20/140 # % FA monomers
        
    elif opt == 'B' or opt == 'b':

        # G:S = 0.78, H = 2 (A); pCA:FA = 6:32

        gmonrat = 0.27


    return frac_mons

# Define linker ratios from literature
def linker_ratios(opt):
# add linker details
    link_rat = dict()

    if opt == 'A' or opt == 'a':

        link_rat['BO4'] = 0.4

    elif opt == 'B' or opt == 'b':

        link_rat['BO4'] = 0

    return link_rat

# Initiate log file
def init_logwrite(flog,casenum,bmtype,M,optv,tfile,pfile,segname,nch):
    flog.write('Creating NAMD file for %s\n' %(bmtype))
    if optv == 'A' | optv == 'a':
        flog.write('Ref: Yan et al., Biomass & Bioener 34, 48-53, 2010\n')
    elif optv == 'B' | optv == 'b':
        flog.write('Ref: Samuel et al., Front. Ener. Res., 1 (14) 2014\n')
    flog.write('Case number: %d\n' %(casenum))
    flog.write('Monomers/Chains: %d\t%d\n' %(M, nch))
    flog.write('Input Topol file/PDB file: %s\t%s\n' %(tfile,pfile))
    flog.write('Segment name: %s\n' %(segname))
    flog.write('Analysis beginning ..\n')

# Create cumulative probability distribution from a dictionary
def cumul_probdist(inpdict):

    dummy_distarr = []
    
    # store first value
    val = list(inpdict.values())[0]
    dummy_distarr.append(val)

    # add rest of the values
    for key in range(1,len(inpdict)-1):#iterate until n-1 elements
        val = dummy_distarr(key) + inpdict[key+1]
        dummy_distarr.append(val)

    # check normalization
    if dummy_distarr[len(dummy_distarr)] != 1:
        print('Warning: data not normalized (', \
              dummy_distarr[len(dummy_distarr)],\
              '). Forcing normalization \n')
        sumval = sum(dummy_distarr)
        
        for cnt in range(len(dummy_distarr))
            dummy_distarr[cnt] = dummy_distarr[cnt]/sumval
            
    return dummy_distarr
    
# Create entire list in one go so that cumulative distribution holds true
def create_segments(flist,nmons,nch,segname,res_dict,cumul_monarr\
                    ,tol,maxattmpt,flog):

    # Write list to a separate file
    flist.write('; Entire segment list\n')
    flist.write('; num_segs\t%d, num_chains\t%d' %(nmons,nch))
    
    for attnum in range(maxattmpt):

        flist.write('; Attempt number \t%d' %(attnum))
        flist.write(' resetpsf\n')
        out_list = [[] for i in range(num_chains)]
   
        for chcnt in range(nch):
            flist.write('; chain number:\t%d\n' %(chcnt+1))
            flist.write(' segment %s {\n' %(segname))

            for moncnt in range(nmons):

                ranval = random.random() #seed is current system time by default
                findflag = 0

                for arrcnt in range(length(distarr)-1):
        
                    if ranval >= distarr[arrcnt] and ranval<distarr[arrcnt+1]:
                
                        flist.write('\tresidue\t%d\t%s\n' \
                                    %(segcnt+1,list(res_dict.key())[arrcnt]))
                        findflag = 1   
                        out_list[chcnt].append(list(res_dict.key())[arrcnt]))
                    
                if findflag != 1:
                print(distarr)
                exit('Error in finding a random residue\n')
            
        
        outdist = []
        for key in inplist_dict:
            outdist.append(sum([i.count(key) for i in cumprobs]))

        #normalize
        sumval = sum(outdist)
        normlist = [x/sumval for x in outdist]
        normval = numpy.linalg.norm(numpy.array(normlist) \
                                    - numpy.array(cumprobs))
    
        if normval <= tol:
            #write to log file
            flog.write('Attempt/avg prob: ')
            break

        else:
            flog.write('Attempt/avg prob: ')

    return out_list
    
# Write segments iteration by iteration
def make_segments(fin,iter_num,nmonsthisiter,segname,res_dict,distarr):

    fin.write('; Iteration number: %d' %(iter_num))
    fin.write('set count %d' %(nmonsthister))
    fin.write(' resetpsf \n')
    fin.write(' segment %s {\n' %(segname))
    
    for segcnt in range(nmonsthisiter):
    
        
                

                
    
