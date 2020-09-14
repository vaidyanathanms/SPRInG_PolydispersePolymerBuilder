#---------------------------------------------------------------------------
# Supporting scripts for making a general psf/pdb file
# Switchgrass variety: Alamo switchgrass

# 'None' is a keyword reserved - DONT USE IT for PDB/PSF filenames.

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

# Define headers for psf files
def psfgen_headers(fin,topname,outname):
    fin.write(';# headers and inputs \n')
    fin.write('package require psfgen \n')
    fin.write('%s\t %s\n' %('topology',topname))
    fin.write('%s\t %s\n' %('set outputname', outname))
#---------------------------------------------------------------------              
# Details for closing input files
def psfgen_postprocess(fin,basic_pdb,writetype,iter_num,segname):
    # outname is already there. no need again
    fin.write(';# Writing output \n')
    fin.write('regenerate angles dihedrals \n')
    if writetype == 'single':
        comnt = 'dimer pdb for reference'
        fin.write('coordpdb %s  ;# %s\n' %(basic_pdb,comnt))
        comnt2 = 'Guesses rest of the coordinates from PDB inp'
    elif writetype == 'multi':
        if iter_num == 1:
            comnt = 'Different in first iteration'
            pdbfyle = basic_pdb
        else:
            comnt = '*.coor is the file generated by NAMD'
            pdbfyle = '$outputname.coor'
        fin.write('coordpdb %s  %s  ;#  %s\n' %(pdbfyle,segname,comnt))
        comnt2 = 'Can create steric clashes and hence the iterations.'
    else:
        exit('ERROR: Unknown option: ' + writetype)
                                                
    fin.write('guesscoord ;#  %s\n' %(comnt2))
    fin.write('writepdb $outputname.pdb \n')
    fin.write('writepsf $outputnmae.psf \n')

#---------------------------------------------------------------------

# Define monomer ratios from literature    
def residue_ratios(opt):
# add monomer details
    frac_res = collections.OrderedDict()
    
    if opt == 'A' or opt == 'a':

        # H:G:S = 26:42:32 (B); pCA:FA = 1
        frac_res['PHP'] = 26/140 # % PHP (H) monomers
        frac_res['GUA'] = 42/140 # % GUA (G) monomers
        frac_res['SYR'] = 32/140 # % SYR (S) monomers
        frac_res['PCA'] = 20/140 # % PCA monomers
        frac_res['FERU'] = 20/140 # % FA monomers
        
    elif opt == 'B' or opt == 'b':

        # G:S = 0.78, H = 2 (A); pCA:FA = 6:32

        gmonrat = 0.27

    return frac_res
#---------------------------------------------------------------------

# Define linker ratios from literature
def linker_ratios(opt):
# add linker details
    frac_link = collections.OrderedDict()

    if opt == 'A' or opt == 'a':

        frac_link['BO4R'] = 0.2
        frac_link['BO4L'] = 0.2
        frac_link['55']  = 0.2
        frac_link['405'] = 0.2
        frac_link['BB']  = 0.2

    elif opt == 'B' or opt == 'b':

        frac_link['BO4'] = 0

    return frac_link
#---------------------------------------------------------------------

# Initiate log file
def init_logwrite(flog,casenum,bmtype,M,optv,tfile,pfile,segname,nch):
    flog.write('Creating NAMD file for %s\n' %(bmtype))
    if optv == 'A' or optv == 'a':
        flog.write('Ref: Yan et al., Biomass & Bioener 34, 48-53, 2010\n')
    elif optv == 'B' or optv == 'b':
        flog.write('Ref: Samuel et al., Front. Ener. Res., 1 (14) 2014\n')
    flog.write('Case number: %d\n' %(casenum))
    flog.write('Monomers/Chains: %d\t%d\n' %(M, nch))
    flog.write('Input Topol file/PDB file: %s\t%s\n' %(tfile,pfile))
    flog.write('Segment name: %s\n' %(segname))
    flog.write('Analysis beginning ..\n')
#---------------------------------------------------------------------

# Create cumulative probability distribution from a dictionary
def cumul_probdist(inpdict,flog):

    dummy_distarr = []

    # store first value
    val = list(inpdict.values())[0]
    dummy_distarr.append(val)

    # add rest of the values
    for key in range(len(inpdict)-1):#iterate until n-1 elements
        val = dummy_distarr[key] + list(inpdict.values())[key+1]
        dummy_distarr.append(val)

    # check normalization
    if dummy_distarr[len(dummy_distarr)-1] != 1:
        print('Warning: data not normalized (', \
              dummy_distarr[len(dummy_distarr)-1],\
              '). Forcing normalization \n')
        flog.write('%s\t%g\t%s\n' %('Warning: data not normalized (', \
                                    dummy_distarr[len(dummy_distarr)-1],\
                                    '). Forcing normalization \n'))
        sumval = sum(dummy_distarr)
        
        # force normalization
        for cnt in range(len(dummy_distarr)):
            dummy_distarr[cnt] = dummy_distarr[cnt]/sumval
            
    else:
        print('Generated target cumulative distribution..')

    return dummy_distarr
#---------------------------------------------------------------------
    
# Create entire list in one go so that cumulative distribution holds true
def create_segments(flist,nmons,nch,segname,inp_dict,cumulprobarr\
                    ,tol,maxattmpt,flog):

    # Write list to a separate file
    flist.write(';#  Entire segment list\n')
    flist.write(';#  num_segs\t%d, num_chains\t%d\n' %(nmons,nch))
    flog.write('Probabilities for each attempt\n')
    flog.write('Attempt#\t')
    for wout in range(len(inp_dict)):
        flog.write('%s\t' %(list(inp_dict.keys())[wout]))
    flog.write('L2norm \n')

    flag_optimal = -1

    for attnum in range(maxattmpt):

        flog.write('%d\t' %(attnum+1))    
        flist.write(';# Attempt number \t%d\n' %(attnum+1))
        flist.write(' resetpsf\n')

        out_list = [[] for i in range(nch)] #reset every attempt
   
        for chcnt in range(nch):
            flist.write(';# chain number:\t%d\n' %(chcnt+1))
            flist.write(' segment %s {\n' %(segname))

            for segcnt in range(nmons):

                ranval = random.random() #seed is current system time by default
                findflag = 0

                for arrcnt in range(len(cumulprobarr)):
        
                    #Only need to check the less than value because
                    #the array is organized in increasing order.
                    #Break the loop once the first point where the
                    #condition is met.
                    if ranval < cumulprobarr[arrcnt]:
                
                        flist.write(' residue\t%d\t%s\n' \
                                    %(segcnt+1,list(inp_dict.keys())[arrcnt]))
                        findflag = 1   
                        out_list[chcnt].append(list(inp_dict.keys())[arrcnt])
                        break

                if findflag != 1:
                    print('Random value/Probarr:', ranval,cumulprobarr)
                    exit('Error in finding a random residue\n')
            
            flist.write(' }')
        # After going through all the chains, count occurence of each res/patch
        outdist = []
        for key in inp_dict:
            outdist.append(sum([i.count(key) for i in out_list]))

        #normalize
        sumval = sum(outdist)
        if sumval != nch*nmons:
            print('Sum from distn,nch*nmons:',sumval,nch*nmons)
            exit('ERROR: Sum not equal to the total # of segments')
        normlist = [x/sumval for x in outdist]

        #extract target probabilities and compare
        targ_probs = list(inp_dict.values())
        normval = numpy.linalg.norm(numpy.array(normlist) \
                                    - numpy.array(targ_probs))
    
        if normval <= tol:
            #write to log file
            for wout in range(len(outdist)):
                flog.write('%g\t' %(outdist[wout]))
            flog.write('%g\n' %(normval))
            flog.write('Found optimal residue configuration\n')
            print('Found optimal residue configuration..')
            flag_optimal = 1
            break

        else:
            flist.write('\n')
            for wout in range(len(outdist)):
                flog.write('%g\t' %(outdist[wout]))
            flog.write('%g\n' %(normval))


    if flag_optimal == -1:
        print('Did not find optimal residue configuration')
        print('Using last residue configuration with L2norm: ', normval)
        flog.write('Did not find optimal residue configuration\n')
        flog.write('Using last configuration with residue L2norm: %g'\
                   %(normval))

    return out_list
#---------------------------------------------------------------------

# Read and check patch constraints -- May not be effective as opposed
# to reading at once and copying to array. Need to think about it
def check_constraints(inpfyle,patchname,resname1,resname2):
    bef_flag = 1; aft_flag = 1 # keep as true
    with open(inpfyle,'r') as fctr: 
        for line in fctr:
            line = line.rstrip('\n')
            all_words = re.split('\W+',line)
            if len(all_words) != 3:
                print('ERR: Line in ctr file does not have 3 entries')
                print(len(all_words),all_words)
                return -2
            if all_words[0] == patchname:
                if all_words[1] == resname1:
                    bef_flag = 0
                elif all_words[2] == resname2:
                    aft_flag = 0

    # Return 0 if any flags are 0, else return 1
    if bef_flag == 0 or aft_flag == 0:
        return 0
    else:
        return 1
#---------------------------------------------------------------------

# Generate patches
def create_patches(flist,nmons,nch,segname,inp_dict,cumulprobarr\
                    ,tol,maxattmpt,flog,ctr_flag,ctrfyle,residlist):

    # Write list to a separate file
    flist.write(';# Entire linker list\n')
    flist.write(';# num_links\t%d, num_chains\t%d\n' %(nmons-1,nch))
    flog.write('Probabilities for each attempt\n')
    flog.write('Attempt#\t')
    for wout in range(len(inp_dict)):
        flog.write('%s\t' %(list(inp_dict.keys())[wout]))
    flog.write('L2norm \n')

    flag_optimal = -1

    for attnum in range(maxattmpt):

        flog.write('%d\t' %(attnum+1))    
        flist.write(';# Attempt number \t%d\n' %(attnum+1))
        out_list = [[] for i in range(nch)] #reset every attempt
   
        for chcnt in range(nch):
            flist.write(';# chain number:\t%d\n' %(chcnt+1))
            flist.write(' segment %s {\n' %(segname))
            segcnt = 0

            while segcnt <= nmons-2: #for checking constraints

                ranval = random.random() #seed is current system time by default
                findflag = 0

                for arrcnt in range(len(cumulprobarr)):
        
                    #Only need to check the less than value because
                    #the array is organized in increasing order.
                    #Break the loop once the first point where the
                    #condition is met.
                    if ranval < cumulprobarr[arrcnt]:
                        patchname = list(inp_dict.keys())[arrcnt]
                        flist.write(' patch\t%d\t%s\t%s:%d\t%s:%d\n' \
                                    %(segcnt+1,patchname,\
                                      segname,segcnt+1,segname,segcnt+2))
                        findflag = 1
                        appendflag = 1#default to 1 so that if
                        #constraints are not there, it will be appended.
                        if ctr_flag:
                            if segcnt == 0:
                                resname1 = 'None'
                            else:
                                resname1 = residlist[chcnt][segcnt]
                            
                            resname2 = residlist[chcnt][segcnt+1]
                            appendflag =check_constraints(ctrfyle,patchname,\
                                                          resname1,resname2)

                        if appendflag == 1: #update while loop if
                            #constraints are met
                            out_list[chcnt].append(patchname)
                            segcnt += 1
                        elif appendflag == -2:
                            return 

                        break

                if findflag != 1:
                    print('Random value/Probarr:', ranval,cumulprobarr)
                    exit('Error in finding a random residue\n')
            
            flist.write(' }')

        # After going through all the chains, count occurence of each res/patch
        outdist = []
        for key in inp_dict:
            outdist.append(sum([i.count(key) for i in out_list]))

        #normalize
        sumval = sum(outdist)
        if sumval != nch*(nmons-1):
            print('Sum from distn,nch*(nmons-1):',sumval,nch*(nmons-1))
            exit('ERROR: Sum not equal to the total # of linkers')
        normlist = [x/sumval for x in outdist]

        #extract target probabilities and compare
        targ_probs = list(inp_dict.values())
        normval = numpy.linalg.norm(numpy.array(normlist) \
                                    - numpy.array(targ_probs))
    
        if normval <= tol:
            #write to log file
            for wout in range(len(outdist)):
                flog.write('%g\t' %(outdist[wout]))
            flog.write('%g\n' %(normval))
            flog.write('Found optimal patch configuration\n')
            print('Found optimal patch configuration..')
            flag_optimal = 1
            break

        else:
            flist.write('\n')
            for wout in range(len(outdist)):
                flog.write('%g\t' %(outdist[wout]))
            flog.write('%g\n' %(normval))


    if flag_optimal == -1:
        print('Did not find optimal patch configuration')
        print('Using last patch configuration with L2norm: ', normval)
        flog.write('Did not find patch optimal configuration\n')
        flog.write('Using last patch configuration with L2norm: %g'\
                   %(normval))

    return out_list
#---------------------------------------------------------------------

# Write residues/patches in one go
def write_segments_onego(fin,nmons,nch,segname,res_list,link_list):

    fin.write(';# Writing % segments' %(nmons))
    fin.write(' resetpsf \n')
    fin.write(' segment %s {\n' %(segname))
    
    #Residues
    for chcnt in range(nch):

        for segcnt in range(nmons):

            fin.write('  residue  %d  %s\n' \
                      %(segcnt+1,res_list[chcnt][segcnt]))

    fin.write('}')        
    fin.write('\n')

    #Patches
    for chcnt in range(nch):

        for segcnt in range(nmons-1):

            fin.write('patch  %s  %s:%d  %s:%d\n' \
                        %(link_list[chcnt][segcnt],\
                          segname,segcnt+1,segname,segcnt+2))

    fin.write('\n')
#---------------------------------------------------------------------

# Write residues/patches iteration by iteration
def write_multi_segments(fin,iter_num,nmonsthisiter,nch,chnum,\
                         segname,res_list,link_list):

    if iter_num == 1:
        fin.write(';# Chain number: %d of %d chains' %(nch,chnum))

    fin.write(';# Iteration number: %d\n' %(iter_num))
    fin.write('set count %d' %(nmonsthisiter))
    fin.write('\n')
    fin.write(' resetpsf \n')
    fin.write(' segment %s {\n' %(segname))

    #Residues -- indices should have -1 for first dimension
    for segcnt in range(nmonsthisiter):

        fin.write('  residue  %d  %s\n' %(segcnt+1,\
                                          res_list[chnum-1][segcnt]))

    fin.write('}')        
    fin.write('\n')
        
    #Patches -- indices should have -1 for first dimension
    for segcnt in range(nmonsthisiter-1):

        fin.write('patch  %s  %s:%d  %s:%d\n' \
                  %(link_list[chnum-1][segcnt],segname,segcnt+1,\
                    segname,segcnt+2))

    fin.write('\n')
#---------------------------------------------------------------------

# Run generic namd script
def run_namd(fin,execfyle,inpfyle,outfyle):
    fin.write(';# Run NAMD\n')
    fin.write('%s  %s  > %s\n' %(execfyle,inpfyle,outfyle))        
    fin.write(';# exit \n')
    fin.write(';# -------------------------------------\n')
    fin.write('\n')
#---------------------------------------------------------------------
