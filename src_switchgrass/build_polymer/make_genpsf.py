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

# Define patch ratios from literature
def patch_ratios(opt,opt_graft,resdict):
# add patch details

    frac_patch = collections.OrderedDict()

    if opt_graft[0] == 1:
        newfrac_patch = collections.OrderedDict() #create new dict
        gr_resname = opt_graft[1]
        gr_patname = opt_graft[2]
        resflag = 0
        for rescnt in range(len(resdict)):
            if list(resdict.keys())[rescnt] == gr_resname:
                resflag = 1
                graft_prob = list(resdict.values())[rescnt]
                frac_patch[gr_patname] = graft_prob
                newfrac_patch[gr_patname] = graft_prob

        if resflag == 0:
            print('ERROR: Could not find ', str(gr_resname))
            return 0

    if opt == 'A' or opt == 'a':

        frac_patch['BO4R'] = 0.2
        frac_patch['BO4L'] = 0.2
        frac_patch['55']  = 0.2
        frac_patch['405'] = 0.2
        frac_patch['BB']  = 0.2

    elif opt == 'B' or opt == 'b':

        frac_patch['BO4'] = 0

    if opt_graft[0] != 1:
        return frac_patch

    # Renormalize if grafts are present and create new dict
   
    if opt_graft[0] == 1:
        sumprob = 0
        for patcnt in range(len(frac_patch)):
            if list(frac_patch.keys())[patcnt] != opt_graft[2]:
                sumprob += list(frac_patch.values())[patcnt]

        normval = sumprob/(1-graft_prob)
        print(sumprob,normval)
        for patcnt in range(len(frac_patch)):
            if list(frac_patch.keys())[patcnt] != opt_graft[2]:

                newprob = list(frac_patch.values())[patcnt]/normval
                keyval = list(frac_patch.keys())[patcnt]
                newfrac_patch[keyval] = newprob

    return newfrac_patch
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
    if abs(dummy_distarr[len(dummy_distarr)-1]-1) > pow(10,-5):
        print('Warning: data not normalized (', \
              dummy_distarr[len(dummy_distarr)-1],\
              '). Forcing normalization')
        flog.write('%s\t%g\t%s\n' %('Warning: data not normalized (', \
                                    dummy_distarr[len(dummy_distarr)-1],\
                                    '). Forcing normalization'))
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
                    ,tol,maxattmpt,flog,graftopt):

    # Write list to a separate file
    flist.write(';#  Entire segment list\n')
    flist.write(';#  num_resds\t%d, num_chains\t%d\n' %(nmons,nch))
    flog.write('Probabilities for each attempt\n')
    flog.write('Attempt#\t')
    for wout in range(len(inp_dict)):
        flog.write('%s (%g)\t' %(list(inp_dict.keys())[wout],\
                                 list(inp_dict.values())[wout]))

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
            rescnt = 0

            while rescnt < nmons:

                ranval = random.random() #seed is current system time by default
                findflag = 0
                consecresflag = 0 #default: consecutive res are NOT found

                for arrcnt in range(len(cumulprobarr)):
        
                    #Only need to check the less than value because
                    #the array is organized in increasing order.
                    #Break the loop once the first point where the
                    #condition is met.
                    if ranval < cumulprobarr[arrcnt]:
                
                        findflag = 1   
                        resname1 = list(inp_dict.keys())[arrcnt]

                        if rescnt == 0: #no restriction. write directly
                            flist.write(' residue\t%d\t%s\n' \
                                        %(rescnt+1,resname1))
                            out_list[chcnt].append(resname1)
                            rescnt = rescnt + 1
                        else:
                            resname2 = out_list[chcnt][rescnt-1]
                            consecresflag = is_res_cons(resname1,resname2\
                                                        ,graftopt)
                            if consecresflag == 0:
                                flist.write(' residue\t%d\t%s\n' \
                                            %(rescnt+1,resname1))
                                out_list[chcnt].append(resname1)
                                rescnt = rescnt + 1

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
            exit('ERROR: Sum not equal to the total # of residues')
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
# to reading at once and copying to array. Need to think about it.
# Check special cases using files
# THIS IS CONSTRAINT FOR RESIDUE1-PATCH-RESIDUE2 combination
def check_constraints(inpfyle,patchname,resname1,resname2):
    
    bef_flag = 1; aft_flag = 1 # keep as true
    with open(inpfyle,'r') as fctr: 
        for line in fctr:
            line = line.rstrip('\n')
            all_words = re.split('\W+',line)
            if len(all_words) != 3:
                print('ERR: Constraint file does not have 3 entries')
                print(len(all_words),all_words)
                return -2 # return -2
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

# check consecutive residues - cannot have graft residue in
# consecutive positions
def is_res_cons(resname1,resname2,graftopt):
    sameflag = 0
    if graftopt[0] == 1:
        if resname1 == graftopt[1] and resname2 == graftopt[1]:
            sameflag = 1
    return sameflag
#---------------------------------------------------------------------

# read all patch incompatibilities
def read_patch_incomp(fname):
    with open(fname,'r') as fin:
        result = [[sval for sval in line.split()] for line in fin]
    return result
#---------------------------------------------------------------------

# check forbidden consecutive patches
# THIS IS FOR RES1-PATCH1-RES2-PATCH2 combination
# Onlu patch1 and patch2 are important. rest is checked in
# check_constraints 
def is_forbid_patch(patchname1,patchname2,patforbid):
    flag = 0 # default not forbidden
    for i in range(len(patforbid)):
        if patforbid[i][0] == patchname1:
            print('Forbidden line', patforbid[i])
            if any(patchname1 in st for st in patforbid[i]):
                flag = 1
        
    return flag
#---------------------------------------------------------------------

# Generate patches
def create_patches(flist,nmons,nch,segname,inp_dict,cumulprobarr\
                    ,tol,maxattmpt,flog,ctr_flag,ctrfyle,residlist,\
                   patforbid,graft_opt):

    # Write list to a separate file
    flist.write(';# Entire patch list\n')
    flist.write(';# num_patches\t%d, num_chains\t%d\n' %(nmons-1,nch))
    flog.write('Probabilities for each attempt\n')
    flog.write('Attempt#\t')
    for wout in range(len(inp_dict)):
        flog.write('%s (%g)\t' %(list(inp_dict.keys())[wout],\
                                 list(inp_dict.values())[wout]))
    flog.write('L2norm \n')

    flag_optimal = -1

    for attnum in range(maxattmpt):
        flog.write('%d\t' %(attnum+1))    
        flist.write(';# Attempt number \t%d\n' %(attnum+1))
        out_list = [[] for i in range(nch)] #reset every attempt
   
        for chcnt in range(nch):

            flist.write(';# chain number:\t%d\n' %(chcnt+1))
            flist.write(';# -- Begin patches for %s ---\n' %(segname))
            patcnt = 0
            branched = 0

            while patcnt <= nmons-2: #for checking constraints

                # Special Case 1: Any residue is PCA (Grafted one)
                # If residue is a grafted one, then the patch is
                # trivial. No need to check for probabilities. Append
                # directly. However, need to reorder patches
                # connecting the n-1 and n+1 residue
                if residlist[patcnt] == graft_opt[1]:
                    patchname = list(graft_opt[2])
                    branched = 1
                    flist.write(' patch\t%d\t%s\t%s:%d\t%s:%d\n' \
                                %(patcnt+1,patchname,\
                                  segname,patcnt+1,segname,patcnt+2))
                    out_list[chcnt].append(patchname)
                    patcnt += 1

                    continue #continue to next residue

                #Special Case 1: Last residue is a pCA.
                #The bond between the last and penultimate residue
                #should be GOG
                if patcnt == nmons-2 and \
                   residlist[patcnt+1] == graft_opt[1]: 
                    
                    if branched == 1:
                        print('ERR: Two consecutive graft resiudes',\
                              residlist[patcnt],residlist[patcnt+1])
                        return -1

                    patchname = list(graft_opt[2])
                    branched = 1
                    flist.write(' patch\t%d\t%s\t%s:%d\t%s:%d\n' \
                                %(patcnt+1,patchname,\
                                  segname,patcnt+1,segname,patcnt+2))
                    out_list[chcnt].append(patchname)
                    patcnt += 1

                    continue #continue to next chain
                    
                
                # If first condition is not met, generate a random
                # number and do the following
                ranval = random.random() #seed is current system time by default
                findflag = 0

                for arrcnt in range(len(cumulprobarr)):
        
                    #Only need to check the less than value because
                    #the array is organized in increasing order.
                    #Break the loop once the first point where the
                    #condition is met.
                    if ranval < cumulprobarr[arrcnt]:
                        patchname = list(inp_dict.keys())[arrcnt]
                        findflag = 1

                        # Add constraint flags: default to 1
                        #so that if constraints are not there, it will
                        #be appended. consec flag has to be 0 for true
                        appendflag = 1; consecpatflag = 0 
                        if ctr_flag:
                            if patcnt == 0:
                                resname1 = 'None'
                                patchname1 = 'None'
                            else:
                                resname1 = residlist[chcnt][patcnt]
                                patchname1 = out_list[chcnt][patcnt-1]
                        
                            # end if patcnt == 0
                            resname2 = residlist[chcnt][patcnt+1]
                            appendflag = check_constraints(ctrfyle,patchname,\
                                                          resname1,resname2)

                            consecpatflag =is_forbid_patch(patchname,\
                                                           patchname1,\
                                                           patforbid)

                        # end if ctr_flag==1

                        if appendflag == 1 and \
                           consecpatflag == 0: #update while loop if
                            #constraints are met
                            out_list[chcnt].append(patchname)

                            if branched == 0: #patch between n and n+1
                                flist.write(' patch\t%d\t%s\t%s:%d\t%s:%d\n' \
                                            %(patcnt+1,patchname,\
                                              segname,patcnt+1,segname,patcnt+2))
                            else: # patch between n-1 and n+1
                                flist.write(' patch\t%d\t%s\t%s:%d\t%s:%d\n' \
                                            %(patcnt+1,patchname,\
                                              segname,patcnt,segname,patcnt+2))
                                branched = 0 # reset branch flag

                            patcnt += 1
                        elif appendflag == -2:
                            return 
                        # end if appendflag == 1

                        break
                    # end ranval < cumulprobarr[]

                # end for arrcnt in range(len(cumulprobarr))

                if findflag != 1:
                    print('Random value/Probarr:', ranval,cumulprobarr)
                    print('Error: Did not find a random residue\n')
                    return -1
                # end if find flag

             # end while loop
            flist.write(';# --End patch list for %d--\n' %(chcnt+1))

        # end for chcnt in range(nch)

        # After going through all the chains, count occurence of each res/patch
        outdist = []
        for key in inp_dict:
            outdist.append(sum([i.count(key) for i in out_list]))

        #normalize
        sumval = sum(outdist)
        if sumval != nch*(nmons-1):
            print('Sum from distn,nch*(nmons-1):',sumval,nch*(nmons-1))
            exit('ERROR: Sum not equal to the total # of patches')
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
def write_segments_onego(fin,nmons,nch,segname,res_list,patch_list):

    fin.write(';# Writing % segments' %(nmons))
    fin.write(' resetpsf \n')
    fin.write(' segment %s {\n' %(segname))
    
    #Residues
    for chcnt in range(nch):

        for rescnt in range(nmons):

            fin.write('  residue  %d  %s\n' \
                      %(rescnt+1,res_list[chcnt][rescnt]))

    fin.write('}')        
    fin.write('\n')

    #Patches
    for chcnt in range(nch):

        for patcnt in range(nmons-1):

            fin.write('patch  %s  %s:%d  %s:%d\n' \
                        %(patch_list[chcnt][patcnt],\
                          segname,patcnt+1,segname,patcnt+2))

    fin.write('\n')
#---------------------------------------------------------------------

# Write residues/patches iteration by iteration
def write_multi_segments(fin,iter_num,nmonsthisiter,nch,chnum,\
                         segname,res_list,patch_list):

    if iter_num == 1:
        fin.write(';# Chain number: %d of %d chains' %(nch,chnum))

    fin.write(';# Iteration number: %d\n' %(iter_num))
    fin.write('set count %d' %(nmonsthisiter))
    fin.write('\n')
    fin.write(' resetpsf \n')
    fin.write(' segment %s {\n' %(segname))

    #Residues -- indices should have -1 for first dimension
    for rescnt in range(nmonsthisiter):

        fin.write('  residue  %d  %s\n' %(rescnt+1,\
                                          res_list[chnum-1][rescnt]))

    fin.write('}')        
    fin.write('\n')
        
    #Patches -- indices should have -1 for first dimension
    for patcnt in range(nmonsthisiter-1):

        fin.write('patch  %s  %s:%d  %s:%d\n' \
                  %(patch_list[chnum-1][patcnt],segname,patcnt+1,\
                    segname,patcnt+2))

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
