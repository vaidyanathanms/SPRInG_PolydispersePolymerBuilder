#---------------------------------------------------------------------------
# Supporting scripts for SPRInG
# This file contains auxiliary function definitions
# Main file: genconf.py
# 'None' is a reserved keyword- DONT USE IT for PDB/PSF filenames
#---------------------------------------------------------------------------

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
import subprocess
import pandas as pd
#---------------------------------------------------------------------------

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
    return 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,3,2.0,50,0.1
#---------------------------------------------------------------------

# Read data for experimental PDI system
def read_expt_pdidata(words):
    mon_mwt = 200; expt_mn = 0; expt_mw = 0; expt_pdi = 0
    npdiatt = 100000; pditol = 0.05
    if len(words) == 3:
        ex_disper_fyle = words[2]
    elif len(words) > 3:
        for wcnt in range(3,len(words),2):
            if words[wcnt] == 'avg_monwt'.lower():
                mon_mwt = float(words[wcnt+1])
            elif words[wcnt] == 'mn'.lower():
                expt_mn = float(words[wcnt+1])
            elif words[wcnt] == 'mw'.lower():
                expt_mw = float(words[wcnt+1])
            elif words[wcnt] == 'pdi'.lower():
                expt_pdi = float(words[wcnt+1])
            elif words[wcnt] == 'ntrials'.lower():
                npdiatt = float(words[wcnt+1])
            elif words[wcnt] == 'pditol'.lower():
                pditol = float(words[wcnt+1])/100
    else:
        raise RuntimeError("Not enough arguments", words,\
                           len(words))

    return ex_disper_fyle,mon_mwt,expt_mn,expt_mw,expt_pdi,\
        npdiatt,pditol
#---------------------------------------------------------------------

def create_new_pdidata(words,line):
    if len(words) < 5:
        raise RuntimeError('Not enough arguments for PDI: '+line)
    inp_pdival = float(words[2])
    if inp_pdival <= 1.0:
        raise RuntimeError('ERR: PDI for polydisperse cases should be > 1.0')
    disper_fyle = words[3]
    npdiatt = int(words[4])
    pditolval = 0; distrange = 0
    if (len(words)-5)%2 != 0:
        raise RuntimeError('ERR: Unknown number of args for disperse')
    for wcnt in range(int(0.5*(len(words)-5))):
        if words[2*wcnt+5].lower() == 'pditol'.lower():
            pditolval = float(words[2*wcnt+6])
        elif words[2*wcnt+5].lower() == 'mwrange'.lower():
            distrange = int(words[2*wcnt+6])
        else:
            exit('ERR: Unknown keyword' + words[2*wcnt+5])
            
    return inp_pdival, disper_fyle, npdiatt, pditolval, distrange
#----------------------------------------------------------------------

# Check all flags 
def check_all_flags(casenum,fresflag,fpatflag,disflag,M,N,\
                    fnamd,fpdbflag,ftopflag):
    outflag = 1
    if casenum < 0:
        print('ERR: Case number not input'); outflag = -1
    elif N == 0:
        print('ERR: No chains found in input'); outflag = -1
    elif disflag == 0 and M == 0:
        print('ERR: Monodisperse systems with no MW'); outflag = -1
    elif fresflag == 0:
        print('ERR: Residue input not entered'); outflag = -1
    elif fpatflag == 0:
        print('ERR: Patch not entered'); outflag = -1
    elif ftopflag == 0:
        print('ERR: Topology file not found'); outflag = -1
    elif fnamd != 0 and fpdbflag == 0:
        print('ERR: To run NAMD, input PDB/top files are required')
        outflag = -1
    return outflag
#---------------------------------------------------------------------

# Define headers for psf files
def psfgen_headers(fin,topname,outname):
    fin.write('; ##*********New Molecule/Segment*********##\n')
    fin.write(';# headers and inputs \n')
    fin.write('package require psfgen \n')
    topinp = topname
    fin.write('%s\t %s\n' %('topology',topinp))
#---------------------------------------------------------------------              
# Details for closing input files
def psfgen_postprocess(fin,writetype,iter_num,segname,fnamdflag,\
                       basic_pdb):
    # outname is already there. no need again
    fin.write(';# Writing output \n')
    fin.write('regenerate angles dihedrals \n')

    if fnamdflag:
        if writetype == 'single':
            comnt = 'Reference PDB'
            comnt2 = 'Guesses rest of the coordinates from PDB inp'
            pdbfyle =  '../' + basic_pdb
        elif writetype == 'multi':
            if iter_num == 1:
                comnt = 'Use reference PDB in the first iteration'
                pdbfyle =  '../' + basic_pdb
            else:
                comnt = '*.coor is the file generated by NAMD'
                pdbfyle = '$outputname.coor'

            comnt2 = 'Can create steric clashes and hence the iterations.'
        else:
            exit('ERROR: Unknown option: ' + writetype)
        fin.write('coordpdb %s  %s  ;#  %s\n' %(pdbfyle,segname,comnt))
        fin.write('guesscoord ;#  %s\n' %(comnt2))        
        fin.write('writepdb $outputname.pdb \n')
        if writetype == 'multi':
            fin.write('writepdb ${outputname}_${count}.pdb \n')#backup


    fin.write('writepsf $outputname.psf \n')
    fin.write('; #--------- End of Molecule/Segment ----------#\n')
    fin.write('\n')
    if writetype == 'multi':
        fin.write('writepdb ${outputname}_${count}.psf \n')
#---------------------------------------------------------------------

# Read monomer ratios from input file
def residue_ratios(inpfyle):

    frac_res = collections.OrderedDict()

    # Check file existence
    if not os.path.exists(inpfyle):
        print('Residue input file not found \n', inpfyle)
        return -1

    # add monomer details
    with open(inpfyle) as fyle_dict:
        for line in fyle_dict:
            line = line.strip()
            (key, val) = line.split()
            frac_res[key] = float(val)

    return frac_res
#---------------------------------------------------------------------

# Define patch ratios from input file
def patch_ratios(opt_branch,resdict,inpfyle):

    frac_patch = collections.OrderedDict()

    # Check file existence
    if not os.path.exists(inpfyle):
        print('Residue input file not found \n', inpfyle)
        return -1

    # add patch details
    with open(inpfyle) as fyle_dict:
        for line in fyle_dict:
            line = line.strip()
            (key, val) = line.split()
            frac_patch[key] = float(val)

    # check for branches and rearrange dictionary
    if opt_branch[0] != 1:
        return frac_patch

    elif opt_branch[0] == 1:
        sumrestot = 0
        for rescnt in range(len(resdict)):
            sumrestot += list(resdict.values())[rescnt]
        newfrac_patch = collections.OrderedDict() #create new dict
        grcnt = 1; branch_prob = 0
        while grcnt < len(opt_branch):
            gr_resname = opt_branch[grcnt]
            gr_patname = opt_branch[grcnt+1]
            resflag = 0
            for rescnt in range(len(resdict)):
                if list(resdict.keys())[rescnt] == gr_resname:
                    resflag = 1
                    branch_prob += list(resdict.values())[rescnt]
                    frac_patch[gr_patname] = branch_prob
                    newfrac_patch[gr_patname] = branch_prob
            grcnt = grcnt + 2 
            if resflag == 0:
                print('ERROR: Could not find ', str(gr_resname))
                return 0
        #sum of norm branch res prob = sum of norm branch patch prob
        branch_prob = branch_prob/sumrestot # normalize branch probs
        # Renormalize if branches are present and create new dict
        sumprob = 0
        for patcnt in range(len(frac_patch)):
            if list(frac_patch.keys())[patcnt] not in opt_branch:
                sumprob += list(frac_patch.values())[patcnt]

        normval = sumprob/(1-branch_prob)

        for patcnt in range(len(frac_patch)):
            if list(frac_patch.keys())[patcnt] not in opt_branch:
                newprob = list(frac_patch.values())[patcnt]/normval
                keyval = list(frac_patch.keys())[patcnt]
                newfrac_patch[keyval] = newprob

        return newfrac_patch

    else:
        print('Unknown option', branch_opt[0])
        return 0
#---------------------------------------------------------------------

# Develop probability distribution for experimental data
def make_expt_pdidata(einp_fyle,nchains,nattempts,pditol):

    expinp_fmt ='NULL'; exptkey = 0
    # Check file existence
    if not os.path.exists(einp_fyle):
        raise RuntimeError('Expt MW distribution file not found \n', inpfyle)

    # Read headers
    df = pd.read_csv(einp_fyle)
    if 'molwt'.lower() not in df.columns:
        raise RuntimeError("Keyword molwt is missing in the header")
    for col in df.columns:
        if 'wlogmw'.lower() == col.lower():
            expinp_fmt = 'WLOGMW'; exptkey +=1
        elif 'wmw'.lower()  == col.lower():
            expinp_fmt = 'WMW'; exptkey += 1
        elif 'pmw'.lower()  == col.lower():
            expinp_fmt = 'PMW'; exptkey += 1
        else:
            raise RuntimeError("Unknown keyword: " + col)
    if exptkey > 1:
        raise RuntimeError("Multiple distribution values cannot be given")

    # Remove duplicates
    df.drop_duplicates(inplace = True)

    # Check for strictly increasing
    xdata = np.array(df['molwt'])
    if not np.all(np.diff(xdata)>0):
        raise RuntimeError('X-data should be monotonically increasing')

    # Convert w(m)/w(logm) to p(m)
    ydata = convert_to_pofm(xdata,np.array(df[expinp_fmt]),\
                            expinp_fmt)

    # Compute Mn, Mw, Mz
    intsum,emn,emw,emz = comp_avgs(xdata,ydata)
    epdi = emw/emn

    # Generate output distributions
    cmn,cmw,cmz,cpdi,eout_file = gen_exptdist(xdata,ydata,intsum,\
                                              nchains,emn,emw,\
                                              emz,epdi,nattempts,\
                                              pditol,einp_file)

    # Output to log file
    expout_log(flog,einp_file,expinp_fmt,emn,emw,emz,epdi,\
               cmn,cmw,cmz,cpdi,eout_file)

    return eout_file
#---------------------------------------------------------------------

# Convert inputs to probability distribution
def convert_to_pofm(xinp,yinp,inpfmt):
    if inpfmt == 'WLOGMW':
        return np.multiply(np.power(xinp,-2),yinp))
    elif inpfmt == 'WMW':
        return np.multiply(np.power(xinp,-1),yinp))
    elif inpfmt == 'PMW':
        return yinp
#---------------------------------------------------------------------

# Compute averages
def comp_avgs(xinp,yinp):
    mu0 = trapz(xinp,yinp)
    mu1 = trapz(xinp,np.multiply(xinp,yinp))
    mu2 = trapz(xinp,np.multiply(xinp**2,yinp))
    mu3 = trapz(xinp,np.multiply(xinp**3,yinp))
    return mu0, mu1/mu0, mu2/mu1, mu3/mu2
#---------------------------------------------------------------------

# Trapezoidal rule
def trapz(xinp,yinp):
    sval = 0
    for indx in range(0,len(xinp)-1)):
        sval += 0.5*(yinp[indx]+yinp[indx+1])*(xinp[indx]-xinp[indx+1])
    return sval
#---------------------------------------------------------------------

# Generate MW distribution and write to output
def gen_exptdist(xinp,pdfy,intsum,emn,emw,emz,epdi,nch,natt,\
                 pditol,inpfyle):

    # Normalize distribution
    normy = pdfy/intsum
    #Generate cumulative distribution
    cdf[0] = normy[0]
    for indx in range(1,len(xinp))):
        cdf.append(trapz(xinp[0:indx],normy[0:indx]))

    if abs(cdf[len(cdf)-1]-1) > 10**-6:
        raise RuntimeError('CDF not adding to one: ',cdf[len(cdf)=1])
    # Generate chains
    print('Generating chains according to expt distribution..')

    for trials in range(natt):
        if trials%1000 == 0:
            print('Trial number: ', trials+1)
        mw_vals = [np.interp(random.random,normy,xinp) \
                   for j in range(nch))]
        comp_mn  = sum(mw_vals)/nch
        comp_mw  = sum(mw_vals**2)/sum(mw_vals)
        comp_pdi = comp_mw/comp_mn
        if abs((comp_pdi-epdi)/epdi) < pditol:
            print('Found chain configuration!..')
            print('Computed Mn: %g, Mw: %g, PDI: %g' \
                  %(comp_mn,comp_mw,comp_pdi))
            comp_mz = sum(mw_vals**3)/sum(mw_vals**2)
            break
    
    if trails > natt:
        raise RuntimeError('Could not find molecular weights of chains correspond to experimental distribution. Try increasing the number of chains or the tolerance')

    # Write output distributions file
    edist_fyle = "pdidistribution_"+inpfyle
    with open(edist_fyle,'w') as fdist:
        fdist.write('# Computed/Experimental Mn: %g, %g\n' %(comp_mn,emn))
        fdist.write('# Computed/Experimental Mw: %g, %g\n' %(comp_mw,emw))
        fdist.write('# Computed/Experimental Mz: %g, %g\n' %(comp_mz,emz))
        fdist.write('# Computed/Experimental PDI: %g, %g\n' %(comp_pdi,epdi))
        fdist.write('%s\t%s\t%s\t%s\t%s\n' %('molwt, pMW, cMW, wMW, wlogMW'))

        for indx in range(0,len(xinp)-1)):
            fdist.write('%g\t%g\t%g\t%g\t%g\n' \
                        %(xinp[indx],pdfy[indx],cdf[indx],\
                          xinp[indx]*pdfy[indx],\
                          xinp[indx]*xinp[indx]*pdfy[indx]))

    # Write output MW file
    eout_fyle = "polydisp_"+inpfyle
    with open(eout_fyle,'w') as fmwvals:
        fmwvals.write('num_chains\t%d\n' %(nch))
        for j in range(nch):
            fmwvals.write('%d\n',mwvals[j]))
    return comp_mn, comp_mw, comp_mz, comp_pdi, eout_fyle
#---------------------------------------------------------------------

def expout_log(flog,einp_file,expinp_fmt,emn,emw,emz,epdi,\
               cmn,cmw,cmz,cpdi,eout_file):
    # Write to log file 
    flog.write('Experimental input data file: %s\n' %(einp_file))
    flog.write('Input prob distribution type: %s\n'%(expinp_fmt))
    flog.write('***********From user input distribution****\n')
    flog.write('Mn: %g; Mw: %g; Mz: %g; PDI: %g\n' %(emn,emw,emz,epdi))
    flog.write('Computed PDI distribution file: %g\n' %(eout_file))
    flog.write('***********From computed distribution******\n'))
    flog.write('Mn: %g; Mw: %g; Mz: %g; PDI: %g\n' %(cmn,cmw,cmz,cpdi))
#---------------------------------------------------------------------

# Initial PDI details if polydisperse chains are to be generated
def init_pdi_write(pdival,avgmw,nch,op_file,npdiatt,pditolval,\
                   min_polysize,distrange):
    pdi_fyl = 'inp_genpdi.txt'
    finit   = open(pdi_fyl,'w')
    finit.write('chain_types\n')
    finit.write('%d\n' %(1)) # for now one chain type
    finit.write('chain_details\n')
    finit.write('%g\t %g\t %g\n' %(pdival,avgmw,nch))
    finit.write('max_attempts\n')
    finit.write('%d\n' %(npdiatt))
    finit.write('min_size\n')
    finit.write('%d\n' %(min_polysize))
    if pditolval != 0:
        finit.write('tolerance\n')
        finit.write('%g\n' %(pditolval))
    if distrange != 0:
        finit.write('distrange\n')
        finit.write('%d\n' %(distrange))
    finit.write('pdi_op_file\n')
    finit.write('%s\n' %(op_file))
    finit.close()
#---------------------------------------------------------------------

# Compile and run PDI. 
def compile_and_run_pdi(destdir):
    
    if not os.path.exists('inp_genpdi.txt'):
        print('inp_genpdi.txt not found')
        return -1
        
    # Generate PDI data
    print("Compiling program to generate polydisperse chains...")
    print("Takes about 10 seconds..")
    if shutil.which("ifort") != None:
        subprocess.call(["ifort","-r8","pdi_dist_params.f90","pdigen.f90",\
                         "-o","pdiinp.o"])
    elif shutil.which("gfortran") != None:
       subprocess.call(["gfortran","-freal-4-real-8",\
                    "pdi_dist_params.f90","pdigen.f90","-o","pdiinp.o"])
    else:
        raise RuntimeError("No Fortran 90 compiler found!")

    print("Compilation successful..")
    subprocess.call(["./pdiinp.o", "inp_genpdi.txt"])
    return 1
#---------------------------------------------------------------------

# Assign MW for polydisperse cases
def make_polydisp_resids(inpfyle, nch, min_polysize):
    if not os.path.exists(inpfyle):
        print('ERR: PDI file: ', inpfyle, 'not found')
        return -1, 0
    chflag = 0
    with open(inpfyle) as fyle_pdi:
        for line in fyle_pdi:
            if line.startswith('#'): # skip lines starting with #
                continue
            if not line: # skip empty lines
                continue
            line = re.split('\W+',line.strip())
            if chflag == 0:
                if len(line) != 2 or line[0] != 'num_chains':
                    print('ERR: Unknown first line in: ', inpfyle)
                    print(line, '\n', len(line), line[0])
                    return -1, 0
                numch = int(line[1])
                chflag = 1
                resmw_data = []
            else:
                if int(line[0]) < min_polysize:
                    print('ERR: Minimum', min_polysize, 'residues should be present')
                    return -1, 0
                resmw_data.append(int(line[0]))

    if len(resmw_data) != nch or nch != numch:
        print('ERR: Mismatch in number of chains')
        print(len(resmw_data), nch, numch)
        return -1, 0

    num_avg_mw = 0; wt_avg_mw = 0
    for mws in range(len(resmw_data)):
        num_avg_mw += resmw_data[mws]
        wt_avg_mw  += pow(resmw_data[mws],2)
    
    wt_avg_mw  = wt_avg_mw/num_avg_mw
    num_avg_mw = num_avg_mw/nch
    pdiout = wt_avg_mw/num_avg_mw

    return resmw_data, pdiout
#---------------------------------------------------------------------

# Initiate log file
def init_logwrite(flog,casenum,bmtype,Marr,tfile,segname,\
                  nch,att,tol,opstyle,fpres_constraint,fpp_constraint,\
                  resfyle,patfyle,disflag,pdiinp):
    flog.write('Case number: %d\n' %(casenum))
    flog.write('Creating TCL file for %s\n' %(bmtype))

    if disflag == 0:
        flog.write('Monodisperse system \n')
        flog.write('Num Chains/num Residues: %d\t%d\n'%(nch,Marr[0]))
    else:
        flog.write('Polydisperse system \n')
        for i in range(len(Marr)):
            flog.write('Chain#/Num Residues: %d\t%d\n' %(i+1,Marr[i]))

    flog.write('PDI: %g\n' %(pdiinp))
    flog.write('Tot res/pat: %d\t%d\n' %(sum(Marr),sum(Marr)-len(Marr)))
    flog.write('Res/patch inputs: %s\t%s\n' %(resfyle,patfyle))
    flog.write('Input Topology file file: %s\n' %(tfile))
    flog.write('Segment name in input (or output prefix): %s\n' \
               %(segname))
    flog.write('#attempts/Tolerance: %d\t%g\n' %(att,tol))
    if fpres_constraint != 0 :
        flog.write('Patch-residue constraints: Yes\n')
    else:
        flog.write('Patch-residue constraints: No\n')
    if fpp_constraint != 0:
        flog.write('Patch-patch constraints: Yes\n')
    else:
        flog.write('Patch-patch constraints: No\n')
    flog.write('Output style: %s\n' %(opstyle))
    
    flog.write('Analysis beginning ..\n')
#---------------------------------------------------------------------

# Check initial files
def find_init_files(fpres_constraint,fpp_constraint,fpdbflag,\
                    fnamdflag,flbdflag,makepdifile,input_top='none',\
                    input_pdb='none',input_pres='none',input_pp='none',\
                    input_lbd='none'):
    # Read defaults and throw exceptions
    if not os.path.exists(input_top):
        print('Topology file not found \n', input_top)
        return -1
    elif not os.path.exists(input_lbd):
        print('Parameter file for LigninBuilder not found \n',\
              input_lbd)
        return -1
    elif fnamdflag:
        if fpdbflag and not os.path.exists(input_pdb):
            print('Initial structure file not found \n', input_pdb)
            return -1
    elif fpres_constraint == 1 and not os.path.exists(input_pres):
        print('Patch-residue constraint file not found', input_pres)
        return -1
    elif fpp_constraint == 1 and not os.path.exists(input_pp): 
        print('Patch-patch constraint file not found', input_pp)
        return -1
    elif makepdifile == 1:
        if not os.path.exists('pdigen.f90') or \
           not os.path.exists('pdi_dist_params.f90'):
            print('Source file to generate PDIs not found')
            return -1
    return 1
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
        sumval = dummy_distarr[len(dummy_distarr)-1]
        
        # force normalization
        for cnt in range(len(dummy_distarr)):
            dummy_distarr[cnt] = dummy_distarr[cnt]/sumval
            
        print('New distribution: ', dummy_distarr)
            
    else:
        print('Generated target cumulative distribution..')

    return dummy_distarr
#---------------------------------------------------------------------

# Check whether the input pdb file is consistent with the inputs given
# for generating the tcl file
def check_pdb_defaults(inpfyle,defa_res,seginp):
    flag = 1 # default true
    resnum = 1
    # Check whether pdb file contains default segment and segment name
    with open(inpfyle) as fpdbin:
        for line in fpdbin:
            line = line.rstrip('\n')
            all_words = re.split('\W+',line)
            if all_words[0] == 'ATOM':
                lenwords = len(all_words)
                if all_words[lenwords-1] != seginp and \
                   all_words[lenwords-2] != seginp:
                    print(all_words[lenwords-1],all_words[lenwords-2])
                    print('Did not find ',seginp,'in ',line)
                    flag = -1
                    break
                if defa_res not in all_words:
                    print('Did not find ',defa_res,'in ',line)
                    flag = -1
                    break
                if all_words[4].isdigit():
                    if int(all_words[4]) > resnum:
                        print('WARNING: More than one res found: ',\
                              resnum)
                        resnum = int(all_words[4])
                else:
                    print('ERR: Unknown value for chain num',\
                          all_words[4])
                    flag = -1


    return flag
#---------------------------------------------------------------------
    
# Create entire list in one go so that cumulative distribution holds true
def create_residues(flist,nresarr,nch,segpref,inp_dict,cumulprobarr\
                    ,tol,maxattmpt,flog,branchopt,defa_res,res_terminator):

    # Write list to a separate file
    flist.write(';#  Entire segment list\n')
    for i in range(nch):
        flist.write(';#  num_resds\t%d, chain#\t%d\n' \
                    %(nresarr[i],i+1))

    sum_of_res = sum(nresarr)
    flist.write('; Total number of residues\t%d\n' %(sum_of_res))
    flog.write('Probabilities for each attempt\n')
    flog.write('Attempt#\t')
    if defa_res != 'none':
        if defa_res not in list(inp_dict.keys()):
            print('FATAL ERR: default residue not in the input')
            print(defa_res, list(inp_dict.keys()))
            return -1

    for wout in range(len(inp_dict)):
        flog.write('%s (%g)\t' %(list(inp_dict.keys())[wout],\
                                 list(inp_dict.values())[wout]))

    flog.write('L2norm \n')

    flag_optimal = -1; oneconfigflag = -1; normold = 999999999

    for attnum in range(maxattmpt):

        flog.write('%d\t' %(attnum+1))    
        flist.write(';# Attempt number \t%d\n' %(attnum+1))
        flist.write(' resetpsf\n')

        out_list = [[] for i in range(nch)] #reset every attempt
   
        for chcnt in range(nch):
            segname = ret_segname(segpref,chcnt+1)
            flist.write(';# chain number:\t%d\n' %(chcnt+1))
            flist.write(' segment %s {\n' %(segname))
            # first is default residue if present
            if defa_res != 'none':
                flist.write(' residue\t%d\t%s\n' %(1,defa_res))
                out_list[chcnt].append(defa_res)
                rescnt = 1
            else:
                rescnt = 0

            deg_poly_chain = nresarr[chcnt]
            while rescnt < deg_poly_chain:

                ranval = random.random() #seed is current system time by default
                findflag = 0
                consecresflag = 0 #default: consecutive res are NOT found
                initres_flag = 0 #def: no initiator res present
                endbranchflag = 0 #def: if last is branch, two prior
                #residues are NOT branch
                for arrcnt in range(len(cumulprobarr)):
        
                    #Only need to check the less than value because
                    #the array is organized in increasing order.
                    #Break the loop once the first point where the
                    #condition is met.

                    if ranval < cumulprobarr[arrcnt]:
                
                        findflag = 1   
                        resname1 = list(inp_dict.keys())[arrcnt]

                        # Conditions for not first and last residues
                        if rescnt != 0: 
                            resname2 = out_list[chcnt][rescnt-1]
                            consecresflag = is_res_cons(resname1,resname2\
                                                        ,branchopt)

                        #terminator condition if present
                        if resname1 == res_terminator:
                            if rescnt != deg_poly_chain-1:
                                initres_flag = 1 #can only be at end
                            elif rescnt == deg_poly_chain-1 and \
                               (resname2 in branchopt):
                                #if end, previous cannot be branch
                                #(this condition is already met by
                                #following if statement. Just to
                                #be extra safe)
                                initres_flag = 1 

                        # Last residue and second last residue cannot
                        # be branches
                        if rescnt == deg_poly_chain-2 or \
                           rescnt == deg_poly_chain-1:
                            if (resname1 in branchopt):
                                endbranchflag = 1                           

#----------------------old version in master branch-----------------------
#+                        # If the last residue is a branch, the previous
#+                        # TWO resiudes cannot be branch
#+                        if rescnt == deg_poly_chain-1 and \
#+                           (resname1 in branchopt):
#+                            resname3 = out_list[chcnt][rescnt-1]
#+                            endbranchflag = is_res_cons(resname1,\
#+                                                       resname3,branchopt)
#---------------------------------------------------------------------------


                        if consecresflag == 0 and initres_flag == 0 \
                           and endbranchflag == 0:

                            flist.write(' residue\t%d\t%s\n' \
                                        %(rescnt+1,resname1))
                            out_list[chcnt].append(resname1)
                            rescnt = rescnt + 1
                            
                        break

                if findflag != 1:
                    print('Random value/Probarr:', ranval,cumulprobarr)
                    exit('Error in finding a random residue\n')
            
            flist.write(' }\n')

        # After going through the chains, count each residue
        outdist = []
        for key in inp_dict:
            outdist.append(sum([i.count(key) for i in out_list]))

        # Sum and normalize
        sumval = sum(outdist)
        if sumval != sum_of_res:
            print('Sum from distn,sum_of_res:',sumval,sum_of_res)
            exit('ERROR: Sum not equal to the total # of residues')
        normlist = [x/sumval for x in outdist]

        # Extract target probabilities and compare
        targ_probs = list(inp_dict.values())
        normval = numpy.linalg.norm(numpy.array(normlist) \
                                    - numpy.array(targ_probs))
    
        # Write to log file
        for wout in range(len(outdist)):
            flog.write('%g\t' %(outdist[wout]/sumval))
        flog.write('%g\n' %(normval))

        if normval <= tol:
            #write to log file
            flog.write('Found optimal residue configuration\n')
            print('Found optimal residue configuration..')
            flag_optimal = 1
            return out_list
        elif normval < normold:
            if oneconfigflag == -1:
                oneconfigflag = 1
                print('Found configuration with res_err: ',normval)
            else:
                print('Updating configuration with res_err: ',normval)
            backup_list = [] #create new_backup list
            backup_list = out_list.copy()
            flist.write('\n')
            normold = normval

    if oneconfigflag == -1:
        print('Could not find a residue list with constraints')
        return -1

    if flag_optimal == -1:
        print('Did not find optimal residue configuration')
        print('Using best residue configuration with L2norm: ',normold)
        flog.write('Did not find optimal residue configuration\n')
        flog.write('Using best configuration with residue L2norm: %g\n'\
                   %(normold))
        return backup_list

#---------------------------------------------------------------------

# Read and check patch constraints -- May not be effective as opposed
# to reading at once and copying to array. Need to think about it.
# Check special cases using files
# THIS IS CONSTRAINT FOR RESIDUE1-PATCH-RESIDUE2 combination
def check_constraints(inpfyle,patchname,resname1,resname2):
    
    bef_flag = 1; aft_flag = 1 # keep as true; default: OK
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
                if all_words[2] == resname2:
                    aft_flag = 0


    if resname1 == 'SYR' and patchname == '55':
        if bef_flag != 0:
            print("ERRRR:", bef_flag)
    if resname2 == 'SYR' and patchname == '55':
        if aft_flag != 0:
            print("ERRRRR:", aft_flag)
    # Return 0 if any flags are 0, else return 1
    if bef_flag == 0 or aft_flag == 0:
        return 0
    else:
        return 1
#---------------------------------------------------------------------

# return segment name according to the length of segpref
def ret_segname(seginp,chval):
    if len(seginp) + len(str(chval)) > 4:
        print('WARNING: Renaming segment names')
        lval = 4-len(str(chval))
        if lval >= 0:
            segout = seginp[0:lval] + str(chval)
        else:
            segout = str(chval)
    else:
        segout = seginp + str(chval)
    return segout
#---------------------------------------------------------------------

# check consecutive residues - cannot have branch residue in
# consecutive positions
def is_res_cons(resname1,resname2,branchopt):
    sameflag = 0
    if branchopt[0] == 1:
        if (resname1 in branchopt) and (resname2 in branchopt):
            sameflag = 1
    return sameflag
#---------------------------------------------------------------------

# read all patch incompatibilities
def read_patch_incomp(fname):
    with open(fname,'r') as fin:
        result = [[sval for sval in line.split()] for line in fin]
    return result
#---------------------------------------------------------------------

# Generate patch rules: patch - m, residues - n
# Structure: RESn-PATm-RESn+1-PATm+1...
# Rule 1: (a)patch "m" between resids n/n+1; check is_forbid_pat(m,m-1)
# (b) check constraints for m between n/n+1
# Rule 2: res_n = branch; (a) branch_patch m between n/(n+1); check
# is_forbid_pat(m-1,m+1); 
# Rule 3: if res_n+1 = branch, patch m between n/n+2; check
# is_forbid_pat(m, m-1) => same as rule 1 (VERY IMP); check constrants
# between n and n+2 (VERY IMP); and is_forbid_pat(m-1,m). Josh
# suggested that GOG-BB (branch-normalpatch) is not compatible.
# Rule 4: if last resiudue is a branch. branch_patch m between n/n-1; no
# checks required
# For left and right residues, reference "n" corresponds to left
# residue. For examples in comments, S/G/H are normal resids, F/PCA
# are branch residue. 
def create_patches(flist,nresarr,nch,segpref,inp_dict,cumulprobarr\
                   ,tol,maxattmpt,flog,fpresflag,fppflag,pres_fyle,\
                   residlist,patforbid,branch_opt):

    # Write list to a separate file
    flist.write(';# Entire patch list\n')
    for i in range(nch):
        flist.write(';#  num_patches\t%d, chain#\t%d\n' \
                    %(nresarr[i]-1,i+1))

    sum_of_res = sum(nresarr)
    sum_of_pat = sum_of_res - nch
    flist.write(';# Total number of patches\t%d\n' %(sum_of_pat))

    flog.write('Probabilities for each attempt\n')
    flog.write('Attempt#\t')
    
    for wout in range(len(inp_dict)):
        flog.write('%s (%g)\t' %(list(inp_dict.keys())[wout],\
                                 list(inp_dict.values())[wout]))
    flog.write('L2norm \n')
    flag_optimal = -1; oneconfigflag = -1; normold = 999999999

    for attnum in range(maxattmpt):

        flog.write('%d\t' %(attnum+1))    
        flist.write(';# Attempt number \t%d\n' %(attnum+1))
        out_list = [[] for i in range(nch)] #reset every attempt
        chcnt = -1 # easier than assigning to 0 and finding end of
        # next loop
        all_patch_flag = 1 # to make sure all patches in all the
        # chains for a given attempt is taken care of. If not move to
        # the next attempt. DEFAULT to TRUE.

        while chcnt < nch-1: #since chcnt is initialized to -1

            if all_patch_flag == -1:
                flog.write('Could not find a sequence \n')
                break #move to next attempt

            chcnt += 1
            segname = ret_segname(segpref,chcnt+1)
            flist.write(';# chain number:\t%d\n' %(chcnt+1))
            flist.write(';# -- Begin patches for %s ---\n' %(segname))

            patcnt = 0
            branched = 0
            patname_L = 'None'
            deg_poly_chain = nresarr[chcnt]

            # aflag: for checking res-pat-res constraints
            # cflag: for checking pat1-pat2 adjancency
            # cflag2/cflag3/cflag4: for checking branchpat-backbonepat
            # Need to check both the monomers a patch connects
            # patch_n between res_n and res_n+1
            innerpatattempt = 0
            while patcnt <= deg_poly_chain-2: #for checking constraints

                innerpatattempt += 1 #update until 100 attempts/res
                if innerpatattempt > deg_poly_chain*100: 
                    all_patch_flag = -1 # uncheck flag 
                    break

                resname1 = residlist[chcnt][patcnt]
                resname2 = residlist[chcnt][patcnt+1]
                
                # Normal case: resname1 and resname2 are "normal" RES
                if (resname1 not in branch_opt) and (resname2 not in\
                   branch_opt):
                    patchname,aflag,cflag = write_normal_patch(cumulprobarr,\
                                                               inp_dict,\
                                                               resname1,\
                                                               resname2,\
                                                               fpresflag,
                                                               fppflag,\
                                                               patcnt,\
                                                               pres_fyle,\
                                                               patforbid,\
                                                               branch_opt,\
                                                               out_list,\
                                                               chcnt,\
                                                               patname_L)

                    if patchname == 'ERR':
                        return -1 

                    # Some of the following conditions maybe redundant
                    # but it is OK.
                    # Two extra conditions to identify: 1) if resname1
                    # is connected to a branch, then the branch-resname1
                    # patch should be compatabile with
                    # resname1-resname2 patch. 2) if resname2 is
                    # connected to a branch, then the branch-resname2
                    # patch should be compatible with
                    # resname1-resname2 patch. Feed patchname from
                    # first condition.
                    #Condition 1 and 2 cannot be simultaneously true
                    cflag4 = 0
                    if patcnt > 0: # Condition 1
                        if (residlist[chcnt][patcnt-1] in branch_opt):
                            resname0  = residlist[chcnt][patcnt-1]
                            resindex  = branch_opt.index(resname0)
                            gr_patname = branch_opt[resindex+1]
                            cflag4 = is_forbid_patch(gr_patname,\
                                                     patchname,patforbid)
                    if (patcnt+2) <= deg_poly_chain - 1: #Condition 2
                        if (residlist[chcnt][patcnt+2] in branch_opt):
                            resname3  = residlist[chcnt][patcnt+2]
                            resindex  = branch_opt.index(resname3)
                            gr_patname = branch_opt[resindex+1]
                            cflag4 = is_forbid_patch(patchname,\
                                                     gr_patname,patforbid)


                    # Update list if conditions are met
                    if aflag == 1 and cflag == 0 and cflag4 == 0:
                        out_list[chcnt].append(patchname)
                        flist.write(' patch\t%d\t%s\t%s:%d\t%s:%d\n' \
                                    %(patcnt+1,patchname,\
                                      segname,patcnt+1,segname,patcnt+2))
                        if patchname == '55' and resname1 == 'SYR':
                            print("ERR:3 respat constraint",\
                                  aflag, c1flag,cflag4) 
                        if patchname == '55' and resname2 == 'SYR':
                            print("ERR:4 respat constraint",\
                                  aflag, c1flag, cflag4)
                        patname_L = patchname # update left "normal" patch
                        patcnt += 1 # update counter


                    elif aflag == -2: #pres_fyle format is wrong
                        return 

                    #end update aflag/cflag

                    continue # continue to while loop
                        
                # Special Case 1: "left RES" of the patch is a branch
                # monomer. Branch patch between left side (res_n) and
                # right side (res_n+1). Patches are assigned to the
                # next residue. But rule 2 written above the function
                # defintion needs to be checked. Therefore don't
                # update patchname_L. Next "normal" patch will be
                # compared alongside patchname_L. Check for cflag
                # between patch m and m-1 to make sure that certain
                # branch-backbone patches are not allowed (for eg: GOG-BB)
                elif resname1 in branch_opt:
                    resindex   = branch_opt.index(resname1)
                    patchname  = branch_opt[resindex+1]
                    flist.write(' patch\t%d\t%s\t%s:%d\t%s:%d\n' \
                                %(patcnt+1,patchname,\
                                  segname,patcnt+1,segname,patcnt+2))
                    out_list[chcnt].append(patchname)
                    patcnt += 1
                    continue #continue to next residue                    


                # Special Case 2: "right RES" of the patch is a branch
                # monomer. 
                elif resname2 in branch_opt:
                    if patcnt > 0:
                        resname0 = residlist[chcnt][patcnt-1]
                    else:
                        resname0 = 'None'
                    # Case 2a: last RES is branch. Patch branch between
                    # n and n+1. branch_at_n is irrelevant here,
                    # because resname1 and resname2 cannot be
                    # simultaneously branches. Again dont update
                    # patchname_L. It is irrelevant
                    # UPDATE: This condition cannot happen for the
                    # current system.
                    if patcnt == deg_poly_chain-2: 
                        resindex  = branch_opt.index(resname2)
                        patchname = branch_opt[resindex+1]
                        flist.write(' patch\t%d\t%s\t%s:%d\t%s:%d\n' \
                                    %(patcnt+1,patchname,\
                                      segname,patcnt+1,segname,patcnt+2))
                        out_list[chcnt].append(patchname)
                        patcnt += 1
                        continue # continue to while loop/next chain

                    #Case 2b: (a) patch normal between n and n+2. But
                    #check patch constraints with the previous
                    #"normal" patch. (b) Also check branch-backbone
                    #patch incompatibility like Josh suggested - for
                    #instance BB-GOG.
                    #Ex ref: S1-P1-G2-P2-F3-P3-G4-P4-F5-P5-G6...;
                    #Consider G4-F5 as reference. P3 and P5 are fixed
                    #because of branches. patname_L will take care of
                    #case (a) whether it has to be P2 or P3 depending
                    #upon residue on the left
                    # (b1) check is_forbid_pat(P4,P5) or (Pn,Pn+1)
                    # if res(n-1) == normal 
                      #(a)  Check is_forbid_pat(P4,P3)  
                    # else res(n-1) == branch
                      #(a) Check is_forbid_pat(P4,P2) 
                      #(b2) Check is_forbid_pat(P4,P3)  or (Pn,Pn-1)
                    else: 
                        resname3 = residlist[chcnt][patcnt+2]
                        # get patchname that satisfies condition (a)
                        patchname,aflag,cflag = write_normal_patch(cumulprobarr,\
                                                                   inp_dict,\
                                                                   resname1,\
                                                                   resname3,\
                                                                   fpresflag,
                                                                   fppflag,\
                                                                   patcnt,\
                                                                   pres_fyle,\
                                                                   patforbid,\
                                                                   branch_opt,\
                                                                   out_list,\
                                                                   chcnt,\
                                                                   patname_L)


                        if patchname == 'ERR':
                            return -1 

                        # feed patchname from condition (a) as left
                        # patch to condition (b1)/(b2)
                        if aflag == 1 and cflag == 0:
                            #check (Pn,Pn+1)
                            resindex  = branch_opt.index(resname2)
                            gr_patname = branch_opt[resindex+1]
                            cflag2 = is_forbid_patch(patchname,\
                                                     gr_patname,patforbid)

                            #check (Pn,Pn-1) if n-1 is branch
                            cflag3 = 0
                            if resname0 in branch_opt:
                                resindex = branch_opt.index(resname0)
                                gr_patname = branch_opt[resindex+1]
                                cflag3 = is_forbid_patch(patchname,\
                                                         gr_patname,patforbid)

                            if cflag2 == 0 and cflag3 == 0: 
                                out_list[chcnt].append(patchname)
                                flist.write(' patch\t%d\t%s\t%s:%d\t%s:%d\n' \
                                            %(patcnt+1,patchname,\
                                              segname,patcnt+1,\
                                              segname,patcnt+3))
                                patname_L = patchname # update "normal" patch
                                patcnt += 1 # update counter


                else: # Unknown condition
                    print('ERROR in sequence')
                    print('ch#/pat#/res1/res2',chcnt+1,patcnt+1,\
                          resname1,resname2)
                    return -1

            # end while loop
            flist.write(';# --End patch list for %d--\n' %(chcnt+1))

        # end for chcnt in range(nch)

        # Sum and update patch list IFF all patches are present
        if all_patch_flag == 1:

            # After going through all the chains, count occurence of 
            #each res/patch
            outdist = []
            for key in inp_dict:
                outdist.append(sum([i.count(key) for i in out_list]))

            #normalize
            sumval = sum(outdist)
            if sumval != sum_of_pat:
                print('Sum from distn,sum_of_pat array:'\
                      ,sumval,sum_of_pat)
                exit('ERROR: Sum not equal to the total # of patches')
            normlist = [x/sumval for x in outdist]

            #extract target probabilities and compare
            targ_probs = list(inp_dict.values())
            normval = numpy.linalg.norm(numpy.array(normlist) \
                                        - numpy.array(targ_probs))

            # Write to log file
            for wout in range(len(outdist)):
                flog.write('%g\t' %(outdist[wout]/sumval))
            flog.write('%g\n' %(normval))

            if normval <= tol:
                #write to log file
                flog.write('Found optimal patch configuration\n')
                print('Found optimal patch configuration..')
                flag_optimal = 1
                return out_list
                break

            elif normval < normold:
                if oneconfigflag == -1:
                    oneconfigflag = 1
                    print('Found configuration with pat_err: ',normval)
                else:
                    print('Updating configuration with pat_err: ',normval)
                backup_pat_list = [] #create new_backup list
                backup_pat_list = out_list.copy()
                flist.write('\n')
                normold = normval

    if oneconfigflag == -1:
        print('Could not find a patch list with constraints')
        return -1
    if flag_optimal == -1:
        print('Did not find optimal patch configuration')
        print('Using best patch configuration with L2norm: ',normold)
        flog.write('Did not find patch optimal configuration\n')
        flog.write('Using last patch configuration with L2norm: %g\n'\
                   %(normold))
        return backup_pat_list

#---------------------------------------------------------------------

# Find patch for Case 1: when RES1 and RES2 are normal residues.
def write_normal_patch(cumulprobarr,pat_dict,resname1,resname2,\
                       fpresflag,fppflag,patincnt,presctrfyle,ppctrlist\
                       ,branch_opt,curpat_list,chcnt,patchname_L):

    ranval = random.random() #seed is current system time by default    
    findflag = 0
    arrcnt = 0

    while arrcnt <= len(cumulprobarr):
        
        #Only need to check the less than value because
        #the array is organized in increasing order.
        #Break the loop once the first point where the
        #condition is met.
        if ranval < cumulprobarr[arrcnt]:

            patchname = list(pat_dict.keys())[arrcnt]
            if patchname in branch_opt: 
                ranval = random.random() #generate new random number
                arrcnt = 0 #reset while loop
                continue # iterate until normal patch

            findflag = 1
            
            # Add constraint flags: default to TRUE
            #so that if constraints are not there, it will
            #be appended. consec flag has to be 0 for true
            appendflag = 1; consecpatflag = 0 
            if fpresflag != 0 or fppflag != 0:
                resname_L = resname1

                # end if patcnt == 0
                resname_R = resname2
                if fpresflag:
                    appendflag = check_constraints(presctrfyle,patchname,\
                                                   resname_L,resname_R)
#                    if resname1 == 'SYR' and patchname == '55':
#                        if appendflag != 0:
#                            print("appendflag", appendflag, patincnt,\
#                                  fpresflag, resname_L)
                if fppflag:
                    consecpatflag =is_forbid_patch(patchname_L,\
                                                   patchname,ppctrlist)
                # patchname cannot follow patchname_L

                # end if fpresflag == 1 or fppflag == 1

            break

        else: # if ranval > cumulprobarr[arrcnt]
            
            arrcnt += 1 # update array counter
                
        # end ranval < cumulprobarr[]

    # end while arrcnt in range(len(cumulprobarr))

    if findflag != 1:
        print('Random value/Probarr:', ranval,cumulprobarr)
        print('Error: Did not find a random residue\n')
        patchname = 'ERR'
    # end if find flag

    if patchname_L == 'B5' and patchname == '55' and \
       consecpatflag == 0:
        print("ERROR",fpresflag,fppflag)
    
    return patchname,appendflag,consecpatflag
#---------------------------------------------------------------------

# check forbidden consecutive patches
# THIS IS FOR RES1-PATCH1-RES2-PATCH2 combination
# Only patch1 and patch2 are important. rest is checked in
# check_constraints: patname1 - leftpatch, patname2 - rightpatch
def is_forbid_patch(patchname1,patchname2,patforbid):
    flag = 0 # default not forbidden
    for i in range(len(patforbid)):
        if patforbid[i][0] == patchname1:
            if any(patchname2 in st \
                   for st in patforbid[i][1:len(patforbid[i])]):
                flag = 1
        
    return flag
#---------------------------------------------------------------------

# Write residues/patches in one go -- OBSOLETE. 
# Added in write_multi_segments
def write_segments_onego(fin,nresarr,nch,chnum,segname,res_list,\
                         patch_list,branch_opt):

    fin.write(';# ------Begin main code -----\n')
    fin.write(';# Writing % segments' %(sum(nresarr)))
    fin.write(';# Writing output for %d' %(chnum))
    fin.write(' resetpsf \n')
    fin.write(' segment %s {\n' %(segname))
    
    #Residues
    for rescnt in range(nresarr):
        fin.write('  residue  %d  %s\n' \
                  %(rescnt+1,res_list[chnum-1][rescnt]))

    fin.write('}')        
    fin.write('\n')

    #Patches
    for patcnt in range(max(nresarr)-2):
            
        resname1 = res_list[chnum-1][patcnt]
        resname2 = res_list[chnum-1][patcnt+1]
        patchname = patch_list[chnum-1][patcnt]

        # Normal Case: (see create_patches)
        if resname1 not in branch_opt and resname2 not in branch_opt:
            fin.write('patch  %s  %s:%d  %s:%d\n' \
                      %(patchname,segname,patcnt+1,segname,patcnt+2))


        # Special Case 1: (see create_patches)
        elif resname1 in branch_opt:
            fin.write('patch  %s  %s:%d  %s:%d\n' \
                          %(patchname,segname,patcnt+1,segname,patcnt+2))

        # Special Case 2: (see create_patches)
        elif resname2 in branch_opt:

            # Case 2a: last RES is branch. Patch branch between
            # n and n+1
            if patcnt == max(nresarr)-2: 
                fin.write('patch  %s  %s:%d  %s:%d\n' \
                          %(patchname,segname,patcnt+1,segname,patcnt+2))

            #Case 2b: patch normal between n and n+2
            else: 
                fin.write('patch  %s  %s:%d  %s:%d\n' \
                          %(patchname,segname,patcnt+1,segname,patcnt+3))
                    
        else: # Error
            print('Unknow res/patch sequence')
            print('ch#/patch#' , chnum, patcnt)
            print(res_list)
            print(patch_list)

 
    fin.write('\n')
#---------------------------------------------------------------------

# Write residues/patches iteration by iteration
def write_multi_segments(fin,iter_num,nresthisiter,nch,chnum,\
                         segpref,res_list,patch_list,branch_opt,\
                         maxnummons):

    # Extra condition to account for the branch monomer happening at
    # the end of a PARTIAL segment. Since mth branch is attached to
    # n+1th residue (except when it is at the end of a FULL segment),
    # the n+1th residue has to be a normal residue. Since two branch
    # residues cannot be adjacent, it suffices to add n+1th residue to
    # that iteration.
    if nresthisiter != maxnummons:
        if res_list[chnum-1][nresthisiter-1] in branch_opt:
            nresthisiter += 1

    if iter_num == -1 or iter_num == 1:
        fin.write(';# Chain number: %d of %d chains\n' %(chnum,nch))
        fin.write(';# ----Begin main code -------\n')
        fin.write('\n')

    if iter_num != -1:
        fin.write(';# Iteration number: %d\n' %(iter_num))
        fin.write('set count %d' %(nresthisiter))
        fin.write('\n')


    segname = ret_segname(segpref,chnum)
    fin.write(' resetpsf \n')
    fin.write(' segment %s {\n' %(segname))

    #Residues -- indices should have -1 for first dimension  
    for rescnt in range(nresthisiter):
        fin.write('  residue  %d  %s\n' %(rescnt+1,\
                                          res_list[chnum-1][rescnt]))

    fin.write('}')        
    fin.write('\n')
    fin.write('\n')

    #Patches -- ch indices should have -1 for first dimension
    for patcnt in range(nresthisiter-1):
        resname1 = res_list[chnum-1][patcnt]
        resname2 = res_list[chnum-1][patcnt+1]
        patchname = patch_list[chnum-1][patcnt]

        # Normal Case: (see create_patches)
        if (resname1 not in branch_opt) and (resname2 not in branch_opt):
            # Add extra letter to patches for consistency with top file
            # Only for B5 residues. May want to write it as an extra
            # wrapper later for generic cases
            if patchname == 'B5' and resname2 == 'GUAI':
                patchname = 'B5G'
            elif patchname == 'B5' and resname2 == 'PHP':
                patchname = 'B5P'
            elif patchname == 'B5' and resname2 == 'CAT':
                patchname = 'B5C'
            fin.write('patch  %s  %s:%d  %s:%d\n' \
                      %(patchname,segname,patcnt+1,segname,patcnt+2))

        # Special Case 1: (see create_patches)
        elif resname1 in branch_opt:
            fin.write('patch  %s  %s:%d  %s:%d\n' \
                      %(patchname,segname,patcnt+1,segname,patcnt+2))
            
        # Special Case 2: (see create_patches)
        elif resname2 in branch_opt:

            # Case 2a: last RES is branch. Patch branch between
            # n and n+1
            if patcnt == nresthisiter-2: 
                fin.write('patch  %s  %s:%d  %s:%d\n' \
                          %(patchname,segname,patcnt+1,segname,patcnt+2))

            #Case 2b: patch normal between n and n+2
            else: 
                # Add extra letter to patches for consistency with top file
                # Only for B5 residues. May want to write it as an extra
                # wrapper later for generic cases
                resname3 = res_list[chnum-1][patcnt+2]
                if patchname == 'B5' and resname3 == 'GUAI':
                    patchname = 'B5G'
                elif patchname == 'B5' and resname3 == 'PHP':
                    patchname = 'B5P'
                elif patchname == 'B5' and resname3 == 'CAT':
                    patchname = 'B5C'

                fin.write('patch  %s  %s:%d  %s:%d\n' \
                          %(patchname,segname,patcnt+1,segname,patcnt+3))
                    
        else: # Error
            print('Unknow res/patch sequence')
            print('ch#/patch#' , chnum, patcnt)
            print(res_list)
            print(patch_list)

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

def initiate_packmol(fpin,inptype,chains,tolval):
    fpin.write('# PACKMOL melt input for %s\n' %(inptype))
    fpin.write('# Contains num chains: %d with tolerance of %g Ang\n'\
               %(chains, tolval))
    fpin.write('\n')
    fpin.write('tolerance %g\n' %(tolval))
    fpin.write('\n')

    fpin.write('# Input filetype\n')
    fpin.write('filetype pdb\n')
    fpin.write('\n')

    outname = 'melt_' + inptype + '_nch_' + str(chains) + '.pdb'
    fpin.write('# Output filename\n')
    fpin.write('output %s\n' %(outname))
    fpin.write('\n')
    fpin.write('# Adding chains\n')
#---------------------------------------------------------------------

# Make packmol input scripts
def make_packmol(fpin,structname,nrepeats,trans_list):
    fpin.write('structure %s\n' %(structname+'.pdb'))
    fpin.write('\t number %d\n' %(nrepeats))
    if trans_list != []:
        fpin.write('\t fixed')
        for k in range(6):
            fpin.write('\t %s' %(trans_list[k]))
        fpin.write('\n')
    fpin.write('end structure\n')
    fpin.write('\n')
#---------------------------------------------------------------------

# Make auxiliary files for NAMD/LigninBuilder/GROMACS
def make_auxiliary_files(tcldir,pref_pdbpsf,nch,topname,flbdflag,\
                         input_lbd):


    # bundle.tcl for generating all psf in one go
    fbund = open(tcldir + '/step1.tcl','w')
    fbund.write('# Combined file to generate psf files for all chains\n')
    fbund.write('# Use source step1.tcl from Tk console to run\n')
    
    # run_ligbuild.tcl to run ligninbuilder
    # flbd = open(tcldir + '/step2.tcl','w')
    # flbd.write('# Run LigninBuilder using this script\n')
    # flbd.write('# Requires psf files in the folder\n')
    # flbd.write('# Use source run_ligbuild.tcl from Tkconsole to run\n')
    # flbd.write('package require ligninbuilder\n')
    # flbd.write('::ligninbuilder::makelignincoordinates . . \n')
    # flbd.close()

    # combine_psf.tcl() to combine psf/pdb files and write GROMACS
    # generator if neeeded
    outname = pref_pdbpsf + '_nch_' + str(nch)
    inpname = pref_pdbpsf + '_chnum_' 
    topinp  = topname
    fcomb = open(tcldir + '/step2.tcl','w')
    fcomb.write('# To generate combined psf/pdb file..\n')    
    fcomb.write('# Use source step2.tcl from Tk console to run\n')
    fcomb.write('package require psfgen\n')
    if flbdflag == 1:
        fcomb.write('package require topotools\n')
    fcomb.write('%s %s %s\n' %('set','name',outname))
    fcomb.write('%s %s\n' %('topology',topinp))
    fcomb.write('resetpsf\n')
    fcomb.write('\n')
    fcomb.write('%s %s %s\n' %('for {set i 1} {$i <= ', str(nch), \
                             ' }  {incr i} {'))
    fcomb.write('%s %s\n' %('readpsf', inpname+str('$i.psf')))
    fcomb.write('%s %s\n' %('coordpdb', inpname+str('$i.pdb')))
    fcomb.write('}\n')
    fcomb.write('writepdb $name.pdb\n')
    fcomb.write('writepsf $name.psf\n')

    if flbdflag == 1:
        gmx_out = outname + '.top'
        fcomb.write('\n')
        fcomb.write("# Generate extraparameters.prm\nexec python3 ../findmissingterms.py\n")
        fcomb.write('# Generate GROMACS *.top file \n')
        fcomb.write('mol new $name.psf\n')
        fcomb.write('mol addfile $name.pdb\n')
        fcomb.write('%s %s %s %s %s %s\n' %('topo','writegmxtop'\
                                            ,gmx_out,'[list ',\
                                            input_lbd,' extraparameters.prm ]'))


    fcomb.write('exit\n')
    fcomb.close()
    
    if flbdflag == 1: #step3.tcl
        fminim = open(tcldir + '/step3.tcl','w')
        fminim.write('# Remove overlapping/pathological structures\n')
        fminim.write('# Use source step3.tcl from Tk console to run\n')
        fminim.write('package require ligninbuilder\n')
        fminim.write('::ligninbuilder::minimizestructure . namd2 +p8'\
                     + '   "parameters extraparameters.prm \n parameters ../par_all36_cgenff.prm \n" \n')
        fminim.write('exit\n')
        fminim.close()
    
    return fbund
#---------------------------------------------------------------------
