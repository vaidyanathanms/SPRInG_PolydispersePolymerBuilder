# Simple Polydisperse Residue Input Generator (SPRInG)

This code will generate polydisperse residues according to Schulz-Zimm
distribution. The generated `.psf` files can be used directly with
[LigninBuilder](https://github.com/jvermaas/LigninBuilder/tree/master/LigninBuilderPlugin)
to generate starting `.psf`, `.pdb` and `.prm` files for 
running with [NAMD](https://www.ks.uiuc.edu/Research/namd/) or
[GROMACS](http://www.gromacs.org/) software. 

Although this code was primarily designed to generate input files for
various types of biomass, it can be used with any protein or polymer
complex. 

## Input Requirements

This code works on a combination of `Python` and `FORTRAN90`
platforms. `Python3.0+` and `ifort`/`gfortran` compilers are required.

## Installation

Dowloading and unzipping the directory can be done using:
```
git clone https://github.com/vaidyanathanms/files_lignin.git
tar xvzf files_lignin.tar.gz
```

The above commands should generate a folder files_lignin with three
sub-folders:

* generic_builder
* src_switchgrass
* src_poplar

`src_switchgrass` and `src_poplar` contains `psf` and `pdb` files
specific to switchgrass and poplar wood. We *recommend* using
`generic_builder` for generating a structure from scratch. Inside
`generic_builder` folder, you should see the following python and
FORTRAN files:

- make_genpsf.py
- genconf.py
- pdi_gen.f90
- pdi_dist_params.f90

If the files are present, you are set to generate a new structure.
There are two ways to generate initial structure. You can either copy
the four files above to a new folder or run from the directory
`generic_builder`. 

## Making Inputs to SPRInG

`inputsforpsfgen.inp` is a sample input file containing all the input
keywords to SPRInG. We will look at the keywords in detail in the *SPRInG
Keywords* section. For running SPRInG, use:

```python
python genconf.py <filename>
```

where `<filename>` is the name of the input file to `genconf.py`.
If this generates a folder with `casenum_ID` and sub-folder
`all_tclfiles` (within `casenum_ID`) which contain a number
of `tcl` files, you are all set. Here, *ID* refers to an integer value
given as input (see SPRInG Keywords).

Several things can go wrong including compiler compatibilities and
incompatible input constraints (see SPRInG Keywords). If you find an
error, please report to [Vaidyanathan M. Sethuraman](v0e@ornl.gov).

## Combining SPRInG with LigninBuilder to Generate Structures

   Outputs from SPRInG can be directly fed into LigninBuilder to
   generate the input structure for `GROMACS` using the following
   three steps. Make sure to follow the order.

   1.  *Step 1*: If SPRInG ran correctly, users should see a folder
       `casenum_ID`, where `ID` is an integer value given as input to
       SPRInG (see SPRInG Keywords). Navigate to this directory using

       ```
       cd <casenum_ID>
       ```

       Within this directory, the inputs given to SPRInG, topology
       inputs and a folder `all_tclfiles` should be present. Navigate
       to `all_tclfiles` using

       ```
       cd all_tclfiles
       ```

       Inside the folder users should see several tcl files *viz.,*

       - bundle.tcl
       - run_ligbuild.tcl
       - combine_all.tcl
       - inpfile_chnum_ID.tcl

       where `inpfile_nch_ID.tcl` are a set of *N* tcl files with *N*
       corresponding to the number of chains (`ID` goes from 1 to *N*)
       in the system and *inpfile* corresponds to the name of the
       input system (see SPRInG Keywords). If all of these files are
       present, Step 1 is complete. 

   2.  *Step 2*: Within `all_tclfiles` directory, execute the
       following:

       ```
       vmd -dispdev text -e bundle.tcl
       ```

       Make sure the path to `vmd` is added to `$BIN` or is given
       correctly. If the command runs smoothly, this should generate
       `psf` files for each chain structure. 

       Following the generation of `psf` files, issue

       ```
       vmd -dispdev text -e run_ligbuild.tcl
       ```

       This should generate `pdb` files corresponding to the `psf`
       files within the directory. This requires `LigninBuilder` to be
       added in `~\.vmdrc` (see `LigninBuilder` on how to do this). 
	
       If the `pdb` file(s) is (are) not generated, please see the
       input constraints. Most likely a particular residue (or patch)
       is incompatible. Please make sure that `LigninBuilder` is added
       to `~\.vmdrc` before executing this command.

   3.  *Step 3*: If all the `psf` and `pdb` files are present, issue

       ```
       vmd -dispdev text -e combine_all.tcl
       ```

       If this command runs smoothly, this should provide an
       output file  of the form `inpfile_nch_N.tcl` where *inpfile* is
       the name of the input system (see SPRInG Keywords) and *N*
       corresponds to the number of chains in the system. Further, it
       should generate an output `top` file of the form
       `inpfile_nch_N.top`. Make sure the `prm` file is in the
       *superfolder* (*casenum_ID*) for this command to run.

   
   Once the `top` file is generated, these can be used as inputs to
   `GROMACS`. It is the user's responsibility to check whether the
   parameters match the atomnames (atomtypes) in the `psf`/`pdb` files.
   



## SPRInG Keywords

In this section, we look at the different keywords that are needed to
generate a polydisperse input structure. 

#### Rules for input file

-  Some keywords are optional and are prefixed with (*Optional*) while
   introducing the keyword. 
-  A space/tab should be present between the keywords and arguments or
   between arguments. 
-  All keywords are case-sensitive. 
-  A line **cannot** be left empty (blank). 
-  A new line can start with an optional `#`. These lines will be
   ignored. However, `#` *cannot* be used in the middle of a line.
-  The input file name (and file names used as arguments) cannot be
   specified as `None` or `none`. These are reserved keywords within
   the program.

#### Keyword list

1.  case_num (*Optional*)

    All input files can start with an *optional* `case_num` keyword. If
    this is used as a keyword, it should be the **first** keyword in the
    input file. Usage:

    ```
    case_num caseID
    ```

    `caseID` should be a positive integer. This will create a folder of
    the name `casenum_caseID` where all the output files will be
    present. Default value for `caseID` is 1.

    Example:

    ```
    case_num 1
    ```


1.  biomass_type

    This is a mandatory keyword and corresponds to the prefix for output
    file. Usage: 

    ```
    biomass_type argname
    ```

    The final tcl files generated will be of the form argname_1_nch.tcl`,
    where `nch` corresponds to number of chains in the system.
    

    Example:

    ```
    biomass_type switchgrass
    ```


1.  num_resids

    This is a mandatory keyword and corresponds to the average number of
    residues per chain (segment). Usage:

    ```
    num_resids nres
    ```

    where `nres` corresponds to the average number of residues per chain
    (segment). Should be an integer value.


    Example:

    ```
    num_resids 20
    ```


1.  num_chains

    Mandatory keyword corresponding to the number of chains in the
    system. Usage:

    ```
    num_chains nch
    ```

    where `nch` corresponds to the number of chains in the system
    (integer value). 

    Example:

    ```
    num_chains 10
    ```

1.  disperse (*Optional*)

    This keyword dictates the polydispersity of the system. If this
    option is not provided, chains are assumed as monodisperse by
    default. If the option is provided, chains will be drawn from a
    Schulz-Zimm distribution. There are two options (and suboptions)
    for this case. Usage:

    ```
    disperse maketype optarg1 optarg2
    ```

    `maketype` can be either `CREATE` or `READ`. `CREATE` corresponds
    to generating a set of polydisperse chains from scratch. `READ`
    corresponds to the molecular weights (degree of polymerization) of
    all the chains from a file. Arguments for each option are
    elaborated below.

    * `CREATE` 
      For this case a new file will be generated according to the
      polydispersity value and the number of chains/number of residues
      per chain. If the option is `CREATE`, it should come *after*
      specifying `num_chains` and `num_resids`. Usage for this option
      is as follows:

      ```
      disperse CREATE PDIval Outputfile ntrials tolerance
      ```

      `PDIval` and `Outputfile` corresponds to the target
      polydispersity value and the output file containing the
      molecular weights (degree of polymerization) of each
      chain. `PDIval` is the target dispersity index and is defined as
      the ratio between the weight average molecular weight and number
      average molecular weight. This number **should* be greater than
      1.0. `Outputfile` is the name of the file that is generated
      where the degree of polymerization of each chain is written out.
      `ntrials` correspond to the number of trials the program attempts
      to generate the polymer chains within the PDI tolerance limit.
      Values between 1000-10000 should be enough for most cases.

      User can also specify an optional tolerance value (0-100). This
      corresponds to the maximum relative error (in %) between the
      target PDI  value and the simulated PDI. Different combinations
      will be tried to obtain either the target PDI value of the
      system exits after 50000 trials. Default value is 10. For all
      practical purposes values between 5 and 15 yield good output
      distribution if the number of chains in the system is less than
      20. This argument is optional.
  
      Examples:

      ```
      disperse CREATE 1.50 polydisp.inp 10000
      disperse CREATE 1.50 polydisp.inp 1000 8.0
      ```

      If this option is used, after running the program, a file with the
      name 'geninp_pdistruct.txt' will be generated and it will contain
      the details of the inputs.

    * `READ`
      Users can also specify a file where the degree of polymerization
      of each chain is specified. In this case, the program will
      directly read this file and create the segments. Usage:

      ```
      disperse READ  inpfilename
      ```

      where `inpfilename` is the name of the file. The `inpfilename`
      should have the following structure. First line **should* have
      the following structure: 

      ```
      num_chains nchains
      ```

      where nchains correspond to the number of chains in the
      system. This **should** be consistent with the `num_chains` in
      the input file used to run `genconf.py`. The next `n` lines
      should correspond to the degree of polymerization of the `n`
      different chains.

1.  top_ipfile

    Mandatory keyword and the argument corresponds to the path to the
    topology file. It is the user's responsibility to check whether the
    residues generated have their monomer structure in the topology
    file. User should also provide the full path to the topology
    file. Default assumption is that the file is present in the path
    from which `genconf.py` is called. Usage:

    ```
    top_ipfile filename
    ```

    Example:
    ```
    top_ipfile	top_lignin.top
    ```

1.  resid_inp

    Mandatory keyword and the argument corresponds to the average
    probability of each residue in the system. It should be provided in
    a file with each line corresponding to the residue name and the
    average probability. Users **should** make sure that the residue
    name matches with the residue name in the topology file. Usage:

    ```
    resid_inp filename
    ```

    Example for formatting filename:

    ```
    SYR	0.4
    TRCN	0.05
    GUAI	0.3
    PCA	0.15
    FERUT 0.1
    ```

    **NOTE**: The sum of the probabilities need not be one. Code
     internally makes the sum to be one. However, a warning will be
     issued if the sum is not one. The inputs should contain the
     details for the branch (graft) monomers or else the code will not
     recongnize any branch monomer.

1.  patch_inp

    Mandatory keyword and the argument corresponds to the average
    probability of each patch in the system. It should be provided in
    a file with each line corresponding to the residue name and the
    average probability. Users **should** make sure that the patch
    name matches with the patch name in the topology file. Usage:

    ```
    resid_inp filename
    ```

    Example for formatting filename:

    ```
    BO4 0.8
    B5 0.1
    BB 0.05
    AO4 0.05
    ```
   
    **NOTES**:

	1.  The sum of the probabilities need not be one. Code
	    internally makes the sum to be one. However, a warning
	    will be issued if the sum is not one. 

	1.  The inputs should **NOT** contain the details for the
	    branch (graft) patches. DO NOT provide patch details for
	    branch monomers here. This is different from inputting
	    residues where the name of the branched residue should be
	    present. 

	1.  If you are using SPRInG with LigninBuilder, please be aware
	    that residues for which there exists equal probability
	    for the tacticities (e.g. BO4R and BO4L for BO4), DO NOT
	    give separate probabilities for each
	    stereoisomer. Provide the overall
	    probability. LigninBuilder will make sure that all the
	    stereoisomers have equal probability. 
   

1.  seg_name

    Mandatory keyword which corresponds to the name of the
    segment. Except for one case (see NOTE below), this will serve as
    the prefix for segment names for different chains in the `psf`
    file output. Usage:

    ```
    seg_name argname
    ```
   
    The output `psf` name will have segment names of the form
    `argname_chainID` where `chainID` is an integer varying from 1 to
    number of chains in the system. 

    **NOTE**: In case, an input PDB file is given to generate the PDB
    file using `genconf.py` and **NOT** LigninBuilder, users must
    make sure that the segment name matches the segment name in the
    input PDB file that is used to generate the initial guesses for
    the initial coordinates (ICs). 

1.  op_style (*Optional*)

    Keyword dictating the output style. There are two argument options
    -- `single` and `multi`. Usage:

    ```
    op_style single
    op_style multi 4
    ```

    For the argument `single`, a single ouput `tcl` file will be
    generated per chain. On running this with `VMD` or
    `LigninBuilder`, you can produce the final structure. However, the
    final structure may have unphysical bonds. It is the user's
    responsibility to check this. If the user is combining this code
    with `LigninBuilder` please use `single` option since
    `LigninBuilder` has capabilities to *untangle* unphysical bonds.

    For the argument `multi`, the program breaks down the `tcl` files
    into *smaller* `tcl` files. This requires an extra integer
    argument. Let us say that we are creating a polymer with degree of
    polymerization 20 and the extra argument is 4, the `tcl` file
    corresponding to this chain will have 5 different builds. First,
    the first four segments of the chain are built. Then NAMD is
    called to minimize the structure. The minimized structure is then
    used as an input to generate the next 4 residues -- so on and so
    forth. This requires `NAMD` path to be added correctly or else
    running `tcl` file may encounter errors.

    **NOTE**: Use the option `single` with `LigninBuilder`. Default is
      `single`.


1.  grafting (*Optional*)

    To define branching of main chain. Branches are single monomer
    long in the current mode. The program can manage multiple types of
    branches. Usage:

    ```
    grafting 1 branch1 patch1 branch2 patch2 ...
    ```

    The keyword `grafting` should be followed by an integer 1 or 0. 1
    corresponds to turning on the branch and 0 corresponds to no
    branch. This gives the user to toggle between branched and
    non-branched system easily. 

    The branch1/patch1 pair corresponds to the name of the branch
    residue and the patch connecting the branch with the backbone. The
    inputs should always be given as pairs. Users must make sure that
    the residue names are already present in the input list for
    residues. 

    Example:
    
    ```
    grafting 1 PCA GOG FERUT GOG
    ```

    **NOTE**: Since, by construction, the number of patches equal to
      the number of residues in the system, the final probabilities
      for the patch values may not reflect the input values.

1.  nattempts (*Optional*)

    Number of attempts to achieve a random configuration that
    corresponds to the input probabilities for residues and
    patches. Each time a better target configuration (smaller residual
    error) is found, the program saves that configuration. In case,
    the target residual error is not met within `nattempts`, the best
    configuration along with the residual error will be generated as
    output. Usage:

    ```
    nattempts intvalue
    ```

    where intvalue is the number of attempts. A value between 50 and
    200 for an average degree of polymerization of 30 works generates
    a target onfiguration in a few minutes. 

    Example:
    
    ```
    nattempts 60
    ```

    Default value for `nattempts` is 50.

1.  tol (*Optional*)

    Relative tolerance between the input probabilities for
    residues/patches and averaged output values for
    residues/patches. L2norm is used to calculate the relative
    error. Program runs until the error is less than tolerance value
    or the number of attempts exceeds `nattempts`. Usage:

    ```
    tol	tolval
    ```

    where `tolval` is a number between 0 and 1. Nominal values are
    between 0.05 and 0.15. Default value is 0.1.

1.  initiator (*Optional*)
    
    Use this option to make sure that certain fraction of the chains
    start (end) with this type of residue. This residue *should* be
    present in the input residue list. The final fraction of this
    residue will have an average probability closer to the value input
    to the program. Usage:

    ```
    initiator resname
    ```

    where `resname` is the name of the residue.

    Example

    ```
    initiator TRCN
    ```

1.  LigninBuilder (*Optional*)

    Use this option to generate output files that can be used in
    conjunction with `LigninBuilder` software to generate `*.top`
    files for `GROMACS` software. Usage:

    ```
    LigninBuilder filename
    ```

    where `filename` corresponds to the file that contains the details
    of the potentials. Usually this will have a `.prm`
    extension. Please make sure to provide the full path to the fil
    unless it is in same directory as `genconf.py`.

    Example:

    ```
    LigninBuilder par_lignin.prm
    ```

1.  clean_directories (*Optional*)

    Use this option to clean the existing output directory
    (`casenum_ID`) and replace with new files. Usage:

    ```
    clean_directories Y
    ```

    Arguments can be `Y` (yes) or `N` (no). Default is `N`.

    **WARNING**: All files will be deleted before the new output files
      are written into the directory.

1.  patch_patch_constraint (*Optional*)

    Optional argument to specify the constraints between adjacent
    patches. This is useful to let the program know that patches
    (linkers) cannot be next to each other. For instance, a patch
    (linker) of type $\beta$-5 cannot be followed by a 55 patch
    (linker) since the 5th position is already occupied. The
    constraints need to be specified in a separate file. Usage:

    ```
    patch_patch_constraint <filename>
    ```

    where `<filename>` contains the patch and the constraints. A
    sample example of this file will be as follows:

    ```
    B5L	55	5BR	5BL
    B5R	55	5BR	5BL
    55	55	B5L	B5R	5BR	5BL
    ```

    The first entry of each row should be the first of the two
    consecutive patches. Next 'n' entries of the row should contain
    all the patches that are incompatible with the first
    entry. 

    **NOTES**

    -  Each row is independent. In other words, if `pat_2` is
       incompatible with `pat_1` **does not** mean that `pat_1` is
       incompatible with `pat_2`. This incompatibility (if necessary)
       should be specified in separate rows.

    -  The patch names should *exactly* match the names of the patches
       in the input patch list given using `patch_inp`. If the patch
       name in the input does not match what is given in the list, it
       will be ignored. However, it is OK for this file to have the
       name of the patches that is not present in the input given
       using `patch_inp`. 

    -  Although this is an optional argument, almost always there will
       be restrictions for real sytems. It is the user's responsiblity
       to make sure that all the constraints are given to the system.
    
1.  patch_res_constraint

    Optional argument to specify the constraints between a patch
    (linker) and a residue. This is useful to let the program know
    that certain patches cannot succeed (or precede) a certain
    residue. For instance, a `55` patch can neither precede nor
    succeed a syringyl residue. However, `405` can precede a syringyl
    residue whereas it cannot succeed a syringyl residue. Similar to
    the `patch_patch_constraint` input, the incompatibility data
    should be input as a file. Usage

    ```
    patch_res_constraint <filename>
    ```

    where `<filename>` corresponds to the name of the file comprising
    the incompatibilities. There are certain keywords that *should* be
    present in this file. Keywords are case-sensitive.

    Example:

    ```
    patch	restrict_before	restrict_after
    55		SYR		SYR
    B5		None		SYR
    55		TRCN		TRCN
    B5		TRCN		TRCN
    ```

    **NOTES**

    -  The first line of this file **should** have the following
       keywords: `patch`, `restrict_before` and `restrict_after`.

    -  Each row should contain only 3 entries. The first entry is the
       name of the patch. The entry under `restrict_before` should
       correspond to the restriction that the patch cannot precede a
       residue and the entry below `restrict_after` should correspond
       to the restriction that the patch cannot succeed a residue. In
       the above example `55` can neither precede nor succeed a `SYR`
       residue. Similarly `B5` can neither precede nor succeed a
       `TRCN` residue. However, `B5` can come before `SYR` where as it
       cannot come after `SYR`.

    -  Use `None` keyword (case-sensitive), if a patch can precede
       (but not succeed) or viceversa a given residue. By default,
       there are no restrictions. In other words, for all the patches
       mentioned in the input will have `None` as default argument for
       incompatibility.

    -  If a patch cannot come before (or after) two residues, it
       should be specified in two separate lines. For instance, in the
       example above, patch `55` should be repeated twice for the
       program to know that it has incompatibility with both `SYR` and
       `TRCN`.

    -  Although this is an optional argument, almost always there will
       be restrictions for real sytems. It is the user's responsiblity
       to make sure that all the constraints are given to the system.


1.  gen_packmol

1.  namd_inp

1.  pdb_ipfile
    

    
