# PolydisperseResidueInputGenerator (PRIG)

This code will generate polydisperse residues according to Schulz-Zimm
distribution. The generated `.psf` files can be used directly with
LigninBuilder to generate starting `.psf`, `.pdb` and `.prm` files for
running with `NAMD` or `GROMACS`. Although this code is primarily used
to generate input files for various types of biomass, it can be used
with any protein or polymer complex. 

## Input Requirements

This code works on a combination of `Python` and `FORTRAN90`
platforms. `Python3.0+` and `ifort` compilers are required.

## Installation

Dowloading and unzipping the directory can be done using:
```
git clone https://github.com/vaidyanathanms/files_lignin.git
tar xvzf files_lignin.tar.gz
```

The above commands should generate a folder files_lignin with three
sub-folders:
*generic_builder
*src_switchgrass
*src_poplar

`src_switchgrass` and `src_poplar` contains `psf` and `pdb` files
specific to switchgrass and poplar wood. We *recommend* using
`generic_builder` for generating a structure from scratch. Inside
`generic_builder` folder, you should see the following python and
FORTRAN files:
-make_genpsf.py
-genconf.py
-pdi_gen.f90
-pdi_dist_params.f90

If the files are present, you are set to generate a new structure.
There are two ways to generate initial structure. You can either copy
the four files above to a new folder or run from the directory
`generic_builder`. 

## Making Inputs to PRIG

`inputsforpsfgen.inp` is a sample input file containing all the input
keywords to PRIG. We will look at the keywords in detail in the next
section (PRIG Keywords). For running PRIG, use:

```python
python genconf.py inputsforpsfgen.inp
```

`inputsforpsfgen.inp` can be substituted with the name of the input
file you are creating. If this generates a folder with `casenum_1` and
sub-folder `all_tclfiles` (within `casenum_1`) which contain a number
of `tcl` files, you are all set. 

Several things can go wrong including compiler compatibilities and
incompatible input constraints (see PRIG Keywords). If you find an
error, please report to [Vaidyanathan M. Sethuraman](v0e@ornl.gov).

## PRIG Keywords

In this section, we look at the different keywords that are needed to
generate a polydisperse input structure. 

#### Rules for input keywords
Some keywords are optional and are prefixed with (*Optional*) while
introducing the keyword. A space/tab should be present between the
keywords and arguments or between arguments. A line **cannot** be left
empty. A new line can start with an optional `#`. These lines will be
ignored. However, `#` *cannot* be used in the middle of a line.

#### Keyword list

1. case_num (*Optional*)

   All input files can start with an *optional* `case_num` keyword. If
   this is used as a keyword, it should be the **first** keyword in the
   input file. Usage:

   ```case_num caseID```

   `caseID` should be a positive integer. This will create a folder of the
   name `casenum_caseID` where all the output files will be
   present. Default value for `caseID` is 1.


   Example:

   ```case_num 1```


2. biomass_type

   This is a mandatory keyword and corresponds to the prefix for output
   file. Usage: 

   ```biomass_type argname```

   The final tcl files generated will be of the form argname_1_nch.tcl`,
   where `nch` corresponds to number of chains in the system.


   Example:

   ```biomass_type switchgrass```


3. num_resids

   This is a mandatory keyword and corresponds to the average number of
   residues per chain (segment). Usage:

   ```num_resids nres```

   where `nres` corresponds to the average number of residues per chain
   (segment). Should be an integer value.


   Example:

   ```num_resids 20```


4. num_chains

   Mandatory keyword corresponding to the number of chains in the
   system. Usage:

   ```num_chains nch```

   where `nch` corresponds to the number of chains in the system (integer
   value). 

   Example:

   ```num_chains 10```

5. disperse (*Optional*)

This keyword dictates the polydispersity of the system. Chains will be
drawn from a Schulz-Zimm distribution. There are two options (and
suboptions) for this case. Usage:

```disperse maketype optarg1 optarg2```

`maketype` can be either `CREATE` or `READ`. `CREATE` corresponds to
generating a set of polydisperse chains from scratch. `READ`
corresponds to the molecular weights (degree of polymerization) of all
the chains from a file. Arguments for each option are elaboratoed
below.

  * `CREATE` 
    For this case a new file will be generated according to the
    polydispersity value and the number of chains/number of residues
    per chain. If the option is `CREATE`, it should come *after*
    specifying `num_chains` and `num_resids`. Usage for this option is
    as follows:

    ```disperse CREATE PDIval Outputfile tolerance```

    `PDIval` and `Outputfile` corresponds to the target polydispersity
    value and the output file containing the molecular weights (degree
    of polymerization) of each chain. `PDIval` is the target dispersity
    index and is defined as the ratio between the weight average
    molecular weight and number average molecular weight. This number
    **should* be greater than 1.0. `Outputfile` is the name of the file
    that is generated where the degree of polymerization of each chain
    is written out.

    User can also specify an optional tolerance value (0-100). This
    corresponds to the maximum relative error (in %) between the target
    PDI  value and the simulated PDI. Different combinations will be
    tried to obtain either the target PDI value of the system exits
    after 50000 trials. Default value is 10. For all practical purposes
    values between 5 and 15 yield good output distribution if the number
    of chains in the system is less than 20. This argument is optional.
  
    Examples:

     ```
     disperse CREATE 1.50 polydisp.inp
     disperse CREATE 1.50 polydisp.inp 8.0
     ```

    If this option is used, after running the program, a file with the
    name 'geninp_pdistruct.txt' will be generated and it will contain
    the details of the inputs.

  * `READ`
    Users can also specify a file where the degree of polymerization of
    each chain is specified. In this case, the program will directly
    read this file and create the segments. Usage:

    ```disperse READ  inpfilename```

    where `inpfilename` is the name of the file. The `inpfilename`
    should have the following structure. First line **should* have the
    structure: 

    ```num_chains nchains```

    where nchains correspond to the number of chains in the system. This
    **should** be consistent with the `num_chains` in the input file
    used to run `genconf.py`. The next `n` lines should correspond to
    the degree of polymerization of the `n` different chains.

* top_ipfile






