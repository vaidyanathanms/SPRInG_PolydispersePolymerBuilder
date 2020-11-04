!------------Module containing parameters for PDI generation---------
!------------Use in conjunction with pdigen.f90----------------------

MODULE PDI_PARAMS

  IMPLICIT NONE

  !--------Input parameters------------------------------------------
  INTEGER,PARAMETER::maxsteps = 100000 ! Integration steps
  INTEGER,PARAMETER::maxiteration = 10000 ! Iteration counts
  REAL,PARAMETER::range = 5 ! Max value of Sigma
  REAL :: tol = 5 ! Tolerance value for PDI(%). Value between 1-100
  INTEGER :: nch_types = 1 ! Number of chain types

  !---------PDI variables-------------------------------------------
  REAL::inv_pdi  ! Related to PDI [k=1/(PDI-1)] 
  REAL::area = 0 ! Area under the curve for chains
  REAL::step ! Step size = range/maxsteps  
  REAL::randnum ! random number from 0 to 1
  REAL::average = 0 ! Sum of all free chain MW, to be divided by num
  REAL::PDIgen ! PDI of generated FREE polymer list
  INTEGER::Mi = 0 ! used to calculate PDI of list
  INTEGER::Mi2 = 0 ! used to calculate PDI of list
  INTEGER::itercnt = 0 !used to count till max iterations
  INTEGER::subtract = 1 ! Conditional: equal to 1 when the program is
  ! generating polymer, 0 when it has decided the M of the chain
  INTEGER::loop = 1 ! Similarly, will be equal to 1 until generated
  ! PDI is within the tolerance, then set to 0

  !---------File names------------------------------------------------
  CHARACTER(len=256) :: dumchar
  CHARACTER(len=256) :: pdi_fname
  INTEGER, PARAMETER :: pdi_fid = 300

  !---------------------Arrays----------------------------------------
  ! Chain list    
  REAL,DIMENSION(1:maxsteps)::P1 ! Probability of SZ distribution 
  REAL,DIMENSION(1:maxsteps)::Normal1 ! Normalized distributation 
  REAL,DIMENSION(1:maxsteps)::IntNormal1 !Integrated area for normdist
  REAL,DIMENSION(1:maxsteps)::S ! Ratios of chainlength to meanlength
  REAL,DIMENSION,ALLOCATABLE(:) :: avg_mw_arr ! Avg MW of chain types
  REAL,DIMENSION,ALLOCATABLE(:) :: MolWt_arr ! Molecular wt o/p list
  REAL,DIMENSION,ALLOCATABLE(:) :: PDI_arr ! PDI of each type

END MODULE PDI_PARAMS
