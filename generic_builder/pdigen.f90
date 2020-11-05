!-------------Main program for generating PDI froM SZ distribution----
!-------------Parameter file: pdigen.f90------------------------------
PROGRAM GENERATE_SZDIST

  USE PDI_PARAMS
  IMPLICIT NONE
  INTEGER :: chid

  CALL READ_PDIINP_FILE()
  CALL INIT_LOGWRITE()
  CALL RANDOM_SEED()
  DO chid = 1,nch_types
     PRINT *, "Generating list for chain type: ", chain_id
     WRITE(log_fid,*) "Generating list for chain type ", chain_id
     CALL GENERATE_MWVALS(chid)
  END DO
  CALL WRITE_STATS()
  CALL CLOSE_FILES()

END PROGRAM GENERATE_SZDIST

!---------------------------------------------------------------------

SUBROUTINE READ_PDIINP_FILE()

  USE PDI_PARAMS
  IMPLICIT NONE

  INTEGER :: nargs, ierr, AllocateStatus,i, chtype_flag = 0
  INTEGER :: chcnt = 0. outflag = 0
  REAL::step ! Step size = range/maxsteps  
  CHARACTER(len=256) :: dumchar

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Incorrect inputs for PDI generation.."

  CALL GETARG(nargs,pdi_fname)

  OPEN(unit = pdi_fid,file=trim(pdi_fname),action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) THEN

     PRINT *, trim(pdi_fname), "not found"
     STOP

  END IF

  ! Read from file

  DO

     READ(pdi_fid,*,iostat=ierr) dumchar
     
     IF(ierr .lt. 0) exit

     IF(trim(adjustl(dumchar)) == 'chain_types') THEN
        READ(pdi_fid,*) nch_types ! number of chain types
        chtype_flag = 1
        ALLOCATE(PDI_arr(1:nch_types),stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "ERR: Did not allocate PDI_arr"
        ALLOCATE(avg_mw_arr(1:nch_types),stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "ERR: Did not allocate avg_mw_arr"
        ALLOCATE(nchain_list(1:nch_types),stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "ERR: Did not allocate nchain_list"
     ELSEIF(trim(adjustl(dumchar)) == 'chain_details') THEN
        IF(chtype_flag == 0) THEN
           PRINT *, "ERR: Input number of chain types"
           STOP
        END IF
        DO i = 1, nch_types
           READ(pdi_fid,*) PDI_arr(i), avg_mw_arr(i), nchain_list(i)
           chcnt = chcnt + 1
        END DO
        IF(chcnt .NE. nch_types) THEN
           PRINT *, "ERR: Not all chain types are present"
           STOP
        END IF
     ELSEIF(trim(adjustl(dumchar)) == 'tolerance') THEN
        READ(pdi_fid,*) tol        
     ELSEIF(trim(adjustl(dumchar)) == 'op_prefix') THEN
        READ(pdi_fid,*) out_fname
        outflag = 1
     ELSE
        PRINT *, trim(pdi_fname), trim(dumchar)
        STOP "ERR: Unknown keyword"
     END IF
  END DO
  CLOSE(pdi_fid)

  tol = real(tol)/100.0 !convert tolerance to fraction
  step = range/maxsteps 

  IF(outflag == 0) out_fname = "polydisp.inp"
  OPEN(unit = out_fid,file=trim(adjustl(out_fname)),action="write"&
       &,status="replace",iostat=ierr)


END SUBROUTINE READ_PDIINP_FILE
!---------------------------------------------------------------------

SUBROUTINE INIT_LOGWRITE()

  USE PDI_PARAMS
  IMPLICIT NONE
  INTEGER :: ierr, i
  OPEN(unit = log_fid,file='log_pdigen.txt',action="write",status="rep&
       &lace",iostat=ierr)
  
  WRITE(log_fid,*) "Log file for Schulz-Zimm distribution"
  WRITE(log_fid,*) "Number of chain types: ", nch_types
  WRITE(log_fid,*) "Tolerance (0,1): ", tol
  WRITE(log_fid,*) "Maximum attempts: " , maxiteration
  WRITE(log_fid,*) "************************************************"
  WRITE(log_fid,*) 

END SUBROUTINE INIT_LOGWRITE
!---------------------------------------------------------------------

SUBROUTINE GENERATE_MWVALS(chain_id)

  USE PDI_PARAMS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: chain_id
  INTEGER :: i,j,chcnt
  REAL::inv_pdi  ! Related to PDI [k=1/(PDI-1)] 
  REAL::randnum ! random number from 0 to 1
  REAL::PDIgen ! PDI of generated FREE polymer list
  REAL::area_val = 0 ! Area under the curve for chains
  INTEGER::Mi = 0 ! used to calculate PDI of list
  INTEGER::Mi2 = 0 ! used to calculate PDI of list
  INTEGER::itercnt = 0 !used to count till max iterations
  INTEGER::subtract = 1 ! Conditional: equal to 1 when the program is
  ! generating polymer, 0 when it has decided the M of the chain
  INTEGER::loop = 1 ! Similarly, will be equal to 1 until generated
  ! PDI is within the tolerance, then set to 0
  CHARACTER(LEN = 256) :: out_fname
  INTEGER:: init_index,nextmw ! for systems with nchains < 8

  inv_pdi = 1.0/real(PDI_arr(chain_id) - 1.0)
  
  ! Initialize all arrays to zero
  DO i = 1, maxsteps
     rat_arr(i) = 0
     Prob_SZ(i) = 0
     NormalDist(i) = 0
     IntNormalDist(i) = 0
  END DO

  ! Allocate MW list for each chain type
  ALLOCATE(MolWt_arr(1:nchain_list(chain_id)), stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "ERR: Did not allocate MolWt_arr"

  ! Sigma (S) is on the range of 0 to infinity but will only count to
  ! "range" where sigma is assumed to be zero         
  DO i=1,maxsteps
     rat_arr(i) = (range/REAL(maxsteps))*i
  END DO

  IF(PDI_arr(chain_id) .GT. 1.0) THEN
     ! Calculate probability function
     DO i=1,maxsteps
        Prob_SZ(i) = (inv_pdi**inv_pdi)*(GAMMA(rat_arr(i))**-1)&
             &*(rat_arr(i)**(inv_pdi-1))*(EXP(-1.0*inv_pdi*rat_arr(i)))
     END DO

     ! Calculate total area for normalization
     area_val = 0.5*(range/REAL(maxsteps))*Prob_SZ(1)

     DO i=2,maxsteps
        area_val = area_val + 0.5*(range/REAL(maxsteps))*(Prob_SZ(i)&
             &+Prob_SZ(i-1))
     END DO

     ! Normalize the probability function
     DO i=1,maxsteps
        NormalDist(i) = Prob_SZ(i)/area_val1 
     END DO

     ! Calculate the normalized area for each slice such that the sum
     ! of normalized areas is one

     IntNormalDist(1) = 0.5*NormalDist(1)*(range/REAL(maxsteps)) ! definte
     ! IntNormal(1) so that it doesnt index out of bounds for first
     ! calculation
     
     DO i=2,maxsteps
        IntNormalDist(i)=0.5*(NormalDist(i)+NormalDist(i-1))*(range&
             &/REAL(maxsteps))
     END DO

  END IF

  ! ====================== Generate Chains  ============================     
  ! Calculate random MW for each polymer using subtraction method
  ! Resetting variables

  Mi = 0
  Mi2 = 0
  itercnt = 0
  MolWt_arr = 0

  IF(ABS(PDI_arr(chain_id) - 1.0) .LT. 10**(-5)) THEN
     
     DO i = 1, nchain_list(chain_id)
        
        MolWt_arr(i) = avg_mw_arr
        
     END DO

     ! Calculate PDI of list
     DO i=1,nchain_list(chain_id)
        Mi2 = Mi2 + (MolWt_arr(chain_id)**2)
        Mi = Mi + Molwt_arr(chain_id)
     END DO
     
     PDIgen = (Mi2*num_free)/(Mi**2)

  ELSE IF(PDI_arr(chain_id) .GT. 1.0 .AND. chain_list(chain_id) .GE.&
       & 8) THEN
     
     DO WHILE(loop == 1 .AND. itercnt .LE. maxiteration)
        
        ! Reset variables
        itercnt = itercnt + 1
        Mi = 0
        Mi2 = 0
        MolWt_arr = 0

        ! ====== Generates polymer list using subtraction method =====
        DO i=1,nchain_list(chain_id)

           CALL RANDOM_NUMBER(randnum)
           subtract = 1
           j = 1
           
           DO WHILE (subtract == 1) 
              randnum = randnum - IntNormalDist(j)
              
              IF (randnum .LE.  0) THEN
                 subtract = 0
                 MolWt_arr(chain_id) = INT(rat_arr(j-1)*avg_mw_arr(chain_id) 
              ELSE  
                 j = j + 1
                 IF(j == maxsteps+1) then
                    PRINT *, 'ERR: Array out of bounds error imminent'
                    loop = 0
                    EXIT
                 END IF

              END IF

           END DO

        END DO
        
        ! Calculate PDI of list
        DO i=chcnt,nchain_list(chain_id)
           Mi2 = Mi2 + (MolWt_arr(chcnt)**2)
           Mi = Mi + Molwt_arr(chcnt)
        END DO
        
        PDIgen = REAL(Mi2*num_free)/real(Mi**2)
        
        ! Checks if generated PDI is within tolerance of desired PDI
        ! and  does  not have an Mi smaller than 2

        IF (ABS(PDIgen - PDI_arr(chain_id)) .LE. (PDI_arr(chain_id)&
             &*tol)) THEN
           IF(MINVAL(MolWt1) .ge. 2) THEN ! not a solvent
              loop = 0
           END IF
        END if

     END DO

     IF(loop == 1) then
        PRINT*, "WARNING: Not converged before maximum iteration"
     END IF

  ELSEIF(nchain_list(chain_id) .LT. 8) THEN
     PRINT *, "WARNING: Number of chains < 8", nchain_list(chain_id)
     PRINT *, "WARNING: Values NOT based on Schulz-Zimm distribution"
     WRITE(log_fid,*) "WARNING: Number of chains < 8",&
          & nchain_list(chain_id)
     WRITE(log_fid,*) "WARNING: Values NOT based on Schulz-Zimm distri&
          &bution"
     
     init_index = 1
     nextmw = INT(avg_mw_arr(chain_id)/8.0)
     IF(nextmw < 1) nextmw = 1
     IF(MOD(nchain_list(chain_id),2) .NE. 0) THEN
        MolWt_arr(1) = avg_mw_arr(chain_id)
        init_index = init_index + 1
     END IF
     DO i = init_index,2,nchain_list(chain_id)
        MolWt_arr(i) = avg_mw_arr(chain_id) + nextmw
        MolWt_arr(i+1) = avg_mw_arr(chain_id) - nextmw
        nextmw = nextmw + INT(avg_mw_arr(chain_id)/8.0)
        IF(nextmw < 1) nextmw = 1
     END DO
     DO i = 1,nchain_list(chain_id)
        IF(MolWt_arr(i) < 2) MolWt_arr = 2
     END DO
     ! Calculate PDI of list
     DO i=chcnt,nchain_list(chain_id)
        Mi2 = Mi2 + (MolWt_arr(chcnt)**2)
        Mi = Mi + Molwt_arr(chcnt)
     END DO
     PDIgen = REAL(Mi2*num_free)/real(Mi**2)

  END IF

  CALL WRITE_STATS(chain_id,Mi,Mi2,PDIgen)

  WRITE(out_fid,*) "num_chains", nchain_list(chain_id)
  DO chcnt=1,nchain_list(chain_id)
     write(out_fid,'(I0)') MolWt_arr(chcnt)
  END DO

END SUBROUTINE GENERATE_MWVALS
!---------------------------------------------------------------------

SUBROUTINE WRITE_STATS(chval,m1,m2,pdival)

  USE PDI_PARAMS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: chval,m1,m2
  REAL, INTENT(IN) :: pdival
  REAL :: mn, mw

  mn = REAL(m1)/REAL(nchain_list(chval))
  mw = mn/REAL(nchain_list(chval))

  WRITE(log_fid,*) "Chain details"
  WRITE(log_fid,*) "-------------"
  WRITE(log_fid,*) "Type","Target PDI_val","Num of Chains","Inp Mn" &
       &,"Out Mn", "Out Mw", "Out PDI"
  WRITE(log_fid,*) i,PDI_arr(chval),,nchain_list(chval)&
       &,avg_mw_arr(chval), mn, mw, pdival
   
END SUBROUTINE WRITE_STATS
!---------------------------------------------------------------------

SUBROUTINE CLOSE_FILES()

  CLOSE(log_fid)
  CLOSE(out_fid)

END SUBROUTINE CLOSE_FILES
!---------------------------------------------------------------------
