!-------------Main program for generating PDI froM SZ distribution----
!-------------Parameter file: pdigen.f90------------------------------
PROGRAM SZdist

  USE PDI_PARAMS
  IMPLICIT NONE

  READ_PDIINP_FILE()
  COMPUTE_INITQUANTS()
  GENERATE_MWVALS()
  WRITE_STATS()
  CLOSE_AND_DEALLOCATE()

END PROGRAM SZdist

!---------------------------------------------------------------------

SUBROUTINE READ_PDIINP_FILE()

  USE PDI_PARAMS
  IMPLICIT NONE

  INTEGER :: nargs, ierr, AllocateStatus,i, chtype_flag = 0

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Incorrect inputs for PDI generation.."

  CALL GETARG(nargs,pdi_fname)

  OPEN(unit = pdi_fid,file=trim(pdi_fname),action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) THEN

     PRINT *, trim(pdi_fname), "not found"
     STOP

  END IF

  ! Retrieve values of PDI and M from Readfile.dat, calculate
  ! corresponding inv_pdi values
  OPEN(unit = 10, file = "init_pdi.txt")

  DO

     READ(pdi_fid,*,iostat=ierr) dumchar
     
     IF(ierr .lt. 0) exit

     IF(trim(dumchar) == 'chain_types') THEN
        READ(pdi_fid,*) nch_types ! number of chain types

        chtype_flag = 1
     ELSEIF(trim(dumchar) == 'new_chain') THEN
        IF(chtype_flag == 0) THEN
           PRINT *, "ERR: Input number of chain types"
           STOP
        END IF

        READ(pdi_fid,*) chid, 
        READ(10,*)PDI1,M1,num_free ! FREE 
             PRINT *, "Unknown keyword", trim(dumchar)
     END IF

  END DO


END SUBROUTINE READ_PDIINP_FILE
!---------------------------------------------------------------------

SUBROUTINE COMPUTE_INITQUANTS()

  USE PDI_PARAMS
  IMPLICIT NONE

  tol = real(tol)/100.0 !convert tolerance to percentage
  step = range/maxsteps 


  allocate(MolWt1(1:num_free)) ! List of MW for N polymers
  allocate(MolWt2(1:num_graft)) ! List of MW for N polymers

END SUBROUTINE COMPUTE_INITQUANTS
!---------------------------------------------------------------------

SUBROUTINE COMPUTE_PDI()


END SUBROUTINE COMPUTE_PDI
!---------------------------------------------------------------------


SUBROUTINE GENERATE_MWVALS()

  USE PDI_PARAMS
  IMPLICIT NONE

  inv_pdi = 1.0/real(PDI1 - 1.0)
  
  ! Sigma (S) is on the range of 0 to infinity but will only count to
  ! "range" where sigma is assumed to be zero         
  do i=1,maxsteps
     S(i) = (range/real(maxsteps))*i
  end do

  if(PDI1 .GT. 1.0) then
     ! Calculate probability function
     do i=1,maxsteps
        P1(i) = (inv_pdi**inv_pdi)*(GAMMA(S(i))**-1)*(S(i)**(inv_pdi&
             &-1))*(EXP(-1*inv_pdi*S(i)))
     end do

     ! Calculate total area for normalization
     area1 = 0.5*(range/real(maxsteps))*P1(1)

     do i=2,maxsteps
        area1 = area1 + 0.5*(range/real(maxsteps))*(P1(i)+P1(i-1))
     end do

     ! Normalize the probability function
     do i=1,maxsteps
        Normal1(i) = P1(i)/area1 
     end do

     ! Calculate the normalized area for each slice such that the sum
     ! of normalized areas is one

     IntNormal1(1) = 0.5*Normal1(1)*(range/real(maxsteps)) ! definte
     ! IntNormal(1) so that it doesnt index out of bounds for first
     ! calculation
     
     do i=2,maxsteps
        IntNormal1(i)=0.5*(Normal1(i)+Normal1(i-1))*(range&
             &/real(maxsteps))
     end do

  end if
     
  ! Calculate random MW for each polymer using subtraction method

  call random_seed()

  ! ====================== Generate Chains  ============================

  ! Resetting variables
  Mi = 0
  Mi2 = 0
  itercnt = 0

  print *, "Generating free polymer list"

  if(PDI1 .GT. 1.0) then
     
     do while(loop == 1 .and. itercnt .le. maxiteration)
        
        itercnt = itercnt + 1
        Mi = 0
        Mi2 = 0

        ! ====== Generates polymer list using subtraction method =====
        do i=1,num_free
           call random_number(randnum)
           subtract = 1
           j = 1
           
           do while (subtract == 1) 
              randnum = randnum - IntNormal1(j)
              
              if (randnum .le.  0) then
                 subtract = 0
                 MolWt1(i) = int(S(j-1)*M1) 
              else      
                 j = j + 1
                 if(j == maxsteps+1) then
                    print *, 'array out of bounds error imminent'
                    loop = 0
                    exit  
                 endif
              endif
           enddo
        end do
        

        !==============================================================
     
        ! Calculate PDI of list
        do i=1,num_free
           Mi2 = Mi2 + (MolWt1(i)**2)
           Mi = Mi + Molwt1(i)
        end do
        
        PDIgen1 = real(Mi2*num_free)/real(Mi**2)
        
        ! Checks if generated PDI is within tolerance of desired PDI
        ! and  does  not have an Mi smaller than 2

        if (ABS(PDIgen1 - PDI1) .le. (PDI1*tol)) then
           if (minval(MolWt1) .ge. 2) then 
              loop = 0
           end if
        endif

     end do

  else

     do i = 1, num_free
        
        MolWt1(i) = M1
        
     end do

     ! Calculate PDI of list
     do i=1,num_free
        Mi2 = Mi2 + (MolWt1(i)**2)
        Mi = Mi + Molwt1(i)
     end do
     
     PDIgen1 = (Mi2*num_free)/(Mi**2)

  end if

  if(PDI1 .GT. 1 .and. loop == 1) then
     print *, "WARNING: Not converged before maximum iteration"
  end if

  open(unit = 20, file = "FreeChains.dat")
  
  write(20,'(3(I0,1X),F16.8)')num_free, Mi,Mi2,PDIgen1
  
  do i=1,num_free
     write(20,'(2(I0,4X))'),i, MolWt1(i)
  end do


END SUBROUTINE GENERATE_MWVALS
!---------------------------------------------------------------------

SUBROUTINE WRITE_STATS()

  !============  Calculates the average MW of the free chains to check
  ! if we got Mn back, reports the PDI and MW ========================
  do i=1,num_free
     average1 = average1 + MolWt1(i)
  end do

  average1 = real(average1)/real(num_free)
  print *, 'the number average molecular wt of the free chains is',&
       & average1 
  print *, 'the PDI of the generated free chain list is', PDIgen1

END SUBROUTINE WRITE_STATS
!---------------------------------------------------------------------

SUBROUTINE CLOSE_AND_DEALLOCATE()

  close(10)
  close(20)
  close(30)

END SUBROUTINE CLOSE_AND_DEALLOCATE
!---------------------------------------------------------------------


