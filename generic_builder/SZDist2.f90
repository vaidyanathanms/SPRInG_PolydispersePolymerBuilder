program SZdist
  implicit none

  !--------Input parameters------------------------------------------
  integer,parameter::maxsteps = 100000 ! Number of steps 
  integer,parameter::maxiteration = 10000 ! sets the maximum amount
  ! of times the program will  try to get a PDI within the tolerance
  ! before exiting
  real:: tol = 5 ! tolerance value for PDI(%). Value between 1-100
  real,parameter::range = 5 ! Max value of Sigma
  
  !--------Other parameters------------------------------------------
  real::PDI1 ! PDI (polydispersive index) of FREE chains
  real::k1 ! Related to PDI [k=1/(PDI-1)]
  real::PDI2 ! PDI of GRAFT chains
  real::k2 ! Similar to k1
  
  real::area1 = 0 ! Area under the curve for FREE chains
  real::area2 = 0 ! same but for GRAFT chains
  real::step ! Step size = range/maxsteps  
  real::randnum ! random number from 0 to 1
  real::average1 = 0 ! Sum of all free chain MW, to be divided by num
  real::average2 = 0 ! Sum of all graft chain MW, to be divided by num
  
  real::PDIgen1 ! PDI of generated FREE polymer list
  real::PDIgen2 ! PDI of generated GRAFT polymer list
  integer::Mi = 0 ! used to calculate PDI of list
  integer::Mi2 = 0 ! used to calculate PDI of list


  integer::i,j,l
  integer::num_free,num_graft
  integer::M1 ! Number avg MW of FREE chains
  integer::M2 ! MW of GRAFT chains

  integer::itercnt = 0 !used to count till max iterations
  integer::subtract = 1 ! Conditional: equal to 1 when the program is
  ! generating polymer, 0 when it has decided the M of the chain
  integer::loop = 1 ! Similarly, will be equal to 1 until generated
  ! PDI is within the tolerance, then set to 0

  character (len=100) :: dumchar
  integer :: ierr
  real,dimension(1:maxsteps)::S ! Ratios of chainlength to meanlength

  ! FREE CHAINS    
  real,dimension(1:maxsteps)::P1 ! Probability of SZ distribution 
  real,dimension(1:maxsteps)::Normal1 ! Normalized distributation 
  real,dimension(1:maxsteps)::IntNormal1 !Integrated area for normdist
  integer,allocatable,dimension(:)::MolWt1 !List of MW for N polymers

  ! GRAFT CHAINS    
  real,dimension(1:maxsteps)::P2 ! Probability of SZ distribution 
  real,dimension(1:maxsteps)::Normal2 ! Normalized distributation 
  real,dimension(1:maxsteps)::IntNormal2 !Integrated Area for normdist
  integer,allocatable,dimension(:)::MolWt2 ! List of MW for N polymers


  !***********ANALYSIS BEGINS HERE************************************

  tol = real(tol)/100.0 !convert tolerance to percentage
  step = range/maxsteps 

  ! Retrieve values of PDI and M from Readfile.dat, calculate
  ! corresponding k values
  open(unit = 10, file = "init_pdi.txt")

  do

     read(10,*,iostat=ierr) dumchar
     
     if(ierr .lt. 0) exit

     if(trim(dumchar) == 'free_data') then
          read(10,*)PDI1,M1,num_free ! FREE 
       else if(trim(dumchar) == 'graft_data') then 
          read(10,*)PDI2,M2,num_graft ! GRAFT
       else
          print *, "unknown keyword", trim(dumchar)
       end if

    end do

  allocate(MolWt1(1:num_free)) ! List of MW for N polymers
  allocate(MolWt2(1:num_graft)) ! List of MW for N polymers

  print *, "PDI1", PDI1
  print *, "PDI2", PDI2

  k1 = 1.0/real(PDI1 - 1.0)
  k2 = 1.0/real(PDI2 - 1.0)


  ! Sigma (S) is on the range of 0 to infinity but will only count to
  ! "range" where sigma is assumed to be zero         
  do i=1,maxsteps
     S(i) = (range/real(maxsteps))*i
  end do

  if(PDI1 .GT. 1.0) then
     ! Calculate probability function
     do i=1,maxsteps
        P1(i) = (k1**k1)*(GAMMA(S(i))**-1)*(S(i)**(k1-1))*(EXP(-1*k1&
             &*S(i)))
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
     
  if(PDI2 .GT. 1.0) then
     do i=1,maxsteps
        P2(i) = (k2**k2)*(GAMMA(S(i))**-1)*(S(i)**(k2-1))*(EXP(-1*k2&
             &*S(i)))
     end do

     ! Calculate total area for normalization
     area2 = 0.5*(range/real(maxsteps))*P2(1)
     
     do i=2,maxsteps
        area2 = area2 + 0.5*(range/real(maxsteps))*(P2(i)+P2(i-1))
     end do

     ! Normalize the probability function
     do i=1,maxsteps
        Normal2(i) = P2(i)/area2 
     end do

     ! Calculate the normalized area for each slice such that the sum
     ! of normalized areas is one

     IntNormal2(1) = 0.5*Normal2(1)*(range/real(maxsteps))
     
     do i=2,maxsteps
        IntNormal2(i)=0.5*(Normal2(i)+Normal2(i-1))*(range&
             &/real(maxsteps))
     end do
     
  end if

  ! Calculate random MW for each polymer using subtraction method

  call random_seed()

  ! ====================== FREE POLYMER CHAINS ======================

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

  ! ====================== GRAFT POLYMER CHAINS ======================

  ! Resetting variables
  Mi = 0
  Mi2 = 0
  loop = 1
  itercnt = 0

  print *, "Generating graft polymer list"


  if(PDI2 .GT. 1.0) then
     
     do while(loop == 1 .and. itercnt .le. maxiteration) 

        Mi = 0
        Mi2 = 0
        itercnt = itercnt + 1

        ! ====== Generates polymer list using subtraction method =====

        do i=1,num_graft
           call random_number(randnum)
           subtract = 1
           j = 1
           
           do while (subtract == 1) 
              randnum = randnum - IntNormal2(j)
              
              if (randnum .le.  0) then
                 subtract = 0
                 MolWt2(i) = int(S(j-1)*M2) 
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
        ! ===========================================================
     
        ! Calculate PDI of list
        do i=1,num_graft
           Mi2 = Mi2 + (MolWt2(i)**2)
           Mi = Mi + Molwt2(i)
        end do
        
        PDIgen2 = real(Mi2*num_graft)/real(Mi**2)
        
        ! Checks if generated PDI is within tolerance of desired PDI
        ! and does  not have an Mi smaller than 4
        
        if (ABS(PDIgen2 - PDI2) .le. (PDI2*tol)) then
           if (minval(MolWt2) .ge. 4) then 
              loop = 0
           end if
        end if
        
     end do

  else
        
     do i = 1, num_graft
        
        MolWt2(i) = M2
        
     end do
     
     ! Calculate PDI of list
     do i=1,num_graft
        Mi2 = Mi2 + (MolWt2(i)**2)
        Mi = Mi + Molwt2(i)
     end do
     
     PDIgen2 = real(Mi2*num_graft)/real(Mi**2)
     
  end if

  if(PDI2 .GT. 1 .and. loop == 1) then
     print *, "WARNING: Not converged before maximum iteration"
  end if

  open(unit = 30, file = "GraftChains.dat")

  write(30,'(3(I0,1X),F16.8)'),num_graft,Mi,Mi2,PDIgen2

  do i=1,num_graft
     write(30,'(2(I0,1X))'),i, MolWt2(i)
  end do

  !============  Calculates the average MW of the free chains to check
  ! if we got Mn back, reports the PDI and MW ========================
  do i=1,num_free
     average1 = average1 + MolWt1(i)
  end do

  average1 = real(average1)/real(num_free)
  print *, 'the number average molecular wt of the free chains is',&
       & average1 
  print *, 'the PDI of the generated free chain list is', PDIgen1

  do i=1,num_graft
     average2 = average2 + MolWt2(i)
  end do

  average2 = average2/real(num_graft)
  print *, 'the number average molecular wt of the graft chains is',&
       & average2 
  print *, 'the PDI of the generated graft chain list is', PDIgen2

  close(10)
  close(20)
  close(30)
end program SZdist
