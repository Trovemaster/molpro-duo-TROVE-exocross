 module FITPES_ch3cl
 
  implicit none
  !
  public do_fitting,lapack_dgesvd,lapack_sdd_pseudo_inverse

  private

 integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
 integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
 integer, parameter :: ark          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
 integer,parameter          :: f_inp=5, f_out=6, f_res=9, f_new =10, f_par  = 11
!, nparam_eq = 5
!  TAKING THIS OUT FOR NOW TO TRY AND MAKE MORE GENERAL (BPM)
  integer(ik), parameter :: verbose = 4                     ! Verbosity level


  contains
  
 subroutine do_fitting
 !
 implicit none
 !
 ! Parameters:  
 ! 
 ! number of parameters and maximum number of ab initio points   
 !
 !
 integer,parameter          ::  enermax =  60000, nmax0=8,taumax0=8, nparam = 10000
 !
 ! where the output goes to 
 ! 
 integer,parameter          :: f_inp=5, f_out=6, f_res=9, f_new =10
 ! , nparam_eq = 5
 ! 
 ! Dinite difference differentiation with 2 or 3 points
 !
 integer,parameter          ::  findif = 3 
 !
 ! the iteration will be stoped when the best standard deviation reached 
 !
 !double precision,parameter :: stadev_best = 0.001d0   
 ! 
 ! parameter to control the correlation when we through out the parameters with large st. error (0.0..1.0)
 !
 double precision,parameter :: fact_corr=0.4, too_big = 1e7
 ! 
 ! Watson's alpha parameter to control the Robust fit
 double precision,parameter :: a_wats = 0.1
 ! parameter for the finite differncies differention 
 !
 double precision,parameter :: factordeltax=0.001
 !
 integer,parameter          ::  ireplace_by = 1
 !
 integer,parameter          ::  Ndim = 12
 !
 ! Universal constants 
 ! 
 double precision,parameter :: bohr = 0.52917715e+00,hartre=2.194746354d+05
 ! Avareging region for the robustfit = 3 
 integer,parameter          :: ishift = 10
  character(len=40)         :: fit_type= 'dgelss'
 ! Variables:  
 double precision  :: maxdeviat  ! st. deviation at the previos iteration 
 !
 integer           :: i,npts,nused,itmax,numpar,ncol,iter,nrow,l,irow,icolumn ! just integers  
 ! 
 character(len=80) :: title(4),longlabel ! input title of the job from the first four input lines 
 character(len=1)  :: mark           ! used at the parameters output to mark the dependent or independent ones 
 character(len=10) :: label          ! Temporary label 
 double precision  :: wtsum,ssq      ! 
 !
 ! some variables 
 !
 double precision  :: stadev,epsil,stadev_old,v,pi
 double precision  :: stability
 double precision  :: corr,ve
 integer           :: ierror,i1,i2,last_i,robustfit
 double precision  :: stab_best,conf_int
 integer           :: maxivar,minivar,tivar
 double precision  :: last_stadev,factor_back,factor_out,fit_scale=1.00d0,gamma_ps,enermax_ps,tol=0.1d0
 double precision  :: rms,ZPE,stadev_best,TEMPE,Ener_min
 integer           :: last_ivar,ndigits,alloc,ipar
 integer           :: i0,l0,parmax,fivar_t,nmax,taumax
 logical           :: yes,still_run,ifopen,find_best_params
 integer           :: ivartmp,imin,imax,ilarge_sterr,ilow_sterr,Nout,parmax_t
 double precision  :: conv_length,conv_energy,paramtmp,relat_sterr,relat_sterr_min,fparam_t,V0,redundancy
 double precision :: theta1,theta2,theta3,theta12,theta13,theta23,tau,Sa,Sb
 character(len=68)  :: fmt1,fmt2
 character(len=2)   :: fmt0
 character(len=80)  ::  char_tmp 
 integer           ::  NDEG,NBONDS,NANGS,nparam_eq,kcol
! MADE NPARAM_EQ A 
 DOUBLE PRECISION, ALLOCATABLE :: local_p(:)
 double precision  :: local(Ndim)
 ! NDEG ADDED BY BPM TO MAKE CODE MORE GENERAL, NOW READ IN NUMBER OF DEGREES
 ! AND ALLLOCATE MATRICES
 ! ALSO ADDED NBONDS
 !
 ! some matrices, to be allocated  
 !
 double precision, allocatable :: rjacob(:,:),last_param(:),param(:)
 double precision, allocatable :: param_tmp(:)
 double precision, allocatable :: crr(:),coord(:,:),coord_tmp(:,:)
 double precision, allocatable :: energy(:),energy_tmp(:)
 double precision, allocatable :: eps(:),wt(:),parold(:),sigma(:),wt_tmp(:),parinitial(:)
 double precision, allocatable :: al(:,:),bl(:),dx(:)
 double precision, allocatable :: ai(:,:),sterr(:),dV(:,:)
 double precision, allocatable :: sens(:),df(:)
 character(len=10), allocatable:: parnam(:),parnam_tmp(:)
 character(len=3)  :: af(12)
 integer, allocatable :: ivar(:),icorr(:),ivar0(:),ivar_tmp(:),ipower_tmp(:,:),ipower(:,:),icol(:),IPOWER_(:),ivar1(:)
 !
 double precision, allocatable :: Tsing(:,:) 
 integer                       :: rank0,info,alloc_p,alloc_p2
 integer                       :: lwork,NumSVD,k
 double precision,allocatable  :: wspace(:)
 !
 ! Common block is to reduce the number of assignments within the potential routine
 !
 ! Here we go!

 !
 !  initial constants 
 !
 last_stadev = 1000000.
 maxivar = 8
 maxdeviat = 90000.0
 tol = -100.0d0
 !
 pi = 4.0d0 * datan2(1.0d0,1.0d0)
 !
 ! Input data from file: 
 !
 call skiplines(f_inp,1)
 !
 !  input the job title 
 !
 do i=1,4
   read  (f_inp,"(a80)") title(i)
 enddo
 !
 !  output the job title 
 !
 write (f_out,"(3(132('*')/' '),16('*'),100x,16('*'))") 
 write (f_out,"(4(' ',16('*'),10x,a80,10x,16('*')/))") (title(i), i=1,4)
 write (f_out,"(' ',16('*'),100x,16('*')/' ',3(132('*')/' '))") 
 !
 call skiplines(f_inp,4)
 !
 read (f_inp,*) NDEG
 read (f_inp,*) NBONDS
 read (f_inp,*) NANGS
! ADDED IN BY BPM TO MAKE MORE GENERAL
 allocate(IPOWER_(NDEG)) 
 
 read (f_inp,*) itmax
 read (f_inp,*) stab_best,stadev_best
 read (f_inp,*) conv_length
 read (f_inp,*) conv_energy
 read (f_inp,*) factor_out
 read (f_inp,*) factor_back
 read (f_inp,*) robustfit
 read (f_inp,*) tol
 read (f_inp,*) fit_scale
 read (f_inp,*) gamma_ps     ! partridge and schwenke weight: gamma 
 read (f_inp,*) enermax_ps   ! partridge and schwenke weight: enermax must be positive to employed 
 !
 write(f_out,"(i10,' <= number of degrees of freedom')")         NDEG
 write(f_out,"(i10,' <= number of bonds')")         NBONDS
 write(f_out,"(i10,' <= number of bond angles')")         NANGS
 write (f_out,"(i10,' <= number of iterations in the fit')")         itmax
 write (f_out,"(f10.6,' <= The fit is over if the stability is reached ')") stab_best
 write (f_out,"(f10.6,' <= The fit is over if the st.dev. best is reached ')") stadev_best
 write (f_out,"(f10.6,' <= sterr parameter out ')")                  factor_out
 write (f_out,"(f10.6,' <= stdev parameter out ')")                  factor_back
 write (f_out,"(i10,  ' <= watson robust on (1) or of (0) ')")       robustfit
 !
 if (Ndim/=Ndeg) then
   write(6,"('Ndim/=Ndeg, update Ndim in the code to ',i8)") Ndeg
   stop 'Ndim/=Ndeg'
 endif
 !
 fit_type = "linur"
 if (tol>0) fit_type = "SDD"
 !
 if (factor_back<0) then 
   find_best_params = .false.
 else
   find_best_params = .true.
 endif 
 !
 factor_back = abs(factor_back)
 !
 ! Converiosn factors:
 ! bond length 
 !
 !if (flag_borh ==1) then 
 !  conv_lenght = bohr 
 !  write (f_out,"(' The bond lengths are in bohr')") 
 !else 
 !  conv_lenght = 1.0
 !  write (f_out,"(' The bond lengths are in Angstrom')") 
 !endif   
 !
 ! energies 
 !if (flag_hartre==1) then 
 !  conv_energy = hartre
 !  write (f_out,"(' The ab initio energies are in hartre')") 
 !else 
 !  conv_energy = 1.0 
 !  write (f_out,"(' The ab initio energies are in 1/cm')") 
 !endif   


 !
 ! Array's allocation:
 !

 parmax_t = nparam
 lwork =  parmax_t*8
 !
 allocate (param_tmp(parmax_t),parnam_tmp(parmax_t),ivar_tmp(parmax_t),ipower_tmp(1:NDEG,parmax_t),stat=alloc)
 if (alloc/=0) then 
  write(6,"('param_tmp - out of memory')")  
    stop 'param_tmp - out of memory'
 endif 

 ipower_tmp = 0


 !
 ! input  parameters of structur
 !
 parnam_tmp = 'xxxxxx'
 param_tmp = 0 
 ivar_tmp = 0 
 !
 call skiplines(f_inp,3)
 READ (f_inp,*) NPARAM_EQ
 WRITE(f_out,"(i10,' <= number of non-linear parameters')")  NPARAM_Eq
 do i = 1,nparam_eq
   read (f_inp,"(a80)") longlabel
   read (longlabel,*) label,ivartmp,paramtmp
     parnam_tmp(i) = label 
     ivar_tmp(i)   = ivartmp
     param_tmp(i)  = paramtmp
 enddo 

 write (f_out,"(' parameters of structure ')")
 do i = 1,nparam_eq
    write (f_out,"(a10,'=',f16.3)") parnam_tmp(i),param_tmp(i)
 enddo 
 !
 Yes = .True. 
 i = nparam_eq
 nmax = 0
 taumax=0
 read_loop: do while ( Yes )
    read (f_inp,"(a80)") char_tmp
    if ( char_tmp(1:3)/='end' ) then
       !
       i = i + 1
       !
       read (char_tmp(1:80),*) parnam_tmp(i),ipower_(1:NDEG),fivar_t,fparam_t
!      THIS PART COULD DO WITH A FORMAT STATEMENT BPM OR JUST DON'T PRINT
       if (verbose >=7) write (f_out,*) trim(parnam_tmp(i)),ipower_(1:NDEG),fparam_t


       param_tmp(i) = fparam_t
       ivar_tmp(i) = fivar_t
       ipower_tmp(1:NDEG,i) = ipower_(1:NDEG)
       !
    else
       Yes = .False.
    endif
 enddo read_loop
 !
 parmax = i 
 !
 allocate (param(parmax),parnam(parmax),ivar(parmax),icorr(parmax),ivar0(parmax),ipower(NDEG,parmax),ivar1(parmax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('param - out of memory')")  
    stop 'param - out of memory'
 endif 
 !
 param(1:parmax) = param_tmp(1:parmax)
 ivar(1:parmax) = ivar_tmp(1:parmax)
 ipower(1:NDEG,1:parmax) = ipower_tmp(1:NDEG,1:parmax)
 parnam(1:parmax) = parnam_tmp(1:parmax)
 ivar0(1:parmax) = ivar_tmp(1:parmax)
 last_i = 0
 !
 deallocate (param_tmp,parnam_tmp,ipower_tmp,ivar_tmp)
 !
 allocate (parold(parmax),parinitial(parmax),last_param(parmax),sterr(parmax),&
            crr(parmax),coord_tmp(NDEG,enermax),energy_tmp(enermax),wt_tmp(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('sens - out of memory')")  
  stop 'sens - out of memory'
 endif 
 !
 if (trim(fit_type)=='dgelss') then 
     allocate (Tsing(parmax,parmax),Wspace(lwork),stat=alloc)
     if (alloc/=0) then 
      write(6,"('Tsing - out of memory')")  
      stop 'Tsing - out of memory'
     endif 
 endif
 !

 parinitial = param 
 last_param = param

! BPM MAKE FIRST LINEAR PARAMATER (ALL 0 COEFFCIENTS ) V0 FOR SUBTRACTING FROM
! ENERGY LATER (CURRENTLY NOT USING)

 V0 = PARAM(NPARAM_EQ+1) 

!---input energy values------!
 npts = 0
 yes = .true. 
 Nout = 1
 do while( yes )
    npts=npts+1
    if (npts>enermax) then
      write(f_out,"('Too many ab initio points, increase enermax:',I6,' vs ',I6)") npts,enermax
      stop 'Too many ab initio points, increase enermax'
    endif
    read (f_inp,*) coord_tmp(1:NDEG,npts),energy_tmp(npts),wt_tmp(npts)
    if ( wt_tmp(npts).lt.0.0 ) wt_tmp(npts) = abs(wt_tmp(npts))*factor_back
    if ( coord_tmp(1,npts).lt.0.0 ) yes = .false.
 enddo 
 !
 ! Number of ab initio points  
 npts=npts-1
 !
 allocate (coord(NDEG,npts),stat=alloc)
 if (alloc/=0) then 
  write(6,"('coord - out of memory')")  
  stop 'coord - out of memory'
 endif 

 coord(1:NDEG,1:npts) = coord_tmp(1:NDEG,1:npts)!
 !
 !coord(4,1:npts)  = coord_tmp(5,1:npts)
 !coord(5,1:npts)  = coord_tmp(4,1:npts)
 !
 !coord(8,1:npts)  = coord_tmp(9,1:npts)
 !coord(9,1:npts)  = coord_tmp(8,1:npts)
 !
 !coord(11,1:npts)  = coord_tmp(12,1:npts)
 !coord(12,1:npts)  = coord_tmp(11,1:npts)
 !
 deallocate(coord_tmp)
 !
 allocate ( energy(npts),wt(npts),stat=alloc)
 if (alloc/=0) then 
  write(6,"('wt,energy - out of memory')")  
  stop 'wt,energy - out of memory'
 endif 
 !
 energy(1:npts) = energy_tmp(1:npts)
 wt(1:npts) = wt_tmp(1:npts)
 !
 deallocate(energy_tmp,wt_tmp)
 !
 allocate (sigma(npts),dV(npts,parmax),wt_tmp(npts),eps(npts),stat=alloc)
 if (alloc/=0) then 
  write(6,"('sigma,vt,esp - out of memory')")  
  stop 'coord,energy,eps - out of memory'
 endif 
 !
 allocate(icol(parmax))
 ncol=0
 do  i=1,parmax
     if (ivar(i) .ne. 0) then
       ncol = ncol + 1
       icol(ncol) = i
     endif
 enddo
 !
 numpar  = 0
 do i=1,parmax
   if (ivar(i) .gt. 0) numpar=numpar+1
 enddo
 !
 !
 ! We introduce a "zero point energy" shift to avoid loosing accuracy 
 ! because of too too large numbers
 Ener_min = minval(energy(1:npts))*conv_energy
 ZPE = 0.0d0 
 !
 write(6,"('Lowest energy = ',g15.9)") minval(energy(1:npts))  
 !
 ! Convert the energies to the internal dimension, e.g. 1/cm 
 energy(1:npts) = ( energy(1:npts)-ZPE )*conv_energy
 !
 ! For the Robust fit we need to define accuracy of the abinitio data "sigma"
 !
 !
 !      ADD IN WEIGHTING FUNCTION HERE. NEED TO CHANGE MAX ENERGIES ETC FOR
 !      DIFFERENT WEIGHTING. COULD ADD AS INPUT PARAMETER
 !
 if (enermax_ps>0) then 
    !
    do i = 1, npts
      !
      tempe = energy(i)- Ener_min ! param(nparam_eq+1)
      !tempe = max(tempe,8000.0d0)
      !
      wt(i) = (tanh(-gamma_ps*( tempe - enermax_ps)) +1.000200020002d0)/2.000200020002d0
      !
    end do
    !
 endif
 ! 
 ! "normalising" the weight factors
 nused=0
 wtsum=0.0d0
 do i=1,npts
   if (wt(i) > 0.0) nused=nused+1
   wtsum=wtsum+wt(i)
 enddo 
 wtsum=wtsum/nused
 wt(:)=wt(:)/wtsum
 !
 do nrow=1,npts
    if (wt(nrow)>0.0) then 
      sigma(nrow) = 1.0d0/wt(nrow)
    else 
      sigma(nrow) = 1000.0 
    endif 
  enddo
 !
 ! a_wats = 1.0d0

 !
 ! Sometimes pot. parameters are dependent. They are nor to be in the fit 
 ! If we meet a dependent parameter, we get rid of it and start the fit again 
 ! This is the outer loop for
 !
 !
 !  make a table header 
 !
 if (itmax.ne.0) then 
   write(f_out,"(/3x,66('-'))") 
   write(f_out,"('   |   iter  |  points |    deviat     |     rms       |  stability  |')")
   write(f_out,"(3x,66('-'))") 
 endif 
 !
 still_run = .true.
 outer_loop: do while (still_run)  
   !
   ! Parameters to control current and previous st.deviations and stability
   !
   stadev_old = 1.e10
   stability =  1.e10
   stadev    =  1.e10
   ! 
   ! Parameters to control excluding dependent parameters   
   ! We need to know the maximum and minimum values of the pot. parameters weights "ivar" 
   !
   ! numpar is to count the number fitted varying pot. parameters  
   !
   numpar  = 0
   maxivar = 0
   minivar = 1000000
   do i=1,parmax
     if (ivar(i) .gt. 0) numpar=numpar+1
     if (ivar(i) .gt. maxivar) maxivar = ivar(i)
     if (ivar(i) .lt. minivar .and. ivar(i).ne.0) minivar = ivar(i)
   enddo 
   !
   if (nused<numpar) then 
       write(f_out,"('warning! the number of parameter >  the number of data points',2I6)") nused,numpar
       !stop 'warning! the number of parameter <  the number of data points'
    endif
   if (nused==numpar) then 
      write(f_out,"('warning! the same number of parameters and data points')")
   endif 
   !
   ncol=0
   do  i=1,parmax
       if (ivar(i) /= 0) then
         ncol = ncol + 1
         icol(ncol) = i
       endif
   enddo
   !
   ! number of actual paramters to vary 
   numpar = ncol 
   !
   if (allocated(al)) then 
     deallocate(al,bl,dx,ai,rjacob,sens)
   endif
   !
   if (allocated(dF)) deallocate(dF)
   !
   ivar1 =0 
   do  i=1,parmax
       if ( ivar(i)/=0 .or. abs(param(i))>1e-12 ) ivar1(i) = 1
   enddo
   !
   !$omp parallel private(dF,alloc_p) shared(dV) 
   allocate (dF(parmax),stat=alloc_p)
   if (alloc_p/=0)  then 
   write(6,"('dF - out of memory')") 
      stop 'dF out of memory'
   endif 
   !
   !$omp do private(nrow,local) schedule(static)
   do nrow=1,npts
     !
     ! Running bond lengths
     !
     local(1:NBONDS)=coord(1:NBONDS,nrow)*conv_length
     !
     ! Running interbond angles
     !
     local(NBONDS+1:NDEG)=coord(NBONDS+1:NBONDS+NANGS,nrow)*pi/180.0d0
     !
     call diff_V_tau(ndeg,nparam_eq,parmax,ipower,ivar,local,param,df)
     !
     dV(nrow,1:parmax) = df(1:parmax) 
     !
   enddo
   !$omp end do
   !
   deallocate(dF,stat=alloc_p2)
   if (alloc_p2/=0) then 
     write (f_out,"('deallocate df - error')")
     stop 'deallocate df - error'
   endif
   !$omp end parallel
   !
   if (verbose>=6) write(6,"('diff_V_tau ...done!')") 
   !
   if (verbose>=6) print*,"Check redundancy ..."
   !
   !$omp parallel do private(ipar,redundancy) shared(param,ivar) schedule(dynamic)
   do ipar = nparam_eq+1,parmax
     !
     redundancy = sum(dv(1:npts,ipar)**2)
     !
     ! exclude from the fit if derivative at all geometetries is zero
     !
     if (redundancy<1e-12) then 
       !
       if (verbose>=6) then
         !
         !omp critical
         if (ivar(ipar)/=0.and.abs(param(ipar))>1e-15) write(f_out,"(i8,'-th is redundant - ',a19,<ndeg>i4)") ipar,parnam(ipar),ipower(1:NDEG,ipar)
         !omp end critical
         !
       endif
       !
       if (ivar(ipar)/=0) then 
         if (ireplace_by==0) param(ipar) = 0 
         ivar(ipar) = 0
       endif 
       !
     endif
     !
   enddo
   !
   ! renumber useful parameters
   !
   ncol=0
   do  i=1,parmax
       if (ivar(i) .ne. 0) then
         ncol = ncol + 1
         icol(ncol) = i
       endif
   enddo
   !
   ! number of actual paramters to vary 
   numpar = ncol
   !
   allocate (al(numpar,numpar),bl(numpar),dx(numpar),ai(numpar,numpar),rjacob(enermax,numpar),&
            sens(numpar),stat=alloc)
   if (alloc/=0) then 
    write(6,"('al,bl - out of memory')")  
    stop 'al,bl - out of memory'
   endif 
   !
   !##################################################c
   !############ start of the fitting ################c
   !##################################################c
   rjacob = 0 
   iter = 0
   sens = 0.d0 
   do while( iter<=itmax .and. stadev>stadev_best .and. (  stability.ge.stab_best) )   
     iter = iter + 1
     !
     if (verbose>=5) write(6,"('iter = ',i9)") iter
     !
     ! if factor_out is negative, the optimisation of the number of paramters won't be performed 
     if (iter>10.and.factor_out<0) still_run = .false.
     !
     ssq=0.0d+00
     rms=0.0d+00
     !
     parold = param
     !
     dV = 0 
     if (verbose>=6) write(6,"('npoints = ',i9)") npts
     !
     if (allocated(dF)) deallocate(dF)
     !
     !$omp parallel private(dF,alloc_p) shared(dV) 
     allocate (dF(parmax),stat=alloc_p)
     if (alloc_p/=0)  then 
     write(6,"('dF - out of memory')") 
        stop 'dF out of memory'
     endif 
     !
     !$omp do private(nrow,local) schedule(static)
     do nrow=1,npts
       !
       !omp critical
       !if (verbose>=6) write(6,"('n = ',i9)") nrow
       !omp end critical
       !
       ! Running bond lengths
       !
       local(1:NBONDS)=coord(1:NBONDS,nrow)*conv_length
       !
       ! Running interbond angles
       !
       local(NBONDS+1:NDEG)=coord(NBONDS+1:NBONDS+NANGS,nrow)*pi/180.0d0
       !
       call diff_V_tau(ndeg,nparam_eq,parmax,ipower,ivar1,local,param,df)
       !
       dV(nrow,1:parmax) = df(1:parmax) 
       !
     enddo
     !$omp end do
     !
     deallocate(dF,stat=alloc_p2)
     if (alloc_p2/=0) then 
       write (f_out,"('deallocate df - error')")
       stop 'deallocate df - error'
     endif
     !$omp end parallel
     !
     if (verbose>=6) write(6,"('diff_V_tau ...done!')") 
     !
     ! renumber useful parameters
     !
     ncol=0
     do  i=1,parmax
         if (ivar(i) .ne. 0) then
           ncol = ncol + 1
           icol(ncol) = i
         endif
     enddo
     !
     ! number of actual paramters to vary 
     numpar = ncol
     !
     if (verbose>=6) then
       write(6,"(' number of actual paramters to vary = ',i8,' out ',i8)") numpar,parmax
     endif
     !
     if (verbose>=6) print*,"...done!"
     !
     if (allocated(dF)) deallocate(dF)
     !
     allocate (dF(parmax),stat=alloc_p)
     if (alloc_p/=0)  then 
     write(6,"('dF2 - out of memory')") 
        stop 'dF2 out of memory'
     endif 
     !
     if (verbose>=6) print*,"...Generate the fitting Jacobi matrix"
     !
     do nrow=1,npts
       !
       if (verbose>=7) write(6,"('nrow = ',i9)") nrow
       !
       ! Running bond lengths
       !
       local(1:NBONDS)=coord(1:NBONDS,nrow)*conv_length
       !
       ! Running interbond angles
       !
       local(NBONDS+1:NBONDS+NANGS)=coord(NBONDS+1:NANGS+NBONDS,nrow)*pi/180.0d0
       ! 
       ! Value of the potential energy function at the current geometry 
       !
       dF = dV(nrow,:)
       call potential(NDEG,NBONDS,NANGS,NPARAM_EQ,parmax,dF,local,param,ipower,ivar1,0,v)
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       if (itmax/=0 .and. wt(nrow)>0.0) then
         !$omp parallel do private(kcol,i) shared(rjacob) schedule(dynamic)
         do kcol=1,ncol                
           !
           i = icol(kcol)
           !
           !if (verbose>=8) write(6,"('k = ',i0,' i= ',i0)") kcol,i
           !
           call potential(NDEG,NBONDS,NANGS,nparam_eq,parmax,dF,local,param,ipower,ivar1,i,rjacob(nrow,kcol))
         enddo ! --- ncol
        !$omp end parallel do
       endif
       !
       eps(nrow) = energy(nrow)-v
       !
       ! Weighted st. square deviation 
       !
       ssq=ssq+eps(nrow)*eps(nrow)*wt(nrow)
       !
       ! mean square 
       !
       rms=rms+eps(nrow)*eps(nrow)
     enddo  ! ---  nrow 
     !
     ! We constract a set of linear equations A x = B
     !

     if (itmax.ne.0) then
        ! form A matrix 
        !$omp parallel do private(irow,icolumn) shared(al) schedule(guided)
        do irow=1,numpar       !==== row-...... ====!
          do icolumn=1,irow    !==== column-....====!
            al(irow,icolumn)=sum(rjacob(1:npts,icolumn)*rjacob(1:npts,irow)*wt(1:npts))
            al(icolumn,irow)=al(irow,icolumn)
          enddo
        enddo
        !$omp end parallel do
        !
        ! form B matrix
        !$omp parallel do private(irow) shared(Bl) schedule(guided) 
        do irow=1,numpar       !==== row-...... ====!
          bl(irow)=sum(eps(1:npts)*rjacob(1:npts,irow)*wt(1:npts))
        enddo   
        !$omp end parallel do
        !
        ! Solve the set of linear equations 
        !
       select case (trim(fit_type))
       !
       case default
         write (6,"('fit_type ',a,' unknown')") trim(fit_type)
         stop 'fit_type unknown'
       case('linur')
         !
         !call linur(numpar,numpar,al,bl,dx,ierror)
         !
         call linur_backwards(numpar,numpar,al,bl,dx,ierror)
         !
          if (ierror.ne.0) then
            ncol=0
            do i=1,parmax
              if (ivar(i) .ne. 0) then
                ncol=ncol+1
                if ( ncol.eq.ierror ) then
                    ivar(i) = 0
                    param(i) = 0 
                    if (ireplace_by/=0) param(i) = parinitial(i)
                    write(f_out,"(i8,'-th is out - ',a7,1x,<ndeg>i4)") i,parnam(i),ipower(1:NDEG,i)
                endif
              endif
            enddo
            cycle outer_loop
          endif
          !
       case ('dgelss')
         !
         !write(6,"('DGELSS: lapack was not presented ')")
         !stop 'dgelss not compiled'
         !
         ai = al
         !
         !
         !call DGELSS(numpar,numpar,1,AL,numpar,BL,numpar,Tsing,1.D-12,RANK0,wspace,lwork,ierror)
         !
         if (info/=0) then
           write(6,"('DGELSS:error',i8)") ierror
           stop 'dgelss'
         endif
         !
         dx = bl
         !
      case ("SVD","svd")
         !
         call lapack_dgesvd(tol,rjacob(1:npts,1:numpar))
         !
         !call dgemv('T',npts,numpar,alpha,rjacob,enermax,eps,1,beta,dx,1) 
         !
         dx(1:numpar) = matmul(eps(1:npts)*wt(1:npts),rjacob(1:npts,1:numpar))
         !
         !
      case ("SDD","sdd")
         !
         call lapack_sdd_pseudo_inverse(tol,rjacob(1:npts,1:numpar),NumSVD)
         !
         !call dgemv('T',npts,numpar,alpha,rjacob,enermax,eps,1,beta,dx,1) 
         !
         dx(1:numpar) = matmul(eps(1:npts)*wt(1:npts),rjacob(1:npts,1:numpar))
         !
         if (iter==1) then 
           write(6,"(/a,1x,i8,1x,a,g12.5)") 'Number of SDD roots = ',NumSVD,'tol = ',tol
         endif
         !
       end select
       !
       !
       if ( nused>numpar ) then 
          stadev=dsqrt(ssq/float(nused-numpar))
       else 
          stadev=dsqrt(ssq/nused)
       endif
       !
        ncol=0
        do i=1,parmax
          if (ivar(i) /= 0) then
             ncol=ncol+1
             param(i)=param(i)+dx(ncol)*fit_scale
          endif
        enddo
        !
        ! Calcualte the inverse to A matrix, needed for the standard error of the fitted pot. parameters
        !
        call invmat(al,ai,numpar,numpar)
        !
        !
        ncol = 0 
        do i=1,parmax
          if (ivar(i) /= 0) then
             ncol=ncol+1
            !  
            ! If nused = numpar, the problem is degenerated, Ax = B is a system of ordinary equations
            ! st. error doesn't make sense  
            !
            if (nused==numpar) then  
              sterr(ncol)=0
            else
              !
              ! Standard definition of the st. error 
              !
              sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev
              !
            endif
          endif
        enddo    
        !
        ! Here we define stability to see if the fit is converged 
        !      
        stability=abs( (stadev-stadev_old)/stadev )
        stadev_old=stadev
        !
        if (RobustFit>0) then 
         !
         !Watson alpha-parameter
         ! 
         !do it0 = 1,20 
         ! da1 = 0
         ! da2 = 0
         ! do nrow=1,npts
         !   if (wt(nrow) .gt. 0.0d0) then 
         !     da1 = da1+eps(nrow)**2/( sigma(nrow)**2+a_wats*eps(nrow)**2 )
         !     da2 = da2+eps(nrow)**4/( sigma(nrow)**2+a_wats*eps(nrow)**2 )**2
         !   endif 
         ! enddo 
         ! !
         ! da =( da1 -  float(nused-numpar) )/da2
         ! a_wats = a_wats + da
         ! !
         ! if (a_wats<0.0000001) a_wats = 0.001+it0*0.01
         !enddo
         !
         ! 
         !  adjusting the weights  
         ! 
         ! We try a_wats = 0.25        
         ! a_wats = 0.25        
         do nrow=1,npts
            if (wt(nrow) .gt. 0.0d0) then 
              wt(nrow) = 1.d0/( sigma(nrow)**2 + a_wats*eps(nrow)**2 )
            endif 
         enddo 
         ! 
         ! "re-normalising" the weight factors
         !
         nused=0
         wtsum=0.0d0
         do i=1,npts
           if (wt(i) > 0.0) nused=nused+1
           wtsum=wtsum+wt(i)
         enddo 
         wtsum=wtsum/nused
         wt(1:npts)=wt(1:npts)/wtsum
       
         ssq = sum(eps(1:npts)**2*wt(1:npts))
         !
         ! We can recalculated the stand. deviation with new weigh factors
         ! otherwise stadev is defined as a weighted RMS 
         !
         if (nused.eq.numpar) then 
           stadev=dsqrt(ssq/float(nused))
         else 
           stadev=dsqrt(ssq/float(nused-numpar))
         endif 
        endif ! --- robust fit
        !    
      else   ! Itermax = 0, so there is no fit
             ! only straightforward calculations 
         !
         stadev_old=stadev
         stadev=dsqrt(ssq/float(nused))
      endif 
    
      !
      ! Here we define the root mean square 
      !
      rms = dsqrt(rms/npts)
      !
      ! Do some output 
      !
      if (itmax/=0 .and. stadev < maxdeviat ) then 
          write (f_out,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |')") &
                 iter,nused,stadev,rms,stability
     endif 
     !
   enddo  ! --- iter
   !################ end of iterations  #######################

   !
   ! update the energies and residuals based on the last change 
   !
   do nrow=1,npts
     !
     ! Running bond lengths
     !
     local(1:NBONDS)=coord(1:NBONDS,nrow)*conv_length
     !
     ! Running interbond angles
     !
     local(NBONDS+1:NBONDS+NANGS)=coord(NBONDS+1:NANGS+NBONDS,nrow)*pi/180.0d0
     ! 
     ! Value of the potential energy function at the current geometry 
     !
     call potential(NDEG,NBONDS,NANGS,NPARAM_EQ,parmax,dV(nrow,:),local,param,ipower,ivar1,0,v)
     !
     eps(nrow) = energy(nrow)-v
     !
     ! Weighted st. square deviation 
     !
     ssq=ssq+eps(nrow)*eps(nrow)*wt(nrow)
     !
     ! mean square 
     !
     rms=rms+eps(nrow)*eps(nrow)
     !
   enddo  ! ---  nrow 
   !
   if (RobustFit>0) write(f_out,"(/'Watson Alpha parameter',f8.4,' fitted with rms2 = ',f8.4)") &
                                        a_wats !,da1/dfloat((nused-numpar))
   !
   ! if RobustFit==2,3 after first iteration it converts to the standard fit RobustFit=0, but with robust weights
   !
   if ( RobustFit==2 ) then
     RobustFit = 0 
   else if ( RobustFit==3 ) then ! Averaging after the first iteration 
     do i=1,npts
       imin = max(   1,i-npts/ishift) 
       imax = min(npts,i+npts/ishift) 
       ssq = sum(wt(imin:imax))/(imax-imin)
       wt_tmp(i) = ssq
     end do
     wtsum=sum(wt_tmp(1:npts))
     wtsum=wtsum/nused
     wt(1:npts)=wt_tmp(1:npts)/wtsum
     RobustFit = 0 
     cycle outer_loop    
   end if 
   !
   ! We need to know the equilibrium value of the pot. function to substract it from the 
   ! output energies 
   !
   ve = param(nparam_eq+1) ! -ZPE*conv_energy

   ! Output some staqtistics and results 
   !
   !  only if we are fitting:  
   !
   if (itmax.ne.0.and.find_best_params) then
     ! 
     ! It's time to do the output of pot. parameters 
     ! We print them out rounded with their st.errors
     ! 
      write (f_out,"(//'Potential parameters rounded in accord. with their standart errors'/)")
      l = 0 
      do i=1,parmax
        if (ivar(i) .ne. 0) then
           l=l+1
           ndigits = 0
           conf_int = sterr(l)
           do while (conf_int.le.10.0.and.ndigits<10)
             ndigits = ndigits +1 
             conf_int = conf_int*10
           enddo
           if (conf_int>1e8) conf_int = 0 
           !
           write(fmt0,"(i2)") ndigits
           fmt0 = adjustl(fmt0)
           !fmt1 = '(a2,I4,2x,f14.'//fmt0//'    ''    ('',i8,'')'')'
           fmt1 = '(a2,12i3,I4,2x,g14.'//fmt0//'    ''    ('',i8,'')'')'
           !
           write(fmt1,"(a,i,a,a,a,a)") '(a2,',ndeg,'i3,I4,2x,f16.',fmt0,'    ''    ('',i8,'')'')'
           !
           write (f_out,FMT=fmt1) trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i),nint(conf_int)
           !write (f_out,FMT=fmt1) trim(parnam(i)),ivar(i),param(i),nint(conf_int)
           !
        else 
           ! 
           ndigits =2
           if (param(i).ne.0.0) ndigits = 8

           write(fmt0,"(i2)") ndigits
           fmt0 = adjustl(fmt0)
           !
           fmt1 = '(a2,<ndeg>i3,I4,2x,g14.'//fmt0//')'
           if (i<=nparam_eq) then 
             fmt1 = '(a10,21x,I4,2x,g14.'//fmt0//')'
           endif
           !
           write(fmt1,"(a,i,a,a,a)") '(a2,',ndeg,'i3,I4,2x,g14.',fmt0,')'
           !
           !fmt1 = '(a2,12i3,I4,2x,g14.'//fmt0//')'
           write (f_out,FMT=fmt1) trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)


           !do i=1,nparam_eq
           !   write (f_out,"(a10,27x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
           !enddo 
           !do i=nparam_eq+1,parmax
           !   write (f_out,"(a10,9i3,i4,2x,f18.8)") parnam(i),ipower(1:9,i),ivar(i),param(i)
           !enddo 

           !
        endif
      enddo  ! --- i


     write (f_out,"('------------------------------')")
    
     write (f_out,"(a80)") (title(i), i=1,4)

     do i=1,nparam_eq
        write (f_out,"(a10,21x,i4,2x,f18.8)") trim(parnam(i)),ivar(i),param(i)
     enddo 
     do i=nparam_eq+1,parmax
        write (f_out,"(a2,<ndeg>i4,i4,2x,f18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
     enddo 
     !
     !--  printing  out the resulted information 
     !
     inquire(f_res,opened=ifopen)
     if ( ifopen ) then
       rewind(f_res)
     else
       open  (f_res,file='00.res',status='replace' )
     endif
 
     if (itmax.ne.0) then 
       write (f_res,"(a80)") (title(i), i=1,4)
       do i=1,nparam_eq
          write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
       enddo 
       do i=nparam_eq+1,parmax
          write (f_res,"(a2,<ndeg>i4,i4,2x,g18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
       enddo 
     endif  ! ---- itmax/=0
     !
     write (f_res,"('------------------------------')")
     !    
     write (f_res,"(a80)") (title(i), i=1,4)
     do i=1,nparam_eq
        write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
     enddo 
     do i=nparam_eq+1,parmax
        write (f_res,"(a2,<ndeg>i4,i4,2x,g18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
     enddo 
     !
     do nrow=1,npts
       !
        v = energy(nrow)-eps(nrow)
        write (f_res,"(<NBONDS>f15.6,<NANGS>f11.4,3f18.7,2x,e9.2)")            & 
              coord(1:NDEG,nrow),                                  &
              energy(nrow)-ve,v-ve,eps(nrow),                   &
              wt(nrow)
     enddo
     !
     if (itmax.ne.0) then 
       write(f_out,"(/3x,66('-'))") 
       write(f_out,"('   |   iter  |  points |    deviat     |     rms       |  stability  |')")
       write(f_out,"(3x,66('-'))") 
     endif 
     !
     7007  format(a10,I4,2x,f14.4,'    (',i8,')')
     7008  format(a10,I4,2x,f14.4)
     !
   endif
   !
   if (itmax.ne.0.and.find_best_params) then
     !
     ! If any of the pot. parameters (last_i) has been removed from the fit at the previous iteration step
     ! at this moment we check whether the st. deviation is worse than it was before we got rid of the pot. parameter 
     ! in this case we take the last_i-pot. parameter back (for one parameter)
     ! this very tricky process is controled by factor_back 
     !
     if (last_i>0) then
       !
       if ( (stadev-last_stadev)/stadev > Nout*factor_back .and.abs(last_param(last_i))<too_big) then 
         ivar = ivar0     
         ivar(last_i)  =  ivar(last_i)+1
         param = last_param
         last_stadev = 10000000.
         write(f_out,"(i6,'-th is back - ',A10)") last_i,parnam(last_i)
         cycle outer_loop    
       endif 
      endif
      !
      ! Another statistics of the fit: 
      ! correlations between pot.parameters 
      !
      write (f_out,"(/'Correlation coefficients larger than 0.9: ')")
      do i1=1,numpar
        do i2=i1+1,numpar
          corr=ai(i1,i2)/sqrt(abs(ai(i1,i1)*ai(i2,i2)))
           if (abs(corr) > 0.9) write (f_out,"(10x,'corr(',i4,',',i4,') = ',f10.7)") i1,i2,abs(corr)
        enddo
      enddo
      ! 
      ! sensitivity
      !
      ncol = 0 
      do i=1,parmax
        if (ivar(i) /=  0) then
          ncol=ncol+1
          if (nused /= numpar) then  
            ssq = sum(rjacob(1:npts,ncol)**2*wt(1:npts))
            sens(ncol)=stadev*0.1d0/dfloat(numpar)/dsqrt( ssq/dfloat(nused) )
          endif
        endif
      enddo    
      ! 
      ! Robust fit secion 
      !
      !
      ! Pot. parameters output with their st. errors 
      !
      write (f_out,"(/15x,80('-')/15x,':     old parm   :     new parm   :  delta parm  :   std. err.  : sensitility'/15x,80('-'))")
      l = 0
      do i=1,parmax
        if (ivar(i) /= 0) then
          l=l+1
          !
          epsil=param(i)-parold(i)
          crr(l)=(sterr(l)*al(l,l)/stadev)**2
          !
          mark  = ' '
          if (abs(sens(l)) >= abs(epsil) ) mark = '!'
          !
          if (abs(param(i)) > abs(sterr(l))) then
             write (f_out,"(i3,2x,a10,1x,2(':',f14.2,2x),':',3(f14.3,':'),a1)") &
                    l,parnam(i),parold(i),param(i),epsil,sterr(l),sens(l),mark
          else
              write (f_out,"(i3,2('#'),a10,1x,2(':',f14.2,2x),':',3(f14.3,':'),a1,'#')") & 
                      l,parnam(i),parold(i),param(i),epsil,sterr(l),sens(l),mark
          endif
        endif
      enddo
      !
      !  Here we check if there is a pot. parameter with too big st. error  
      !  We consider a parameter to be unneccesary and remove from the fit 
      !  keeping in mind that we might have to take it back in case the st.deviation will get worse
      !
      !
      ! Search for the largest relative sterr
      !
      tivar = minivar
      ivar0 = ivar
      icorr = 0
      last_param = param
      Nout = 0
      do while ( tivar.lt.maxivar .and.factor_out>0)
        if (verbose>=5) write (f_out,"(' tivar, factor_out =',i5,f18.8)") tivar,factor_out
        relat_sterr_min = 1.0e20
        ilow_sterr = -1

        relat_sterr = 0.0
        ilarge_sterr = -1

        l = numpar+1
        do i=parmax,1,-1
          if (ivar(i) /= 0) then
            l=l-1
            !
            ! Get read of parameters with large rel.st.error and their dependencies 
            !
            if (verbose>=5) write (f_out,"(' param,number, ivar, st.error ',f18.7,i5,i5,f18.7)") &
                 param(i),i,ivar(i),abs(sterr(l)/param(i))
            !
            if ( i>nparam_eq.and.(abs(param(i))>too_big.or.param(i)<-too_big*10.0).and.ivar(i)/=maxivar .and. last_i/=i) then
              last_param = param
              last_stadev = stadev
              last_i = i
              ivar(i) = 0
              !
              ! The parameter gets zero value  
              !
              param(i) = 0.d0
              write(f_out,"(' parameter ',i6,' is  out - ',A10)") i,parnam(i) 
              cycle outer_loop    
            endif 
            !
            if (param(i)/=0.0d0.and.ivar(i) == tivar) then
              if (icorr(i)==0.and.factor_out<abs(sterr(l)/param(i))) then
                if (verbose>=5) write (f_out,"(' param, number, st.err., ivar',f18.7,i5,f18.7,i5)") &
                 param(i),i,abs(sterr(l)/param(i)),ivar(i)
                ! relat_sterr=abs(sterr(l)/param(i))
                ! ilarge_sterr = i
                icorr(i) = 1
                l0=0
                relat_sterr = abs(sterr(l)/param(i))
                ilarge_sterr = i
                do i0=1,parmax,1
                  if (ivar(i0) /= 0) then
                    l0=l0+1
                    corr=abs(ai(l,l0))/sqrt(abs(ai(l,l)*ai(l0,l0)))
                    !   if (verbose>=4) then 
                    !     write (f_out,"('corr(',i,',',i,') = ',f14.7)") l,l0,corr
                    !   endif 
                    if (i0/=i.and.param(i0)/=0.0d0.and.ivar(i0) == tivar) then
                      if (icorr(i0)==0.and.corr>fact_corr ) then
                        if (verbose>=4) write (f_out,"('corr, param, (st. err.)  number',2f18.7,f18.7,i4)") & 
                                      corr,param(i0),sterr(l0),i0
                        icorr(i0) = 1
                        !
                        ! Here one finds the (lowest) largest st. error and which par. it belongs to
                        if ( relat_sterr<abs(sterr(l0)/param(i0))) then
                            relat_sterr=abs(sterr(l0)/param(i0))
                            ilarge_sterr = i0
                        endif
                      endif
                    endif
                  endif
                enddo
                ivar(ilarge_sterr)  = 0
                param(ilarge_sterr) = 0.0
                if (ireplace_by/=0) param(ilarge_sterr) = parinitial(ilarge_sterr)
                if (verbose>2) then 
                  write(6,"('relat_sterr,ilarge_sterr, tivar:',f18.8,2i5)") relat_sterr,ilarge_sterr, tivar
                endif 
                Nout = Nout+1
                !
                ! Here one finds the (lowest) largest st. error and which par. it belongs to
                if ( relat_sterr_min>relat_sterr) then
                    relat_sterr_min=relat_sterr
                   ilow_sterr = ilarge_sterr
                endif
              endif
            endif
          endif
        enddo
         if (verbose>2) then 
           write(6,"('relat_sterr_min,ilow_sterr, tivar:',d18.8,2i5)") relat_sterr_min,ilow_sterr, tivar
        endif 
        if (verbose>=4) then 
          do i=1,nparam_eq
             write (f_out,"(a10,22x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
          enddo 
          do i=nparam_eq+1,parmax
             write (f_out,"(a2,<ndeg>i4,i4,2x,f18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
          enddo 
        endif 
        !
        if (relat_sterr_min > factor_out .and. ilow_sterr/=-1) then
          last_ivar = ilow_sterr
          ! last_param = param
          last_stadev = stadev
          last_i = ilow_sterr
          ! ivar(ilarge_sterr) = 0
          !
          ! The parameter gets zero value  
          !
          ! param(ilarge_sterr) = 0.d0
          write(f_out,"(I6,' parameters are out - ',A10)") Nout
          cycle outer_loop    
        endif 
        tivar = tivar + 1
      enddo

      still_run = .false.


!     tivar = minivar
!     !
!     do while ( tivar.lt.maxivar )
!       l = numpar+1
!       !
!       do i=parmax,1,-1
!         if (ivar(i) .ne. 0) then
!           l=l-1
!           epsil=param(i)-parold(i)
!           crr(l)=(sterr(l)*al(l,l)/stadev)**2
!           if (factor_out*abs(param(i)) < abs(sterr(l)) .and. ivar(i) == tivar) then
!              last_ivar = ivar(i)
!              last_param = param
!              last_stadev = stadev
!              last_i = i
!              ivar(i) = 0
!              !
!              ! The parameter gets zero value  
!              !
!              param(i) = 0.d0
!              write(f_out,"(I6,'-th is out - ',A10)") i,parnam(i)
!              cycle outer_loop    
!           endif
!         endif
!       enddo
!       tivar = tivar + 1
!     enddo

    endif 
    still_run = .false.


!--  printing  out the resulted information 

   inquire(f_res,opened=ifopen)
   if ( ifopen ) then
     rewind(f_res)
   else
     open  (f_res,file='00.res',status='replace' )
   endif

  if (itmax.ne.0) then 
   write (f_res,"(a80)") (title(i), i=1,4)
   do i=1,nparam_eq
      write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
   enddo 
   do i=nparam_eq+1,parmax
      write (f_res,"(a2,<ndeg>i4,i4,2x,g18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
   enddo 
  endif  ! ---- itmax/=0
  
  write (f_res,"('------------------------------')")

    
   write (f_res,"(a80)") (title(i), i=1,4)
   do i=1,nparam_eq
      write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
   enddo 
   do i=nparam_eq+1,parmax
      write (f_res,"(a2,<ndeg>i4,i4,2x,g18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
   enddo 
   !
   do nrow=1,npts
     !
      v = energy(nrow)-eps(nrow)
      write (f_res,"(<NBONDS>f15.6,<NANGS>f11.4,3f18.7,2x,e9.2)")            & 
            coord(1:NDEG,nrow),                                  &
            energy(nrow)-ve,v-ve,eps(nrow),                   &
            wt(nrow)
   enddo
   !
   if (itmax.ne.0) then 
     write(f_out,"(/3x,66('-'))") 
     write(f_out,"('   |   iter  |  points |    deviat     |     rms       |  stability  |')")
     write(f_out,"(3x,66('-'))") 
   endif 
   !
   ! Here we can stop the fit 
   !
 enddo outer_loop   


!--  printing  out the resulted information 

 inquire(f_res,opened=ifopen)
 if ( ifopen ) then
    rewind(f_res)
 else
    open  (f_res,file='00.res',status='replace' )
 endif

 if (itmax.ne.0) then 
   write (f_res,"(a80)") (title(i), i=1,4)
   do i=1,nparam_eq
      write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
   enddo 
   do i=nparam_eq+1,parmax
      write (f_res,"(a2,<ndeg>i4,i4,2x,g18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
   enddo 
 endif  ! ---- itmax/=0
  
 write (f_res,"('------------------------------')")
 !
 write (f_res,"(a80)") (title(i), i=1,4)
 do i=1,nparam_eq
    write (f_res,"(a10,27x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
 enddo 
 do i=nparam_eq+1,parmax
    if (ivar(i)/=0.or.abs(param(i))>1e-10) &
       write (f_res,"(a2,<ndeg>i4,i4,2x,g18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
 enddo 
 

do nrow=1,npts
   !
   theta1 = coord(10,nrow)*pi/180.d0
   theta1 = mod(theta1+2.0d0*pi,2.0d0*pi)
   !
   ! subtract equilbrium theta values to make a1/a2 zero at equilibrium
   ! and ensure consistent transfroms
   !
   theta2 = coord(11,nrow)*pi/180.d0
   theta3 = coord(12,nrow)*pi/180.d0
   !
   ! Swap dihedrals as in kinetic energy operator
   theta12 = mod(theta2-theta1+2.0d0*pi,2.0d0*pi) ! 1 and 2 swapped
   theta23 = mod(theta3-theta2+2.0d0*pi,2.0d0*pi) ! 
   theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
   !
   Sa   = ( 2.0d0*theta23 - theta13 - theta12 )/sqrt(6.0d0)*180.d0/pi
   Sb   = (                 theta13 - theta12 )/sqrt(2.0d0)*180.d0/pi
   !
   if (theta1<-1e-9) then
     theta1 = theta1 + 2.0_ark*pi
   endif
   if (theta1>2.0_ark*pi+1e-12) then
     theta1 = theta1 - 2.0_ark*pi
   endif
   !
   tau = (theta1+theta1+theta12+theta1+theta12+theta23-2.0d0*pi)/3.0d0*180.d0/pi
   !
    v = energy(nrow)-eps(nrow)
    write (f_res,"(<NBONDS>f15.6,<NANGS>f15.8,2x,3f15.8,2x,3f18.7,2x,e9.2)")     &
          coord(1:NDEG,nrow),tau,Sa,Sb,                                  &
          energy(nrow)-ve,v-ve,eps(nrow),                   &
          wt(nrow)
   !
 enddo



 open  (f_par,file='00.par',status='replace' )

 do i=nparam_eq+1,parmax
    !
    if (abs(param(i))<1e-12) cycle
    !
    af(:) = 'i1,'
    do k =1,12
      if ( ipower(k,i) >9) af(k) = 'i2,'
    enddo
    !
    write(fmt1,"(a,<ndeg>a,a)") '(''f'',',af(1),af(2),af(3),af(4),af(5),af(6),af(7),af(8),af(9),af(10),af(11),af(12),'    ''  := '',g17.8'':'')'
    write(fmt2,"(a,<ndeg>a,a)") '(''h'',',af(1),af(2),af(3),af(4),af(5),af(6),af(7),af(8),af(9),af(10),af(11),af(12),'    ''  := '',g17.8'':'')'
    !
    if (ipower(12,i)>=0) then 
      write (f_par,fmt1) ipower(1:NDEG,i),param(i)
    else
      write (f_par,fmt2) abs(ipower(1:NDEG,i)),param(i)
    endif
 enddo 



!--  printing  out the fitted PES in the form of the original ab initio
!--  for refitting with the lineariyed coordinates later 

 !open  (f_new,file='new.en',status='replace' )


 !ve = param(10)
 !write (f_new,"('ve=',e20.8)") ve


 !do nrow=1,npts
 !  !
 !  v = (energy(nrow)-eps(nrow))+ZPE
 !  !energy(:) = ( energy(:)-ZPE )*conv_energy
 !  if (v<20000.) &
 !  write (f_res,"(4f9.6,5f10.4,3g18.5,2x,e8.2)")    &
 !       coord(1:9,nrow),                           &
 !       v/conv_energy,-wt(nrow)
 !enddo


 write (f_res,"(/3x,59('-')/'   |  param  |  points |   deviat    |    rms      | stability |')")
 write (f_res,"(3x,59('-'))")
 write (f_res,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e11.3,' |')") &
               numpar,nused,stadev,rms,stability

 write (f_res,"(/92('#'))")



 write (f_out,"(/3x,59('-')/'   |  param  |  points |   deviat    |    rms      | stability |')")
 write (f_out,"(3x,59('-'))")
 write (f_out,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e11.3,' |')") &
               numpar,nused,stadev,rms,stability

 write (f_out,"(/92('#'))")

 end subroutine do_fitting
 !
 !

 !
 recursive subroutine potential(NDEG,NBONDS,NANGS,NPARAM_EQ,parmax,dF,local,param,ipower,ivar,ipar,V)
 !
 implicit none 
 !
 integer,intent(in) ::           parmax,NDEG,NBONDS,NANGS,nparam_eq
 double precision,intent(in) ::  dF(parmax),local(NDEG),param(parmax)
 integer,intent(in) ::           ipar,ipower(NDEG,parmax),ivar(parmax)
 
 double precision,allocatable         ::  param_(:),dF_(:)
 !
 double precision         ::  deltax,potright,potleft,V
 integer                  ::  i,alloc
 integer                  ::  iverbose = 6
 !
 ! parameter for the finite differncies differention 
 !
 double precision,parameter :: factordeltax=0.001
  !
  V = 0
  !
  if (ipar==0) then
    !
    do i = nparam_eq+1,parmax
      !
      !k = ind(i)
      !
      V = V + param(i)*dF(i)
      !
    enddo
    !
  elseif (ipar>nparam_eq) then
    !
    !k = ind(ipar)
    !
    V = dF(ipar)
    !
  else
    !
    !if (iverbose>=6) write(6,"('Allocation ...')") 
    !
    allocate (param_(parmax),dF_(parmax),stat=alloc)
    if (alloc/=0) then 
     write(6,"('param_ - out of memory')")  
     stop 'param_ - out of memory'
    endif 
    !
    deltax=factordeltax*abs(param(ipar))
    if (deltax .le. 1e-15) deltax=1e-6
    !
    param_ = param
    !
    param_(ipar) = param(ipar)+deltax
    !
    !call diff_V_tau_MEP(nmax,local,param_(1:2),dF_)
    !
    call diff_V_tau(NDEG,NPARAM_EQ,PARMAX,ipower,ivar,local,param_(1:parmax),dF_)
    !
    call potential(NDEG,NBONDS,NANGS,NPARAM_EQ,parmax,dF_,local,param_,ipower,ivar,0,potright)
    !
    param_(ipar) = param(ipar)-deltax
    !
    !call diff_V_tau(nmax,local,param_(1:2),dF_)
    !
    call diff_V_tau(NDEG,NPARAM_EQ,parmax,ipower,ivar,local,param_,dF_)
    !
    call potential(NDEG,NBONDS,NANGS,NPARAM_EQ,parmax,dF_,local,param_,ipower,ivar,0,potleft)
    !
    V=(potright-potleft)/(2.d0*deltax)
    !
    deallocate(param_,dF_)
    !
  endif
  !
 end subroutine potential
 !
 function aacos(x,ierror,txt) result (f)
    !
    double precision,intent(in)   :: x
    character(len=40),intent(in),optional :: txt
    integer,intent(out)    :: ierror
    double precision              :: f,small_,pi
    !
    small_ = epsilon(1.0d0)
    pi = 4.0d0 * datan2(1.0d0,1.0d0)
    !
    ierror = 0
    !
    if ( abs(x)>1.0d0+10.0*sqrt(small_) ) then 
       !
       write (6,"('|cos(x)|>1: ',f18.8,2x,a)") x,txt
       !
       ierror = 1
       !
       !stop 'aacos - bad cosalpha'
       !
    elseif ( x>=1.0d0) then 
       f = 0.0d0
    elseif ( x<=-1.0d0) then 
       f = pi
    else 
       f = acos(x)
    endif
    !
 end function aacos
 !

  recursive subroutine diff_V_tau(NDEG,nparam_eq,parmax,ipower,ivar,local,param,dV)
  !
  implicit none
  !
  !integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
  !integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
  !
  integer, parameter :: nsym = 6
  integer, intent(in) :: parmax,NDEG,nparam_eq
  integer,intent(in)  :: ipower(NDEG,parmax),ivar(parmax)
  double precision,intent(in)  :: local(NDEG),param(parmax)
  double precision,intent(out) :: dV(parmax)
  !
  double precision         ::  deg
  double precision         ::  xi(NDEG),pi,term
  DOUBLE PRECISION         ::  TAU,a1,a2
  integer                  ::  ioper,iparam,itau
  character(len=40)        ::  txt

  double precision :: rCO,ROH,RCOe,ROHe,alpha_COHe,alpha_OCHe,deltae,RH1,RH2,RH3,RH1E
  double precision :: b1,b3,b4,b5,theta1,theta2,theta3,y(NDEG),theta12,theta23,theta13,theta1e,theta2e,theta3e
  double precision :: chi(NDEG,nsym)
  logical :: do_cos


!c 
!o , 1, rco 
!h , 2, rch0, 1, acoh
!h1, 1, rch1, 2, aoch1  3  ahh1
!h2, 1, rch2, 2, aoch2, 3, ahh2
!h3, 1, rch3, 2, aoch3, 3, ahh3}
!
!Variables:
!re1= 1.4296
!re2= 0.95887
!re3= 1.092294
!re4= 1.092294
!re5= 1.092294
!ae1= 107.9812
!ae2= 110.6646
!ae3= 110.6646
!ae4= 110.6646
!rco   =re1+ s1
!rch0  =re2+ s2
!rch1  =re3+ s3
!rch2  =re4+ s4
!rch3  =re5+ s5
!acoh  =ae1+ s6
!aoch1 =ae2+ s7
!aoch2 =ae2+ s8
!aoch3 =ae2+ s9

!A1   = s10
!A2   = s11
!tbar = s12
!
!ahh1 = tbar+1.0/3.0*sqrt(2.0)*A2                             !10
!ahh2 = 120.0+tbar-1.0/6.0*sqrt(2.0)*A2-1.0/6.0*sqrt(6.0)*A1  !11
!ahh3 = 240.0+tbar+1.0/6.0*sqrt(6.0)*A1-1.0/6.0*sqrt(2.0)*A2  !12

  !
  if (verbose>=6) print*,"diff_V_tau"
  !
  txt = 'diff_V_tau'
  !
  pi = 4.0d0 * datan2(1.0d0,1.0d0)
  !
  deg=pi/180.0d0
  !
  !
  rCO      = local(1)
  rOH      = local(2)
  rH1      = local(3)
  rH2      = local(4)
  rH3      = local(5)
  !
  rCOe      = param(1)
  rOHe      = param(2)
  rH1e      = param(3)
  !
  alpha_COHe = param(4)*deg
  alpha_OCHe = param(5)*deg
  deltae = 2.0d0*pi/3.0d0
  !
  theta1e = 0.0d0
  theta2e = 120.0d0*deg
  theta3e = 240.0d0*deg
  !
  a1       = param(6)
  a2       = param(7)
  b1       = param(8)
  ! C-O
  xi(1)=1.0d+00-exp(-a1*(rCO-rCOe))
  xi(2)=1.0d+00-exp(-a2*(rOH-rOHe))
  !CH3
  xi(3)=1.0d+00-exp(-b1*(rH1-rH1e))
  xi(4)=1.0d+00-exp(-b1*(rH2-rH1e))
  xi(5)=1.0d+00-exp(-b1*(rH3-rH1e))
  !
  !OCH
  xi(6) = (local(6))- (alpha_COHe)
  !
  ! alphas
  !H2-C-H3,H1-C-H2,H1-C-H3
  xi(7) = local(7)- alpha_OCHe
  xi(8) = local(8)- alpha_OCHe
  xi(9) = local(9)- alpha_OCHe
  !
  theta1 = local(10)
  theta1 = mod(theta1+2.0d0*pi,2.0d0*pi)
  !
  ! subtract equilbrium theta values to make a1/a2 zero at equilibrium
  ! and ensure consistent transfroms
  !
  theta2 = local(11)
  theta3 = local(12)
  !
  ! Swap dihedrals as in kinetic energy operator
  theta12 = mod(theta2-theta1+2.0d0*pi,2.0d0*pi) ! 1 and 2 swapped
  theta23 = mod(theta3-theta2+2.0d0*pi,2.0d0*pi) ! 
  !theta13 = mod(theta1-theta3+2.0d0*pi,2.0d0*pi)
  theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
  !
 !xi(10)   = ( 2.0d0*theta23 - theta13 - theta12 )/sqrt(6.0d0)
 !xi(11)   = (                 theta13 - theta12 )/sqrt(2.0d0)

  xi(10)   = ( 2.0d0*theta23 - theta13 - theta12 )/sqrt(6.0d0)
  xi(11)   = (                 theta13 - theta12 )/sqrt(2.0d0)
  !
  if (theta1<-1e-9) then
    theta1 = theta1 + 2.0_ark*pi
  endif
  if (theta1>2.0_ark*pi+1e-12) then
    theta1 = theta1 - 2.0_ark*pi
  endif
  !
  tau = (theta1+theta1+theta12+theta1+theta12+theta23-2.0d0*pi)/3.0d0
  !
  xi(12) = mod(tau,2.0_ark*pi)
  !
  call ML_symmetry_transformation_XY3_VI(nsym,xi,chi,ndeg)
  !
  dV = 0
  !
  do iparam = nparam_eq+1,parmax
    !
    if (ivar(iparam)==0) cycle
    !
    do_cos = .true.
    itau = ipower(12,iparam)
    !    
    if (itau<0) then
      do_cos = .false.
      itau = abs(itau)
    endif
    !
    term = 0
    !
    do ioper =1,Nsym
      !
      y(1:11) = chi(1:11,ioper)**ipower(1:11,iparam)
      if (do_cos) then
        !print*,itau,chi(21,ioper)
        y(12) = cos(real(itau,8)*chi(12,ioper))
      else
        y(12) = sin(real(itau,8)*chi(12,ioper))
      endif
      !
      term = term + product(y(:))
      !
    enddo
    !
    dV(iparam) = term/6.0d0
    !
  enddo
  !
end subroutine diff_V_tau



!
!   Skip n lines  in the input file 
!
  subroutine skiplines( inpunit,n )
  integer,intent(in) :: n,inpunit
  character(len=80)  :: label
  integer            :: i0

    do i0=1,n
       read  (inpunit,"(a80)") label 
    enddo

  end subroutine skiplines


 subroutine linur_backwards(dimen,npar,coeff,constant,solution,error)
   !
   implicit none

  integer,intent(in)  :: dimen,npar
  integer,intent(out) :: error 
  double precision,intent(in)  :: coeff(npar,npar),constant(npar)
  double precision,intent(out) :: solution(npar)
  double precision          :: a0(npar,npar)
  double precision          :: c
  integer                   :: i1,i2,i,k8,k,l,k9

  !----- begin ----!
  
    do i1=1,dimen
    do i2=1,dimen 
       a0(i1,i2)=coeff(i1,i2)
    enddo
    enddo

    do i=1,dimen
      solution(i)=constant(i)
    enddo
    !
    error=0
    do i=dimen,1,-1
      c=0
      do k=dimen,i+1,-1
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
        !      write(6,*) '(',i,'-th adj. parameter is wrong )'
        error=i
        return
      endif

      a0(i,i)=sqrt(a0(i,i)-c)
      !
      if (a0(i,i).eq.0) then
         !
         !  write(6,*) '(',i,'-th adj. parameter is wrong )'
         error=i
         return
      endif
      !
      do l=i-1,1,-1
         c=0.0
         do k=dimen,i+1,-1
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    !
    do i=dimen,1,-1
      c=0.0
      do k=dimen,i+1,-1
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    !
    do i1=dimen,1,-1
      i=1+dimen-i1
      c=0.0
      do k=i-1,1,-1
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    !
  return
  end subroutine linur_backwards


 subroutine linur(dimen,npar,coeff,constant,solution,error)
   !
   implicit none

  integer,intent(in)  :: dimen,npar
  integer,intent(out) :: error 
  double precision,intent(in)  :: coeff(npar,npar),constant(npar)
  double precision,intent(out) :: solution(npar)
  double precision          :: a0(npar,npar)
  double precision          :: c
  integer                   :: i1,i2,i,k8,k,l,k9

  !----- begin ----!
  
    do i1=1,dimen
    do i2=1,dimen 
       a0(i1,i2)=coeff(i1,i2)
    enddo
    enddo

    do i=1,dimen
      solution(i)=constant(i)
    enddo
    error=0
    do i=1,dimen
      c=0
      k8=i-1
      do k=1,k8
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
       error=i
       return
      endif

      a0(i,i)=sqrt(a0(i,i)-c)
      if (a0(i,i).eq.0) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
         error=i
         return
      endif
      k8=i+1
      do l=k8,dimen
         k9=i-1
         c=0.0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0.0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0.0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
  return
  end subroutine linur

!------------------------------------------!
  subroutine invmat(al,ai,dimen,npar)
  implicit none
  integer,intent(in)           :: npar,dimen
  double precision,intent(in)  :: al(npar,npar)
  double precision,intent(out) :: ai(npar,npar)
  double precision             :: h(npar),p,q
  integer                      :: i1,k,i,j,k8,k9
      

    ai(1:dimen,1:dimen)=al(1:dimen,1:dimen)
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
        h(i)=q/p
        if(i.le.k) h(i)=-q/p
        do j=2,i
          k8=i-1
          k9=j-1
          ai(k8,k9)=ai(i,j)+q*h(j)
        enddo 
      enddo 
      ai(dimen,dimen)=1.0/p
      do i=2,dimen
        k8=i-1
        ai(dimen,k8)=h(i)
      enddo 
   end do 
   do i=1,dimen
     k8=i-1
     do j=1,k8
       ai(j,i)=ai(i,j)
     enddo 
   enddo 
   return
 end subroutine invmat
!------------------------------------------!


  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  recursive subroutine ML_symmetry_transformation_XY3_III(nsym,src,dst,ndeg)
    implicit none 
    !
    !integer, parameter :: ark         = selected_real_kind(12,25)  ! "Accurate" reals
    integer,intent(in)    :: nsym  ! n of irreps
    integer,intent(in)    :: ndeg  ! n degrees of freedom 
    double precision,intent(in)      :: src(1:ndeg)
    double precision,intent(out)     :: dst(1:ndeg,nsym)
    !
    integer :: ioper
    !
    double precision         :: repres(nsym,ndeg,ndeg),a,b,e,o,pi
    !
    a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
    pi = 4.0d0 * datan2(1.0d0,1.0d0)
    !
    if (nsym>6) then
      write (6,"('symmetry_transformation_local: illegal nsym = ',i8)") nsym
      stop 'symmetry_transformation_local: illegal nsym'
    endif
    !
    repres = 0 
    !
    repres = 0
    !
    repres(:,1,1) = 1.0_ark
    repres(:,2,2) = 1.0_ark
    !
    repres(:,6,6) = 1.0_ark
    !
    repres(:,12,12) = 1.0_ark
    !
    ! E
    ! r123
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    ! a123
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    !d9
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    !
    !C3+/(132)
    repres(2,3,5) = 1.0_ark
    repres(2,4,3) = 1.0_ark
    repres(2,5,4) = 1.0_ark
    !
    repres(2,7,9) = 1.0_ark
    repres(2,8,7) = 1.0_ark
    repres(2,9,8) = 1.0_ark
    !
    repres(2,10,10) = -a
    repres(2,10,11) = -b
    repres(2,11,10) =  b
    repres(2,11,11) = -a
    !
    !C3-/(93)
    !
    repres(3,3,4) = 1.0_ark
    repres(3,4,5) = 1.0_ark
    repres(3,5,3) = 1.0_ark
    !
    repres(3,7,8) = 1.0_ark
    repres(3,8,9) = 1.0_ark
    repres(3,9,7) = 1.0_ark
    !
    repres(3,10,10) = -a
    repres(3,10,11) =  b
    repres(3,11,10) = -b
    repres(3,11,11) = -a
    !
    !C2/(23)->(45)
    !
    repres(4,3,3) = 1.0_ark
    repres(4,4,5) = 1.0_ark
    repres(4,5,4) = 1.0_ark
    !
    repres(4,7,7) = 1.0_ark
    repres(4,8,9) = 1.0_ark
    repres(4,9,8) = 1.0_ark
    !
    repres(4,10,10) =  1.0_ark
    repres(4,11,11) = -1.0_ark
    !
    !C2'/(9)->(34)
    repres(5,3,4) = 1.0_ark
    repres(5,4,3) = 1.0_ark
    repres(5,5,5) = 1.0_ark
    !
    repres(5,7,8)  = 1.0_ark
    repres(5,8,7)  = 1.0_ark
    repres(5,9,9)  = 1.0_ark
    !
    repres(5,10,10) = -a
    repres(5,10,11) =  b
    repres(5,11,10) =  b
    repres(5,11,11) =  a
    !
    !(13)->(35)
    repres(6,3,5) = 1.0_ark
    repres(6,4,4) = 1.0_ark
    repres(6,5,3) = 1.0_ark
    !
    repres(6,7,9) = 1.0_ark
    repres(6,8,8) = 1.0_ark
    repres(6,9,7) = 1.0_ark
    !
    repres(6,10,10) = -a
    repres(6,10,11) = -b
    repres(6,11,10) = -b
    repres(6,11,11) =  a
    !
    do ioper = 1,nsym
      dst(:,ioper) = matmul(repres(ioper,:,:),src) 
    enddo
    !
    dst(12,1) = src(12)
    dst(12,2) = src(12)+2.0d0*pi/3.0d0
    dst(12,3) = src(12)-2.0d0*pi/3.0d0
    dst(12,4) =-src(12)
    dst(12,5) =-src(12)-2.0d0*pi/3.0d0
    dst(12,6) =-src(12)+2.0d0*pi/3.0d0
    !
  end subroutine ML_symmetry_transformation_XY3_III      
  


  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  recursive subroutine ML_symmetry_transformation_XY3_IV(nsym,src,dst,ndeg)
    implicit none 
    !
    !integer, parameter :: ark         = selected_real_kind(12,25)  ! "Accurate" reals
    integer,intent(in)    :: nsym  ! n of irreps
    integer,intent(in)    :: ndeg  ! n degrees of freedom 
    double precision,intent(in)      :: src(1:ndeg)
    double precision,intent(out)     :: dst(1:ndeg,nsym)
    !
    integer :: ioper
    !
    double precision         :: repres(nsym,ndeg,ndeg),a,b,e,o,pi
    !
    a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
    pi = 4.0d0 * datan2(1.0d0,1.0d0)
    !
    if (verbose>=6) print*,"ML_symmetry_transformation_XY3_IV"
    !
    if (nsym>6) then
      write (6,"('symmetry_transformation_local: illegal nsym = ',i8)") nsym
      stop 'symmetry_transformation_local: illegal nsym'
    endif
    !
    repres = 0 
    !
    repres = 0
    !
    repres(:,1,1) = 1.0_ark
    repres(:,2,2) = 1.0_ark
    !
    repres(:,6,6) = 1.0_ark
    !
    repres(:,12,12) = 1.0_ark
    !
    ! E
    ! r123
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    ! a123
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    !d9
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    !
    !C3+/(132)
    repres(2,3,5) = 1.0_ark
    repres(2,4,3) = 1.0_ark
    repres(2,5,4) = 1.0_ark
    !
    repres(2,7,9) = 1.0_ark
    repres(2,8,7) = 1.0_ark
    repres(2,9,8) = 1.0_ark
    !
    repres(2,10,10) = -a
    repres(2,10,11) = -b
    repres(2,11,10) =  b
    repres(2,11,11) = -a
    !
    !C3-/(93)
    !
    repres(3,3,4) = 1.0_ark
    repres(3,4,5) = 1.0_ark
    repres(3,5,3) = 1.0_ark
    !
    repres(3,7,8) = 1.0_ark
    repres(3,8,9) = 1.0_ark
    repres(3,9,7) = 1.0_ark
    !
    repres(3,10,10) = -a
    repres(3,10,11) =  b
    repres(3,11,10) = -b
    repres(3,11,11) = -a
    !
    !C2/(23)->(45)
    !
    repres(4,3,3) = 1.0_ark
    repres(4,4,5) = 1.0_ark
    repres(4,5,4) = 1.0_ark
    !
    repres(4,7,7) = 1.0_ark
    repres(4,8,9) = 1.0_ark
    repres(4,9,8) = 1.0_ark
    !
    repres(4,10,10) =  1.0_ark
    repres(4,11,11) = -1.0_ark
    !
    !C2'/(9)->(34)
    repres(6,3,4) = 1.0_ark
    repres(6,4,3) = 1.0_ark
    repres(6,5,5) = 1.0_ark
    !
    repres(6,7,8)  = 1.0_ark
    repres(6,8,7)  = 1.0_ark
    repres(6,9,9)  = 1.0_ark
    !
    repres(6,10,10) = -a
    repres(6,10,11) =  b
    repres(6,11,10) =  b
    repres(6,11,11) =  a
    !
    !(13)->(35)
    repres(5,3,5) = 1.0_ark
    repres(5,4,4) = 1.0_ark
    repres(5,5,3) = 1.0_ark
    !
    repres(5,7,9) = 1.0_ark
    repres(5,8,8) = 1.0_ark
    repres(5,9,7) = 1.0_ark
    !
    repres(5,10,10) = -a
    repres(5,10,11) = -b
    repres(5,11,10) = -b
    repres(5,11,11) =  a
    !
    do ioper = 1,nsym
      dst(:,ioper) = matmul(repres(ioper,:,:),src) 
    enddo
    !
    dst(12,1) = src(12)
    dst(12,2) = src(12)-2.0d0*pi/3.0d0
    dst(12,3) = src(12)+2.0d0*pi/3.0d0
    dst(12,4) =-src(12)
    dst(12,5) =-src(12)-2.0d0*pi/3.0d0
    dst(12,6) =-src(12)+2.0d0*pi/3.0d0
    !
  end subroutine ML_symmetry_transformation_XY3_IV  


   recursive subroutine ML_symmetry_transformation_XY3_V(nsym,src,dst,ndeg)
    implicit none 
    !
    !integer, parameter :: ark         = selected_real_kind(12,25)  ! "Accurate" reals
    integer,intent(in)    :: nsym  ! n of irreps
    integer,intent(in)    :: ndeg  ! n degrees of freedom 
    double precision,intent(in)      :: src(1:ndeg)
    double precision,intent(out)     :: dst(1:ndeg,nsym)
    !
    integer :: ioper
    !
    double precision         :: repres(nsym,ndeg,ndeg),a,b,e,o,pi
    !
    a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
    pi = 4.0d0 * datan2(1.0d0,1.0d0)
    !
    if (verbose>=6) print*,"ML_symmetry_transformation_XY3_V"
    !
    if (nsym>6) then
      write (6,"('symmetry_transformation_local: illegal nsym = ',i8)") nsym
      stop 'symmetry_transformation_local: illegal nsym'
    endif
    !
    repres = 0 
    !
    ! repres = 0
    !
    repres(:,1,1) = 1.0_ark
    repres(:,2,2) = 1.0_ark
    !
    repres(:,6,6) = 1.0_ark
    !
    repres(:,12,12) = 1.0_ark
    !
    ! E
    ! r123
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    ! a123
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    !d9
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    !
    !C3+/(132)
    repres(3,3,5) = 1.0_ark
    repres(3,4,3) = 1.0_ark
    repres(3,5,4) = 1.0_ark
    !
    repres(3,7,9) = 1.0_ark
    repres(3,8,7) = 1.0_ark
    repres(3,9,8) = 1.0_ark
    !
    repres(3,10,10) = -a
    repres(3,10,11) = -b
    repres(3,11,10) =  b
    repres(3,11,11) = -a
    !
    !C3-/(93)
    !
    repres(2,3,4) = 1.0_ark
    repres(2,4,5) = 1.0_ark
    repres(2,5,3) = 1.0_ark
    !
    repres(2,7,8) = 1.0_ark
    repres(2,8,9) = 1.0_ark
    repres(2,9,7) = 1.0_ark
    !
    repres(2,10,10) = -a
    repres(2,10,11) =  b
    repres(2,11,10) = -b
    repres(2,11,11) = -a
    !
    !C2/(23)->(45)
    !
    repres(4,3,3) = 1.0_ark
    repres(4,4,5) = 1.0_ark
    repres(4,5,4) = 1.0_ark
    !
    repres(4,7,7) = 1.0_ark
    repres(4,8,9) = 1.0_ark
    repres(4,9,8) = 1.0_ark
    !
    repres(4,10,10) =  1.0_ark
    repres(4,11,11) = -1.0_ark
    !
    !C2'/(9)->(34)
    repres(5,3,4) = 1.0_ark
    repres(5,4,3) = 1.0_ark
    repres(5,5,5) = 1.0_ark
    !
    repres(5,7,8)  = 1.0_ark
    repres(5,8,7)  = 1.0_ark
    repres(5,9,9)  = 1.0_ark
    !
    repres(5,10,10) = -a
    repres(5,10,11) =  b
    repres(5,11,10) =  b
    repres(5,11,11) =  a
    !
    !(13)->(35)
    repres(6,3,5) = 1.0_ark
    repres(6,4,4) = 1.0_ark
    repres(6,5,3) = 1.0_ark
    !
    repres(6,7,9) = 1.0_ark
    repres(6,8,8) = 1.0_ark
    repres(6,9,7) = 1.0_ark
    !
    repres(6,10,10) = -a
    repres(6,10,11) = -b
    repres(6,11,10) = -b
    repres(6,11,11) =  a
    !
    do ioper = 1,nsym
      dst(:,ioper) = matmul(repres(ioper,:,:),src) 
    enddo
    !
    dst(12,1) = src(12)
    dst(12,2) = mod(src(12)-2.0d0*pi/3.0d0,2.0_ark*pi)
    dst(12,3) = mod(src(12)+2.0d0*pi/3.0d0,2.0_ark*pi)
    dst(12,4) =-src(12)+2.0_ark*pi
    dst(12,5) =2.0_ark*pi-mod(src(12)+4.0d0*pi/3.0d0,2.0_ark*pi)
    dst(12,6) =2.0_ark*pi-mod(src(12)+2.0d0*pi/3.0d0,2.0_ark*pi)
    !
  end subroutine ML_symmetry_transformation_XY3_V  



 recursive subroutine ML_symmetry_transformation_XY3_VI(nsym,src,dst,ndeg)
    implicit none 
    !
    !integer, parameter :: ark         = selected_real_kind(12,25)  ! "Accurate" reals
    integer,intent(in)    :: nsym  ! n of irreps
    integer,intent(in)    :: ndeg  ! n degrees of freedom 
    double precision,intent(in)      :: src(1:ndeg)
    double precision,intent(out)     :: dst(1:ndeg,nsym)
    !
    integer :: ioper
    !
    double precision         :: repres(nsym,ndeg,ndeg),a,b,e,o,pi
    !
    a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
    pi = 4.0d0 * datan2(1.0d0,1.0d0)
    !
    if (verbose>=6) print*,"ML_symmetry_transformation_XY3_V"
    !
    if (nsym>6) then
      write (6,"('symmetry_transformation_local: illegal nsym = ',i8)") nsym
      stop 'symmetry_transformation_local: illegal nsym'
    endif
    !
    repres = 0 
    !
    ! repres = 0
    !
    repres(:,1,1) = 1.0_ark
    repres(:,2,2) = 1.0_ark
    !
    repres(:,6,6) = 1.0_ark
    !
    repres(:,12,12) = 1.0_ark
    !
    ! E
    ! r123
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    ! a123
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    !d9
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    !
    !C3+/(132)
    repres(3,3,5) = 1.0_ark
    repres(3,4,3) = 1.0_ark
    repres(3,5,4) = 1.0_ark
    !
    repres(3,7,9) = 1.0_ark
    repres(3,8,7) = 1.0_ark
    repres(3,9,8) = 1.0_ark
    !
    repres(3,10,10) = -a
    repres(3,10,11) = -b
    repres(3,11,10) =  b
    repres(3,11,11) = -a
    !
    !C3-/(93)
    !
    repres(2,3,4) = 1.0_ark
    repres(2,4,5) = 1.0_ark
    repres(2,5,3) = 1.0_ark
    !
    repres(2,7,8) = 1.0_ark
    repres(2,8,9) = 1.0_ark
    repres(2,9,7) = 1.0_ark
    !
    repres(2,10,10) = -a
    repres(2,10,11) =  b
    repres(2,11,10) = -b
    repres(2,11,11) = -a
    !
    !C2/(23)->(45)
    !
    repres(4,3,3) = 1.0_ark
    repres(4,4,5) = 1.0_ark
    repres(4,5,4) = 1.0_ark
    !
    repres(4,7,7) = 1.0_ark
    repres(4,8,9) = 1.0_ark
    repres(4,9,8) = 1.0_ark
    !
    repres(4,10,10) =  1.0_ark
    repres(4,11,11) = -1.0_ark
    !
    !C2'/(9)->(34)
    repres(5,3,4) = 1.0_ark
    repres(5,4,3) = 1.0_ark
    repres(5,5,5) = 1.0_ark
    !
    repres(5,7,8)  = 1.0_ark
    repres(5,8,7)  = 1.0_ark
    repres(5,9,9)  = 1.0_ark
    !
    repres(5,10,10) = -a
    repres(5,10,11) =  b
    repres(5,11,10) =  b
    repres(5,11,11) =  a
    !
    !(13)->(35)
    repres(6,3,5) = 1.0_ark
    repres(6,4,4) = 1.0_ark
    repres(6,5,3) = 1.0_ark
    !
    repres(6,7,9) = 1.0_ark
    repres(6,8,8) = 1.0_ark
    repres(6,9,7) = 1.0_ark
    !
    repres(6,10,10) = -a
    repres(6,10,11) = -b
    repres(6,11,10) = -b
    repres(6,11,11) =  a
    !
    do ioper = 1,nsym
      dst(:,ioper) = matmul(repres(ioper,:,:),src) 
    enddo
    !
    dst(12,1) = src(12)
    dst(12,2) = mod(src(12)+2.0d0*pi/3.0d0,2.0_ark*pi)
    dst(12,3) = mod(src(12)+4.0d0*pi/3.0d0,2.0_ark*pi)
    dst(12,4) =-src(12)+2.0_ark*pi
    dst(12,5) =2.0_ark*pi-mod(src(12)+2.0d0*pi/3.0d0,2.0_ark*pi)
    dst(12,6) =2.0_ark*pi-mod(src(12)+4.0d0*pi/3.0d0,2.0_ark*pi)
    !
  end subroutine ML_symmetry_transformation_XY3_VI


  
  recursive subroutine ML_symmetry_transformation_XY3_II(nsym,src,dst,ndeg)
    implicit none 
    !
    !integer, parameter :: ark         = selected_real_kind(12,25)  ! "Accurate" reals
    integer,intent(in)    :: nsym  ! n of irreps
    integer,intent(in)    :: ndeg  ! n degrees of freedom 
    double precision,intent(in)      :: src(1:ndeg)
    double precision,intent(out)     :: dst(1:ndeg,nsym)
    !
    integer :: ioper
    !
    double precision         :: repres(nsym,ndeg,ndeg),a,b,e,o,pi
    !
    a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
    pi = 4.0d0 * datan2(1.0d0,1.0d0)
    !
    if (nsym>6) then
      write (6,"('symmetry_transformation_local: illegal nsym = ',i8)") nsym
      stop 'symmetry_transformation_local: illegal nsym'
    endif
    !
    repres = 0 
    !
    repres = 0
    !
    repres(:,1,1) = 1.0_ark
    repres(:,2,2) = 1.0_ark
    !
    repres(:,6,6) = 1.0_ark
    !
    repres(:,12,12) = 1.0_ark
    !
    ! E
    ! r123
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    ! a123
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    !d9
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    !
    !C3+/(132)
    repres(2,3,5) = 1.0_ark
    repres(2,4,3) = 1.0_ark
    repres(2,5,4) = 1.0_ark
    !
    repres(2,7,9) = 1.0_ark
    repres(2,8,7) = 1.0_ark
    repres(2,9,8) = 1.0_ark
    !
    repres(2,10,10) = -a
    repres(2,10,11) = -b
    repres(2,11,10) =  b
    repres(2,11,11) = -a
    !
    !C3-/(93)
    !
    repres(3,3,4) = 1.0_ark
    repres(3,4,5) = 1.0_ark
    repres(3,5,3) = 1.0_ark
    !
    repres(3,7,8) = 1.0_ark
    repres(3,8,9) = 1.0_ark
    repres(3,9,7) = 1.0_ark
    !
    repres(3,10,10) = -a
    repres(3,10,11) =  b
    repres(3,11,10) = -b
    repres(3,11,11) = -a
    !
    !C2/(23)->(45)
    !
    repres(4,3,3) = 1.0_ark
    repres(4,4,5) = 1.0_ark
    repres(4,5,4) = 1.0_ark
    !
    repres(4,7,7) = 1.0_ark
    repres(4,8,9) = 1.0_ark
    repres(4,9,8) = 1.0_ark
    !
    repres(4,10,10) =  1.0_ark
    repres(4,11,11) = -1.0_ark
    !
    !C2'/(9)->(34)
    repres(5,3,4) = 1.0_ark
    repres(5,4,3) = 1.0_ark
    repres(5,5,5) = 1.0_ark
    !
    repres(5,7,8)  = 1.0_ark
    repres(5,8,7)  = 1.0_ark
    repres(5,9,9)  = 1.0_ark
    !
    repres(5,10,10) = -a
    repres(5,10,11) =  b
    repres(5,11,10) =  b
    repres(5,11,11) =  a
    !
    !(13)->(35)
    repres(6,3,5) = 1.0_ark
    repres(6,4,4) = 1.0_ark
    repres(6,5,3) = 1.0_ark
    !
    repres(6,7,9) = 1.0_ark
    repres(6,8,8) = 1.0_ark
    repres(6,9,7) = 1.0_ark
    !
    repres(6,10,10) = -a
    repres(6,10,11) = -b
    repres(6,11,10) = -b
    repres(6,11,11) =  a
    !
    do ioper = 1,nsym
      dst(:,ioper) = matmul(repres(ioper,:,:),src) 
    enddo
    !
    dst(12,1) = src(12)
    dst(12,2) = src(12)+2.0d0*pi/3.0d0
    dst(12,3) = src(12)-2.0d0*pi/3.0d0
    dst(12,4) =-src(12)
    dst(12,5) =-src(12)-2.0d0*pi/3.0d0
    dst(12,6) =-src(12)+2.0d0*pi/3.0d0
    !
  end subroutine ML_symmetry_transformation_XY3_II      
  
  subroutine lapack_dgesvd(tol,h)

    double precision, intent(in)  :: tol 
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    character(len=1) :: jobu,jobvt

    double precision,allocatable    :: work(:),u(:,:),vt(:,:),s(:),v(:,:),ut(:,:)
    integer           :: info
    integer           :: nh1, nh2,i,Nkeep,j,nu1,nu2,nvt1,nvt2,LDVT,LDU
    integer           :: lwork
    double precision  :: tol_
    double precision  :: alpha = 1.0d0,beta=0
    !
    jobu = 'A'
    !
    jobvt = 'A'
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = 50*nh1
    !
    select case (jobu)
       case('A')
         LDU = nh1
         nu1 = nh1
         nu2 = nh1
       case('S')
         LDU = nh1
         nu1 = nh1
         nu2 = min(nh1,nh2)
    end select
    !
    select case (jobvt)
       case('A')
         LDVT = nh2
         nvt1 = nh2
         nvt2 = nh2
       case('S')
         LDVT = min(nh1,nh2)
         nvt1 = min(nh1,nh2)
         nvt2 = nh2
    end select
    !
    allocate(work(lwork),u(LDU,nu2),ut(nu2,LDU),vt(LDVT,nvt2),v(nvt2,LDVT),s(min(nh1,nh2)),stat=info)
    !
    if (info/=0) then
      write (6,"(' lapack_dgesvd cannot allocate',i8)") info
      stop 'lapack_dgesvd cannot allocate'
    end if
    !
    call dgesvd(jobu,jobvt,nh1,nh2,h,nh1,s,u,LDU,vt,LDVT,work,-1,info)
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
    endif
    !
    call dgesvd(jobu,jobvt,nh1,nh2,h,nh1,s,u,LDU,vt,LDVT,work,lwork,info)
    !
    if (info/=0) then
      write (6,"(' dgesvd returned ',i8)") info
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    !
    ! pseudo-inverse 
    !
    Nkeep = 0
    !
    tol_  = tol*maxval(s,dim=1)
    !
    !do  i = 1,nh2
    !  !
    !  if (s(i)<tol_) then
    !      !
    !      Nkeep  = i - 1
    !      exit
    !      !
    !  endif 
    !  !
    !  s_inv(i) = 1.0d0/s(i)
    !  !
    !enddo
    !
    Nkeep = minloc(s,dim=1,mask=s.ge.tol_)
    !
    Nkeep = min(Nkeep,nh1)
    !
    ut = 0 
    !
    do  i = 1,nu1
      do  j = 1,Nkeep
         ut(j,i) = U(i,j)/s(j)
      enddo
    enddo
    !
    !v = transpose(vt)
    u = transpose(ut)
    !
    !ut(1:nh2,:) = matmul(v(1:nh2,1:nkeep),ut(1:Nkeep,:))
    !
    !h = 0
    !
    !h = transpose(ut)
    !
    h(1:LDU,1:nvt2) = matmul(u(1:LDU,1:Nkeep),vt(1:nkeep,1:nvt2))
    !
    !call dgemm('T','N',nh2,nh1,nkeep,alpha,vt,nh2,ut,nh2,beta,ut,nh2)
    !h = transpose(ut)
    !
    !call dgemm('T','N',nh1,nkeep,nh2,alpha,ut,nh2,vt,nh2,beta,h,nh1)
    !
    deallocate(work,u,vt,ut,v,s,stat=info)
    !
    if (info/=0) then
      write (6,"(' lapack_dgesvd cannot allocate',i8)") info
      stop 'lapack_dgesvd cannot allocate'
    end if
    !
  end subroutine lapack_dgesvd  
  
  
  
    subroutine lapack_sdd_pseudo_inverse(tol,h,Nkeep)

    double precision, intent(in)    :: tol 
    double precision, intent(inout) :: h(:,:)  ! In:  matrix
    !                                          ! Out: pseudo-inverse matrix V Sigma- VT
    integer, intent(out)  :: Nkeep
    character(len=1) :: jobu,jobvt,jobz
    !
    double precision,allocatable    :: work(:),u(:,:),vt(:,:),s(:),v(:,:),ut(:,:),A(:,:)
    integer,allocatable    :: iwork(:)
    integer           :: info
    integer           :: nh1, nh2,i,j,nu1,nu2,nvt1,nvt2,LDVT,LDU
    integer           :: lwork
    double precision  :: tol_
    double precision  :: alpha = 1.0d0,beta=0
    !
    jobu  = 'A'
    jobvt = 'A'
    jobz  = 'A'
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = 50*nh1
    !
    select case (jobz)
       case('A')
         LDU = nh1
         nu1 = nh1
         nu2 = nh1
         LDVT = nh2
         nvt1 = nh2
         nvt2 = nh2
       case('S')
         LDU = nh1
         nu1 = nh1
         nu2 = min(nh1,nh2)
         LDVT = min(nh1,nh2)
         nvt1 = min(nh1,nh2)
         nvt2 = nh2
    end select
    !
    allocate(A(nh1,nh2),stat=info)
    !
    A = h
    !
    allocate(work(lwork),iwork(8*min(nh1,nh2)),u(LDU,nu2),ut(nu2,LDU),vt(LDVT,nvt2),v(nvt2,LDVT),s(min(nh1,nh2)),stat=info)
    !
    call dgesdd(jobz,nh1,nh2,A,nh1,s,u,LDU,vt,LDVT,work,-1,iwork,info)
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
    endif
    !
    call dgesdd(jobz,nh1,nh2,A,nh1,s,u,LDU,vt,LDVT,work,lwork,iwork,info)
    !
    if (info/=0) then
      write (6,"(' dgesdd returned ',i8)") info
      stop 'lapack_dgesdd - dgesvd failed'
    end if
    !
    ! pseudo-inverse 
    !
    Nkeep = 0
    !
    tol_  = tol*maxval(s,dim=1)
    !
    Nkeep = minloc(s,dim=1,mask=s.ge.tol_)
    !
    Nkeep = min(Nkeep,nh1)
    !
    ut = 0 
    !
    do  i = 1,nu1
      do  j = 1,Nkeep
         ut(j,i) = U(i,j)/s(j)
      enddo
    enddo
    !
    !v = transpose(vt)
    u = transpose(ut)
    !
    !ut(1:nh2,:) = matmul(v(1:nh2,1:nkeep),ut(1:Nkeep,:))
    !
    !h = 0
    !
    !h = transpose(ut)
    !
    h(1:LDU,1:nvt2) = matmul(u(1:LDU,1:Nkeep),vt(1:nkeep,1:nvt2))
    !
    !call dgemm('T','N',nh2,nh1,nkeep,alpha,vt,nh2,ut,nh2,beta,ut,nh2)
    !h = transpose(ut)
    !
    !call dgemm('T','N',nh1,nkeep,nh2,alpha,ut,nh2,vt,nh2,beta,h,nh1)
    !
    deallocate(work,iwork,u,vt,ut,v,s,A)
    !
  end subroutine lapack_sdd_pseudo_inverse


  
  
  
  
  end module FITPES_ch3cl
  
  
  
  
  program driver
  use FITPES_ch3cl
  
  call do_fitting
  
  end program driver
