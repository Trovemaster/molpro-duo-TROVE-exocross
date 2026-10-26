 module FITDMS_ch3cl
  !
  public ik,rk,do_fitting,coordtype

  private
 

 integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
 integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
 integer,parameter          :: f_inp=5, f_out=6, f_res=9, f_new =10, nparam_eq = 3
 !character(len=40),parameter   ::  coordtype = 'abinitio'
 !character(len=40),parameter   ::  coordtype = '1-2-y'
 character(len=40),parameter   ::  coordtype = 'sym'
 !character(len=40),parameter   ::  coordtype = 'sym-x'



  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains

  
 subroutine do_fitting
 !
 implicit none
 !
 ! Parameters:  
 ! 
 ! number of parameters and maximum number of ab initio points   
 !
 integer,parameter          ::  enermax =  150000, nmax0=8,taumax0=8, nparam = 2000
 !
 ! where the output goes to 
 ! 
 integer,parameter          :: f_inp=5, f_out=6, f_res=9, f_new =10
 ! 
 ! Dinite difference differentiation with 2 or 3 points
 !
 integer,parameter          ::  findif = 3 
 !
 ! the iteration will be stoped when the best standard deviation reached 
 !
 !real(rk),parameter :: stadev_best = 0.001_rk
 !
 integer,parameter          ::  verbose = 4
 ! 
 ! parameter to control the correlation when we through out the parameters with large st. error (0.0..1.0)
 !
 real(rk),parameter :: fact_corr=0.4, too_big = 1e9
 ! 
 ! Watson's alpha parameter to control the Robust fit
 !real(rk),parameter :: a_wats = 0.1
 ! parameter for the finite differncies differention 
 !
 real(rk),parameter :: factordeltax=0.001
 !
 integer,parameter          ::  ireplace_by = 0
 !
 ! Universal constants 
 ! 
 real(rk),parameter :: bohr = 0.52917715e+00,hartre=2.194746354d+05
 ! Avareging region for the robustfit = 3 
 integer,parameter          :: ishift = 10
 !character(len=40),parameter  :: fit_type= 'dgelss'
 character(len=40),parameter  :: fit_type= 'linur'
 !
 ! Variables:  
 real(rk)  :: maxdeviat  ! st. deviation at the previos iteration 
 !
 integer           :: i,npts,nused,itmax,numpar,ncol,kcol,iter,nrow,l,itertypemax,irow,icolumn ! just integers  
 ! 
 character(len=80) :: title(4),longlabel ! input title of the job from the first four input lines 
 character(len=1)  :: mark           ! used at the parameters output to mark the dependent or independent ones 
 character(len=10) :: label          ! Temporary label 
 real(rk)  :: wtsum,ssq,ssq3(3)      ! 
 !
 ! some variables 
 !
 real(rk)  :: stadev,epsil,stadev_old,v,potential,dd,pi
 real(rk)  :: tempx,deltax,potright,potleft,stability
 real(rk)  :: r14,r24,r34,alpha1,alpha2,alpha3
 real(rk)  :: corr,ve
 integer           :: ierror,n0,i1,i2,i3,i4,i5,i6,itertype,last_i,robustfit
 real(rk)  :: stab_best,acc_rms,conf_int,rms0
 integer           :: maxivar,minivar,tivar
 real(rk)  :: last_stadev,factor_back,factor_out
 real(rk)  :: rms,s,da1,da2,da,ZPE,local(10),stadev_best,rms3(3)
 integer           :: last_ivar,ndigits,npoints_delta,rand_itmin,alloc
 integer           :: max_pnts,ipar,i0,l0,parmax,it(9),fivar_t,nmax,taumax
 logical           :: yes,still_run,ifopen,find_best_params
 integer           :: flag_borh, flag_hartre,ivartmp,it0,imin,imax,ilarge_sterr,ilow_sterr,Nout,parmax_t,i_t
 real(rk)  :: conv_lenght,conv_dip,paramtmp,relat_sterr,relat_sterr_min,fparam_t,e_t,ecm_t,edip(3,1),tsing3(3)

 real(rk)  :: dip0,x(4,3),cm(3),y(4,3),r(3),v12(3),v23(3),v31(3),n3(3)
 real(rk)  :: tmat(4,3),rmat(3,3),rmat0(3,3),dip(3),mu(3),cdip(3),crr
 real(rk)  :: a_wats = 0.1,small_,redundancy,tau1,tau2,tau3,phi,tauA,tauB
 
 character(len=68)  :: fmt1,fmt2
 character(len=2)   :: fmt0
 character(len=80)  ::  char_tmp 
 !
 ! some matrices, to be allocated  
 !
 real(rk), allocatable :: rjacob(:,:),last_param(:),param(:)
 real(rk), allocatable :: param_tmp(:)
 real(rk), allocatable :: cart(:,:,:),dipol(:,:),cart_(:,:,:),dipol_(:,:)
 real(rk), allocatable :: eps(:),wt(:),parold(:),sigma(:),wt_tmp(:),parinitial(:)
 real(rk), allocatable :: al(:,:),bl(:),dx(:)
 real(rk), allocatable :: ai(:,:),sterr(:),dmu(:,:,:)
 real(rk), allocatable :: sens(:)
 character(len=10), allocatable:: parnam(:),parnam_tmp(:)
 integer, allocatable :: ivar(:),icorr(:),ivar0(:),ivar_tmp(:),ipower_tmp(:,:),ipower(:,:),icol(:)
 !
 real(rk), allocatable :: Tsing(:) 
 integer                       :: rank0,info,icomponent=0
 integer                       :: lwork
 real(rk),allocatable  :: wspace(:)
 integer           :: ipower_(0:9),iparam
 !
 double precision,allocatable :: amat(:,:),bmat(:,:)
 !
 ! Common block is to reduce the number of assignments within the potential routine
 !
 ! Here we go!

 !
 !  initial constants 
 !
 last_stadev = 1000000.
 maxivar = 8
 maxdeviat = 200.0
 small_ = epsilon(1.0_rk)
 !
 pi = 4.0_rk * atan2(1.0_rk,1.0_rk)
 !
 open(f_res,file='00.res',status='replace' )
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
 read (f_inp,*) itmax
 read (f_inp,*) stab_best,stadev_best
 read (f_inp,*) conv_lenght
 read (f_inp,*) conv_dip
 read (f_inp,*) factor_out
 read (f_inp,*) factor_back
 read (f_inp,*) robustfit
 !
 read (f_inp,*) icomponent ! component to fit
 !
 if (icomponent<0.or.icomponent>4) stop 'illegal icomponent'
 ! 
 write (f_out,"(i10,' <= number of iterations in the fit')")         itmax
 write (f_out,"(f10.6,' <= The fit is over if the stability is reached ')") stab_best
 write (f_out,"(f10.6,' <= sterr parameter out ')")                  factor_out
 write (f_out,"(f10.6,' <= stdev parameter out ')")                  factor_back
 write (f_out,"(i10,  ' <= watson robust on (1) or of (0) ')")       robustfit
 write (f_out,"(i10,  ' <= component to fit (0: fit all three) ')")       icomponent
 !
 call skiplines(f_inp,3) 
 !
 if (robustfit<0) then 
   find_best_params = .false.
 else
   find_best_params = .true.
 endif 
 !
 ! Array's allocation:
 !
 parmax_t = nparam
 lwork =  parmax_t*8
 !
 allocate (param_tmp(parmax_t),parnam_tmp(parmax_t),ivar_tmp(parmax_t),ipower_tmp(0:9,parmax_t),stat=alloc)
 if (alloc/=0) then 
  write(6,"('param_tmp - out of memory')")  
    stop 'param_tmp - out of memory'
 endif 
 !
 ipower_tmp = 0
 !
 ! input  parameters of structur
 !
 parnam_tmp = 'xxxxxx'
 param_tmp = 0 
 ivar_tmp = 0 
 !
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
       read (char_tmp(1:80),*) parnam_tmp(i),ipower_(0:9),fivar_t,fparam_t
       !
       write (f_out,"(a2,10i3,' =',f18.8)") trim(parnam_tmp(i)),ipower_(0:9),fparam_t
       param_tmp(i) = fparam_t
       ivar_tmp(i) = fivar_t
       ipower_tmp(0:9,i) = ipower_(0:9)
       !
    else
       Yes = .False.
    endif
 enddo read_loop
 !
 parmax = i 
 !
 allocate (param(parmax),parnam(parmax),ivar(parmax),icorr(parmax),ivar0(parmax),ipower(0:9,parmax),icol(parmax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('param - out of memory')")  
    stop 'param - out of memory'
 endif 
 !
 if (parmax>nparam) then 
  write(6,"('nparam is too small')")  
    stop 'nparam is too small'
 endif 
 !
 param(1:parmax) = param_tmp(1:parmax)
 ivar(1:parmax) = ivar_tmp(1:parmax)
 ipower(0:9,1:parmax) = ipower_tmp(0:9,1:parmax)
 parnam(1:parmax) = parnam_tmp(1:parmax)
 !
 deallocate (param_tmp,parnam_tmp,ipower_tmp,ivar_tmp)
 !
 allocate (parold(parmax),parinitial(parmax),last_param(parmax),sterr(parmax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('sterr - out of memory')")  
  stop
 endif 
 !
 allocate (Tsing(parmax),Wspace(lwork),al(parmax,parmax),bl(parmax),&
            dx(parmax),ai(parmax,parmax),cart_(4,3,enermax),dipol_(3,enermax),wt_tmp(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('sens - out of memory')")  
  stop
 endif
 !
 parinitial = param 
 !
 !---input dipxy values------!
 npts = 0
 yes = .true. 
 Nout = 1
 !
 cart_ = 0 
 !
 do
    npts=npts+1
    if (npts>enermax) then
         write(f_out,"('Too many ab initio points, increase enermax:',I9,' vs ',I9)") npts,enermax
    endif
    !
    read (f_inp,*) cart_(1,1,npts),cart_(1,2,npts),cart_(1,3,npts),&
                   cart_(2,1,npts),cart_(2,2,npts),cart_(2,3,npts),&
                   cart_(3,1,npts),cart_(3,2,npts),cart_(3,3,npts),&
                   cart_(4,1,npts),cart_(4,2,npts),cart_(4,3,npts),&
                   edip(1,1),edip(2,1),edip(3,1),wt_tmp(npts)!,ecm_t,e_t
    !
    if ( cart_(1,1,npts).lt.-99.0 ) exit
    if ( npts>1e9 ) stop 'npts>1e9'
    !
    if ( wt_tmp(npts).lt.0.0 ) wt_tmp(npts) = abs(wt_tmp(npts))*factor_back
    forall(i=1:4) x(i,1:3) = cart_(i,1:3,npts)
    call cartesian_to_local(x,local,rmat0)
    !
    !if (abs(rmat0(1,1))>small_.and.abs(rmat0(2,2))>small_) then
    !   print*,'non standard rmat',npts
    !endif
    !
    !call linur(3,3,rmat0,edip(1:3,1),mu,ierror)
    !
    !rmat = transpose(rmat0)
    !
    rmat = rmat0
    !
    !if (ierror==0) then 
    !  !
    !  dipol_(:,npts) = mu(1:3)
    !  !
    !else
      !
      ierror = 0
      !
      !call dgelss(3,3,1,rmat,3,edip(1:3,1),3,tsing3,10.0d0*spacing(1.0d0),rank0,wspace,lwork,ierror)
      !
      dip0 = sqrt(sum(edip(1:3,1)**2))
      !
      dipol_(:,npts) = matmul(rmat,edip(1:3,1))
      !
      dip0 = sqrt(sum(dipol_(:,npts)**2))
      !
      phi = atan2(-edip(2,1),edip(1,1))
      !
      rmat0(1,1) = cos(phi)
      rmat0(1,2) =-sin(phi)
      rmat0(2,1) = sin(phi)
      rmat0(2,2) = cos(phi)
      !
      !dipol_(1:2,npts) = matmul(rmat0(1:2,1:2),edip(1:2,1))
      !
      !dipol_(:,npts) = edip(1:3,1)
      !
      !dipol_(1,npts) = edip(2,1)
      !dipol_(2,npts) = edip(1,1)
      !dipol_(3,npts) = edip(3,1)
      !
      if (ierror/=0) then 
        !
        print*,'dgelss 3x3: error',ierror
        stop 'dgelss 3x3'
        !
      endif 
      !
    !endif
    !
    !dipol_(3,npts) = edip(3,1)
    !dipol_(1,npts) = edip(2,1)
    !dipol_(2,npts) = edip(1,1)
    !
 enddo 
 !
 ! Number of ab initio points  
 npts=npts-1
 ! 
 ! "normalising" the weight factors
 nused=0
 wtsum=0
 do i=1,npts
   if (wt_tmp(i) > 0.0) nused=nused+1
   wtsum=wtsum+wt_tmp(i)
 enddo
 !
 numpar  = 0
 do i=1,parmax
   if (ivar(i) .gt. 0) numpar=numpar+1
 enddo 
 !
 allocate ( rjacob(3*npts,numpar),&
            eps(3*npts),wt(3*npts),sigma(3*npts),&
            dmu(3,parmax,npts),cart(4,3,npts),dipol(3,npts),&
            sens(parmax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('rjacob,eps,wt,sigma,dmu,cart,dipol,sens - out of memory')")  
  stop
 endif
 !
 cart(:,:,1:npts) = cart_(:,:,1:npts)*conv_lenght
 dipol(:,1:npts)  = dipol_(:,1:npts)*conv_dip
 !
 deallocate(cart_,dipol_)
 !
 nused = 3*nused
 wtsum=wtsum*3.0_rk/real(nused,8)
 wt(1:npts)=wt_tmp(1:npts)/wtsum
 wt(npts+1:2*npts)= wt(1:npts)
 wt(2*npts+1:3*npts)= wt(1:npts)
 !
 if (icomponent/=0) then
  !
  wtsum=0
  do i=1,npts
    if (wt_tmp(i) > 0.0) nused=nused+1
    wtsum=wtsum+wt_tmp(i)
  enddo
  wtsum=wtsum/real(nused,8)
  !
  wt = 0
  wt((icomponent-1)*npts+1:(icomponent)*npts)= wt_tmp(1:npts)
  !
 endif
 !
 deallocate(wt_tmp)
 !
 ! For the Robust fit we need to dedine accuracy of the abinitio data "sigma"
 ! 

 do nrow=1,3*npts
   if (wt(nrow)>0.0) then 
     sigma(nrow) = 0.0001_rk/wt(nrow)
   else 
     sigma(nrow) = 1000.0 
   endif 
 enddo
 !
 a_wats = 1.0_rk
 !
 ! Sometimes pot. parameters are dependent. They are nor to be in the fit 
 ! If we meet a dependent parameter, we get rid of it and start the fit again 
 ! This is the outer loop for
 !
 !
 !  make a table header 
 !
 if (itmax.ne.0) then 
   write(f_out,"(/3x,67('-'))") 
   write(f_out,"('   |   iter  |  points |    deviat     |     rms       |  stability  |')")
   write(f_out,"(3x,67('-'))") 
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
   dmu = 0
   !
   !$omp parallel do private(nrow,dip,x,local,rmat,iparam) shared(dmu) schedule(dynamic)
   do nrow=1,npts
     !
     if (verbose>=6) write(6,"('n = ',i9)") nrow
     !
     dip(1:3) = dipol(1:3,nrow)
     forall(i=1:4) x(i,1:3) = cart(i,1:3,nrow)
     call cartesian_to_local(x,local,rmat)
     !
     continue
     !
     do iparam = nparam_eq+1,parmax
        !
        dmu(:,iparam,nrow) = diff_mu(ipower(0:9,iparam),local,param(1:nparam_eq))
        !
     enddo
     !
   enddo
   !$omp end parallel do
   ! 
   ! check for redundant parameteres
   !
   do iparam = nparam_eq+1,parmax
     !
     ! dipole component
     !
     i = ipower(0,iparam)
     !
     redundancy = sum(dmu(i,iparam,:)**2)
     !
     ! exclude from the fit if derivative at all geometetries is zero
     !
     if (redundancy<sqrt(small_)) then 
       !
       if (ivar(iparam)/=0.or.abs(param(iparam))>small_) write(f_out,"(i3,'-th component,',i8,'-th is out - ',a10,10i2)") i,iparam,parnam(iparam),ipower(0:9,iparam)
       !
       param(iparam) = 0 
       ivar(iparam) = 0 
       !
     endif
     !
   enddo
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
     if (ivar(i) .gt. 0) then 
       numpar=numpar+1
       icol(numpar) = i
     endif 
     !
     if (ivar(i) .gt. maxivar) maxivar = ivar(i)
     if (ivar(i) .lt. minivar .and. ivar(i).ne.0) minivar = ivar(i)
   enddo 
   !
   if (nused*3<numpar) then 
       write(f_out,"('worning! the number of data points <  the number of parameter',2I6)") nused*3,numpar
       stop 'worning! the number of data points <  the number of parameter'
    endif
   if (nused*3==numpar) then 
      write(f_out,"('worning! the same number of parameters and data points')")
   endif 
   !
   if (allocated(amat)) deallocate(amat)
   if (allocated(bmat)) deallocate(bmat)
   !
   allocate (amat(numpar,numpar),bmat(numpar,1),stat=alloc)
   if (alloc/=0) then 
     write(6,"('amat,bmat - out of memory')")  
     stop
   endif
   !
   !##################################################c
   !############ start fitting #######################c
   !##################################################c
   rjacob = 0 
   iter = 0
   sens = 0 
   do while( iter<=itmax .and. stadev>stadev_best .and. (  stability.ge.stab_best) )   
     iter = iter + 1
     ssq=0.0
     rms=0.0
     ssq3 = 0 
     rms3 = 0
     !
     parold = param
     !
     if (verbose>=5) write(6,"('iter = ',i9)") iter

     if (verbose>=5) write(6,"('dmu - done!')")
     !
     do nrow=1,npts
       !
       if (verbose>=6) write(6,"('n = ',i9)") nrow
       !
       ! Running geometry
       !
       forall(i=1:4) x(i,1:3) = cart(i,1:3,nrow)
       !call cartesian_to_local(x,local,rmat)
       !
       call dipole(parmax,dmu(:,:,nrow),mu,param,0)
       !
       !if (verbose>=6) write(6,"('dipole0 - done! ')")
       !
       dip(1:3) = dipol(1:3,nrow)
       !
       eps(nrow       ) = dip(1)-mu(1)
       eps(nrow+npts  ) = dip(2)-mu(2)
       eps(nrow+2*npts) = dip(3)-mu(3)
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       if (itmax/=0 .and. wt(nrow)+wt(nrow+npts)+wt(nrow+2*npts)>0.0) then
         !
         !$omp parallel do private(kcol,i,mu) shared(rjacob) schedule(dynamic)
         do  kcol=1,numpar
             !
             i  = icol(kcol)
             !
             call dipole(parmax,dmu(:,:,nrow),mu,param,i)
             !
             rjacob(nrow       ,kcol)=mu(1)
             rjacob(nrow+  npts,kcol)=mu(2)
             rjacob(nrow+2*npts,kcol)=mu(3)
             !
         enddo ! --- kcol
         !$omp end parallel do
         !
       endif
       !
       ! Weighted st. square deviation 
       !
       ! root mean square 
       !
       rms0=eps(nrow)**2+eps(nrow+npts)**2+eps(nrow+2*npts)**2
       rms=rms+rms0
       !
       rms3(1)=rms3(1)+eps(nrow)**2
       rms3(2)=rms3(2)+eps(nrow+npts)**2
       rms3(3)=rms3(3)+eps(nrow+2*npts)**2
       !
       ! Weighted st. square deviation 
       !
       ssq3(1)=ssq3(1)+eps(nrow)**2*wt(nrow)
       ssq3(2)=ssq3(2)+eps(nrow+npts)**2*wt(nrow+npts)
       ssq3(3)=ssq3(3)+eps(nrow+2*npts)**2*wt(nrow+2*npts)
       !
       ssq=ssq+eps(nrow)**2*wt(nrow)+eps(nrow+npts)**2*wt(nrow+npts)+ &
               eps(nrow+2*npts)**2*wt(nrow+2*npts)
       !
       !
     enddo  ! ---  nrow
     !
     if (verbose>=5) write(6,"('rjacob - done!')")
     !
     ! We constract a set of linear equations A x = B
     !
     if (itmax.ne.0) then
        ! form A matrix 
        do irow=1,numpar       !==== row-...... ====!
          do icolumn=1,irow    !==== column-....====!
            al(irow,icolumn)=sum(rjacob(:,icolumn)*rjacob(:,irow)*wt(:))
            al(icolumn,irow)=al(irow,icolumn)
          enddo
          !
          if (abs(al(irow,irow))<small_) then 
            i  = icol(irow)
            ivar(i) = 0
            param(i) = 0 
            if (ireplace_by/=0) param(i) = parinitial(i)
            !
            call print_res
            !
            write(f_res,"(i,'-th is out - ',a7,10i2)") i,parnam(i),ipower(0:9,i)
            !
            cycle outer_loop
          endif
          !
        enddo
        ! form B matrix 
        do irow=1,numpar       !==== row-...... ====!
          bl(irow)=sum(eps(:)*rjacob(:,irow)*wt(:))
        enddo   
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
         call linur(numpar,parmax,al,bl,dx,ierror)
         !
         if (verbose>=5) write(6,"('linur - done!')")
         !
         if (ierror.ne.0) then
            i  = icol(ierror)
            ivar(i) = 0
            param(i) = 0 
            if (ireplace_by/=0) param(i) = parinitial(i)
            write(f_res,"(i,'-th is out - ',a7,10i2)") i,parnam(i),ipower(0:9,i)
            if (itmax/=0 .and. stadev < maxdeviat ) then 
                write (f_out,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |')") &
                       iter,nused,stadev,rms,stability
            endif 
            !
            call print_res
            !
            cycle outer_loop
         endif
         !
       case ('dgelss')
         !
         !write(6,"('DGELSS: lapack was not presented ')")
         !stop 'dgelss not compiled'
         !
         amat(1:numpar,1:numpar) = al(1:numpar,1:numpar)
         bmat(1:numpar,1) = bl(1:numpar)
         !
         call DGELSS(numpar,numpar,1,amat,numpar,bmat,numpar,Tsing,1.D-12,rank0,wspace,lwork,info)
         !
         if (info/=0) then
           write(6,"('DGELSS:error',i)") info
           stop 'dgelss'
         endif
         !
         dx(1:numpar) = bmat(1:numpar,1)
         !
       end select
       !
       if ( nused.ne.numpar ) then 
          stadev=dsqrt(ssq/float(nused-numpar))
       else 
          stadev=dsqrt(ssq/nused)
       endif
       !
       do  kcol=1,numpar
          i  = icol(kcol)
          param(i)=param(i)+dx(kcol)
       enddo
       !
       ! Calcualte the inverse to A matrix, needed for the standard error of the fitted pot. parameters
       !
       call invmat(al,ai,numpar,parmax,ierror)
       !
       !if (ierror.ne.0) then
       !   i  = icol(ierror)
       !   ivar(i) = 0
       !   param(i) = 0 
       !   if (ireplace_by/=0) param(i) = parinitial(i)
       !   write(f_res,"(i,'-th is out - ',a7,10i2)") i,parnam(i),ipower(0:9,i)
       !   cycle outer_loop
       !endif
       !
       do  kcol=1,numpar
           !  
           ! If nused = numpar, the problem is degenerated, Ax = B is a system of ordinary equations
           ! st. error doesn't make sense  
           !
           if (nused==numpar) then  
             sterr(kcol)=0
           else
             !
             ! Standard definition of the st. error 
             !
             sterr(kcol)=sqrt(abs(ai(kcol,kcol)))*stadev
             !
           endif
           !
       enddo    
       !
       ! Here we define stability to see if the fit is converged 
       !      
       stability=abs( (stadev-stadev_old)/stadev )
       stadev_old=stadev
       !

       if (RobustFit/=0) then 
        !
        !Watson alpha-parameter
        ! 
        do it0 = 1,20 
          da1 = 0
          da2 = 0
          do nrow=1,3*npts
            if (wt(nrow) .gt. 0.0_rk) then 
              da1 = da1+eps(nrow)**2/( sigma(nrow)**2+a_wats*eps(nrow)**2 )
              da2 = da2+eps(nrow)**4/( sigma(nrow)**2+a_wats*eps(nrow)**2 )**2
            endif 
          enddo 
          !
          da =( da1 -  float(nused-numpar) )/da2
          a_wats = a_wats + da
          !
          if (a_wats<0.0000001) a_wats = 0.001+it0*0.01
        enddo
        !
        ! 
        !  adjusting the weights  
        ! 
        do nrow=1,3*npts
           if (wt(nrow) .gt. 0.0_rk) then 
             wt(nrow) = 1._rk/( sigma(nrow)**2 + a_wats*eps(nrow)**2 )
           endif 
        enddo 

        ! 
        ! "re-normalising" the weight factors
        !
        nused=0
        wtsum=0.0_rk
        do i=1,3*npts
          if (wt(i) > 0.0) nused=nused+1
          wtsum=wtsum+wt(i)
        enddo 
        wtsum=wtsum/nused
        wt(:)=wt(:)/wtsum

        ssq = sum(eps(:)**2*wt(:))
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
      ! Here we define the root mean square 
      !            
     else   ! Itermax = 0, so there is no fit
            ! only straightforward calculations 
      !
      stadev_old=stadev
      stadev=dsqrt(ssq/float(nused))
      ssq3(:)=dsqrt(ssq3(:)/float(nused/3))
     endif 
     !
     ! Here we define the root mean square 
     !
     rms = dsqrt(rms/real(npts*3,8))
     rms3(:) = dsqrt(rms3(:)/real(npts))
     !
     ! Do some output 
     !
     if (itmax/=0 .and. stadev < maxdeviat ) then 
         write (f_out,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |',6e12.5)") &
                iter,nused,stadev,rms,stability,ssq3(:),rms3(:)
     endif 
     !
   enddo  ! --- iter
   !################ end of iterations  #######################
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
       wt(i) = ssq
     end do
     wtsum=sum(wt(1:npts))
     wtsum=wtsum/nused
     wt(1:npts)=wt(1:npts)/wtsum
     wt(npts+1:2*npts)=wt(1:npts)
     wt(2*npts+1:3*npts)=wt(1:npts)
     RobustFit = 0 
     cycle outer_loop    
   end if 
   !
   ! We need to know the equilibrium value of the pot. function to substract it from the 
   ! output energies 
   !
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
           do while (conf_int.le.10.0)
             ndigits = ndigits +1 
             conf_int = conf_int*10
           enddo
           write(fmt0,"(i2)") ndigits
           fmt0 = adjustl(fmt0)
           fmt1 = '(a10,I4,2x,f14.'//fmt0//'    ''    ('',i8,'')'')'
           write (f_out,FMT=fmt1) parnam(i),ivar(i),param(i),nint(conf_int)
           !
        else 
           ndigits =2
           if (param(i).ne.0.0) ndigits = 8

           write(fmt0,"(i2)") ndigits
           fmt0 = adjustl(fmt0)
           fmt1 = '(a10,I4,2x,f14.'//fmt0//')'

           fmt1 = '(a10,I4,2x,f14.'//fmt0//')'
           write (f_out,FMT=fmt1) parnam(i),ivar(i),param(i)
           !
        endif
      enddo  ! --- i


     write (f_out,"('------------------------------')")
    
     write (f_out,"(a80)") (title(i), i=1,4)
     do i=1,parmax
        write (f_out,"(a10,1x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
     enddo 
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
     if ( (stadev-last_stadev)/stadev > Nout*factor_back .and.last_i/=0) then 
       ivar = ivar0     
       ivar(last_i)  =  ivar(last_i)+1
       param = last_param
       last_stadev = 10000000.
       write(f_out,"(I6,'-th is back - ',A10)") last_i,parnam(last_i)
       !
       call print_res
       !
       cycle outer_loop    
     endif 

      !
      ! Another statistics of the fit: 
      ! correlations between pot.parameters 
      !
      write (f_out,"(/'Correlation coefficients larger than 0.9: ')")
      do i1=1,numpar
        do i2=i1+1,numpar
           corr=ai(i1,i2)/sqrt(abs(ai(i1,i1)*ai(i2,i2)))
           if (abs(corr) .gt. 0.9) write (f_out,"(10x,'corr(',i2,',',i2,') = ',f10.7)") i1,i2,abs(corr)
        enddo
      enddo
      ! 
      ! sensitivity
      !
      do  kcol=1,numpar
          if (nused /= numpar) then  
            ssq = sum(rjacob(:,kcol)**2*wt(:))
            sens(kcol)=stadev*0.1_rk/real(numpar,8)/sqrt( ssq/real(nused,8) )
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
          crr=(sterr(l)*al(l,l)/stadev)**2
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
            if (verbose>=5) write (f_out,"(' param, number, ivar, st.error ',f18.7,i5,i5,f18.7)") param(i),i,ivar(i),abs(sterr(l)/param(i))
            !
            if ( i>4.and.abs(param(i))>too_big ) then
              last_param = param
              last_stadev = stadev
              last_i = i
              ivar(i) = 0
              !
              ! The parameter gets zero value  
              !
              param(i) = 0._rk
              write(f_out,"(' parameter ',i6,' is  out - ',A10)") i,parnam(i) 
              call print_res
              cycle outer_loop    
            endif 
            !
            if (param(i)/=0.0_rk.and.ivar(i) == tivar) then
              if (icorr(i)==0.and.factor_out<abs(sterr(l)/param(i))) then
                if (verbose>=5) write (f_out,"(' param, number, st.err., ivar',f18.7,i5,f18.7,i5)") param(i),i,abs(sterr(l)/param(i)),ivar(i)
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
                    if (i0/=i.and.param(i0)/=0.0_rk.and.ivar(i0) == tivar) then
                      if (icorr(i0)==0.and.corr>fact_corr ) then
                        if (verbose>=4) write (f_out,"('corr, param, (st. err.)  number',2f18.7,f18.7,i4)") corr,param(i0),sterr(l0),i0
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
          do i=1,parmax
             write (f_out,"(a10,1x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
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
          ! param(ilarge_sterr) = 0._rk
          write(f_out,"(I6,' parameters are out - ',A10)") Nout
          call print_res
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
!              param(i) = 0._rk
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

  if (itmax.ne.0) then
    !
    !inquire(f_res,opened=ifopen)
    !if ( ifopen ) then
    !
    rewind(f_res)
    !
    !else
    !  open  (f_res,file='00.res',status='replace' )
    !endif
    !
    write (f_res,"(a80)") (title(i), i=1,4)
    do i=1,parmax
       write (f_res,"(a10,10i3,i4,2x,f18.8)") parnam(i),ipower(0:9,i),ivar(i),param(i)
    enddo
    !
  endif  ! ---- itmax/=0
  
  !write (f_res,"('------------------------------')")
  write (f_res,"(a80)") (title(i), i=1,4)
  do i=1,parmax
     write (f_res,"(a10,10i3,i4,2x,f18.8)") parnam(i),ipower(0:9,i),ivar(i),param(i)
  enddo
  !
  do nrow=1,npts
     !
     forall(i=1:4) x(i,1:3) = cart(i,1:3,nrow)
     !
     call cartesian_to_local(x,local,rmat,tau1,tau2,tau3)
     !
     dip(:) = dipol(:,nrow)
     !
     cdip(1) = dip(1)-eps(nrow)
     cdip(2) = dip(2)-eps(nrow+npts)
     cdip(3) = dip(3)-eps(nrow+npts*2)
     !
     tauA=(2.0_rk*tau1-tau2-tau3)/sqrt(6.0_rk)
     tauB= (tau2-tau3)/sqrt(2.0_rk)
     !
     write (f_res,"(4f11.4,5f11.4,14f14.8)") & 
            local(1:4),local(5:10)*180.0_rk/pi,&
            dip(1),cdip(1),eps(nrow),        &
            dip(2),cdip(2),eps(nrow+npts),   &
            dip(3),cdip(3),eps(nrow+npts*2), &
            wt(nrow),tau1,tau2,tau3,tauA,tauB
      !
  enddo
  !
  ! Here we can stop the fit 
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
   do i=1,parmax
      write (f_res,"(a10,10i3,i4,2x,f18.8)") parnam(i),ipower(0:9,i),ivar(i),param(i)
   enddo 
  endif  ! ---- itmax/=0
  
  write (f_res,"('------------------------------')")

    
   write (f_res,"(a80)") (title(i), i=1,4)
   do i=1,parmax
      write (f_res,"(a10,10i3,i4,2x,f18.8)") parnam(i),ipower(0:9,i),ivar(i),param(i)
   enddo 
   !
   do nrow=1,npts
      !
      forall(i=1:4) x(i,1:3) = cart(i,1:3,nrow)
      !
     call cartesian_to_local(x,local,rmat,tau1,tau2,tau3)
     !
     dip(:) = dipol(:,nrow)
     !
     cdip(1) = dip(1)-eps(nrow)
     cdip(2) = dip(2)-eps(nrow+npts)
     cdip(3) = dip(3)-eps(nrow+npts*2)
     !
     tauA=(2.0_rk*tau1-tau2-tau3)/sqrt(6.0_rk)
     tauB= (tau2-tau3)/sqrt(2.0_rk)
     !
     write (f_res,"(4f11.4,5f11.4,10f14.8,2x,d13.6,3f12.4,3x,2f12.4)") & 
            local(1:4),local(5:10)*180.0_rk/pi,&
            dip(1),cdip(1),eps(nrow),        &
            dip(2),cdip(2),eps(nrow+npts),   &
            dip(3),cdip(3),eps(nrow+npts*2), &
            wt(nrow),tau1*180.0_rk/pi,tau2*180.0_rk/pi,tau3*180.0_rk/pi,tauA*180.0_rk/pi,tauB*180.0_rk/pi
       !
   enddo


!--  printing  out the fitted PES in the form of the original ab initio
!--  for refitting with the lineariyed coordinates later 

 open  (f_new,file='new.en',status='replace' )
 !
 do nrow=1,npts
    !
    forall(i=1:4) x(i,1:3) = cart(i,1:3,nrow)
    !
    call cartesian_to_local(x,local,rmat)
    !
    dip(:) = dipol(:,nrow)
    !
    cdip(1) = dip(1)-eps(nrow)
    cdip(2) = dip(2)-eps(nrow+npts)
    cdip(3) = dip(3)-eps(nrow+npts*2)
    !
    write (f_new,"(4f11.4,5f11.4,10f14.8,2x,d13.6)") & 
           local(1:4),local(5:10)*180.0_rk/pi,&
           dip(1),cdip(1),eps(nrow),        &
           dip(2),cdip(2),eps(nrow+npts),   &
           dip(3),cdip(3),eps(nrow+npts*2), &
           wt(nrow)
     !
 enddo
 !
 write (f_res,"(/3x,67('-')/'   |  param  |  points |    deviat     |     rms       |  stability  |')")
 write (f_res,"(3x,67('-'))")
 write (f_res,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |',6e12.5)") &
               numpar,nused,stadev,rms,stability,ssq3(:),rms3(:)


 write (f_res,"(/92('#'))")



 write (f_out,"(/3x,67('-')/'   |  param  |  points |    deviat     |     rms       |  stability  |')")
 write (f_out,"(3x,67('-'))")
 write (f_out,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |',6e12.5)") &
               numpar,nused,stadev,rms,stability,ssq3(:),rms3(:)

 write (f_out,"(/92('#'))")
 !
 contains 
    !
    subroutine print_res
     !
     integer :: i, nrow
     !
     rewind(f_res)
     !
     write (f_res,"(a80)") (title(i), i=1,4)
     do i=1,parmax
        write (f_res,"(a10,10i3,i4,2x,f18.8)") parnam(i),ipower(0:9,i),ivar(i),param(i)
     enddo
     !
     do nrow=1,npts
        !
        forall(i=1:4) x(i,1:3) = cart(i,1:3,nrow)
        !
        call cartesian_to_local(x,local,rmat)
        !
        dip(:) = dipol(:,nrow)
        !
        cdip(1) = dip(1)-eps(nrow)
        cdip(2) = dip(2)-eps(nrow+npts)
        cdip(3) = dip(3)-eps(nrow+npts*2)
        !
        write (f_res,"(4f11.4,5f11.4,10f14.8,2x,d13.6)") & 
               local(1:4),local(5:10)*180.0_rk/pi,&
               dip(1),cdip(1),eps(nrow),        &
               dip(2),cdip(2),eps(nrow+npts),   &
               dip(3),cdip(3),eps(nrow+npts*2), &
               wt(nrow)
         !
     enddo


 end subroutine print_res
 !

 end subroutine do_fitting

 subroutine cartesian_to_local(y,local,rmat,tau1_,tau2_,tau3_)
    !
    implicit none 
    !
    real(rk),intent(in)  :: y(4,3)
    real(rk),intent(out) :: local(10),rmat(3,3)
    real(rk),intent(out),optional :: tau1_,tau2_,tau3_
    real(rk)             :: tmat(4,3),x(4,3),cosbeta,r1,r2,r3,r4,cosa12,cosa13,cosa14,cosa23,cosa24,cosa34,alpha12,alpha13,alpha14,alpha23,alpha24,alpha34
    real(rk)             :: beta312,beta412,betae,pi,phi1,phi2,phi3,beta(3),beta314,sinbeta,tau0,beta_,tau1,tau2,tau3,p1,p2,p3,p4
    integer             :: i,ierror,ibeta
    character(len=40)   ::  txt
       !
       txt = 'cartesian_to_local'
       !
       pi = 4.0_rk * datan2(1.0_rk,1.0_rk)
       !
       x = y
       !
       ! transform XYZ system into the molecular system
       !
       r1=sqrt(sum(x(1,:)**2))
       r2=sqrt(sum(x(2,:)**2))
       r3=sqrt(sum(x(3,:)**2))
       r4=sqrt(sum(x(4,:)**2))
       !
       p1=sqrt(sum(x(1,:2)**2))
       p2=sqrt(sum(x(2,:2)**2))
       p3=sqrt(sum(x(3,:2)**2))
       p4=sqrt(sum(x(4,:2)**2))
       !
       cosa12=(sum(x(1,:)*x(2,:)))/(r1*r2)
       cosa13=(sum(x(1,:)*x(3,:)))/(r1*r3)
       cosa14=(sum(x(1,:)*x(4,:)))/(r1*r4)
       cosa23=(sum(x(2,:)*x(3,:)))/(r2*r3)
       cosa24=(sum(x(2,:)*x(4,:)))/(r2*r4)
       cosa34=(sum(x(3,:)*x(4,:)))/(r3*r4)
       !
       alpha12=acos(cosa12)
       alpha13=acos(cosa13)
       alpha14=acos(cosa14)
       alpha23=acos(cosa23)
       alpha24=acos(cosa24)
       alpha34=acos(cosa34)

       local(1) = r1
       local(2) = r2
       local(3) = r3
       local(4) = r4
       !
       local(5)=alpha12
       local(6)=alpha13
       local(7)=alpha14

       local( 8)=alpha23
       local( 9)=alpha24
       local(10)=alpha34
       !
       tmat(1, :) = x(1,:) / r1
       tmat(2, :) = x(2,:) / r2
       tmat(3, :) = x(3,:) / r3
       tmat(4, :) = x(4,:) / r4
       !
       !tmat(1, :) = x(1,:) / r1
       !tmat(2, :) = x(2,:) / p2
       !tmat(3, :) = x(3,:) / p3
       !tmat(4, :) = x(4,:) / p4
       !
       if ( trim(coordtype)=='1-2-y' ) then
         !
         rmat(3,:) = tmat(1, :)
         !
         rmat(2,:) = vector_product(rmat(3,:),tmat(2, :))
         !
         rmat(1,:) = vector_product(-rmat(2, :),rmat(3, :))
         !
         forall(i=1:3) rmat(i,:) = rmat(i,:)/sqrt(sum(rmat(i,:)**2))
         !
       elseif ( trim(coordtype)=='abinitio' ) then
         !
         rmat = 0
         !
         rmat(1,2) = 1.0_rk
         rmat(2,1) = 1.0_rk
         rmat(3,3) = 1.0_rk
         !
       elseif (trim(coordtype)=='sym') then
         !
         !tmat(2:4,3) = 0
         !
         rmat(1,:) = 1._rk/sqrt(6._rk)*(2.0_rk*tmat(2,:)-tmat(3,:)-tmat(4,:))
         rmat(2,:) = 1._rk/sqrt(2._rk)*(                 tmat(3,:)-tmat(4,:))
         rmat(3,:) = tmat(1, :)
         !
         !rmat(1,:) = 1._rk/sqrt(6._rk)*(2.0_rk*tmat(2,:)-tmat(4,:)-tmat(3,:))
         !rmat(2,:) = 1._rk/sqrt(2._rk)*(                 tmat(4,:)-tmat(3,:))
         rmat(3,:) = tmat(1, :)
         !
         !forall(i=1:3) rmat(i,:) = rmat(i,:)/sqrt(sum(rmat(i,:)**2))
         !
       elseif (trim(coordtype)=='sym-x') then
         !
         !tmat(2:4,3) = 0
         !
         rmat(3,:) = tmat(1, :)
         rmat(1,:) = 1._rk/sqrt(6._rk)*(2.0_rk*tmat(2,:)-tmat(3,:)-tmat(4,:))
         !
         rmat(2,:) = vector_product(rmat(3,:),rmat(1, :))
         !
         !rmat(2,:) = 1._rk/sqrt(2._rk)*(                 tmat(3,:)-tmat(4,:))
         !rmat(3,:) = tmat(1, :)
         !
         !rmat(1,:) = 1._rk/sqrt(6._rk)*(2.0_rk*tmat(2,:)-tmat(4,:)-tmat(3,:))
         !rmat(2,:) = 1._rk/sqrt(2._rk)*(                 tmat(4,:)-tmat(3,:))
         !rmat(3,:) = tmat(1, :)
         !
         forall(i=1:3) rmat(i,:) = rmat(i,:)/sqrt(sum(rmat(i,:)**2))
         !
       else
         !
         stop 'illegal coordtype'
         !
       endif
       !
       tmat(1:3,1:3) = matmul(transpose(rmat),rmat)
       !
       cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
       beta312 = aacos(cosbeta,ierror,txt)
       !
       cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
       beta412 = aacos(cosbeta,ierror,txt)
       !
       cosbeta = (cos(alpha34)-cos(alpha13)*cos(alpha14) )/(sin(alpha13)*sin(alpha14))
       beta314 = aacos(cosbeta,ierror,txt)
       !
       cosbeta = cos(alpha34)-cos(alpha13)*cos(alpha14)
       sinbeta = sin(alpha13)*sin(alpha14)
       !
       tau0 = 2.0_rk*pi-(beta312+beta412)
       !
       beta(1) = beta314
       beta(2) = beta312
       beta(3) = beta412
       !
       if (abs(sum(beta)-2.0*pi)>1e-4) then 
         !
         ibeta = maxloc(beta,dim=1)
         beta_ = beta(ibeta)
         !
         tau0 = 2.0_rk*pi-(sum(beta)-beta_)
         !
         if (tau0>pi) beta(ibeta) = 2.0_rk*pi-beta_
         !
       endif 
       !
       tau1 = beta(1)  !; if (tau0>pi) tau1 = 2.0_rk*pi-beta314 
       tau2 = beta(2)
       tau3 = beta(3)
       !
       if (present(tau1_).and.present(tau2_).and.present(tau3_)) then
         tau1_ = tau1
         tau2_ = tau2
         tau3_ = tau3
       endif


       !
 end subroutine cartesian_to_local


   function vector_product(v1,v2) result (v)
     !
     real(8),intent(in) :: v1(3),v2(3)
     real(8) :: v(3)
     !
     v(1) = v1(2)*v2(3)-v1(3)*v2(2)
     !
     v(2) = v1(3)*v2(1)-v1(1)*v2(3)
     !
     v(3) = v1(1)*v2(2)-v1(2)*v2(1)
     !
   end function vector_product






 function aacos(x,ierror,txt) result (f)
    !
    implicit none 
    !
    real(rk),intent(in)   :: x
    character(len=40),intent(in),optional :: txt
    integer,intent(out)    :: ierror
    real(rk)              :: f,small_,pi
    !
    small_ = epsilon(1.0_rk)
    pi = 4.0_rk * atan2(1.0_rk,1.0_rk)
    !
    ierror = 0
    !
    if ( abs(x)>1.0_rk+10.0*sqrt(small_) ) then 
       !
       write (6,"('|cos(x)|>1: ',f18.8,2x,a)") x,txt
       !
       ierror = 1
       !
       !stop 'aacos - bad cosalpha'
       !
    elseif ( x>=1.0_rk) then 
       f = 0.0_rk
    elseif ( x<=-1.0_rk) then 
       f = pi
    else 
       f = acos(x)
    endif
    !
 end function aacos


  function diff_mu(ipower,local,param)
  !
  implicit none
  !
  integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
  integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
  !
  integer,intent(in)  :: ipower(0:9)
  real(rk),intent(in)  :: local(10),param(nparam_eq)
  real(rk)             :: s1,s2,s3,s4 
  real(rk)             :: v,r1,r2,r3,r4,alpha12,alpha13,alpha14,alpha23,alpha24,alpha34,deg,r0,r12,r13,r14,r23,r24,r34,t12,t13,t14,t23,t24,t34
  real(rk)             :: cosa12,cosa13,cosa14,cosa23,cosa24,cosa34
  real(rk)             :: beta312,beta412,betae,pi,phi1,phi2,phi3
  real(rk)             :: r1e,r2e,a,y1,y2,y3,y4,y5,y6,y7,y8,y9,xi(9),chi(9),b,tau1,tau2,tau3,cosbeta,chi0(9)
  real(rk)             :: mu(3),diff_mu(3),term(3),nu(3),beta314,tau0,alpha34_,sinbeta,beta(3),beta_
  integer              ::  ierror,ioper,imu,ibeta
  character(len=40)    ::  txt
  !
  txt = 'diff_mu'
  !
  pi = 4.0_rk * atan2(1.0_rk,1.0_rk)
  !
  deg=pi/180.0_rk
  !
  r1      = local(1)
  r2      = local(2)
  r3      = local(3)
  r4      = local(4)
  !
  alpha12 =local( 5)
  alpha13 =local( 6)
  alpha14 =local( 7)
  !
  alpha23 =local( 8)
  alpha24 =local( 9)
  alpha34 =local(10)
  !
  alpha34_ =  calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
  !
  cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
  beta312 = aacos(cosbeta,ierror,txt)
  !
  cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
  beta412 = aacos(cosbeta,ierror,txt)
  !
  cosbeta = (cos(alpha34)-cos(alpha13)*cos(alpha14) )/(sin(alpha13)*sin(alpha14))
  beta314 = aacos(cosbeta,ierror,txt)
  !
  cosbeta = cos(alpha34)-cos(alpha13)*cos(alpha14)
  sinbeta = sin(alpha13)*sin(alpha14)
  !
  tau0 = 2.0_rk*pi-(beta312+beta412)
  !
  beta(1) = beta314
  beta(2) = beta412
  beta(3) = beta312
  !
  beta(1) = beta314
  beta(2) = beta312
  beta(3) = beta412
  !
  if (abs(sum(beta)-2.0*pi)>1e-4) then 
    !
    ibeta = maxloc(beta,dim=1)
    beta_ = beta(ibeta)
    !
    tau0 = 2.0_rk*pi-(sum(beta)-beta_)
    !
    if (tau0>pi) beta(ibeta) = 2.0_rk*pi-beta_
    !
  endif 
  !
  !tau1 = beta(1)  !; if (tau0>pi) tau1 = 2.0_rk*pi-beta314 
  !tau3 = beta(2)
  !tau2 = beta(3)
  !
  tau1 = beta(1)  !; if (tau0>pi) tau1 = 2.0_rk*pi-beta314 
  tau2 = beta(2)
  tau3 = beta(3)
  !
  !
  if (abs(tau1+tau2+tau3-2.0*pi)>1e-4 )  then
    !
    write(6,"('tau1+tau2+tau3-pi)>1e-4',4f12.4)") tau1,tau2,tau3,tau1+tau2+tau3
    stop 'tau1+tau2+tau3-pi)>1e-4'
    !
  endif 
  !
  !tau1 = cos(tau1)  !; if (tau0>pi) tau1 = 2.0_rk*pi-beta314 
  !tau2 = cos(tau2)
  !tau3 = cos(tau3)
  !
  r1e      = param(1)
  r2e      = param(2)
  betae    = param(3)*deg
  !
  xi(1)=(r1-r1e)
  xi(2)=(r2-r2e)
  xi(3)=(r3-r2e)
  xi(4)=(r4-r2e)
  !
  xi(5) = alpha12-betae
  xi(6) = alpha13-betae
  xi(7) = alpha14-betae
  !
  !xi(5) = cos(alpha12)-cos(betae)
  !xi(6) = cos(alpha13)-cos(betae)
  
  !xi(7) = cos(alpha14)-cos(betae)
  !
  xi(8)=(2.0_rk*tau1-tau2-tau3)/sqrt(6.0_rk)
  xi(9)= (tau2-tau3)/sqrt(2.0_rk)
  !
  !xi(8)=(2.0_rk*cos(tau1)-cos(tau2)-cos(tau3))/sqrt(6.0_rk)
  !xi(9)= (cos(tau2)-cos(tau3))/sqrt(2.0_rk)
  !
  imu = ipower(0)
  !
  mu = 0
  term = 0 
  !
  do ioper =1,6
    !
    !call ML_symmetry_transformation_XY3(ioper,xi,chi0)
    !
    call ML_symmetry_transformation_XY3_II(ioper,xi,chi)
    !
    !if (any( (chi(:)-chi0(:))**2>1e-12 ) ) then
    !  continue
    !endif 
    !
    mu(imu) = product(chi(:)**ipower(1:))
    !
    call ML_symmetry_transformation_E(ioper,mu,nu)
    !
    term(:) = term(:) + nu(:)
    !
  enddo
  !
  diff_mu(:) = term(:)
  !
end function diff_mu

 !
 recursive subroutine dipole(parmax,dmu,mu,param,ipar)
 !
 implicit none 
 !
 integer,intent(in)           ::  parmax
 real(rk),intent(in)  ::  dmu(3,parmax),param(parmax)
 real(rk),intent(out) ::  mu(3)
 integer,intent(in) ::           ipar
 
 !real(rk)         ::  param_(parmax),local(10)
 !
 !real(rk)         ::  deltax,potright(3),potleft(3)
 !real(rk)         ::  dF_(nmax,3),rmat(3,3),dmu(3)
 integer                  :: i,k
 !character(len=40)        ::  txt
 !
 ! parameter for the finite differncies differention 
 !
 real(rk),parameter :: factordeltax=0.001
  !
  !call cartesian_to_local(x,local,rmat)
  !
  mu = 0
  !
  if (ipar==0) then
    !
    do i = nparam_eq+1,parmax
      !
      mu(:) = mu(:) + param(i)*dmu(:,i)
      !
    enddo
    !
  elseif (ipar>nparam_eq) then
    !
    mu(:) = dmu(:,ipar)
    !
  else
    !
    mu(:) = 0
    !
  endif
  !
 end subroutine dipole
 !
  !
  function find_alpha34(alpha12,alpha13,alpha23,alpha14,alpha24,ierror) result (s10_f)
    !
    implicit none
    !
    real(8),intent(in)  ::alpha12,alpha13,alpha14,alpha23,alpha24
    integer,intent(out) :: ierror
    real(8) s10_f
    real(8) xtmp(4),ytmp(4),ztmp(4),temp

     ierror = 0
     !
     s10_f = 0
     !
     xtmp(1)=1.0_rk
     ytmp(1)=0
     ztmp(1)=0
     xtmp(2)=cos(alpha12)
     ytmp(2)=sin(alpha12)
     ztmp(2)=0


     xtmp(3)=cos(alpha13)
     ytmp(3)=(cos(alpha23)-xtmp(2)*xtmp(3))/ytmp(2)
     !
     temp = 1.0_rk-xtmp(3)**2-ytmp(3)**2
     !
     if (temp<0) then 
       !
       ierror = 1 ; return 
       !
     endif
     !
     ztmp(3)=sqrt(temp)
     !
     xtmp(4)=cos(alpha14)
     ytmp(4)=(cos(alpha24)-xtmp(2)*xtmp(4))/ytmp(2)
     !
     temp = 1.0_rk-xtmp(4)**2-ytmp(4)**2
     !
     if (temp<0) then 
       !
       ierror = 1 ; return 
       !
     endif
     !
     ztmp(4)=-sqrt(temp)
     !
     temp  = xtmp(3)*xtmp(4)+ytmp(3)*ytmp(4)+ztmp(3)*ztmp(4)
     !
     if (abs(temp)>1._rk) then 
       !
       ierror = 1 ; return 
       !
     endif

     s10_f=acos(temp)

 end function find_alpha34


 function calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24) result (alpha34)
   !
   implicit none
   !
   real(8),intent(in)  ::alpha12,alpha13,alpha14,alpha23,alpha24
   real(8) :: alpha34,ca12,ca13,ca23,ca14,ca24
   real(8) :: cosalpha,small_

       small_ = epsilon(1.0_rk)
       ca12 = cos(alpha12)
       ca13 = cos(alpha13)
       ca14 = cos(alpha14)
       ca23 = cos(alpha23)
       ca24 = cos(alpha24)

       cosalpha = 0.5_rk*( &
                  2.0_rk*ca12*ca23*ca14-     &
                  2.0_rk*ca13*ca14+     &
                  2.0_rk*ca12*ca24*ca13-     &
                  2.0_rk*ca23*ca24+     &
                  2.0_rk*sqrt(     &
                  1.0_rk-ca24**2+     &
                  2.0_rk*ca12*ca24*ca14+     &
                  2.0_rk*ca12*ca23*ca13-     &
                  2.0_rk*ca12**3*ca23*ca13-     &
                  2.0_rk*ca12**2-     &
                  ca14**2-     &
                  ca13**2+     &
                  ca14**2*ca23**2+     &
                  ca13**2*ca24**2+     &
                  4.0_rk*ca12**2*ca23*ca14*ca24*ca13+     &
                  ca12**4-     &
                  2.0_rk*ca12*ca23*ca14**2*ca13-     &
                  2.0_rk*ca12*ca23**2*ca14*ca24-     &
                  2.0_rk*ca13**2*ca14*ca12*ca24-     &
                  2.0_rk*ca12*ca24**2*ca13*ca23+     &
                  ca13**2*ca14**2+     &
                  ca23**2*ca24**2+     &
                  ca12**2*ca13**2+     &
                  ca12**2*ca14**2+     &
                  ca12**2*ca23**2+     &
                  ca12**2*ca24**2-     &
                  2.0_rk*ca12**3*ca24*ca14-     &
                  ca23**2))/( -sin(alpha12)**2 )

       if ( abs(cosalpha)>1.0_rk+sqrt(small_) ) then 
          !
          write (6,"('calc_alpha34: cosalpha>1: ',f18.8)") cosalpha
          stop 'calc_alpha34 - bad cosalpha'
          !
       elseif ( cosalpha>=1.0_rk) then 
          alpha34 = 0.0_rk
       else 
          alpha34 = acos(cosalpha)
       endif
       !

 end function calc_alpha34


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


 subroutine linur(dimen,npar,coeff,constant,solution,error)
   !
   implicit none

  integer,intent(in)  :: dimen,npar
  integer,intent(out) :: error 
  real(rk),intent(in)  :: coeff(npar,npar),constant(npar)
  real(rk),intent(out) :: solution(npar)
  real(rk)          :: a0(npar,npar)
  real(rk)          :: c
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
  subroutine invmat(al,ai,dimen,npar,info)
  implicit none
  integer,intent(in)           :: npar,dimen
  real(rk),intent(in)  :: al(npar,npar)
  real(rk),intent(out) :: ai(npar,npar)
  integer(ik),intent(out) :: info
  real(rk)             :: h(npar),p,q,small_
  integer                      :: i1,i2,k,i,j,k8,k9
  
    !
    small_ = epsilon(1.0_rk)  
    !
    info = 0
    !

    ai(1:dimen,1:dimen)=al(1:dimen,1:dimen)
 
   do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      !
      do i=2,dimen
        q=ai(i,1)
        !
        if (abs(q)<100.0*small_) then 
          h(i) = 0
        elseif (abs(p)<small_) then
          info = i-1
          return
        else 
          h(i)=q/p
        endif
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
  subroutine ML_symmetry_transformation_XY3(ioper,src,dst)
    implicit none 
    !
    integer,intent(in)    :: ioper  ! group operation  
    real(rk),intent(in)      :: src(1:9)
    real(rk),intent(out)     :: dst(1:9)
    !
    real(rk)         :: repres(12,2,2),a,b,e,o
    !
    if (verbose>=5) write(f_out,"('ML_symmetry_transformation_XY3/start')") 
    !
    dst(1) = src(1)
    !
    select case(ioper)
    !
    case (1) ! identity 

      dst = src

    case (3) ! (123)

      dst(2) = src(3)
      dst(3) = src(4)
      dst(4) = src(2)
      !
      dst(5) = src(6)
      dst(6) = src(7)
      dst(7) = src(5)

    case (2) ! (321)

      dst(2) = src(4)
      dst(3) = src(2)
      dst(4) = src(3)
      !
      dst(5) = src(7)
      dst(6) = src(5)
      dst(7) = src(6)
 
    case (6) ! (12)

      dst(2) = src(3)
      dst(3) = src(2)
      dst(4) = src(4)
      !
      dst(5) = src(6)
      dst(6) = src(5)
      dst(7) = src(7)

    case (5) ! (13)

      dst(2) = src(4)
      dst(3) = src(3)
      dst(4) = src(2)
      !
      dst(5) = src(7)
      dst(6) = src(6)
      dst(7) = src(5)

    case (4) ! (23)

      dst(2) = src(2)
      dst(3) = src(4)
      dst(4) = src(3)
      !
      dst(5) = src(5)
      dst(6) = src(7)
      dst(7) = src(6)

    case default

      write (f_out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
      stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
    end select 
    !
    a = 0.5_rk ; b = 0.5_rk*sqrt(3.0_rk) ; e = 1.0_rk ; o = 0.0_rk
    !
    repres ( 1,:,:)= reshape((/ e, o,  & 
                                o, e/),(/2,2/))
    !
    repres ( 3,:,:)= reshape((/-a,-b,  &
                                b,-a/),(/2,2/))
    !
    repres ( 2,:,:)= reshape((/-a, b,  &
                               -b,-a/),(/2,2/))
    !
    repres ( 4,:,:)= reshape((/ e, o,  &
                                o,-e/),(/2,2/))
    !
    repres ( 6,:,:)= reshape((/-a, b,  &
                                b, a/),(/2,2/))
    !
    repres ( 5,:,:)= reshape((/-a,-b,  &
                               -b, a/),(/2,2/))
    !
    dst(8) = repres(ioper,1,1)*src(8)+repres(ioper,1,2)*src(9)
    dst(9) = repres(ioper,2,1)*src(8)+repres(ioper,2,2)*src(9)
           
  end subroutine ML_symmetry_transformation_XY3        
  

  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  subroutine ML_symmetry_transformation_XY3_II(ioper,src,dst)
    implicit none 
    !
    integer,intent(in)    :: ioper  ! group operation  
    real(rk),intent(in)      :: src(1:9)
    real(rk),intent(out)     :: dst(1:9)
    !
    real(rk)         :: repres(6,9,9),a,b,e,o
    !
    if (verbose>=5) write(f_out,"('ML_symmetry_transformation_XY3/start')") 
    !
    repres = 0 
    a = 0.5_rk ; b = 0.5_rk*sqrt(3.0_rk) ; e = 1.0_rk ; o = 0.0_rk
    !
    repres(1,1,1) = 1.0_rk
    repres(1,2,2) = 1.0_rk
    repres(1,3,3) = 1.0_rk
    repres(1,4,4) = 1.0_rk
    repres(1,5,5) = 1.0_rk
    repres(1,6,6) = 1.0_rk
    repres(1,7,7) = 1.0_rk
    repres(1,8,8) = 1.0_rk
    repres(1,9,9) = 1.0_rk
    !
    ! 321
    !
    repres(2,1,1) = 1.0_rk
    repres(2,2,4) = 1.0_rk
    repres(2,3,2) = 1.0_rk
    repres(2,4,3) = 1.0_rk
    repres(2,5,7) = 1.0_rk
    repres(2,6,5) = 1.0_rk
    repres(2,7,6) = 1.0_rk
    repres(2,8,8) = -a
    repres(2,8,9) =  b
    repres(2,9,8) = -b
    repres(2,9,9) = -a
    !
    ! 123
    !
    repres(3,1,1) = 1.0_rk
    repres(3,2,3) = 1.0_rk
    repres(3,3,4) = 1.0_rk
    repres(3,4,2) = 1.0_rk
    repres(3,5,6) = 1.0_rk
    repres(3,6,7) = 1.0_rk
    repres(3,7,5) = 1.0_rk
    repres(3,8,8) = -a
    repres(3,8,9) = -b
    repres(3,9,8) =  b
    repres(3,9,9) = -a
    !
    ! 23
    !
    repres(4,1,1) = 1.0_rk
    repres(4,2,2) = 1.0_rk
    repres(4,3,4) = 1.0_rk
    repres(4,4,3) = 1.0_rk
    repres(4,5,5) = 1.0_rk
    repres(4,6,7) = 1.0_rk
    repres(4,7,6) = 1.0_rk
    repres(4,8,8) = 1.0_rk
    repres(4,8,9) = 0.0_rk
    repres(4,9,8) = 0.0_rk
    repres(4,9,9) = -1.0_rk
    !
    ! 13
    repres(5,1,1) = 1.0_rk
    repres(5,2,4) = 1.0_rk
    repres(5,3,3) = 1.0_rk
    repres(5,4,2) = 1.0_rk
    repres(5,5,7) = 1.0_rk
    repres(5,6,6) = 1.0_rk
    repres(5,7,5) = 1.0_rk
    repres(5,8,8) = -a
    repres(5,8,9) =  b
    repres(5,9,8) =  b
    repres(5,9,9) =  a
    !
    ! 12
    repres(6,1,1) = 1.0_rk
    repres(6,2,3) = 1.0_rk
    repres(6,3,2) = 1.0_rk
    repres(6,4,4) = 1.0_rk
    repres(6,5,6) = 1.0_rk
    repres(6,6,5) = 1.0_rk
    repres(6,7,7) = 1.0_rk
    repres(6,8,8) = -a
    repres(6,8,9) = -b
    repres(6,9,8) = -b
    repres(6,9,9) = a
    !
    if (ioper<0.or.ioper>6) then
      write (6,"('symmetry_transformation_local: operation ',i8,' unknown')") ioper
      stop 'symmetry_transformation_local - bad operation. type'
    endif
    !
    !repres(ioper,8:9,8:9) = transpose(repres(ioper,8:9,8:9))
    !
    dst = matmul(repres(ioper,:,:),src) 
    !
    !dst = matmul(transpose(repres(ioper,:,:)),src) 
           
  end subroutine ML_symmetry_transformation_XY3_II      
  
  
  subroutine ML_symmetry_transformation_E(ioper,src,dst)
    implicit none 
    !
    integer,intent(in)    :: ioper  ! group operation  
    real(rk),intent(in)      :: src(1:3)
    real(rk),intent(out)     :: dst(1:3)
    !
    real(rk)         :: repres(6,3,3),a,b,e,o
    !
    if (verbose>=5) write(f_out,"('ML_symmetry_transformation_XY3/start')")
    !
    a = 0.5_rk ; b = 0.5_rk*sqrt(3.0_rk) ; e = 1.0_rk ; o = 0.0_rk
    !
    repres = 0
    !
    repres(:,3,3) = 1.0_rk 
    !
    repres(1,1,1) = 1.0_rk
    repres(1,2,2) = 1.0_rk
    !
    ! 321
    repres(2,1,1) = -a
    repres(2,1,2) = -b
    repres(2,2,1) =  b
    repres(2,2,2) = -a
    !
    ! 123
    repres(3,1,1) = -a
    repres(3,1,2) =  b
    repres(3,2,1) = -b
    repres(3,2,2) = -a
    !
    ! 23
    repres(4,1,1) = 1.0_rk
    repres(4,1,2) = 0.0_rk
    repres(4,2,1) = 0.0_rk
    repres(4,2,2) = -1.0_rk
    !
    ! 13
    repres(5,1,1) = -a
    repres(5,1,2) = -b
    repres(5,2,1) = -b
    repres(5,2,2) =  a
    !
    ! 12
    repres(6,1,1) = -a
    repres(6,1,2) = b
    repres(6,2,1) = b
    repres(6,2,2) = a
    !
    if (ioper<0.or.ioper>6) then
      write (6,"('symmetry_transformation_local: operation ',i8,' unknown')") ioper
      stop 'symmetry_transformation_local - bad operation. type'
    endif
    !
    !dst = matmul(repres(ioper,:,:),src) 
    !
    dst = matmul(transpose(repres(ioper,:,:)),src) 
           
  end subroutine ML_symmetry_transformation_E       
  


  
  end module FITDMS_ch3cl
  
  
  
  
  program driver
  use FITDMS_ch3cl
  
  call do_fitting
  
  end program driver
  
  