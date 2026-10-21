 module FITPES_ch3cl
  !
  public ik,rk,do_fitting

  private
 

 integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
 integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
 integer,parameter          :: f_inp=5, f_out=6, f_res=9, f_new =10, nparam_eq = 5

  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains
  
 subroutine do_fitting
 !
 ! Parameters:  
 ! 
 ! number of parameters and maximum number of ab initio points   
 !
 !
 integer,parameter          ::  enermax =  150000, nmax0=8,taumax0=8, nparam = 400
 !
 ! where the output goes to 
 ! 
 integer,parameter          :: f_inp=5, f_out=6, f_res=9, f_new =10, nparam_eq = 5
 ! 
 ! Dinite difference differentiation with 2 or 3 points
 !
 integer,parameter          ::  findif = 3 
 !
 ! the iteration will be stoped when the best standard deviation reached 
 !
 !double precision,parameter :: stadev_best = 0.001d0   
 !
 integer,parameter          ::  verbose = 5
 ! 
 ! parameter to control the correlation when we through out the parameters with large st. error (0.0..1.0)
 !
 double precision,parameter :: fact_corr=0.4, too_big = 1e9
 ! 
 ! Watson's alpha parameter to control the Robust fit
 double precision,parameter :: a_wats = 0.1
 ! parameter for the finite differncies differention 
 !
 double precision,parameter :: factordeltax=0.001
 !
 integer,parameter          ::  ireplace_by = 0
 !
 ! Universal constants 
 ! 
 double precision,parameter :: bohr = 0.52917715e+00,hartre=2.194746354d+05
 ! Avareging region for the robustfit = 3 
 integer,parameter          :: ishift = 10
  character(len=40),parameter  :: fit_type= 'linur'  !'dgelss'
 ! Variables:  
 double precision  :: maxdeviat  ! st. deviation at the previos iteration 
 !
 integer           :: i,npts,nused,itmax,numpar,ncol,iter,nrow,l,itertypemax,irow,icolumn ! just integers  
 ! 
 character(len=80) :: title(4),longlabel ! input title of the job from the first four input lines 
 character(len=1)  :: mark           ! used at the parameters output to mark the dependent or independent ones 
 character(len=10) :: label          ! Temporary label 
 double precision  :: wtsum,ssq      ! 
 !
 ! some variables 
 !
 double precision  :: stadev,epsil,stadev_old,v,dd,pi
 double precision  :: tempx,deltax,potright,potleft,stability
 double precision  :: r14,r24,r34,alpha1,alpha2,alpha3
 double precision  :: corr,ve
 integer           :: ierror,n0,i1,i2,i3,i4,i5,i6,itertype,last_i,robustfit
 double precision  :: stab_best,acc_rms,conf_int
 integer           :: maxivar,minivar,tivar
 double precision  :: last_stadev,factor_back,factor_out
 double precision  :: rms,s,da1,da2,da,ZPE,local(18),stadev_best
 integer           :: last_ivar,ndigits,npoints_delta,rand_itmin,alloc
 integer           :: max_pnts,ipar,i0,l0,parmax,it(18),fivar_t,nmax,taumax
 logical           :: yes,still_run,ifopen,find_best_params
 integer           :: flag_borh, flag_hartre,ivartmp,it0,imin,imax,ilarge_sterr,ilow_sterr,Nout,parmax_t,i_t
 double precision  :: conv_length,conv_energy,paramtmp,relat_sterr,relat_sterr_min,fparam_t,beta312,alpha34
 character(len=68)  :: fmt1,fmt2
 character(len=2)   :: fmt0
 character(len=80)  ::  char_tmp 
 integer           :: ipower_(18),iparam
 !
 ! some matrices, to be allocated  
 !
 double precision, allocatable :: rjacob(:,:),last_param(:),param(:)
 double precision, allocatable :: param_tmp(:)
 double precision, allocatable :: crr(:),coord(:,:)
 double precision, allocatable :: energy(:)
 double precision, allocatable :: eps(:),wt(:),parold(:),sigma(:),wt_tmp(:),parinitial(:)
 double precision, allocatable :: al(:,:),bl(:),dx(:)
 double precision, allocatable :: ai(:,:),sterr(:),dV(:,:)
 double precision, allocatable :: sens(:)
 character(len=10), allocatable:: parnam(:),parnam_tmp(:)
 integer, allocatable :: ivar(:),icorr(:),ivar0(:),ivar_tmp(:),ind_tmp(:),ipower_tmp(:,:),ipower(:,:)
 !
 double precision, allocatable :: Tsing(:,:) 
 integer                       :: rank0,info
 integer                       :: lwork
 double precision,allocatable  :: wspace(:)
 double precision              :: alpha13,alpha12,alpha14,alpha23,alpha24,cosbeta,aacos,beta412,cosa34
 character(len=40)        ::  txt
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
 read (f_inp,*) itmax
 read (f_inp,*) stab_best,stadev_best
 read (f_inp,*) conv_length
 read (f_inp,*) conv_energy
 read (f_inp,*) factor_out
 read (f_inp,*) factor_back
 read (f_inp,*) robustfit
 ! 
 write (f_out,"(i10,' <= number of iterations in the fit')")         itmax
 write (f_out,"(f10.6,' <= The fit is over if the stability is reached ')") stab_best
 write (f_out,"(f10.6,' <= The fit is over if the st.dev. best is reached ')") stadev_best
 write (f_out,"(f10.6,' <= sterr parameter out ')")                  factor_out
 write (f_out,"(f10.6,' <= stdev parameter out ')")                  factor_back
 write (f_out,"(i10,  ' <= watson robust on (1) or of (0) ')")       robustfit
 !
 if (robustfit<0) then 
   find_best_params = .false.
 else
   find_best_params = .true.
 endif 
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
 allocate (param_tmp(parmax_t),parnam_tmp(parmax_t),ivar_tmp(parmax_t),ipower_tmp(1:9,parmax_t),stat=alloc)
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
       read (char_tmp(1:80),*) parnam_tmp(i),ipower_(1:9),fivar_t,fparam_t
       !
       write (f_out,"(a2,9i3,' =',f18.8)") trim(parnam_tmp(i)),ipower_(1:9),fparam_t

       param_tmp(i) = fparam_t
       ivar_tmp(i) = fivar_t
       ipower_tmp(1:9,i) = ipower_(1:9)
       !
    else
       Yes = .False.
    endif
 enddo read_loop
 !
 parmax = i 
 !
 allocate (param(parmax),parnam(parmax),ivar(parmax),icorr(parmax),ivar0(parmax),ipower(1:9,parmax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('param - out of memory')")  
    stop 'param - out of memory'
 endif 
 !
 param(1:parmax) = param_tmp(1:parmax)
 ivar(1:parmax) = ivar_tmp(1:parmax)
 ipower(1:9,1:parmax) = ipower_tmp(1:9,1:parmax)
 parnam(1:parmax) = parnam_tmp(1:parmax)
 !
 deallocate (param_tmp,parnam_tmp,ipower_tmp,ivar_tmp)
 !
 allocate (parold(parmax),parinitial(parmax),last_param(parmax),sterr(parmax),&
            Tsing(parmax,parmax),Wspace(lwork),crr(parmax),al(parmax,parmax),bl(parmax),&
            dx(parmax),ai(parmax,parmax),rjacob(enermax,parmax),coord(10,enermax),energy(enermax),&
            eps(enermax),wt(enermax),wt_tmp(enermax),sigma(enermax),sens(parmax),&
            dV(enermax,nparam),stat=alloc)
 if (alloc/=0) then 
  write(6,"('sens - out of memory')")  
  stop
 endif 

 parinitial = param 


!---input energy values------!
 npts = 0
 yes = .true. 
 Nout = 1
 do while( yes )
    npts=npts+1
    if (npts>enermax) then
      write(f_out,"('Too many ab initio points, increase enermax:',I6,' vs ',I6)") npts,enermax
      endif
    read (f_inp,*) coord(1:9,npts),energy(npts),wt(npts)
    if ( wt(npts).lt.0.0 ) wt(npts) = abs(wt(npts))*factor_back
    if ( coord(1,npts).lt.0.0 ) yes = .false.
 enddo 
 !
 ! Number of ab initio points  
 npts=npts-1
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
 ! We introduce a "zero point energy" shift to avoid loosin accuracy 
 ! because of too too large numbers
 ZPE = 0 !energy(1)
 ! Convert the energies to the internal dimension, e.g. 1/cm 
 energy(:) = ( energy(:)-ZPE )*conv_energy
 !
 ! For the Robust fit we need to dedine accuracy of the abinitio data "sigma"
 ! 

  do nrow=1,npts
    if (wt(nrow)>0.0) then 
      sigma(nrow) = 1.0d0/wt(nrow)
    else 
      sigma(nrow) = 1000.0 
    endif 
  enddo

 ! a_wats = 1.0d0

 !
 ! Sometimes pot. parameters are dependent. They are nor to be in the fit 
 ! If we meet a dependent parameter, we get rid of it and start the fit again 
 ! This is the outer loop for
 !

 still_run = .true.
 outer_loop: do while (still_run)  
   !
   !  make a table header 
   !
   if (itmax.ne.0) then 
     write(f_out,"(/3x,66('-'))") 
     write(f_out,"('   |   iter  |  points |    deviat     |     rms       |  stability  |')")
     write(f_out,"(3x,66('-'))") 
   endif 
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
       write(f_out,"('worning! the number of parameter <  the number of data points',2I6)") nused,numpar
       stop 'worning! the number of parameter <  the number of data points'
    endif
   if (nused==numpar) then 
      write(f_out,"('worning! the same number of parameters and data points')")
   endif 
   !
   !
   !##################################################c
   !############ start of the fitting ################c
   !##################################################c
   rjacob = 0 
   iter = 0
   sens = 0.d0 
   do while( iter<=itmax .and. stadev>stadev_best .and. (  stability.ge.stab_best) )   
     iter = iter + 1
     ssq=0.0d+00
     rms=0.0d+00
     !
     parold = param
     !
     dV = 0 
     !
     do nrow=1,npts
       !
       if (verbose>=6) write(6,"('n = ',i9)") nrow
       !
       ! Running bond lenght 
       !
       local(1:4)=coord(1:4,nrow)*conv_length
       !
       ! Running interbond angles
       !
       local(5:9)=coord(5:9,nrow)*pi/180.0d0
       !
       do iparam = nparam_eq+1,parmax
         !
         dV(nrow,iparam) = diff_V(ipower(1:9,iparam),local,param(1:nparam_eq))
         !
       enddo
       !
     enddo
     !
     do nrow=1,npts
       !
       ! Running bond lenght 
       !
       local(1:4)=coord(1:4,nrow)*conv_length
       !
       ! Running interbond angles
       !
       local(5:9)=coord(5:9,nrow)*pi/180.0d0
       ! 
       ! Value of the potential energy fucntion at the current geometry 
       !
       v =  potential(parmax,dV(nrow,:),local,param,ipower,0)
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       if (itmax/=0 .and. wt(nrow)>0.0) then
         ncol=0
         do  i=1,parmax
           if (ivar(i) .ne. 0) then
             !
             ncol=ncol+1
             rjacob(nrow,ncol)=potential(parmax,dV(nrow,:),local,param,ipower,i)
             !
           endif
         enddo ! --- ncol
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
        do irow=1,numpar       !==== row-...... ====!
          do icolumn=1,irow    !==== column-....====!
            al(irow,icolumn)=sum(rjacob(:,icolumn)*rjacob(:,irow)*wt(:))
            al(icolumn,irow)=al(irow,icolumn)
          enddo
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
          if (ierror.ne.0) then
            ncol=0
            do i=1,parmax
              if (ivar(i) .ne. 0) then
                ncol=ncol+1
                if ( ncol.eq.ierror ) then
                    ivar(i) = 0
                    param(i) = 0 
                    if (ireplace_by/=0) param(i) = parinitial(i)
                    write(f_out,"(i,'-th is out - ',a7)") i,parnam(i)
                endif
              endif
            enddo
            cycle outer_loop
          endif
       case ('dgelss')
         !
         write(6,"('DGELSS: lapack was not presented ')")
         stop 'dgelss not compiled'
         !
         ai = al
         !
         !
         !call DGELSS(numpar,numpar,1,AL,numpar,BL,numpar,Tsing,1.D-12,RANK0,wspace,lwork,ierror)
         !
         if (info/=0) then
           write(6,"('DGELSS:error',i)") ierror
           stop 'dgelss'
         endif
         !
         dx = bl
         !
       end select
       !
       !
       if ( nused.ne.numpar ) then 
          stadev=dsqrt(ssq/float(nused-numpar))
       else 
          stadev=dsqrt(ssq/nused)
       endif
       !
        ncol=0
        do i=1,parmax
          if (ivar(i) /= 0) then
             ncol=ncol+1
             param(i)=param(i)+dx(ncol)
          endif
        enddo
        !
        ! Calcualte the inverse to A matrix, needed for the standard error of the fitted pot. parameters
        !
        call invmat(al,ai,numpar,parmax)
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
     wtsum=sum(wt_tmp(:))
     wtsum=wtsum/nused
     wt(:)=wt_tmp(:)/wtsum
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
           do while (conf_int.le.10.0)
             ndigits = ndigits +1 
             conf_int = conf_int*10
           enddo
           write(fmt0,"(i2)") ndigits
           fmt0 = adjustl(fmt0)
           fmt1 = '(a2,I4,2x,f14.'//fmt0//'    ''    ('',i8,'')'')'
           write (f_out,FMT=fmt1) trim(parnam(i)),ivar(i),param(i),nint(conf_int)
           !
        else 
           ndigits =2
           if (param(i).ne.0.0) ndigits = 8

           write(fmt0,"(i2)") ndigits
           fmt0 = adjustl(fmt0)
           !
           fmt1 = '(a2,9i3,I4,2x,f14.'//fmt0//')'
           if (i<=nparam_eq) then 
             fmt1 = '(a10,21x,I4,2x,f14.'//fmt0//')'
           endif
           !
           fmt1 = '(a2,9i3,I4,2x,f14.'//fmt0//')'
           write (f_out,FMT=fmt1) trim(parnam(i)),ipower(1:9,i),ivar(i),param(i)


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
        write (f_out,"(a2,9i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:9,i),ivar(i),param(i)
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
      ncol = 0 
      do i=1,parmax
        if (ivar(i) /=  0) then
          ncol=ncol+1
          if (nused /= numpar) then  
            ssq = sum(rjacob(:,ncol)**2*wt(:))
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
              param(i) = 0.d0
              write(f_out,"(' parameter ',i6,' is  out - ',A10)") i,parnam(i) 
              cycle outer_loop    
            endif 
            !
            if (param(i)/=0.0d0.and.ivar(i) == tivar) then
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
                    if (i0/=i.and.param(i0)/=0.0d0.and.ivar(i0) == tivar) then
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
          do i=1,nparam_eq
             write (f_out,"(a10,22x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
          enddo 
          do i=nparam_eq+1,parmax
             write (f_out,"(a2,9i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:9,i),ivar(i),param(i)
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
      write (f_res,"(a2,9i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:9,i),ivar(i),param(i)
   enddo 
  endif  ! ---- itmax/=0
  
  write (f_res,"('------------------------------')")

    
   write (f_res,"(a80)") (title(i), i=1,4)
   do i=1,nparam_eq
      write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
   enddo 
   do i=nparam_eq+1,parmax
      write (f_res,"(a2,9i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:9,i),ivar(i),param(i)
   enddo 

  

 do nrow=1,npts
   !
    v = energy(nrow)-eps(nrow)
    write (f_res,"(4f15.6,5f11.4,3f18.7,2x,e9.2)")            & 
          coord(1:9,nrow),                                  &
          energy(nrow)-ve,v-ve,eps(nrow),                   &
          wt(nrow)

 enddo


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
    do i=1,nparam_eq
       write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
    enddo 
    do i=nparam_eq+1,parmax
       write (f_res,"(a2,9i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:9,i),ivar(i),param(i)
    enddo 
  endif  ! ---- itmax/=0
  
  write (f_res,"('------------------------------')")
  !
  write (f_res,"(a80)") (title(i), i=1,4)
  do i=1,nparam_eq
     write (f_res,"(a10,27x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
  enddo 
  do i=nparam_eq+1,parmax
     write (f_res,"(a2,9i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:9,i),ivar(i),param(i)
  enddo 
  
  

do nrow=1,npts
   !
    v = energy(nrow)-eps(nrow)
    write (f_res,"(4f15.6,1x,5f11.4,3f18.7,2x,e9.2)")            &
          coord(1:9,nrow),                                  &
          energy(nrow)-ve,v-ve,eps(nrow),                   &
          wt(nrow)
   !
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
 recursive function potential(parmax,dF,local,param,ipower,ipar) result (V)

 integer,intent(in) ::           parmax
 double precision,intent(in) ::  dF(parmax),local(18),param(parmax)
 integer,intent(in) ::           ipar,ipower(18,parmax)
 
 double precision         ::  param_(parmax)
 !
 double precision         ::  deltax,potright,potleft,V
 double precision         ::  dF_(parmax)
 integer                  ::  ierror,i,k
 character(len=40)        ::  txt
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
    deltax=factordeltax*abs(param(ipar))
    if (deltax .le. 1e-15) deltax=1e-6
    !
    param_ = param
    !
    param_(ipar) = param(ipar)+deltax
    !
    !call diff_V(nmax,local,param_(1:2),dF_)
    !
    do iparam = nparam_eq+1,parmax
      !
      dF_(iparam) = diff_V(ipower(1:9,iparam),local,param_(1:nparam_eq))
      !
    enddo
    !
    potright=potential(parmax,dF_,local,param_,ipower,0)
    !
    param_(ipar) = param(ipar)-deltax
    !
    !call diff_V(nmax,local,param_(1:2),dF_)
    !
    do iparam = 5,parmax
      !
      dF_(iparam) = diff_V(ipower(1:9,iparam),local,param_(1:nparam_eq))
      !
    enddo
    !
    potleft=potential(parmax,dF_,local,param_,ipower,0)
    !
    V=(potright-potleft)/(2.d0*deltax)
    !
  endif
  !
 end function potential
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

  double precision function diff_V(ipower,local,param)
  !
  implicit none
  !
  integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
  integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
  !
  integer,intent(in)  :: ipower(18)
  double precision,intent(in)  :: local(18),param(5)
  double precision             :: s1,s2,s3,s4 
  double precision         ::  v,r1,r2,r3,r4,alpha12,alpha13,alpha14,alpha23,alpha24,alpha34,deg,r0,r12,r13,r14,r23,r24,r34,t12,t13,t14,t23,t24,t34
  double precision         ::  r1e,r2e,a,y1,y2,y3,y4,y5,y6,y7,y8,y9,xi(18),chi(18),betae,pi,b,tau12,tau13,term,tau23
  integer                  ::  ierror,ioper
  character(len=40)        ::  txt
  !
  txt = 'diff_V'
  !
  pi = 4.0d0 * datan2(1.0d0,1.0d0)
  !
  deg=pi/180.0d0
  !
  r1      = local(1)
  r2      = local(2)
  r3      = local(3)
  r4      = local(4)
  !
  r1e      = param(1)
  r2e      = param(2)
  betae    = param(3)*deg
  !
  a       = param(4)
  b       = param(5)
  !
  xi(1)=1.0d+00-exp(-a*(r1-r1e))
  xi(2)=1.0d+00-exp(-b*(r2-r2e))
  xi(3)=1.0d+00-exp(-b*(r3-r2e))
  xi(4)=1.0d+00-exp(-b*(r4-r2e))
  !
  xi(5) = local(5)-betae
  xi(6) = local(6)-betae
  xi(7) = local(7)-betae
  tau12 = local(8)
  tau13 = local(9)
  tau23 = 2.0d0*pi-tau12-tau13
  !
  xi(8)=(2.0_rk*tau23-tau13-tau12)/sqrt(6.0_rk)
  xi(9)=(tau13-tau12)/sqrt(2.0_rk)
  !
  !cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
  !beta312 = aacos(cosbeta,ierror,txt)
  !
  !cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
  !
  !alpha34 = aacos(cosa34,ierror,txt)
  !
  term = 0
  !
  do ioper =1,6
    !
    call ML_symmetry_transformation_XY3_II(ioper,xi,chi)
    !
    term = term + product(chi(:)**ipower(:))
    !
  enddo
  !
  diff_V = term 
  !
!vdump = 1.0d0
!
!dF(10:289) = vdump*dF(10:289)
!dF(7:8) =       dF(7:8)*va
!dF(5:6) = va*dF(5:6)

end function diff_V

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
     xtmp(1)=1.0d0
     ytmp(1)=0
     ztmp(1)=0
     xtmp(2)=cos(alpha12)
     ytmp(2)=sin(alpha12)
     ztmp(2)=0


     xtmp(3)=cos(alpha13)
     ytmp(3)=(cos(alpha23)-xtmp(2)*xtmp(3))/ytmp(2)
     !
     temp = 1.0d0-xtmp(3)**2-ytmp(3)**2
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
     temp = 1.0d0-xtmp(4)**2-ytmp(4)**2
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
     if (abs(temp)>1.d00) then 
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

       small_ = epsilon(1.0d0)
       ca12 = cos(alpha12)
       ca13 = cos(alpha13)
       ca14 = cos(alpha14)
       ca23 = cos(alpha23)
       ca24 = cos(alpha24)

       cosalpha = 0.5d0*( &
                  2.0d0*ca12*ca23*ca14-     &
                  2.0d0*ca13*ca14+     &
                  2.0d0*ca12*ca24*ca13-     &
                  2.0d0*ca23*ca24+     &
                  2.0d0*sqrt(     &
                  1.0d0-ca24**2+     &
                  2.0d0*ca12*ca24*ca14+     &
                  2.0d0*ca12*ca23*ca13-     &
                  2.0d0*ca12**3*ca23*ca13-     &
                  2.0d0*ca12**2-     &
                  ca14**2-     &
                  ca13**2+     &
                  ca14**2*ca23**2+     &
                  ca13**2*ca24**2+     &
                  4.0d0*ca12**2*ca23*ca14*ca24*ca13+     &
                  ca12**4-     &
                  2.0d0*ca12*ca23*ca14**2*ca13-     &
                  2.0d0*ca12*ca23**2*ca14*ca24-     &
                  2.0d0*ca13**2*ca14*ca12*ca24-     &
                  2.0d0*ca12*ca24**2*ca13*ca23+     &
                  ca13**2*ca14**2+     &
                  ca23**2*ca24**2+     &
                  ca12**2*ca13**2+     &
                  ca12**2*ca14**2+     &
                  ca12**2*ca23**2+     &
                  ca12**2*ca24**2-     &
                  2.0d0*ca12**3*ca24*ca14-     &
                  ca23**2))/( -sin(alpha12)**2 )

       if ( abs(cosalpha)>1.0d0+sqrt(small_) ) then 
          !
          write (6,"('calc_alpha34: cosalpha>1: ',f18.8)") cosalpha
          stop 'calc_alpha34 - bad cosalpha'
          !
       elseif ( cosalpha>=1.0d0) then 
          alpha34 = 0.0d0
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
  integer                      :: i1,i2,k,i,j,k8,k9
      

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
  subroutine ML_symmetry_transformation_XY3(ioper,src,dst)
    implicit none 
    !
    integer,intent(in)    :: ioper  ! group operation  
    double precision,intent(in)      :: src(1:9)
    double precision,intent(out)     :: dst(1:9)
    !
    double precision         :: repres(12,2,2),a,b,e,o
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
    a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
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
    double precision,intent(in)      :: src(1:9)
    double precision,intent(out)     :: dst(1:9)
    !
    double precision         :: repres(6,9,9),a,b,e,o
    !
    if (verbose>=5) write(f_out,"('ML_symmetry_transformation_XY3/start')") 
    !
    a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
    !
    repres = 0 
    !
    repres(1,1,1) = 1.0d0
    repres(1,2,2) = 1.0d0
    repres(1,3,3) = 1.0d0
    repres(1,4,4) = 1.0d0
    repres(1,5,5) = 1.0d0
    repres(1,6,6) = 1.0d0
    repres(1,7,7) = 1.0d0
    repres(1,8,8) = 1.0d0
    repres(1,9,9) = 1.0d0
    !
    repres(2,1,1) = 1.0d0
    repres(2,2,4) = 1.0d0
    repres(2,3,2) = 1.0d0
    repres(2,4,3) = 1.0d0
    repres(2,5,7) = 1.0d0
    repres(2,6,5) = 1.0d0
    repres(2,7,6) = 1.0d0
    repres(2,8,8) = -a
    repres(2,8,9) = -b
    repres(2,9,8) =  b
    repres(2,9,9) = -a
    !
    repres(3,1,1) = 1.0d0
    repres(3,2,3) = 1.0d0
    repres(3,3,4) = 1.0d0
    repres(3,4,2) = 1.0d0
    repres(3,5,6) = 1.0d0
    repres(3,6,7) = 1.0d0
    repres(3,7,5) = 1.0d0
    repres(3,8,8) = -a
    repres(3,8,9) =  b
    repres(3,9,8) = -b
    repres(3,9,9) = -a
    !
    repres(4,1,1) = 1.0d0
    repres(4,2,2) = 1.0d0
    repres(4,3,4) = 1.0d0
    repres(4,4,3) = 1.0d0
    repres(4,5,5) = 1.0d0
    repres(4,6,7) = 1.0d0
    repres(4,7,6) = 1.0d0
    repres(4,8,8) = 1.0d0
    repres(4,8,9) = 0.0d0
    repres(4,9,8) = 0.0d0
    repres(4,9,9) = -1.0d0
    !
    repres(5,1,1) = 1.0d0
    repres(5,2,4) = 1.0d0
    repres(5,3,3) = 1.0d0
    repres(5,4,2) = 1.0d0
    repres(5,5,7) = 1.0d0
    repres(5,6,6) = 1.0d0
    repres(5,7,5) = 1.0d0
    repres(5,8,8) = -a
    repres(5,8,9) = -b
    repres(5,9,8) = -b
    repres(5,9,9) =  a
    !
    repres(6,1,1) = 1.0d0
    repres(6,2,3) = 1.0d0
    repres(6,3,2) = 1.0d0
    repres(6,4,4) = 1.0d0
    repres(6,5,6) = 1.0d0
    repres(6,6,5) = 1.0d0
    repres(6,7,7) = 1.0d0
    repres(6,8,8) = -a
    repres(6,8,9) = b
    repres(6,9,8) = b
    repres(6,9,9) = a
    !
    if (ioper<0.or.ioper>6) then
      write (6,"('symmetry_transformation_local: operation ',i8,' unknown')") ioper
      stop 'symmetry_transformation_local - bad operation. type'
    endif
    !
    dst = matmul(repres(ioper,:,:),src) 
           
  end subroutine ML_symmetry_transformation_XY3_II      
  
  
  
  
  end module FITPES_ch3cl
  
  
  
  
  program driver
  use FITPES_ch3cl
  
  call do_fitting
  
  end program driver
