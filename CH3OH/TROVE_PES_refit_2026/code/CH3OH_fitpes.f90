#ifdef OPENBLAS64
module ob64_lapack
  ! Optional multithreaded OpenBLAS (ILP64, scipy-openblas symbol names)
  ! interfaces; compile with -cpp -DOPENBLAS64 and link libscipy_openblas64_.
  use iso_c_binding
  implicit none
  interface
    subroutine ob_dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc, l1, l2) &
        bind(C, name='scipy_dsyrk_64_')
      import :: c_char, c_int64_t, c_double, c_size_t
      character(kind=c_char) :: uplo, trans
      integer(c_int64_t) :: n, k, lda, ldc
      real(c_double) :: alpha, beta, a(lda, *), c(ldc, *)
      integer(c_size_t), value :: l1, l2
    end subroutine ob_dsyrk
    subroutine ob_dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy, l1) &
        bind(C, name='scipy_dgemv_64_')
      import :: c_char, c_int64_t, c_double, c_size_t
      character(kind=c_char) :: trans
      integer(c_int64_t) :: m, n, lda, incx, incy
      real(c_double) :: alpha, beta, a(lda, *), x(*), y(*)
      integer(c_size_t), value :: l1
    end subroutine ob_dgemv
    subroutine ob_dpotrf(uplo, n, a, lda, info, l1) bind(C, name='scipy_dpotrf_64_')
      import :: c_char, c_int64_t, c_double, c_size_t
      character(kind=c_char) :: uplo
      integer(c_int64_t) :: n, lda, info
      real(c_double) :: a(lda, *)
      integer(c_size_t), value :: l1
    end subroutine ob_dpotrf
    subroutine ob_dpotrs(uplo, n, nrhs, a, lda, b, ldb, info, l1) bind(C, name='scipy_dpotrs_64_')
      import :: c_char, c_int64_t, c_double, c_size_t
      character(kind=c_char) :: uplo
      integer(c_int64_t) :: n, nrhs, lda, ldb, info
      real(c_double) :: a(lda, *), b(ldb, *)
      integer(c_size_t), value :: l1
    end subroutine ob_dpotrs
    subroutine ob_dpotri(uplo, n, a, lda, info, l1) bind(C, name='scipy_dpotri_64_')
      import :: c_char, c_int64_t, c_double, c_size_t
      character(kind=c_char) :: uplo
      integer(c_int64_t) :: n, lda, info
      real(c_double) :: a(lda, *)
      integer(c_size_t), value :: l1
    end subroutine ob_dpotri
  end interface
end module ob64_lapack
#endif

module FITPES_ch3oh
  !
  ! Least-squares fitting of the 12D C3v(M) symmetry-adapted potential energy
  ! function (PEF) of CH3OH to ab initio energies.
  !
  ! Built on the structure of CH3Br_fitpes.f90 (FITPES_ch3br: compact
  ! one-pass evaluation of all symmetrised terms, OpenMP, analytic Jacobian
  ! for the linear parameters, TROVE output) and on the CH3OH specific
  ! coordinates/symmetry of CH3OH_v10_theta12_14-MAPLESym-13-09.f90
  ! (FITPES_ch3cl: diff_V_tau + ML_symmetry_transformation_XY3_IV,
  ! Partridge-Schwenke weights, Watson robust weights). The input file
  ! format is the same as for CH3OH_v10_theta12_14-MAPLESym-13-09.f90.
  !
  ! Coordinates (identical to TROVE TRANSFORM R-ALPHA132-THETA-S-TAU, the
  ! AB3CD.nb exact KEO and the working potential.chk):
  !   input (Molpro z-matrix): rCO, rOH, rCH1, rCH2, rCH3, aCOH, aOCH1..3,
  !                            dHOCH1..3 (dihedrals H-O-C-Hi)
  !   xi(1)  = 1-exp(-a1*(rCO-rCOe))         xi(2) = 1-exp(-a2*(rOH-rOHe))
  !   xi(3)  = 1-exp(-b1*(rCH1-rCHe))
  !   xi(4)  = 1-exp(-b1*(rCH3-rCHe))        xi(5) = 1-exp(-b1*(rCH2-rCHe))
  !   xi(6)  = aCOH-aCOHe
  !   xi(7)  = aOCH1-aOCHe  xi(8) = aOCH3-aOCHe  xi(9) = aOCH2-aOCHe
  !   th1 = d1, th2 = d3, th3 = d2 (H2<->H3 relabelling of the z-matrix)
  !   th12 = (th1-th2) mod 2pi, th23 = (th2-th3) mod 2pi, th13 = 2pi-th12-th23
  !   xi(10) = S_a = (2*th23-th12-th13)/sqrt(6)
  !   xi(11) = S_b = (th13-th12)/sqrt(2)
  !   xi(12) = tau = (th1 + sqrt(2)/3*S_b) mod 2pi   [chart-consistent]
  ! The chart-consistent tau is the value for which the TROVE/AB3CD.nb inverse
  ! transform  t1 = tau - sqrt(2)/3 S_b, t2 = 4pi/3 + tau + S_a/sqrt(6) +
  ! S_b/(3 sqrt(2)), t3 = 2pi/3 + tau - S_a/sqrt(6) + S_b/(3 sqrt(2))
  ! reproduces the input dihedrals. It coincides with the CH3OH_v10 value
  ! (th1+th2+th3-2pi)/3 whenever the latter is unambiguous (true for the whole
  ! CH3OH-CBS grid), but, unlike it, it is invariant under relabelling of the
  ! methyl protons (option 'tau raw' restores the old definition).
  !
  ! PEF: V = sum_i f_i * (1/6) sum_{g in C3v(M)} prod_{m=1}^{11} chi_g(m)^p_im
  !                                                * T_{k_i}(chi_g(12))
  !      T_k = cos(k tau) (k>=0), sin(|k| tau) (k<0); chi_g = g.xi with the
  !      C3v(M)-2 operations of TROVE (ML_symmetry_transformation_XY3_IV).
  !
  implicit none
  public :: do_fitting
  private

  integer, parameter :: rk = selected_real_kind(12, 25)
  integer, parameter :: f_inp = 5, f_out = 6, f_res = 9, f_trove = 11
  integer, parameter :: ndeg = 12, nsym = 6, nparam_eq_ = 8
  integer, parameter :: maxterm_nz = 11
  real(rk), parameter :: pi = 3.14159265358979323846264338327950288_rk
  real(rk), parameter :: a_wats = 0.1_rk
  real(rk), parameter :: fact_corr = 0.4_rk, too_big = 1.0e7_rk

  ! compressed representation of the linear terms
  integer              :: nlin = 0, pmax = 0, kmax = 0
  integer, allocatable :: tnz(:), tmode(:, :), tpow(:, :), tk(:)
  logical              :: tau_chart = .true.

contains

  !===========================================================================
  subroutine do_fitting
    integer  :: NDEG_in, NBONDS, NANGS, itmax, robustfit, nparam_eq
    real(rk) :: stab_best, stadev_best, conv_length, conv_energy, factor_out
    real(rk) :: factor_back, tol, fit_scale, gamma_ps, enermax_ps, ridge
    character(len=80)  :: title(4)
    character(len=2048) :: line, key
    character(len=10)  :: label
    character(len=256) :: trove_file
    integer  :: i, j, k, n, ios, npts, parmax, ivartmp, ipw(ndeg), alloc
    real(rk) :: paramtmp
    character(len=10), allocatable :: parnam(:)
    integer, allocatable  :: ivar(:), ipower(:, :), ivar0(:)
    real(rk), allocatable :: param(:), parold(:), last_param(:)
    real(rk), allocatable :: coord(:, :), energy(:), wt(:), wt0(:), sigma(:)
    real(rk), allocatable :: tmp2(:, :), tmp1(:), tmpw(:)
    integer, allocatable  :: tmpi(:)
    logical, allocatable  :: real_pt(:)
    integer :: idxtmp, nreal
    real(rk) :: wref
    real(rk), allocatable :: eps(:), vcalc(:), rjacob(:, :), al(:, :), ai(:, :)
    real(rk), allocatable :: bl(:), dx(:), sterr(:)
    integer, allocatable  :: icol(:)
    integer  :: numpar, nused, iter, l, l0, last_i, Nout, tivar
    integer  :: minivar, maxivar, ilarge, ilow, ierr
    integer, allocatable :: icorr(:)
    real(rk) :: ssq, rms, stadev, stadev_old, stability, Emin, ve
    real(rk) :: relat, relat_min, corr, last_stadev
    logical  :: still_run, find_best
    integer  :: nband
    real(rk) :: bands(7)
    data bands / 2000.0_rk, 5000.0_rk, 10000.0_rk, 15000.0_rk, 20000.0_rk, &
                 30000.0_rk, 1.0e30_rk /

    ridge = 0.0_rk
    trove_file = 'ch3oh_poten.trove'

    call skiplines(f_inp, 1)
    do i = 1, 4
      read (f_inp, "(a80)") title(i)
    enddo
    write (f_out, "(3(132('*')/' '),16('*'),100x,16('*'))")
    write (f_out, "(4(' ',16('*'),10x,a80,10x,16('*')/))") (title(i), i = 1, 4)
    write (f_out, "(' ',16('*'),100x,16('*')/' ',3(132('*')/' '))")
    call skiplines(f_inp, 4)
    read (f_inp, *) NDEG_in
    read (f_inp, *) NBONDS
    read (f_inp, *) NANGS
    read (f_inp, *) itmax
    read (f_inp, *) stab_best, stadev_best
    read (f_inp, *) conv_length
    read (f_inp, *) conv_energy
    read (f_inp, *) factor_out
    read (f_inp, *) factor_back
    read (f_inp, *) robustfit
    read (f_inp, *) tol
    read (f_inp, *) fit_scale
    read (f_inp, *) gamma_ps
    read (f_inp, *) enermax_ps
    !
    ! optional "keyword value" lines (until a line starting with '*')
    do
      read (f_inp, "(a)") line
      line = adjustl(line)
      if (line(1:1) == '*') exit
      read (line, *) key
      select case (trim(key))
      case ('ridge', 'RIDGE')
        read (line, *) key, ridge
      case ('tau', 'TAU')
        read (line, *) key, key
        tau_chart = (trim(key) /= 'raw' .and. trim(key) /= 'RAW')
      case ('trove', 'TROVE')
        read (line, *) key, trove_file
      case default
        write (f_out, "('unknown option line: ',a)") trim(line)
        stop 'unknown option'
      end select
    enddo
    ! the '*' line already consumed is the first of the 3 header lines
    call skiplines(f_inp, 2)
    !
    if (NDEG_in /= ndeg .or. NBONDS /= 5 .or. NANGS /= 7) &
      stop 'CH3OH_fitpes: NDEG/NBONDS/NANGS must be 12/5/7'
    write (f_out, "(i10,' <= number of degrees of freedom')") NDEG_in
    write (f_out, "(i10,' <= number of iterations in the fit')") itmax
    write (f_out, "(es10.2,' <= stability / ',es10.2,' st.dev. convergence')") &
      stab_best, stadev_best
    write (f_out, "(f10.4,' <= factor_out ; ',es10.2,' <= factor_back')") &
      factor_out, factor_back
    write (f_out, "(i10,' <= Watson robust fit: 0 off, 1/2 Watson (CH3OH_v10), 4 scaled Watson')") robustfit
    write (f_out, "(es10.2,' <= tol (>0: eigen-cutoff pseudo-inverse, <=0: Cholesky)')") tol
    write (f_out, "(es10.2,' <= ridge (Tikhonov, relative to diag(A))')") ridge
    write (f_out, "(f10.6,' <= Partridge-Schwenke gamma; Emax = ',f12.2)") gamma_ps, enermax_ps
    if (tau_chart) then
      write (f_out, "(' tau definition: chart-consistent (TROVE inverse transform)')")
    else
      write (f_out, "(' tau definition: raw (th1+th2+th3-2pi)/3 as in CH3OH_v10')")
    endif
    find_best = (robustfit >= 0)
    !
    read (f_inp, *) nparam_eq
    if (nparam_eq /= nparam_eq_) stop 'CH3OH_fitpes: 8 non-linear parameters expected'
    !
    ! --- read parameters (dynamic) ---
    allocate (parnam(20000), ivar(20000), param(20000), ipower(ndeg, 20000))
    ipower = 0
    do i = 1, nparam_eq
      read (f_inp, "(a)") line
      read (line, *) label, ivartmp, paramtmp
      parnam(i) = label; ivar(i) = ivartmp; param(i) = paramtmp
    enddo
    parmax = nparam_eq
    do
      read (f_inp, "(a)") line
      if (len_trim(adjustl(line)) == 0) cycle
      line = adjustl(line)
      if (line(1:3) == 'end' .or. line(1:3) == 'END') exit
      parmax = parmax + 1
      if (parmax > size(param)) stop 'too many parameters (20000)'
      read (line, *) parnam(parmax), ipw(1:ndeg), ivar(parmax), param(parmax)
      ipower(:, parmax) = ipw
    enddo
    nlin = parmax - nparam_eq
    write (f_out, "(/i8,' linear parameters read')") nlin
    call compress_terms(ipower(:, nparam_eq + 1:parmax))
    !
    ! --- read ab initio points (dynamic) ---
    allocate (tmp2(ndeg, 200000), tmp1(200000), tmpw(200000), tmpi(200000))
    npts = 0
    do
      read (f_inp, "(a)", iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      npts = npts + 1
      if (npts > size(tmp1)) stop 'too many points (200000)'
      ! optional 15th column: point index (negative = pseudo-point, excluded
      ! from the statistics)
      read (line, *, iostat=ios) tmp2(:, npts), tmp1(npts), tmpw(npts), idxtmp
      if (ios /= 0) then
        read (line, *) tmp2(:, npts), tmp1(npts), tmpw(npts)
        idxtmp = npts
      endif
      tmpi(npts) = idxtmp
      if (tmp2(1, npts) < 0.0_rk) then
        npts = npts - 1
        exit
      endif
    enddo
    allocate (coord(ndeg, npts), energy(npts), wt(npts), wt0(npts), sigma(npts), &
              eps(npts), vcalc(npts), real_pt(npts), stat=alloc)
    if (alloc /= 0) stop 'points - out of memory'
    do n = 1, npts
      coord(1:5, n) = tmp2(1:5, n) * conv_length
      coord(6:12, n) = tmp2(6:12, n) * pi / 180.0_rk
      energy(n) = tmp1(n) * conv_energy
      wt(n) = tmpw(n)
      if (wt(n) < 0.0_rk) wt(n) = abs(wt(n)) * factor_back
      real_pt(n) = (tmpi(n) >= 0)
    enddo
    deallocate (tmp2, tmp1, tmpw, tmpi)
    nreal = count(real_pt)
    Emin = minval(energy, mask=real_pt)
    write (f_out, "(i8,' points read (',i8,' ab initio, ',i8,' pseudo); lowest energy = ',f20.6)") &
      npts, nreal, npts - nreal, Emin
    !
    ! Partridge-Schwenke energy weights
    if (enermax_ps > 0.0_rk) then
      do n = 1, npts
        wt(n) = wt(n) * (tanh(-gamma_ps * (energy(n) - Emin - enermax_ps)) + &
                1.000200020002_rk) / 2.000200020002_rk
      enddo
    endif
    call normalise_weights(npts, wt, nused)
    do n = 1, npts
      if (wt(n) > 0.0_rk) then
        sigma(n) = 1.0_rk / wt(n)
      else
        sigma(n) = 1000.0_rk
      endif
    enddo
    wt0 = wt
    ! robust weights of real points with small residuals become ~wt0^2; to keep
    ! the intended pseudo/real weight ratio pseudo-points get wt0*wref, with
    ! wref = the largest initial weight of an ab initio point
    wref = maxval(wt0, mask=real_pt)
    !
    allocate (parold(parmax), last_param(parmax), ivar0(parmax), icorr(parmax), &
              sterr(parmax), icol(parmax))
    sterr = 0.0_rk
    ivar0 = ivar
    last_param = param
    last_i = 0
    last_stadev = 1.0e10_rk
    Nout = 1
    !
    still_run = .true.
    outer_loop: do while (still_run)
      !
      numpar = 0
      maxivar = 0
      minivar = huge(1)
      do i = 1, parmax
        if (ivar(i) /= 0) then
          numpar = numpar + 1
          icol(numpar) = i
          maxivar = max(maxivar, ivar(i))
          minivar = min(minivar, ivar(i))
        endif
      enddo
      if (allocated(rjacob)) deallocate (rjacob, al, ai, bl, dx)
      allocate (rjacob(npts, max(numpar, 1)), al(numpar, numpar), ai(numpar, numpar), &
                bl(numpar), dx(numpar), stat=alloc)
      if (alloc /= 0) stop 'jacobian - out of memory'
      !
      if (itmax /= 0) then
        write (f_out, "(/3x,66('-'))")
        write (f_out, "('   |   iter  |  points |    deviat     |     rms       |  stability  |')")
        write (f_out, "(3x,66('-'))")
      endif
      stadev_old = 1.0e10_rk
      stability = 1.0e10_rk
      stadev = 1.0e10_rk
      iter = 0
      !
      iter_loop: do while (iter <= itmax .and. stadev > stadev_best .and. stability >= stab_best)
        iter = iter + 1
        parold = param
        !
        call build_jacobian(npts, coord, param, ivar, numpar, icol, nparam_eq, &
                            vcalc, rjacob, itmax /= 0)
        eps = energy - vcalc
        ssq = sum(eps**2 * wt)
        rms = sqrt(sum(eps**2) / npts)
        if (itmax == 0) then
          stadev = sqrt(ssq / nused)
          exit iter_loop
        endif
        !
        call normal_equations(npts, numpar, rjacob, wt, eps, al, bl)
        call solve_normal(numpar, al, bl, dx, ai, tol, ridge, ierr)
        if (ierr > 0) then
          ! singular in Cholesky: drop the offending parameter and restart
          i = icol(ierr)
          write (f_out, "(i8,'-th is out (singular normal matrix) - ',a10,12i3)") &
            i, parnam(i), ipower(:, i)
          ivar(i) = 0
          param(i) = 0.0_rk
          cycle outer_loop
        endif
        !
        do l = 1, numpar
          param(icol(l)) = param(icol(l)) + dx(l) * fit_scale
        enddo
        ! residuals after the update (exact for the linear parameters,
        ! first order for the non-linear ones)
        !$omp parallel do schedule(static)
        do n = 1, npts
          eps(n) = eps(n) - fit_scale * dot_product(rjacob(n, 1:numpar), dx(1:numpar))
        enddo
        !$omp end parallel do
        ssq = sum(eps**2 * wt)
        rms = sqrt(sum(eps**2) / npts)
        if (nused > numpar) then
          stadev = sqrt(ssq / real(nused - numpar, rk))
        else
          stadev = sqrt(ssq / real(nused, rk))
        endif
        do l = 1, numpar
          sterr(l) = sqrt(abs(ai(l, l))) * stadev
        enddo
        stability = abs((stadev - stadev_old) / stadev)
        stadev_old = stadev
        !
        if (robustfit > 0) then
          do n = 1, npts
            ! Watson robust weights for the ab initio points only; pseudo-points
            ! keep their (pre-normalisation) relative weights
            if (wt(n) > 0.0_rk .and. real_pt(n)) then
              if (robustfit == 4) then
                ! scaled Watson: initial weight times 1/(1 + a (w0 eps)^2), i.e. the
                ! energy weighting is kept and only outliers are down-weighted
                wt(n) = wt0(n) / (1.0_rk + a_wats * (wt0(n) * eps(n))**2)
              else
                ! original Watson (CH3OH_v10): 1/(sigma^2 + a eps^2), sigma = 1/w0
                wt(n) = 1.0_rk / (sigma(n)**2 + a_wats * eps(n)**2)
              endif
            endif
            if (.not. real_pt(n)) then
              if (robustfit == 4) then
                wt(n) = wt0(n)
              else
                wt(n) = wt0(n) * wref
              endif
            endif
          enddo
          call normalise_weights(npts, wt, nused)
          ssq = sum(eps**2 * wt)
          if (nused > numpar) then
            stadev = sqrt(ssq / real(nused - numpar, rk))
          else
            stadev = sqrt(ssq / real(nused, rk))
          endif
        endif
        write (f_out, "('   |  ',i4,'   | ',i6,'  |  ',es12.5,' | ',es12.5,'  |  ',es10.3,' |')") &
          iter, nused, stadev, rms, stability
        flush (f_out)
      enddo iter_loop
      !
      ! residuals for the final parameters
      call build_jacobian(npts, coord, param, ivar, numpar, icol, nparam_eq, &
                          vcalc, rjacob, .false.)
      eps = energy - vcalc
      rms = sqrt(sum(eps**2) / npts)
      if (robustfit == 2) robustfit = 0
      !
      if (itmax /= 0 .and. find_best) then
        ! parameter taken back if the st.dev. got too much worse
        if (last_i > 0) then
          if ((stadev - last_stadev) / stadev > Nout * factor_back .and. &
              abs(last_param(last_i)) < too_big) then
            ivar = ivar0
            ivar(last_i) = ivar(last_i) + 1
            param = last_param
            last_stadev = 1.0e10_rk
            write (f_out, "(i8,'-th is back - ',a10)") last_i, parnam(last_i)
            last_i = 0
            cycle outer_loop
          endif
        endif
        call print_params(.true.)
        !
        ! elimination of parameters with large relative st. errors
        tivar = minivar
        ivar0 = ivar
        icorr = 0
        last_param = param
        Nout = 0
        do while (tivar < maxivar .and. factor_out > 0.0_rk)
          relat_min = 1.0e20_rk
          ilow = -1
          do l = numpar, 1, -1
            i = icol(l)
            if (i <= nparam_eq) cycle
            if (abs(param(i)) > too_big .and. ivar(i) /= maxivar .and. last_i /= i) then
              last_param = param; last_stadev = stadev; last_i = i
              ivar(i) = 0; param(i) = 0.0_rk
              write (f_out, "(' parameter ',i6,' is out (too big) - ',a10)") i, parnam(i)
              cycle outer_loop
            endif
            if (param(i) /= 0.0_rk .and. ivar(i) == tivar) then
              if (icorr(i) == 0 .and. factor_out < abs(sterr(l) / param(i))) then
                icorr(i) = 1
                relat = abs(sterr(l) / param(i))
                ilarge = i
                do l0 = 1, numpar
                  j = icol(l0)
                  corr = abs(ai(l, l0)) / sqrt(abs(ai(l, l) * ai(l0, l0)))
                  if (j /= i .and. param(j) /= 0.0_rk .and. ivar(j) == tivar .and. j > nparam_eq) then
                    if (icorr(j) == 0 .and. corr > fact_corr) then
                      icorr(j) = 1
                      if (relat < abs(sterr(l0) / param(j))) then
                        relat = abs(sterr(l0) / param(j))
                        ilarge = j
                      endif
                    endif
                  endif
                enddo
                ivar(ilarge) = 0
                param(ilarge) = 0.0_rk
                Nout = Nout + 1
                if (relat_min > relat) then
                  relat_min = relat
                  ilow = ilarge
                endif
              endif
            endif
          enddo
          if (relat_min > factor_out .and. ilow /= -1) then
            last_stadev = stadev
            last_i = ilow
            write (f_out, "(i8,' parameters are out')") Nout
            cycle outer_loop
          endif
          tivar = tivar + 1
        enddo
      endif
      still_run = .false.
    enddo outer_loop
    !
    call print_params(.false.)
    call write_res
    call write_trove
    !
    ! RMS by energy bands (relative to the lowest ab initio energy)
    write (f_out, "(/' Unweighted rms deviations by energy band (E-Emin < Emax):')")
    do nband = 1, size(bands)
      k = count(energy - Emin < bands(nband) .and. real_pt)
      if (k > 0) then
        rms = sqrt(sum(eps**2, mask=(energy - Emin < bands(nband) .and. real_pt)) / k)
        write (f_out, "('   Emax = ',f12.0,'  npts = ',i7,'  rms = ',f12.4,'  max|dE| = ',f12.4)") &
          min(bands(nband), 1.0e9_rk), k, rms, &
          maxval(abs(eps), mask=(energy - Emin < bands(nband) .and. real_pt))
      endif
    enddo
    if (nreal < npts) then
      k = npts - nreal
      write (f_out, "('   pseudo-points: n = ',i7,'  rms = ',f12.4,'  min(V-E) = ',f12.4,'  max(V-E) = ',f12.4)") &
        k, sqrt(sum(eps**2, mask=.not. real_pt) / k), minval(-eps, mask=.not. real_pt), &
        maxval(-eps, mask=.not. real_pt)
    endif
    rms = sqrt(sum(eps**2) / npts)
    write (f_out, "(/3x,59('-')/'   |  param  |  points |   deviat    |    rms      | stability |')")
    write (f_out, "(3x,59('-'))")
    write (f_out, "('   |  ',i5,'  | ',i6,'  |  ',es12.5,' | ',es12.5,'  |  ',es11.3,' |')") &
      numpar, nused, stadev, rms, stability
    write (f_out, "(/92('#'))")

  contains

    subroutine print_params(with_errors)
      logical, intent(in) :: with_errors
      integer :: ii, ll
      if (with_errors) then
        write (f_out, "(//'Potential parameters (with st. errors):'/)")
      else
        write (f_out, "(//'Final potential parameters:'/)")
      endif
      ll = 0
      do ii = 1, parmax
        if (ii <= nparam_eq) then
          if (ivar(ii) /= 0) then
            ll = ll + 1
            write (f_out, "(a10,4x,i4,2x,f20.10,'   (',es10.2,')')") parnam(ii), ivar(ii), param(ii), &
              merge(sterr(ll), 0.0_rk, with_errors)
          else
            write (f_out, "(a10,4x,i4,2x,f20.10)") parnam(ii), ivar(ii), param(ii)
          endif
        else
          if (ivar(ii) /= 0) then
            ll = ll + 1
            write (f_out, "(1x,a2,12i4,i4,2x,es22.14,'   (',es10.2,')')") trim(parnam(ii)), &
              ipower(:, ii), ivar(ii), param(ii), merge(sterr(ll), 0.0_rk, with_errors)
          elseif (with_errors .eqv. .false.) then
            if (param(ii) /= 0.0_rk) &
              write (f_out, "(1x,a2,12i4,i4,2x,es22.14)") trim(parnam(ii)), ipower(:, ii), &
                ivar(ii), param(ii)
          endif
        endif
      enddo
    end subroutine print_params

    subroutine write_res
      integer :: ii, nn
      open (f_res, file='00.res', status='replace')
      write (f_res, "(a80)") (title(ii), ii = 1, 4)
      do ii = 1, nparam_eq
        write (f_res, "(a10,24x,i4,2x,f20.10)") parnam(ii), ivar(ii), param(ii)
      enddo
      do ii = nparam_eq + 1, parmax
        if (ivar(ii) /= 0 .or. param(ii) /= 0.0_rk) &
          write (f_res, "(1x,a2,12i4,i4,2x,es22.14)") trim(parnam(ii)), ipower(:, ii), &
            ivar(ii), param(ii)
      enddo
      write (f_res, "('end')")
      ve = 0.0_rk
      do ii = nparam_eq + 1, parmax
        if (all(ipower(:, ii) == 0)) ve = param(ii)
      enddo
      do nn = 1, npts
        write (f_res, "(5f12.7,7f12.6,3f18.6,2x,es9.2)") coord(1:5, nn) / conv_length, &
          coord(6:12, nn) * 180.0_rk / pi, energy(nn) - Emin, vcalc(nn) - Emin, eps(nn), wt(nn)
      enddo
      write (f_res, "(/3x,59('-')/'   |  param  |  points |   deviat    |    rms      | stability |')")
      write (f_res, "(3x,59('-'))")
      write (f_res, "('   |  ',i5,'  | ',i6,'  |  ',es12.5,' | ',es12.5,'  |  ',es11.3,' |')") &
        numpar, nused, stadev, sqrt(sum(eps**2) / npts), stability
      close (f_res)
    end subroutine write_res

    subroutine write_trove
      integer :: ii, nt, nn
      real(rk) :: v0
      nt = 0
      v0 = 0.0_rk
      do ii = nparam_eq + 1, parmax
        if (all(ipower(:, ii) == 0)) v0 = param(ii)
        if (param(ii) /= 0.0_rk) nt = nt + 1
      enddo
      open (f_trove, file=trim(trove_file), status='replace')
      write (f_trove, "('! CH3OH PEF fitted with CH3OH_fitpes (TROVE R-ALPHA132-THETA-S-TAU, C3v(M)-2)')")
      write (f_trove, "('! Constant term V0 = ',es24.16,' cm-1 (energies relative to the lowest')") v0
      write (f_trove, "('! ab initio point); it is dropped when potential.chk is generated.')")
      write (f_trove, "('! The TROVE potential coordinates must use the SAME equilibrium and Morse')")
      write (f_trove, "('! parameters as the fit:')")
      write (f_trove, "('EQUILIBRIUM')")
      write (f_trove, "('re1          0  ',f16.10)") param(1)
      write (f_trove, "('re2          0  ',f16.10)") param(2)
      write (f_trove, "('re3          0  ',f16.10)") param(3)
      write (f_trove, "('re4          0  ',f16.10)") param(3)
      write (f_trove, "('re5          0  ',f16.10)") param(3)
      write (f_trove, "('alphae1      0  ',f16.10,' DEG')") param(4)
      write (f_trove, "('alphae2      0  ',f16.10,' DEG')") param(5)
      write (f_trove, "('alphae3      0  ',f16.10,' DEG')") param(5)
      write (f_trove, "('alphae4      0  ',f16.10,' DEG')") param(5)
      write (f_trove, "('taue1        0       0.0000      DEG')")
      write (f_trove, "('taue2        0     120.0000      DEG')")
      write (f_trove, "('taue3        0     240.0000      DEG')")
      write (f_trove, "('end'/)")
      write (f_trove, "('SPECPARAM')")
      write (f_trove, "('beta1        0  ',f16.10)") param(6)
      write (f_trove, "('beta2        0  ',f16.10)") param(7)
      write (f_trove, "('beta3        0  ',f16.10)") param(8)
      write (f_trove, "('beta4        0  ',f16.10)") param(8)
      write (f_trove, "('beta5        0  ',f16.10)") param(8)
      write (f_trove, "('END'/)")
      write (f_trove, "('POTEN')")
      write (f_trove, "('NPARAM  ',i0)") nparam_eq + nt
      write (f_trove, "('POT_TYPE  poten_ch3oh_sym_IV')")
      write (f_trove, "('COEFF  powers  (powers or list)')")
      do ii = 1, nparam_eq
        write (f_trove, "(a4,1x,12i4,i4,2x,f22.12)") parnam(ii), (0, nn = 1, 12), 1, param(ii)
      enddo
      do ii = nparam_eq + 1, parmax
        if (param(ii) == 0.0_rk) cycle
        write (f_trove, "(1x,a2,1x,12i4,i4,2x,es24.16)") trim(parnam(ii)), ipower(:, ii), 1, param(ii)
      enddo
      write (f_trove, "('end')")
      close (f_trove)
    end subroutine write_trove

  end subroutine do_fitting

  !===========================================================================
  subroutine compress_terms(ip)
    integer, intent(in) :: ip(:, :)
    integer :: i, m, n
    nlin = size(ip, 2)
    allocate (tnz(nlin), tmode(maxterm_nz, nlin), tpow(maxterm_nz, nlin), tk(nlin))
    tmode = 1; tpow = 0
    pmax = 0; kmax = 0
    do i = 1, nlin
      n = 0
      do m = 1, 11
        if (ip(m, i) > 0) then
          n = n + 1
          tmode(n, i) = m
          tpow(n, i) = ip(m, i)
          pmax = max(pmax, ip(m, i))
        elseif (ip(m, i) < 0) then
          stop 'negative power for a rigid mode'
        endif
      enddo
      tnz(i) = n
      tk(i) = ip(12, i)
      kmax = max(kmax, abs(ip(12, i)))
    enddo
    pmax = max(pmax, 1)
  end subroutine compress_terms

  !===========================================================================
  subroutine compute_xi(c, pnl, xi)
    real(rk), intent(in)  :: c(ndeg), pnl(nparam_eq_)
    real(rk), intent(out) :: xi(ndeg)
    real(rk) :: th1, th2, th3, t12, t23, t13, twopi, deg
    twopi = 2.0_rk * pi
    deg = pi / 180.0_rk
    xi(1) = 1.0_rk - exp(-pnl(6) * (c(1) - pnl(1)))
    xi(2) = 1.0_rk - exp(-pnl(7) * (c(2) - pnl(2)))
    xi(3) = 1.0_rk - exp(-pnl(8) * (c(3) - pnl(3)))
    xi(4) = 1.0_rk - exp(-pnl(8) * (c(5) - pnl(3)))
    xi(5) = 1.0_rk - exp(-pnl(8) * (c(4) - pnl(3)))
    xi(6) = c(6) - pnl(4) * deg
    xi(7) = c(7) - pnl(5) * deg
    xi(8) = c(9) - pnl(5) * deg
    xi(9) = c(8) - pnl(5) * deg
    th1 = modulo(c(10), twopi)
    th2 = c(12)
    th3 = c(11)
    t12 = modulo(th1 - th2, twopi)
    t23 = modulo(th2 - th3, twopi)
    t13 = modulo(twopi - t12 - t23, twopi)
    xi(10) = (2.0_rk * t23 - t12 - t13) / sqrt(6.0_rk)
    xi(11) = (t13 - t12) / sqrt(2.0_rk)
    if (tau_chart) then
      xi(12) = modulo(th1 + sqrt(2.0_rk) / 3.0_rk * xi(11), twopi)
    else
      xi(12) = (th1 + th2 + th3 - twopi) / 3.0_rk
    endif
  end subroutine compute_xi

  !===========================================================================
  subroutine sym_ops(xi, chi)
    ! C3v(M)-2 operations (ML_symmetry_transformation_XY3_IV of CH3OH_v10,
    ! identical as a set to TROVE's C3v(M)-2 for R-ALPHA132-THETA-S-TAU)
    real(rk), intent(in)  :: xi(ndeg)
    real(rk), intent(out) :: chi(ndeg, nsym)
    real(rk) :: a, b, p
    a = 0.5_rk
    b = 0.5_rk * sqrt(3.0_rk)
    p = 2.0_rk * pi / 3.0_rk
    chi(:, 1) = xi
    ! (132): r1'=r3, r2'=r1, r3'=r2 ; tau+2pi/3
    chi(1:2, 2) = xi(1:2); chi(6, 2) = xi(6)
    chi(3, 2) = xi(5); chi(4, 2) = xi(3); chi(5, 2) = xi(4)
    chi(7, 2) = xi(9); chi(8, 2) = xi(7); chi(9, 2) = xi(8)
    chi(10, 2) = -a * xi(10) - b * xi(11)
    chi(11, 2) = b * xi(10) - a * xi(11)
    chi(12, 2) = xi(12) + p
    ! (123): r1'=r2, r2'=r3, r3'=r1 ; tau-2pi/3
    chi(1:2, 3) = xi(1:2); chi(6, 3) = xi(6)
    chi(3, 3) = xi(4); chi(4, 3) = xi(5); chi(5, 3) = xi(3)
    chi(7, 3) = xi(8); chi(8, 3) = xi(9); chi(9, 3) = xi(7)
    chi(10, 3) = -a * xi(10) + b * xi(11)
    chi(11, 3) = -b * xi(10) - a * xi(11)
    chi(12, 3) = xi(12) - p
    ! (23)*
    chi(1:2, 4) = xi(1:2); chi(6, 4) = xi(6)
    chi(3, 4) = xi(3); chi(4, 4) = xi(5); chi(5, 4) = xi(4)
    chi(7, 4) = xi(7); chi(8, 4) = xi(9); chi(9, 4) = xi(8)
    chi(10, 4) = xi(10)
    chi(11, 4) = -xi(11)
    chi(12, 4) = -xi(12)
    ! (13)*
    chi(1:2, 5) = xi(1:2); chi(6, 5) = xi(6)
    chi(3, 5) = xi(5); chi(4, 5) = xi(4); chi(5, 5) = xi(3)
    chi(7, 5) = xi(9); chi(8, 5) = xi(8); chi(9, 5) = xi(7)
    chi(10, 5) = -a * xi(10) - b * xi(11)
    chi(11, 5) = -b * xi(10) + a * xi(11)
    chi(12, 5) = -xi(12) - p
    ! (12)*
    chi(1:2, 6) = xi(1:2); chi(6, 6) = xi(6)
    chi(3, 6) = xi(4); chi(4, 6) = xi(3); chi(5, 6) = xi(5)
    chi(7, 6) = xi(8); chi(8, 6) = xi(7); chi(9, 6) = xi(9)
    chi(10, 6) = -a * xi(10) + b * xi(11)
    chi(11, 6) = b * xi(10) + a * xi(11)
    chi(12, 6) = -xi(12) + p
  end subroutine sym_ops

  !===========================================================================
  subroutine basis_functions(c, pnl, dF)
    ! all symmetrised linear basis functions at one geometry
    real(rk), intent(in)  :: c(ndeg), pnl(nparam_eq_)
    real(rk), intent(out) :: dF(nlin)
    real(rk) :: xi(ndeg), chi(ndeg, nsym)
    real(rk) :: pw(0:pmax, 11, nsym), cs(0:kmax, nsym), sn(0:kmax, nsym)
    real(rk) :: s, t, c1, s1
    integer  :: g, p, k, i, j
    call compute_xi(c, pnl, xi)
    call sym_ops(xi, chi)
    do g = 1, nsym
      pw(0, :, g) = 1.0_rk
      do p = 1, pmax
        pw(p, :, g) = pw(p - 1, :, g) * chi(1:11, g)
      enddo
      c1 = cos(chi(12, g)); s1 = sin(chi(12, g))
      cs(0, g) = 1.0_rk; sn(0, g) = 0.0_rk
      do k = 1, kmax
        cs(k, g) = cs(k - 1, g) * c1 - sn(k - 1, g) * s1
        sn(k, g) = sn(k - 1, g) * c1 + cs(k - 1, g) * s1
      enddo
    enddo
    do i = 1, nlin
      s = 0.0_rk
      do g = 1, nsym
        if (tk(i) >= 0) then
          t = cs(tk(i), g)
        else
          t = sn(-tk(i), g)
        endif
        do j = 1, tnz(i)
          t = t * pw(tpow(j, i), tmode(j, i), g)
        enddo
        s = s + t
      enddo
      dF(i) = s / real(nsym, rk)
    enddo
  end subroutine basis_functions

  !===========================================================================
  subroutine build_jacobian(npts, coord, param, ivar, numpar, icol, nparam_eq, &
                            vcalc, rjacob, want_jac)
    integer, intent(in)   :: npts, numpar, icol(:), nparam_eq, ivar(:)
    real(rk), intent(in)  :: coord(:, :), param(:)
    real(rk), intent(out) :: vcalc(npts)
    real(rk), intent(inout) :: rjacob(:, :)
    logical, intent(in)   :: want_jac
    real(rk), allocatable :: dF(:), dFp(:), dFm(:)
    real(rk) :: pnl(nparam_eq_), pp(nparam_eq_), delta
    integer  :: n, l, i
    integer  :: lmap(size(param))
    lmap = 0
    do l = 1, numpar
      lmap(icol(l)) = l
    enddo
    pnl = param(1:nparam_eq_)
    !$omp parallel default(shared) private(n, l, i, dF, dFp, dFm, pp, delta)
    allocate (dF(nlin), dFp(nlin), dFm(nlin))
    !$omp do schedule(static)
    do n = 1, npts
      call basis_functions(coord(:, n), pnl, dF)
      vcalc(n) = dot_product(param(nparam_eq + 1:), dF)
      if (want_jac) then
        do i = nparam_eq + 1, size(param)
          if (lmap(i) > 0) rjacob(n, lmap(i)) = dF(i - nparam_eq)
        enddo
        do i = 1, nparam_eq
          if (lmap(i) > 0) then
            delta = 1.0e-3_rk * abs(param(i))
            if (delta < 1.0e-6_rk) delta = 1.0e-6_rk
            pp = pnl; pp(i) = pnl(i) + delta
            call basis_functions(coord(:, n), pp, dFp)
            pp(i) = pnl(i) - delta
            call basis_functions(coord(:, n), pp, dFm)
            rjacob(n, lmap(i)) = dot_product(param(nparam_eq + 1:), dFp - dFm) / (2.0_rk * delta)
          endif
        enddo
      endif
    enddo
    !$omp end do
    deallocate (dF, dFp, dFm)
    !$omp end parallel
  end subroutine build_jacobian

  !===========================================================================
  subroutine normal_equations(npts, numpar, rjacob, wt, eps, al, bl)
    ! A = J^T W J, b = J^T W eps   (rjacob is left unscaled)
#ifdef OPENBLAS64
    use ob64_lapack
#endif
    integer, intent(in)   :: npts, numpar
    real(rk), intent(in)  :: rjacob(:, :)
    real(rk), intent(in)  :: wt(npts), eps(npts)
    real(rk), intent(out) :: al(numpar, numpar), bl(numpar)
    real(rk), allocatable :: tmp(:)
    integer :: i, j
#ifdef OPENBLAS64
    real(rk), allocatable :: B(:, :), sw(:)
    integer(c_int64_t) :: n8, k8, one8
    allocate (B(npts, numpar), sw(npts))
    sw = sqrt(wt)
    !$omp parallel do schedule(static)
    do i = 1, numpar
      B(:, i) = rjacob(1:npts, i) * sw
    enddo
    !$omp end parallel do
    n8 = numpar; k8 = npts; one8 = 1
    call ob_dsyrk('L', 'T', n8, k8, 1.0_rk, B, k8, 0.0_rk, al, n8, 1_c_size_t, 1_c_size_t)
    sw = sw * eps
    call ob_dgemv('T', k8, n8, 1.0_rk, B, k8, sw, one8, 0.0_rk, bl, one8, 1_c_size_t)
    deallocate (B, sw)
    if (.false.) allocate (tmp(1))
#else
    !$omp parallel default(shared) private(i, j, tmp)
    allocate (tmp(npts))
    !$omp do schedule(dynamic, 2)
    do i = 1, numpar
      tmp = rjacob(1:npts, i) * wt
      do j = 1, i
        al(i, j) = dot_product(tmp, rjacob(1:npts, j))
      enddo
      bl(i) = dot_product(tmp, eps)
    enddo
    !$omp end do
    deallocate (tmp)
    !$omp end parallel
#endif
    do i = 1, numpar
      do j = i + 1, numpar
        al(i, j) = al(j, i)
      enddo
    enddo
  end subroutine normal_equations

  !===========================================================================
  subroutine solve_normal(n, al, bl, dx, ai, tol, ridge, ierr)
#ifdef OPENBLAS64
    use ob64_lapack
#endif
    ! Solve (A + ridge*diag(A)) dx = b on the diagonally scaled system.
    ! tol <= 0 : Cholesky (dpotrf/dpotrs/dpotri); ierr = index of the first
    !            non-positive pivot (singular parameter) on failure.
    ! tol  > 0 : eigen-decomposition pseudo-inverse (dsyevd), eigenvalues
    !            below tol*lambda_max are discarded.
    integer, intent(in)   :: n
    real(rk), intent(in)  :: al(n, n), bl(n), tol, ridge
    real(rk), intent(out) :: dx(n), ai(n, n)
    integer, intent(out)  :: ierr
    real(rk), allocatable :: s(:), a(:, :), w(:), work(:), b(:)
    integer, allocatable  :: iwork(:)
    integer :: i, j, info, lwork, liwork, nkeep
    allocate (s(n), a(n, n), b(n))
    ierr = 0
    do i = 1, n
      s(i) = sqrt(max(al(i, i), tiny(1.0_rk)))
    enddo
    do j = 1, n
      a(:, j) = al(:, j) / (s * s(j))
    enddo
    do i = 1, n
      a(i, i) = a(i, i) * (1.0_rk + ridge)
    enddo
    b = bl / s
    if (tol <= 0.0_rk) then
#ifdef OPENBLAS64
      block
        integer(c_int64_t) :: n8, info8, one8
        n8 = n; one8 = 1
        call ob_dpotrf('L', n8, a, n8, info8, 1_c_size_t)
        if (info8 > 0) then
          ierr = int(info8)
          return
        endif
        dx = b
        call ob_dpotrs('L', n8, one8, a, n8, dx, n8, info8, 1_c_size_t)
        call ob_dpotri('L', n8, a, n8, info8, 1_c_size_t)
      end block
#else
      call dpotrf('L', n, a, n, info)
      if (info > 0) then
        ierr = info
        return
      endif
      dx = b
      call dpotrs('L', n, 1, a, n, dx, n, info)
      call dpotri('L', n, a, n, info)
#endif
      do j = 1, n
        do i = 1, j - 1
          a(i, j) = a(j, i)
        enddo
      enddo
    else
      allocate (w(n))
      lwork = 1 + 6 * n + 2 * n * n
      liwork = 3 + 5 * n
      allocate (work(lwork), iwork(liwork))
      call dsyevd('V', 'L', n, a, n, w, work, lwork, iwork, liwork, info)
      if (info /= 0) stop 'dsyevd failed'
      nkeep = 0
      do i = 1, n
        if (w(i) > tol * w(n)) then
          w(i) = 1.0_rk / w(i)
          nkeep = nkeep + 1
        else
          w(i) = 0.0_rk
        endif
      enddo
      ! pseudo-inverse: V diag(w) V^T
      dx = matmul(a, w * matmul(transpose(a), b))
      do j = 1, n
        ai(:, j) = a(:, j) * w(j)
      enddo
      ai = matmul(ai, transpose(a))
      a = ai
      write (f_out, "('   eigen-cutoff: kept ',i6,' of ',i6,' directions')") nkeep, n
      deallocate (w, work, iwork)
    endif
    do j = 1, n
      ai(:, j) = a(:, j) / (s * s(j))
    enddo
    dx = dx / s
    deallocate (s, a, b)
  end subroutine solve_normal

  !===========================================================================
  subroutine normalise_weights(npts, wt, nused)
    integer, intent(in)     :: npts
    real(rk), intent(inout) :: wt(npts)
    integer, intent(out)    :: nused
    real(rk) :: wtsum
    nused = count(wt > 0.0_rk)
    wtsum = sum(wt) / real(max(nused, 1), rk)
    wt = wt / wtsum
  end subroutine normalise_weights

  subroutine skiplines(inpunit, n)
    integer, intent(in) :: n, inpunit
    character(len=80)   :: label
    integer             :: i0
    do i0 = 1, n
      read (inpunit, "(a80)") label
    enddo
  end subroutine skiplines

end module FITPES_ch3oh

program driver
  use FITPES_ch3oh
  implicit none
  call do_fitting
end program driver
