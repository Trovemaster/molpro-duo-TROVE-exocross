program pes_eval
  !
  ! Fast OpenMP evaluator of a fitted CH3OH PEF (CH3OH_fitpes .trove POTEN
  ! block) at points given in the fit/TROVE coordinates xi(1..12).
  !
  ! usage: pes_eval.x params.trove xi.bin v.bin
  !   xi.bin : stream file, int64 N followed by real64 xi(12,N)
  !   v.bin  : stream file, real64 V(N) (constant term V0 included)
  !
  implicit none
  integer, parameter :: rk = selected_real_kind(12, 25), ndeg = 12, nsym = 6
  real(rk), parameter :: pi = 3.14159265358979323846264338327950288_rk
  character(len=512) :: fpar, fxi, fv, line
  character(len=16)  :: lab
  integer :: nterm, i, m, n, ios, nnz, pmax, kmax, ipw(12), iv
  integer(8) :: npts8
  integer, allocatable :: tnz(:), tmode(:, :), tpow(:, :), tk(:)
  real(rk), allocatable :: c(:), xi(:, :), v(:)
  real(rk) :: val
  logical :: inpot

  call get_command_argument(1, fpar)
  call get_command_argument(2, fxi)
  call get_command_argument(3, fv)
  ! --- parameters -----------------------------------------------------------
  allocate (tnz(20000), tmode(11, 20000), tpow(11, 20000), tk(20000), c(20000))
  open (10, file=trim(fpar), status='old')
  nterm = 0; inpot = .false.; pmax = 1; kmax = 0
  do
    read (10, "(a)", iostat=ios) line
    if (ios /= 0) exit
    line = adjustl(line)
    if (line(1:5) == 'POTEN') then
      inpot = .true.; cycle
    endif
    if (.not. inpot) cycle
    if (line(1:3) == 'end' .or. line(1:3) == 'END') exit
    read (line, *) lab
    if (trim(lab) /= 'f' .and. trim(lab) /= 'h') cycle
    read (line, *) lab, ipw, iv, val
    nterm = nterm + 1
    nnz = 0
    do m = 1, 11
      if (ipw(m) > 0) then
        nnz = nnz + 1; tmode(nnz, nterm) = m; tpow(nnz, nterm) = ipw(m)
        pmax = max(pmax, ipw(m))
      endif
    enddo
    tnz(nterm) = nnz; tk(nterm) = ipw(12); kmax = max(kmax, abs(ipw(12)))
    c(nterm) = val
  enddo
  close (10)
  ! --- points ----------------------------------------------------------------
  open (11, file=trim(fxi), access='stream', form='unformatted', status='old')
  read (11) npts8
  allocate (xi(ndeg, npts8), v(npts8))
  read (11) xi
  close (11)
  !$omp parallel do schedule(dynamic, 64)
  do n = 1, int(npts8)
    v(n) = value_at(xi(:, n))
  enddo
  !$omp end parallel do
  open (12, file=trim(fv), access='stream', form='unformatted', status='replace')
  write (12) v
  close (12)

contains

  function value_at(x) result(s)
    real(rk), intent(in) :: x(ndeg)
    real(rk) :: s, chi(ndeg, nsym), pw(0:pmax, 11, nsym), cs(0:kmax, nsym), sn(0:kmax, nsym)
    real(rk) :: t, tt, a, b, p, c1, s1
    integer :: g, q, k, ii, j
    a = 0.5_rk; b = 0.5_rk * sqrt(3.0_rk); p = 2.0_rk * pi / 3.0_rk
    chi(:, 1) = x
    chi(1:2, 2) = x(1:2); chi(6, 2) = x(6)
    chi(3, 2) = x(5); chi(4, 2) = x(3); chi(5, 2) = x(4)
    chi(7, 2) = x(9); chi(8, 2) = x(7); chi(9, 2) = x(8)
    chi(10, 2) = -a * x(10) - b * x(11); chi(11, 2) = b * x(10) - a * x(11)
    chi(12, 2) = x(12) + p
    chi(1:2, 3) = x(1:2); chi(6, 3) = x(6)
    chi(3, 3) = x(4); chi(4, 3) = x(5); chi(5, 3) = x(3)
    chi(7, 3) = x(8); chi(8, 3) = x(9); chi(9, 3) = x(7)
    chi(10, 3) = -a * x(10) + b * x(11); chi(11, 3) = -b * x(10) - a * x(11)
    chi(12, 3) = x(12) - p
    chi(1:2, 4) = x(1:2); chi(6, 4) = x(6)
    chi(3, 4) = x(3); chi(4, 4) = x(5); chi(5, 4) = x(4)
    chi(7, 4) = x(7); chi(8, 4) = x(9); chi(9, 4) = x(8)
    chi(10, 4) = x(10); chi(11, 4) = -x(11); chi(12, 4) = -x(12)
    chi(1:2, 5) = x(1:2); chi(6, 5) = x(6)
    chi(3, 5) = x(5); chi(4, 5) = x(4); chi(5, 5) = x(3)
    chi(7, 5) = x(9); chi(8, 5) = x(8); chi(9, 5) = x(7)
    chi(10, 5) = -a * x(10) - b * x(11); chi(11, 5) = -b * x(10) + a * x(11)
    chi(12, 5) = -x(12) - p
    chi(1:2, 6) = x(1:2); chi(6, 6) = x(6)
    chi(3, 6) = x(4); chi(4, 6) = x(3); chi(5, 6) = x(5)
    chi(7, 6) = x(8); chi(8, 6) = x(7); chi(9, 6) = x(9)
    chi(10, 6) = -a * x(10) + b * x(11); chi(11, 6) = b * x(10) + a * x(11)
    chi(12, 6) = -x(12) + p
    do g = 1, nsym
      pw(0, :, g) = 1.0_rk
      do q = 1, pmax
        pw(q, :, g) = pw(q - 1, :, g) * chi(1:11, g)
      enddo
      c1 = cos(chi(12, g)); s1 = sin(chi(12, g))
      cs(0, g) = 1.0_rk; sn(0, g) = 0.0_rk
      do k = 1, kmax
        cs(k, g) = cs(k - 1, g) * c1 - sn(k - 1, g) * s1
        sn(k, g) = sn(k - 1, g) * c1 + cs(k - 1, g) * s1
      enddo
    enddo
    s = 0.0_rk
    do ii = 1, nterm
      tt = 0.0_rk
      do g = 1, nsym
        if (tk(ii) >= 0) then
          t = cs(tk(ii), g)
        else
          t = sn(-tk(ii), g)
        endif
        do j = 1, tnz(ii)
          t = t * pw(tpow(j, ii), tmode(j, ii), g)
        enddo
        tt = tt + t
      enddo
      s = s + c(ii) * tt / 6.0_rk
    enddo
  end function value_at

end program pes_eval
