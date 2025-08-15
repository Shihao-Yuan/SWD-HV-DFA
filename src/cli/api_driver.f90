! ======================================================================================
! CLI for HV-DFA
! Reads 'model.txt', computes H/V via hv_compute.
! Original algorithm by HV-INV team; API driver by Shihao Yuan (syuan@mines.edu).
! ======================================================================================
program api_driver
  use hv_dfa_api
  use TYPES
  implicit none

  character(len=256) :: model_path
  integer :: n_layers, i, ios
  real(long_float), allocatable :: vp(:), vs(:), rho(:), h(:)
  real(long_float), allocatable :: f(:), hv(:)
  
  real(long_float) :: t, vp_i, vs_i, rho_i
  integer :: nf
  integer :: status

  ! Settings
  model_path = 'model.txt'
  nf = 100

  ! Read model file: first line n_layers; then n_layers-1 with h,vp,vs,rho; last line 0,vp,vs,rho
  open(unit=10, file=trim(model_path), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Failed to open ', trim(model_path)
    stop 1
  end if
  read(10,*,iostat=ios) n_layers
  if (ios /= 0 .or. n_layers < 1) then
    print *, 'Invalid model header'
    stop 1
  end if
  allocate(vp(n_layers), vs(n_layers), rho(n_layers), h(n_layers-1))
  do i = 1, n_layers-1
    read(10,*,iostat=ios) t, vp_i, vs_i, rho_i
    if (ios /= 0) then
      print *, 'Invalid layer line at ', i
      stop 1
    end if
    h(i)   = t
    vp(i)  = vp_i
    vs(i)  = vs_i
    rho(i) = rho_i
  end do
  read(10,*,iostat=ios) t, vp_i, vs_i, rho_i
  if (ios /= 0) then
    print *, 'Invalid halfspace line'
    stop 1
  end if
  vp(n_layers)  = vp_i
  vs(n_layers)  = vs_i
  rho(n_layers) = rho_i
  close(10)

  allocate(f(nf), hv(nf))
  call build_logspace(0.1_long_float, 10.0_long_float, nf, f)

  ! Final H/V only
  call hv_compute(f, vp, vs, rho, h, hv, n_rayleigh_modes=20, n_love_modes=20, precision_percent=1.0e-4, nks=200, status=status)
  call ensure_results_dir()
  open(unit=13, file='results/HV_api.dat', status='replace', action='write')
  do i = 1, nf
    write(13,'(F12.6,1X,F12.6)') f(i), hv(i)
  end do
  close(13)
  print *, 'Wrote results/HV_api.dat'

contains

  subroutine build_logspace(a, b, n, out)
    real(long_float), intent(in) :: a, b
    integer, intent(in) :: n
    real(long_float), intent(out) :: out(n)
    integer :: j
    real(long_float) :: la, lb
    if (n == 1) then
      out(1) = a
      return
    end if
    la = log10(a)
    lb = log10(b)
    do j = 1, n
      out(j) = 10.0_long_float**(la + (j-1)*(lb-la)/real(n-1,kind=long_float))
    end do
  end subroutine build_logspace

  subroutine ensure_results_dir()
    implicit none
    call execute_command_line('mkdir -p results')
  end subroutine ensure_results_dir

end program api_driver


