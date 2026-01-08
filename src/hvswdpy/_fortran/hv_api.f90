! ======================================================================================
! HV-DFA API module
! Fortran API (hv_compute, hv_components) wrapping the original HV-DFA routines.
! Algorithms follow the original implementation (Sanchez-Sesma et al., HV-INV team).
! 
! Modifications and API by: Shihao Yuan (syuan@mines.edu)
! ======================================================================================
module hv_dfa_api
  ! Public API for computing H/V using the legacy core algorithms without CLI or file I/O
  use TYPES
  use Globales
  use Marc
  implicit none
  private
  public :: hv_compute
  public :: hv_compute_components
  public :: hv_compute_f2py
  public :: hv_compute_components_f2py
  public :: hv_compute_dispersion_f2py
  public :: hv_components
  public :: hv_from_terms

  ! Concise generic alias for components API
  interface hv_components
    module procedure hv_compute_components
  end interface hv_components

  ! Interface to legacy dispersion routine
  interface
    logical function DISPERSION(VALUES, VALID) result(RETORNO)
      real,    intent(inout), dimension(:), target :: VALUES
      logical, intent(inout), dimension(:), target :: VALID
    end function DISPERSION
  end interface

contains

  subroutine hv_compute(frequencies_hz, vp_in, vs_in, rho_in, thickness_in, &
                        hv_out, &
                        n_rayleigh_modes, n_love_modes, precision_percent, &
                        nks, sh_damp, psv_damp, status)
    ! Computes the H/V spectral ratio for a layered halfspace model
    ! Inputs
    real(long_float), intent(in)  :: frequencies_hz(:)          ! size nf
    real(long_float), intent(in)  :: vp_in(:)                    ! size n_layers
    real(long_float), intent(in)  :: vs_in(:)                    ! size n_layers
    real(long_float), intent(in)  :: rho_in(:)                   ! size n_layers
    real(long_float), intent(in)  :: thickness_in(:)             ! size n_layers-1 (last is halfspace)
    ! Outputs
    real(long_float), intent(out) :: hv_out(size(frequencies_hz))
    ! Optional controls
    integer, intent(in),  optional :: n_rayleigh_modes
    integer, intent(in),  optional :: n_love_modes
    real,    intent(in),  optional :: precision_percent          ! default 1e-4 (% units, legacy semantics)
    integer, intent(in),  optional :: nks                        ! number of k samples for BW integrals (0 to skip)
    real(long_float), intent(in), optional :: sh_damp, psv_damp  ! dimensionless attenuation factors
    integer, intent(out), optional :: status                     ! 0=ok, non-zero means a failure in dispersion

    ! Locals matching legacy globals
    integer :: nf, nlayers
    integer :: nm_rayleigh, nm_love, nm_rl
    integer :: i, imode

    integer :: nks_eff
    real    :: prec_eff
    real(long_float), pointer :: G1(:,:), IMG2(:,:), IMG3(:,:)
    real(long_float), allocatable :: IMG11_pihalf(:), IMG33(:)
    real(long_float), allocatable :: IMVV(:), IMHPSV(:), IMHSH(:)
    integer, allocatable :: OFFSET_R(:), OFFSET_L(:)
    logical :: ok

    ! Pointing slices
    ! Defaults
    nf = size(frequencies_hz)
    nlayers = size(vp_in)
    if (present(n_rayleigh_modes)) then
      nm_rayleigh = n_rayleigh_modes
    else
      nm_rayleigh = 0
    end if
    if (present(n_love_modes)) then
      nm_love = n_love_modes
    else
      nm_love = 0
    end if
    if (present(precision_percent)) then
      prec_eff = precision_percent
    else
      prec_eff = 1.0e-4
    end if
    if (present(nks)) then
      nks_eff = nks
    else
      nks_eff = 0
    end if


    ! Short-circuit: if no Rayleigh, no Love, and no body waves requested, return zeros
    if (nm_rayleigh == 0 .and. nm_love == 0 .and. nks_eff == 0) then
      hv_out = 0.0_long_float
      if (present(status)) status = 0
      return
    end if

    ! Initialize legacy modules (safest: reset everything we touch)
    UNIT_RDPGS = -1
    if (present(sh_damp)) then
      SHDAMP = sh_damp
    else
      SHDAMP = 1.0d-5
    end if
    if (present(psv_damp)) then
      PSVDAMP = psv_damp
    else
      PSVDAMP = 1.0d-5
    end if

    ! Build circular frequency vector X
    G_NX = nf
    if (allocated(X)) deallocate(X)
    allocate(X(G_NX))
    do i = 1, G_NX
      X(i) = TWOPI * frequencies_hz(i)
    end do

    ! Load model into legacy globals
    NCAPAS = nlayers
    if (allocated(ALFA)) deallocate(ALFA)
    if (allocated(BTA))  deallocate(BTA)
    if (allocated(RHO))  deallocate(RHO)
    if (allocated(H))    deallocate(H)
    if (allocated(MU))   deallocate(MU)
    if (allocated(G_SLOWS)) deallocate(G_SLOWS)
    if (allocated(G_SLOWP)) deallocate(G_SLOWP)
    allocate(ALFA(NCAPAS), BTA(NCAPAS), RHO(NCAPAS), MU(NCAPAS))
    allocate(G_SLOWS(NCAPAS), G_SLOWP(NCAPAS))
    allocate(H(NCAPAS-1))
    ALFA = vp_in
    BTA  = vs_in
    RHO  = rho_in
    H    = thickness_in
    MU   = BTA**2.0_long_float * RHO

    ! Compute limits and inversion flags
    call set_model_parameters_api()

    ! Prepare dispersion search controls
    G_PRECISION = prec_eff * 1.0e-2         ! convert % to fraction as in legacy
    G_DX = 0.1
    G_DXTYPE = .false.

    ! Rayleigh dispersion
    if (nm_rayleigh > 0) then
      G_NMODES = nm_rayleigh
      if (allocated(VALUES_R)) deallocate(VALUES_R)
      if (allocated(VALID_R))  deallocate(VALID_R)
      allocate(VALUES_R(G_NX*G_NMODES), VALID_R(G_NX*G_NMODES))
      VALUES_R = 0.0
      VALID_R  = .false.
      ISRAYLEIGH = .true.
      ok = DISPERSION(VALUES_R, VALID_R)
      if (.not. ok) then
        if (present(status)) status = 2
        hv_out = 0.0_long_float
        return
      end if
    else
      ok = .true.
    end if

    ! Love dispersion
    if (ok .and. nm_love > 0) then
      G_NMODES = nm_love
      if (allocated(VALUES_L)) deallocate(VALUES_L)
      if (allocated(VALID_L))  deallocate(VALID_L)
      allocate(VALUES_L(G_NX*G_NMODES), VALID_L(G_NX*G_NMODES))
      VALUES_L = 0.0
      VALID_L  = .false.
      ISRAYLEIGH = .false.
      ok = DISPERSION(VALUES_L, VALID_L) .and. ok
      if (.not. ok) then
        if (present(status)) status = 3
        hv_out = 0.0_long_float
        return
      end if
    end if

    if (present(status)) status = merge(0, 1, ok)

    ! Prepare Green function accumulators
    nm_rl = max(nm_rayleigh, nm_love)
    if (nm_rl > 0) then
      allocate(G1(G_NX, nm_rl), IMG2(G_NX, nm_rl), IMG3(G_NX, nm_rl))
      G1 = 0.0_long_float; IMG2 = 0.0_long_float; IMG3 = 0.0_long_float
    else
      nullify(G1); nullify(IMG2); nullify(IMG3)
    end if
    allocate(IMG11_pihalf(G_NX), IMG33(G_NX))
    IMG11_pihalf = 0.0_long_float; IMG33 = 0.0_long_float

    ! Rayleigh modal contributions
    if (nm_rayleigh > 0) then
      allocate(OFFSET_R(nm_rayleigh))
      call offsets_from_valid(VALID_R, nm_rayleigh, OFFSET_R)
      do imode = 1, nm_rayleigh
        if (OFFSET_R(imode) >= 0 .and. OFFSET_R(imode) < G_NX) then
          if (G_NX - OFFSET_R(imode) <= 0) cycle
          call GR(G1(:, imode), IMG3(:, imode), imode, OFFSET_R(imode))
          IMG11_pihalf(OFFSET_R(imode)+1:G_NX) = IMG11_pihalf(OFFSET_R(imode)+1:G_NX) + &
               0.5_long_float * real(G1(OFFSET_R(imode)+1:G_NX, imode), kind=long_float)
          IMG33(OFFSET_R(imode)+1:G_NX) = IMG33(OFFSET_R(imode)+1:G_NX) + &
               IMG3(OFFSET_R(imode)+1:G_NX, imode)
        end if
      end do
      deallocate(OFFSET_R)
    end if

    ! Love modal contributions
    if (nm_love > 0) then
      allocate(OFFSET_L(nm_love))
      call offsets_from_valid(VALID_L, nm_love, OFFSET_L)
      do imode = 1, nm_love
        if (OFFSET_L(imode) >= 0 .and. OFFSET_L(imode) < G_NX) then
          call GL(IMG2(:, imode), imode, OFFSET_L(imode))
          IMG11_pihalf(OFFSET_L(imode)+1:G_NX) = IMG11_pihalf(OFFSET_L(imode)+1:G_NX) - &
               0.5_long_float * IMG2(OFFSET_L(imode)+1:G_NX, imode)
        end if
      end do
      deallocate(OFFSET_L)
    end if

    ! Body-wave integrals (BW)
    allocate(IMVV(G_NX), IMHPSV(G_NX), IMHSH(G_NX))
    if (nks_eff > 0) then
      do i = 1, G_NX
        call BWR(IMVV(i), IMHPSV(i), IMHSH(i), nks_eff, real(X(i), kind=kind(0.0)))
      end do
    else
      IMVV  = 0.0_long_float
      IMHPSV= 0.0_long_float
      IMHSH = 0.0_long_float
    end if

    ! Final H/V
    do i = 1, G_NX
      hv_out(i) = hv_from_terms(IMG11_pihalf(i), IMHPSV(i), IMHSH(i), IMG33(i), IMVV(i))
    end do

    ! Cleanup local arrays (module arrays kept allocated for potential reuse by caller if multiple calls)
    if (associated(G1))   deallocate(G1)
    if (associated(IMG2)) deallocate(IMG2)
    if (associated(IMG3)) deallocate(IMG3)
    deallocate(IMG11_pihalf, IMG33, IMVV, IMHPSV, IMHSH)

  end subroutine hv_compute


  subroutine hv_compute_components(frequencies_hz, vp_in, vs_in, rho_in, thickness_in, &
                                   img11_pihalf_total, img33_total, &
                                   img11_pihalf_rayleigh, img11_pihalf_love, &
                                   imvv, imhpsv, imhsh, &
                                   n_rayleigh_modes, n_love_modes, precision_percent, nks, status)
    ! Returns per-frequency totals of each contribution used in H/V
    real(long_float), intent(in)  :: frequencies_hz(:)
    real(long_float), intent(in)  :: vp_in(:), vs_in(:), rho_in(:), thickness_in(:)
    real(long_float), intent(out) :: img11_pihalf_total(:)  ! sum over Rayleigh (+) and Love (-)
    real(long_float), intent(out) :: img33_total(:)         ! Rayleigh only
    real(long_float), intent(out) :: img11_pihalf_rayleigh(:)
    real(long_float), intent(out) :: img11_pihalf_love(:)
    real(long_float), intent(out) :: imvv(:), imhpsv(:), imhsh(:)
    integer, intent(in),  optional :: n_rayleigh_modes, n_love_modes, nks
    real,    intent(in),  optional :: precision_percent
    integer, intent(out), optional :: status

    integer :: nf, nlayers
    integer :: nm_rayleigh, nm_love
    integer :: i, imode
    integer :: nks_eff
    real    :: prec_eff
    real(long_float), pointer :: G1(:,:), IMG2(:,:), IMG3(:,:)
    integer, allocatable :: OFFSET_R(:), OFFSET_L(:)
    logical :: ok

    nf = size(frequencies_hz)
    nlayers = size(vp_in)
    if (present(n_rayleigh_modes))then; nm_rayleigh = n_rayleigh_modes; else; nm_rayleigh = 0; end if
    if (present(n_love_modes))     then; nm_love     = n_love_modes;     else; nm_love     = 0; end if
    if (present(precision_percent))then; prec_eff    = precision_percent; else; prec_eff    = 1.0e-4; end if
    if (present(nks))              then; nks_eff     = nks;              else; nks_eff     = 0; end if

    UNIT_RDPGS = -1
    SHDAMP = 1.0d-5
    PSVDAMP = 1.0d-5

    G_NX = nf
    if (allocated(X)) deallocate(X)
    allocate(X(G_NX))
    do i=1,G_NX; X(i)=TWOPI*frequencies_hz(i); end do

    NCAPAS = nlayers
    if (allocated(ALFA)) deallocate(ALFA)
    if (allocated(BTA))  deallocate(BTA)
    if (allocated(RHO))  deallocate(RHO)
    if (allocated(H))    deallocate(H)
    if (allocated(MU))   deallocate(MU)
    if (allocated(G_SLOWS)) deallocate(G_SLOWS)
    if (allocated(G_SLOWP)) deallocate(G_SLOWP)
    allocate(ALFA(NCAPAS), BTA(NCAPAS), RHO(NCAPAS), MU(NCAPAS))
    allocate(G_SLOWS(NCAPAS), G_SLOWP(NCAPAS))
    allocate(H(NCAPAS-1))
    ALFA = vp_in; BTA = vs_in; RHO = rho_in; H = thickness_in
    MU = BTA**2.0_long_float * RHO

    call set_model_parameters_api()
    G_PRECISION = prec_eff * 1.0e-2
    G_DX = 0.1
    G_DXTYPE = .false.

    img11_pihalf_total = 0.0_long_float
    img33_total = 0.0_long_float
    img11_pihalf_rayleigh = 0.0_long_float
    img11_pihalf_love = 0.0_long_float
    imvv = 0.0_long_float; imhpsv = 0.0_long_float; imhsh = 0.0_long_float

    ! Rayleigh
    if (nm_rayleigh > 0) then
      G_NMODES = nm_rayleigh
      if (allocated(VALUES_R)) deallocate(VALUES_R)
      if (allocated(VALID_R))  deallocate(VALID_R)
      allocate(VALUES_R(G_NX*G_NMODES), VALID_R(G_NX*G_NMODES))
      VALUES_R = 0.0; VALID_R = .false.
      ISRAYLEIGH = .true.
      ok = DISPERSION(VALUES_R, VALID_R)
      allocate(G1(G_NX,nm_rayleigh), IMG3(G_NX,nm_rayleigh))
      G1=0.0; IMG3=0.0
      allocate(OFFSET_R(nm_rayleigh))
      call offsets_from_valid(VALID_R, nm_rayleigh, OFFSET_R)
      do imode=1,nm_rayleigh
        if (OFFSET_R(imode) >= 0) then
          call GR(G1(:,imode), IMG3(:,imode), imode, OFFSET_R(imode))
          img11_pihalf_rayleigh(OFFSET_R(imode)+1:G_NX) = img11_pihalf_rayleigh(OFFSET_R(imode)+1:G_NX) + &
            0.5_long_float * real(G1(OFFSET_R(imode)+1:G_NX, imode), kind=long_float)
          img33_total(OFFSET_R(imode)+1:G_NX) = img33_total(OFFSET_R(imode)+1:G_NX) + &
            IMG3(OFFSET_R(imode)+1:G_NX, imode)
        end if
      end do
      deallocate(OFFSET_R)
    else
      ok = .true.
    end if

    ! Love
    if (ok .and. nm_love > 0) then
      G_NMODES = nm_love
      if (allocated(VALUES_L)) deallocate(VALUES_L)
      if (allocated(VALID_L))  deallocate(VALID_L)
      allocate(VALUES_L(G_NX*G_NMODES), VALID_L(G_NX*G_NMODES))
      VALUES_L = 0.0; VALID_L = .false.
      ISRAYLEIGH = .false.
      ok = DISPERSION(VALUES_L, VALID_L) .and. ok
      allocate(IMG2(G_NX,nm_love))
      IMG2=0.0
      allocate(OFFSET_L(nm_love))
      call offsets_from_valid(VALID_L, nm_love, OFFSET_L)
      do imode=1,nm_love
        if (OFFSET_L(imode) >= 0) then
          call GL(IMG2(:,imode), imode, OFFSET_L(imode))
          img11_pihalf_love(OFFSET_L(imode)+1:G_NX) = img11_pihalf_love(OFFSET_L(imode)+1:G_NX) + &
            0.5_long_float * IMG2(OFFSET_L(imode)+1:G_NX, imode)
        end if
      end do
      deallocate(OFFSET_L)
    end if

    img11_pihalf_total = img11_pihalf_rayleigh - img11_pihalf_love

    ! Body waves
    if (nks_eff > 0) then
      do i=1,G_NX
        call BWR(imvv(i), imhpsv(i), imhsh(i), nks_eff, real(X(i),kind=kind(0.0)))
      end do
    end if

    if (present(status)) status = merge(0, 1, ok)

    if (associated(G1))   deallocate(G1)
    if (associated(IMG2)) deallocate(IMG2)
    if (associated(IMG3)) deallocate(IMG3)

  end subroutine hv_compute_components


  subroutine hv_compute_f2py(nf, nl, frequencies_hz, vp_in, vs_in, rho_in, thickness_in, &
                             hv_out, &
                             n_rayleigh_modes, n_love_modes, precision_percent, &
                             nks, sh_damp, psv_damp, status)
    ! f2py-friendly wrapper with explicit sizes to avoid assumed-shape ABI issues
    integer, intent(in) :: nf, nl
    real(long_float), intent(in)  :: frequencies_hz(nf)
    real(long_float), intent(in)  :: vp_in(nl), vs_in(nl), rho_in(nl)
    real(long_float), intent(in)  :: thickness_in(nl-1)
    real(long_float), intent(out) :: hv_out(nf)
    integer, intent(in),  optional :: n_rayleigh_modes, n_love_modes, nks
    real,    intent(in),  optional :: precision_percent
    real(long_float), intent(in), optional :: sh_damp, psv_damp
    integer, intent(out), optional :: status
    call hv_compute(frequencies_hz, vp_in, vs_in, rho_in, thickness_in, &
                    hv_out, n_rayleigh_modes, n_love_modes, precision_percent, &
                    nks, sh_damp, psv_damp, status)
  end subroutine hv_compute_f2py


  subroutine hv_compute_components_f2py(nf, nl, frequencies_hz, vp_in, vs_in, rho_in, thickness_in, &
                                        img11_pihalf_total, img33_total, &
                                        img11_pihalf_rayleigh, img11_pihalf_love, &
                                        imvv, imhpsv, imhsh, &
                                        n_rayleigh_modes, n_love_modes, precision_percent, nks, status)
    ! f2py-friendly wrapper with explicit sizes to avoid assumed-shape ABI issues
    integer, intent(in) :: nf, nl
    real(long_float), intent(in)  :: frequencies_hz(nf)
    real(long_float), intent(in)  :: vp_in(nl), vs_in(nl), rho_in(nl)
    real(long_float), intent(in)  :: thickness_in(nl-1)
    real(long_float), intent(out) :: img11_pihalf_total(nf), img33_total(nf)
    real(long_float), intent(out) :: img11_pihalf_rayleigh(nf), img11_pihalf_love(nf)
    real(long_float), intent(out) :: imvv(nf), imhpsv(nf), imhsh(nf)
    integer, intent(in),  optional :: n_rayleigh_modes, n_love_modes, nks
    real,    intent(in),  optional :: precision_percent
    integer, intent(out), optional :: status
    call hv_compute_components(frequencies_hz, vp_in, vs_in, rho_in, thickness_in, &
                               img11_pihalf_total, img33_total, &
                               img11_pihalf_rayleigh, img11_pihalf_love, &
                               imvv, imhpsv, imhsh, &
                               n_rayleigh_modes, n_love_modes, precision_percent, nks, status)
  end subroutine hv_compute_components_f2py


  subroutine hv_compute_dispersion_f2py(nf, nl, frequencies_hz, vp_in, vs_in, rho_in, thickness_in, &
                                        rayleigh_slowness, rayleigh_valid, love_slowness, love_valid, &
                                        n_rayleigh_modes, n_love_modes, precision_percent, status)
    !f2py intent(in) nf, nl, frequencies_hz, vp_in, vs_in, rho_in, thickness_in, n_rayleigh_modes, n_love_modes, precision_percent
    !f2py intent(out) rayleigh_slowness, rayleigh_valid, love_slowness, love_valid, status
    ! f2py-friendly wrapper for dispersion curve calculation
    integer, intent(in) :: nf, nl, n_rayleigh_modes, n_love_modes
    real(long_float), intent(in)  :: frequencies_hz(nf)
    real(long_float), intent(in)  :: vp_in(nl), vs_in(nl), rho_in(nl)
    real(long_float), intent(in)  :: thickness_in(nl-1)
    real(long_float), intent(out) :: rayleigh_slowness(nf, n_rayleigh_modes), love_slowness(nf, n_love_modes)
    logical, intent(out) :: rayleigh_valid(nf, n_rayleigh_modes), love_valid(nf, n_love_modes)
    real,    intent(in),  optional :: precision_percent
    integer, intent(out), optional :: status
    
    ! Local variables
    integer :: nm_rayleigh_eff, nm_love_eff, prec_eff
    logical :: ok
    integer :: i, j, imode
    real(long_float) :: NORMT, NORMH, NORMV, NORMD, NORMG, NORMM
    
    ! Set default values
    nm_rayleigh_eff = n_rayleigh_modes
    nm_love_eff = n_love_modes
    prec_eff = merge(int(precision_percent), 1, present(precision_percent))
    
    ! Initialize output arrays
    rayleigh_slowness = 0.0_long_float
    love_slowness = 0.0_long_float
    rayleigh_valid = .false.
    love_valid = .false.
    
    ! Set up global variables
    UNIT_RDPGS = -1
    SHDAMP = 1.0d-5
    PSVDAMP = 1.0d-5
    G_NX = nf
    NCAPAS = nl
    if (allocated(ALFA)) deallocate(ALFA)
    if (allocated(BTA))  deallocate(BTA)
    if (allocated(RHO))  deallocate(RHO)
    if (allocated(H))    deallocate(H)
    if (allocated(MU))   deallocate(MU)
    if (allocated(G_SLOWS)) deallocate(G_SLOWS)
    if (allocated(G_SLOWP)) deallocate(G_SLOWP)
    allocate(ALFA(NCAPAS), BTA(NCAPAS), RHO(NCAPAS), MU(NCAPAS))
    allocate(G_SLOWS(NCAPAS), G_SLOWP(NCAPAS))
    allocate(H(NCAPAS-1))
    ALFA = vp_in; BTA = vs_in; RHO = rho_in; H = thickness_in
    MU = BTA**2.0_long_float * RHO
    
    ! Do NOT normalize for dispersion API to match original CLI default
    ! (HV.f90 sets NORMALIZE = .FALSE. by default.)
    NORMT = 1.0_long_float
    NORMH = 1.0_long_float
    NORMV = 1.0_long_float
    NORMD = 1.0_long_float
    NORMG = 1.0_long_float
    NORMM = 1.0_long_float

    call set_model_parameters_api()
    G_PRECISION = prec_eff * 1.0e-2
    G_DX = 0.1
    G_DXTYPE = .false.
    
    ! Set up frequency array (no normalization)
    if (allocated(X)) deallocate(X)
    allocate(X(G_NX))
    do i=1,G_NX; X(i)=TWOPI*frequencies_hz(i); end do
    
    ! Initialize additional global variables that might be needed
    G_OMEGA = 0.05_long_float  ! Initial value for polarity calculation
    G_POLARITY = 0.0_long_float
    G_X2 = 0.0_long_float
    
    ok = .true.
    
    ! Rayleigh dispersion curves
    if (nm_rayleigh_eff > 0) then
      G_NMODES = nm_rayleigh_eff
      if (allocated(VALUES_R)) deallocate(VALUES_R)
      if (allocated(VALID_R))  deallocate(VALID_R)
      allocate(VALUES_R(G_NX*G_NMODES), VALID_R(G_NX*G_NMODES))
      VALUES_R = 0.0; VALID_R = .false.
      ISRAYLEIGH = .true.
      ok = DISPERSION(VALUES_R, VALID_R)
      
      ! Copy results to output arrays (with normalization)
      do imode = 1, nm_rayleigh_eff
        do i = 1, G_NX
          j = (imode-1)*G_NX + i
          rayleigh_slowness(i, imode) = VALUES_R(j)
          rayleigh_valid(i, imode) = VALID_R(j)
        end do
      end do
      
    end if
    
    ! Love dispersion curves
    if (ok .and. nm_love_eff > 0) then
      G_NMODES = nm_love_eff
      if (allocated(VALUES_L)) deallocate(VALUES_L)
      if (allocated(VALID_L))  deallocate(VALID_L)
      allocate(VALUES_L(G_NX*G_NMODES), VALID_L(G_NX*G_NMODES))
      VALUES_L = 0.0; VALID_L = .false.
      ISRAYLEIGH = .false.
      ok = DISPERSION(VALUES_L, VALID_L) .and. ok
      
      ! Copy results to output arrays (with normalization)
      do imode = 1, nm_love_eff
        do i = 1, G_NX
          j = (imode-1)*G_NX + i
          love_slowness(i, imode) = VALUES_L(j)
          love_valid(i, imode) = VALID_L(j)
        end do
      end do
    end if
    
    if (present(status)) status = merge(0, 1, ok)
    
  end subroutine hv_compute_dispersion_f2py


  pure function hv_from_terms(img11_pihalf, imhpsv, imhsh, img33, imvv) result(hv)
    real(long_float), intent(in) :: img11_pihalf, imhpsv, imhsh, img33, imvv
    real(long_float) :: hv
    real(long_float) :: num, den
    num = 2.0_long_float * ( img11_pihalf + imhpsv + imhsh )
    den = img33 + imvv
    if (abs(den) <= tiny(den)) then
      hv = 0.0_long_float
    else
      if (num < 0.0_long_float .and. den < 0.0_long_float) then
        hv = sqrt( abs(num/den) )
      else
        hv = sqrt( max(0.0_long_float, num/den) )
      end if
    end if
  end function hv_from_terms

  subroutine offsets_from_valid(valid_flat, nmodes, offsets)
    ! Computes the first valid index (offset) per mode along the frequency axis
    logical, intent(in) :: valid_flat(:)
    integer, intent(in) :: nmodes
    integer, intent(out):: offsets(nmodes)
    integer :: i_mode, i, start_idx, end_idx
    offsets = -1
    do i_mode = 1, nmodes
      start_idx = (i_mode-1)*G_NX + 1
      end_idx   = i_mode*G_NX
      do i = start_idx, end_idx
        if (valid_flat(i)) then
          offsets(i_mode) = mod(i-1, G_NX)
          exit
        end if
      end do
    end do
  end subroutine offsets_from_valid


  subroutine set_model_parameters_api()
    ! Copy of SET_MODEL_PARAMETERS from HV.f90, kept local to avoid depending on the main program unit
    use Marc, only : NCAPAS, G_MAXRAYLEIGHSLOWNESS, G_SLOWS, G_SLOWP, ALFA, BTA, &
                     G_VELOCITYINVERSION, G_SLOWSMIN, G_SLOWSMAX, HALFSPACE_RAYLEIGH
    implicit none
    integer :: i, imaxslow
    G_SLOWS = 1.0/real(BTA)
    G_SLOWP = 1.0/real(ALFA)
    imaxslow = maxloc(G_SLOWS, 1)
    G_MAXRAYLEIGHSLOWNESS = real(HALFSPACE_RAYLEIGH(imaxslow))
    G_SLOWSMAX = G_SLOWS(imaxslow)
    G_SLOWSMIN = G_SLOWS(NCAPAS)
    do i = 2, NCAPAS
      if (G_SLOWS(i) > G_SLOWS(i-1) .or. G_SLOWP(i) > G_SLOWP(i-1)) then
        G_VELOCITYINVERSION = i
        return
      end if
    end do
    G_VELOCITYINVERSION = -1
  end subroutine set_model_parameters_api

end module hv_dfa_api


