module dftbp_common_accuracy

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: mc = 50

end module dftbp_common_accuracy

!> Connects to an external energy/hamiltonian model or provides a dummy external model if no
!> external library is linked
module dftbp_externalmodel
  use iso_c_binding, only : c_int, c_char, c_bool, c_null_char, c_ptr, c_double, c_loc,&
      & c_f_pointer, c_intptr_t
  use dftbp_common_accuracy, only : dp, mc
  implicit none

  private
  public :: TExtModelProvides, TExtModel, getExtModelCapabilities, externalModel_init

  !> Return type for capabilities of external electronic structure/energy model
  type TExtModelProvides

    !> Name of external model
    character(:), allocatable :: modelName

    !> Is a hamiltonian provided
    logical :: gives_ham = .false.

    !> Is an overlap provdided, or should a diagonal matrix be substituted
    logical :: gives_overlap = .false.

    !> Are double-counting energies provided (i.e. non-bandstructure)
    logical :: gives_dc_energy = .false.

    !> Order of derivatives returned by the model
    integer :: model_derivative_order = 0

    !> Is the model self-consistently evaluated
    logical :: is_self_consistent = .false.

    !> Number of spin channels (0 for spin free)
    integer :: nSpin = 0

  end type TExtModelProvides


  !> Type for instance of model
  type TExtModel

    !> External model's capabilities
    type(TExtModelProvides) :: capabilities

    !> Distance cutoff for interactions and local environment
    real(dp) :: interactionCutoff

    !> Distance cutoff for surrounding environment
    real(dp) :: environmentCutoff

    !> Distance atom atom sites such that environmentCutoff can be checked around bonds
    real(dp) :: maxCutoff

    !> Number of shells for each atom
    integer, allocatable :: nShellsOnSpecies(:)

    !> Angular momentum of each atomic shell
    integer, allocatable :: shellLValues(:,:)

    !> Reference (neutral) occupations for atomic shells
    real(dp), allocatable :: shellOccs(:,:)

    !> C pointer to internal state of the model (assume that model manages this)
    integer(c_intptr_t) :: modelState

  contains

    !> Update external model for geometric changes
    procedure :: update

  end type TExtModel

  !> ctype for receiving external capabilities
  type, bind(C) :: extModelAbilities

    logical(c_bool) :: gives_ham
    logical(c_bool) :: gives_overlap
    logical(c_bool) :: gives_dc_energy
    integer(c_int) :: model_derivative_order
    logical(c_bool) :: requires_self_consistency
    integer(c_int) :: spinchannels

  end type extModelAbilities


  !> C code interface
  interface

    !> External model declared API/ABI level (for any checks)
    subroutine APBI(major, minor, patch) bind(C, name='dftbp_model_apbi')
      import :: c_int
      implicit none
      !> Major version of release (potentially breaking changes)
      integer(c_int), intent(out) :: major
      !> Minor version (revised on non-breaking extensions)
      integer(c_int), intent(out) :: minor
      !> Patch level (transparent changes)
      integer(c_int), intent(out) :: patch
    end subroutine APBI

    !> External model declared capabilities
    subroutine model_provides(modelname, capabilities)&
        & bind(C, name='dftbp_provided_with')
      import :: extModelAbilities
      import :: c_char
      implicit none
      character(c_char), intent(out) :: modelname(*)
      type(extModelAbilities), intent(out) :: capabilities
    end subroutine model_provides


    !> Initialise external model for calculation
    function model_init(nspecies, speciesnames, interactionCutoff, environmentCutoff,&
        & nShellsOnSpecies, shellLValues, shellOccs, modelstate, errString)&
        & bind(C, name='initialise_model_for_dftbp')
      import :: c_int, c_ptr, c_char, c_double, c_intptr_t
      implicit none
      integer(c_int), intent(in) :: nspecies
      type(c_ptr), target, intent(in) :: speciesnames(nspecies)
      real(c_double), intent(out) :: interactionCutoff
      real(c_double), intent(out) :: environmentCutoff
      type(c_ptr), target, intent(out) :: nShellsOnSpecies
      type(c_ptr), target, intent(out) :: shellLValues
      type(c_ptr), target, intent(out) :: shellOccs
      integer(c_intptr_t) :: modelstate
      character(c_char), intent(out) :: errString(*)
      integer(c_int) :: model_init
    end function model_init


    function model_update(modelstate, errString) bind(C, name='update_model_for_dftbp')
      import :: c_int, c_ptr, c_char, c_intptr_t
      implicit none
      !> Internal state of model, opaque to us
      integer(c_intptr_t) :: modelstate
      !> Any returned error string
      character(c_char), intent(out) :: errString(*)
      !> Model state after operation
      integer(c_int) :: model_update
    end function model_update


    !> Clean up memory attached to a C pointer
    subroutine c_free(ptr) bind(c, name="free")
      import :: c_ptr
      implicit none
      !> Pointer to nullify
      type(c_ptr), value :: ptr
    end subroutine c_free

  end interface

  !> Buffer size on the Fortran side
  integer, parameter :: nBufChar = 256

contains

  !> What are the capabilities of the attached external model
  subroutine getExtModelCapabilities(modelProvides)

    !> Capabilities of externally provided hamiltonian/energy model
    type(TExtModelProvides), intent(out) :: modelProvides

    !> Structure for returned model capabilities
    type(extModelAbilities) :: capabilities

    !> Buffer on Fortran side, expecting a null termination somewhere inside of this, or throws an
    !> error
    character(nBufChar) :: modelname = " "

    integer :: major, minor, patch

    call apbi(major, minor, patch)

    write(*,'(1X,A,I0,A,I0,A,I0)')'External model API/ABI : ',major,'.',minor,'.',patch
    if (major /= 0 .and. minor /= 1) then
      write(*,*)"External model API/ABI non compliant"
      stop
    end if

    call model_provides(modelname, capabilities)
    if (.not.isCStringOK(modelname, nBufChar)) then
      write(*,*)"External model name string does not fit in buffer"
      stop
    end if
    modelProvides%modelName = trim(modelname)
    write(*,'(1X,A,A,A)')'External model used : "', modelProvides%modelName,'"'
    modelProvides%gives_ham = capabilities%gives_ham
    modelProvides%gives_overlap = capabilities%gives_overlap
    modelProvides%gives_dc_energy = capabilities%gives_dc_energy
    modelProvides%model_derivative_order = capabilities%model_derivative_order
    modelProvides%is_self_consistent = capabilities%requires_self_consistency
    modelProvides%nSpin = capabilities%spinchannels

  end subroutine getExtModelCapabilities


  !> Initialise external model for current calculation
  subroutine externalModel_init(this, speciesNames)

    !> Instance of external model
    type(TExtModel), intent(out) :: this

    !> labels of atomic species
    character(mc), intent(in) :: speciesNames(:)

    integer :: iStatus, ii, iSp, iSh
    character(nBufChar) :: errString = " "
    integer(c_int) :: nspecies
    type(c_ptr), dimension(size(speciesnames)) :: speciesPtrs
    character(mc+1), allocatable, target :: stringArray(:)
    real(c_double) :: interactionCutoff
    real(c_double) :: environmentCutoff
    type(c_ptr) :: cptr_nshells, cptr_shells, cptr_shellOccs
    integer, pointer :: fptr_nShells(:), fptr_shells(:)
    real(dp), pointer :: fptr_shellOccs(:)

    call getExtModelCapabilities(this%capabilities)

    nspecies = size(speciesNames)
    allocate(stringArray(nspecies))

    do ii = 1, nspecies
      stringArray(ii) = trim(speciesNames(ii)) // c_null_char
      speciesPtrs(ii) = c_loc(stringArray(ii))
    end do

    iStatus = model_init(nspecies, speciesPtrs, interactionCutoff, environmentCutoff, cptr_nShells,&
        & cptr_shells, cptr_shellOccs, this%modelState, errString)

    if (iStatus /= 0) then
      if (.not.isCStringOK(errString, nBufChar)) then
        write(*,*)"External model name string does not fit in buffer"
        stop
      end if
      write(*,*)trim(errString)
      stop
    end if

    this%interactionCutoff = interactionCutoff
    ! scale to allow for a cylinder around the bond between interacting atoms
    this%environmentCutoff = environmentCutoff
    this%maxCutoff = sqrt(this%environmentCutoff**2 + (this%interactionCutoff**2)/4.0_dp)

    allocate(this%nShellsOnSpecies(nspecies))
    call c_f_pointer(cptr_nShells, fptr_nShells, [nSpecies])
    this%nShellsOnSpecies(:) = fptr_nShells

    call c_f_pointer(cptr_shells, fptr_shells, [sum(this%nShellsOnSpecies)])
    allocate(this%shellLValues(maxval(this%nShellsOnSpecies), nspecies))
    this%shellLValues(:,:) = 0
    ii = 1
    do iSp = 1, nSpecies
      do iSh = 1, this%nShellsOnSpecies(iSp)
        this%shellLValues(iSh,iSp) = fptr_shells(ii)
        ii = ii + 1
      end do
    end do

    call c_f_pointer(cptr_shellOccs, fptr_shellOccs, [sum(this%nShellsOnSpecies)])
    allocate(this%shellOccs(maxval(this%nShellsOnSpecies), nspecies))
    this%shellOccs(:,:) = 0.0_dp
    ii = 1
    do iSp = 1, nSpecies
      do iSh = 1, this%nShellsOnSpecies(iSp)
        this%shellOccs(iSh,iSp) = fptr_shellOccs(ii)
        ii = ii + 1
      end do
    end do

    ! clean up
    fptr_nShells => null()
    call c_free(cptr_nShells)
    fptr_shells => null()
    call c_free(cptr_shells)
    fptr_shellOccs => null()
    call c_free(cptr_shellOccs)

  end subroutine externalModel_init


  !> Updates instance of external model for change in the geometry, and hence local environments
  !> around atoms/bonds
  subroutine update(this)

    !> Instance of external model
    class(TExtModel), intent(inout) :: this

    integer :: iStatus
    character(nBufChar) :: errString = " "

    iStatus = model_update(this%modelState, errString)
    if (iStatus /= 0) then
      if (.not.isCStringOK(errString, nBufChar)) then
        write(*,*)"External model error string does not fit in buffer"
        stop
      end if
      write(*,*)trim(errString)
      stop
    end if

  end subroutine update


  !> Checks if string has a null termination within the expected length
  function isCStringOK(string, nBufferChar)

    !> String to check
    character(c_char), intent(in) :: string(*)

    !> Expected max length of string
    integer, intent(in) :: nBufferChar

    !> Is the string within the length and null terminated
    logical isCStringOK

    integer :: ii

    isCStringOK = .false.
    lpBufCheck: do ii = 1, nBufferChar
      if (string(ii) == c_null_char) then
        isCStringOK = .true.
        exit lpBufCheck
      end if
    end do lpBufCheck

  end function isCStringOK

end module dftbp_externalmodel


program main

  use dftbp_common_accuracy
  use dftbp_externalmodel
  implicit none

  type(TExtModelProvides) :: modelAbilties
  type(TExtModel) :: model
  character(mc) :: speciesNames(2) = ["C", "H"]

  call getExtModelCapabilities(modelAbilties)
  call externalModel_init(model, speciesNames)

  call model%update()


end program main
