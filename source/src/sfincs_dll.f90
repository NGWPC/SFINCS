module sfincs_dll
  !!
  !! XMI-style shim around sfincs_bmi, exposing a simple C API:
  !!   int initialize(const char* cfg);
  !!   int update();
  !!   int finalize();
  !!   int get_value_ptr(const char* name, void** ptr);
  !!   int get_var_size(const char* name);
  !!
  !! This version is aligned with the newer sfincs_bmi2.f90, which no longer
  !! exposes nx/ny or direct BMI pointer-return arrays.
  !!
  use, intrinsic :: iso_c_binding, only: &
       c_int, c_ptr, c_char, c_null_ptr, c_null_char, &
       c_f_pointer, c_loc, c_associated
  use, intrinsic :: iso_fortran_env, only: real32

  use bmif_2_0,    only: BMI_SUCCESS, BMI_FAILURE
  use sfincs_bmi2, only: sfincs_bmi

  implicit none

  ! Single global BMI instance
  type(sfincs_bmi), save :: M

  ! DLL-owned export buffers used by get_value_ptr
  real(real32), allocatable, target, save :: export_f32(:)

contains

  !===================================================================
  !  C-callable initialize(config_file)
  !===================================================================
  function initialize(cfg_c) bind(C, name="initialize") result(ierr)
    type(c_ptr), value :: cfg_c
    integer(c_int)     :: ierr

    character(len=:), allocatable :: cfg
    integer :: status

    call cstring_to_fortran(cfg_c, cfg)

    if (.not. allocated(cfg)) then
      status = M%initialize('')
    else
      status = M%initialize(trim(cfg))
    end if

    ierr = int(status, c_int)
  end function initialize

  !===================================================================
  !  C-callable update()
  !===================================================================
  function update() bind(C, name="update") result(ierr)
    integer(c_int) :: ierr
    integer :: status

    status = M%update()
    ierr   = int(status, c_int)
  end function update

  !===================================================================
  !  C-callable finalize()
  !===================================================================
  function finalize() bind(C, name="finalize") result(ierr)
    integer(c_int) :: ierr
    integer :: status

    status = M%finalize()

    if (allocated(export_f32)) deallocate(export_f32)

    ierr = int(status, c_int)
  end function finalize

  !===================================================================
  !  C-callable get_value_ptr(name, ptr)
  !
  !  Since sfincs_bmi2 no longer supports direct BMI pointer-return access,
  !  this shim maintains DLL-owned export buffers and returns a pointer to
  !  those buffers.
  !
  !  Supported here:
  !    zs, eta2, troute_eta2, q, uv
  !
  !  Returned pointer is valid until the next get_value_ptr call or finalize().
  !===================================================================
  function get_value_ptr(name_c, ptr_c) bind(C, name="get_value_ptr") result(ierr)
    type(c_ptr), value :: name_c
    type(c_ptr)        :: ptr_c
    integer(c_int)     :: ierr

    character(len=:), allocatable :: name
    integer :: stat
    integer :: n

    call cstring_to_fortran(name_c, name)

    ptr_c = c_null_ptr
    ierr  = int(BMI_FAILURE, c_int)

    if (.not. allocated(name)) return

    n = get_var_size_from_bmi(trim(name))
    if (n <= 0) return

    if (allocated(export_f32)) then
      if (size(export_f32) /= n) then
        deallocate(export_f32)
        allocate(export_f32(n))
      end if
    else
      allocate(export_f32(n))
    end if

    export_f32 = 0.0_real32

    select case (trim(name))
    case ('zs', 'eta2', 'troute_eta2', 'q', 'uv')
      stat = M%get_value_float(trim(name), export_f32)
      if (stat /= BMI_SUCCESS) then
        ptr_c = c_null_ptr
        ierr  = int(BMI_FAILURE, c_int)
        return
      end if

      if (n > 0) then
        ptr_c = c_loc(export_f32(1))
        ierr  = int(BMI_SUCCESS, c_int)
      end if

    case default
      ptr_c = c_null_ptr
      ierr  = int(BMI_FAILURE, c_int)
    end select
  end function get_value_ptr

  !===================================================================
  !  C-callable get_var_size(name)
  !
  !  Returns number of elements, not bytes.
  !===================================================================
  function get_var_size(name_c) bind(C, name="get_var_size") result(out)
    type(c_ptr), value :: name_c
    integer(c_int)     :: out

    character(len=:), allocatable :: name
    integer :: n

    call cstring_to_fortran(name_c, name)

    if (.not. allocated(name)) then
      out = 0_c_int
      return
    end if

    n = get_var_size_from_bmi(trim(name))
    out = int(n, c_int)
  end function get_var_size

  !===================================================================
  !  Helper: return element count using BMI metadata
  !===================================================================
  integer function get_var_size_from_bmi(name) result(n)
    character(len=*), intent(in) :: name
    integer :: stat
    integer :: nbytes
    integer :: itemsize

    n = 0

    stat = M%get_var_nbytes(trim(name), nbytes)
    if (stat /= BMI_SUCCESS) return

    stat = M%get_var_itemsize(trim(name), itemsize)
    if (stat /= BMI_SUCCESS) return

    if (itemsize <= 0) return
    n = nbytes / itemsize
  end function get_var_size_from_bmi

  !===================================================================
  !  Helper: Convert C string (char*) -> allocatable Fortran string
  !===================================================================
  subroutine cstring_to_fortran(cstr, fstr)
    type(c_ptr), value :: cstr
    character(len=:), allocatable :: fstr

    character(kind=c_char), pointer :: p_chars(:)
    integer :: maxlen, n, i
    character(len=1, kind=c_char) :: ch_c
    character(len=1)              :: ch_f

    if (.not. c_associated(cstr)) then
      return
    end if

    maxlen = 10000
    call c_f_pointer(cstr, p_chars, [maxlen])

    n = 0
    do i = 1, maxlen
      if (p_chars(i) == c_null_char) exit
      n = n + 1
    end do

    if (n <= 0) then
      allocate(character(len=0) :: fstr)
      return
    end if

    allocate(character(len=n) :: fstr)
    do i = 1, n
      ch_c = p_chars(i)
      ch_f = transfer(ch_c, ch_f)
      fstr(i:i) = ch_f
    end do
  end subroutine cstring_to_fortran

end module sfincs_dll
