module sfincs_bmi2
  !! BMI 2.0 wrapper for SFINCS
  !! Quadtree-aware / unstructured-grid implementation
  !!
  !! Direct wrapper around SFINCS globals.
  !! No mirror arrays, no changes required in sfincs_data.f90.

  use, intrinsic :: iso_fortran_env, only: real32, real64
  use bmif_2_0_iso, only: bmi, BMI_SUCCESS, BMI_FAILURE

  use sfincs_data, only: &
      np, npuv, t0, t1, &
      zs, q, uv, zsmax, z_volume, &
      qext, prcp, windu, windv, patm, uorb, &
      z_xz, z_yz, &
      uv_index_z_nm, uv_index_z_nmu, &
      kcs, z_flags_iref, uv_flags_dir, uv_flags_type

  use sfincs_lib, only: sfincs_initialize, sfincs_update, sfincs_finalize, t, dt

  implicit none
  private

  integer, parameter :: GRID_CELL = 0
  integer, parameter :: GRID_EDGE = 1

  character(len=*), parameter :: VAR_ZS      = 'zs'
  character(len=*), parameter :: VAR_Q       = 'q'
  character(len=*), parameter :: VAR_UV      = 'uv'
  character(len=*), parameter :: VAR_ZSMAX   = 'zsmax'
  character(len=*), parameter :: VAR_ZVOL    = 'z_volume'

  character(len=*), parameter :: VAR_QEXT    = 'qext'
  character(len=*), parameter :: VAR_PRCP    = 'prcp'
  character(len=*), parameter :: VAR_WINDU   = 'windu'
  character(len=*), parameter :: VAR_WINDV   = 'windv'
  character(len=*), parameter :: VAR_PATM    = 'patm'
  character(len=*), parameter :: VAR_UORB    = 'uorb'

  character(len=*), parameter :: VAR_Z_XZ    = 'z_xz'
  character(len=*), parameter :: VAR_Z_YZ    = 'z_yz'
  character(len=*), parameter :: VAR_KCS     = 'kcs'
  character(len=*), parameter :: VAR_UV_NM   = 'uv_index_z_nm'
  character(len=*), parameter :: VAR_UV_NMU  = 'uv_index_z_nmu'
  character(len=*), parameter :: VAR_UV_DIR  = 'uv_flags_dir'
  character(len=*), parameter :: VAR_UV_TYPE = 'uv_flags_type'
  character(len=*), parameter :: VAR_Z_IREF  = 'z_flags_iref'

  ! Legacy aliases
  character(len=*), parameter :: VAR_ETA2        = 'eta2'
  character(len=*), parameter :: VAR_TROUTE_ETA2 = 'troute_eta2'

  character(len=:), allocatable, target, save :: g_component_name
  character(len=:), allocatable, target, save :: g_time_units

  character(len=16), target, save :: g_input_names(6) = [ &
      'qext            ', &
      'prcp            ', &
      'windu           ', &
      'windv           ', &
      'patm            ', &
      'uorb            '  ]

  character(len=16), target, save :: g_output_names(12) = [ &
      'zs              ', &
      'eta2            ', &
      'troute_eta2     ', &
      'q               ', &
      'uv              ', &
      'zsmax           ', &
      'z_volume        ', &
      'z_xz            ', &
      'z_yz            ', &
      'kcs             ', &
      'uv_index_z_nm   ', &
      'uv_index_z_nmu  '  ]

  type, extends(bmi) :: sfincs_bmi
    real(real64) :: t       = 0.0d0
    real(real64) :: dt      = 0.0d0
    real(real64) :: t_start = 0.0d0
    real(real64) :: t_end   = 0.0d0

    integer :: n_cells = 0
    integer :: n_edges = 0

    character(len=:), pointer :: component_name => null()

    logical :: is_initialized = .false.

  contains
    procedure :: initialize                 => sfincs_bmi_initialize
    procedure :: update                     => sfincs_bmi_update
    procedure :: update_until               => sfincs_bmi_update_until
    procedure :: finalize                   => sfincs_bmi_finalize

    procedure :: get_component_name         => sfincs_bmi_get_component_name
    procedure :: get_input_item_count       => sfincs_bmi_get_input_item_count
    procedure :: get_output_item_count      => sfincs_bmi_get_output_item_count
    procedure :: get_input_var_names        => sfincs_bmi_get_input_var_names
    procedure :: get_output_var_names       => sfincs_bmi_get_output_var_names

    procedure :: get_start_time             => sfincs_bmi_get_start_time
    procedure :: get_end_time               => sfincs_bmi_get_end_time
    procedure :: get_current_time           => sfincs_bmi_get_current_time
    procedure :: get_time_step              => sfincs_bmi_get_time_step
    procedure :: get_time_units             => sfincs_bmi_get_time_units

    procedure :: get_var_grid               => sfincs_bmi_get_var_grid
    procedure :: get_var_type               => sfincs_bmi_get_var_type
    procedure :: get_var_units              => sfincs_bmi_get_var_units
    procedure :: get_var_itemsize           => sfincs_bmi_get_var_itemsize
    procedure :: get_var_nbytes             => sfincs_bmi_get_var_nbytes
    procedure :: get_var_location           => sfincs_bmi_get_var_location

    procedure :: get_grid_rank              => sfincs_bmi_get_grid_rank
    procedure :: get_grid_size              => sfincs_bmi_get_grid_size
    procedure :: get_grid_type              => sfincs_bmi_get_grid_type
    procedure :: get_grid_shape             => sfincs_bmi_get_grid_shape
    procedure :: get_grid_spacing           => sfincs_bmi_get_grid_spacing
    procedure :: get_grid_origin            => sfincs_bmi_get_grid_origin
    procedure :: get_grid_x                 => sfincs_bmi_get_grid_x
    procedure :: get_grid_y                 => sfincs_bmi_get_grid_y
    procedure :: get_grid_z                 => sfincs_bmi_get_grid_z

    procedure :: get_grid_edge_count        => sfincs_bmi_get_grid_edge_count
    procedure :: get_grid_face_count        => sfincs_bmi_get_grid_face_count
    procedure :: get_grid_node_count        => sfincs_bmi_get_grid_node_count
    procedure :: get_grid_edge_nodes        => sfincs_bmi_get_grid_edge_nodes
    procedure :: get_grid_face_nodes        => sfincs_bmi_get_grid_face_nodes
    procedure :: get_grid_nodes_per_face    => sfincs_bmi_get_grid_nodes_per_face
    procedure :: get_grid_face_edges        => sfincs_bmi_get_grid_face_edges

    procedure :: get_value_float            => sfincs_bmi_get_value_float
    procedure :: set_value_float            => sfincs_bmi_set_value_float
    procedure :: get_value_ptr_float        => sfincs_bmi_get_value_ptr_float
    procedure :: get_value_at_indices_float => sfincs_bmi_get_value_at_indices_float
    procedure :: set_value_at_indices_float => sfincs_bmi_set_value_at_indices_float

    procedure :: get_value_double           => sfincs_bmi_get_value_double
    procedure :: set_value_double           => sfincs_bmi_set_value_double
    procedure :: get_value_ptr_double       => sfincs_bmi_get_value_ptr_double
    procedure :: get_value_at_indices_double=> sfincs_bmi_get_value_at_indices_double
    procedure :: set_value_at_indices_double=> sfincs_bmi_set_value_at_indices_double

    procedure :: get_value_int              => sfincs_bmi_get_value_int
    procedure :: set_value_int              => sfincs_bmi_set_value_int
    procedure :: get_value_ptr_int          => sfincs_bmi_get_value_ptr_int
    procedure :: get_value_at_indices_int   => sfincs_bmi_get_value_at_indices_int
    procedure :: set_value_at_indices_int   => sfincs_bmi_set_value_at_indices_int
  end type sfincs_bmi

  public :: sfincs_bmi

contains

  function sfincs_bmi_initialize(this, config_file) result(status)
    class(sfincs_bmi), intent(out) :: this
    character(len=*),  intent(in)  :: config_file
    integer :: status
    integer :: ierr

    write(*,*) 'sfincs_bmi_initialize config_file = ', trim(config_file)

    ierr = sfincs_initialize()
    if (ierr /= 0) then
      status = BMI_FAILURE
      return
    end if

    this%t_start = t0
    this%t_end   = t1
    this%t       = t

    if (this%t_end > this%t_start) then
      this%dt = max(1.0d0, (this%t_end - this%t_start) / 100.d0)
    else
      this%dt = 1.0d0
    end if

    this%n_cells = np
    this%n_edges = npuv

    if (.not. allocated(g_component_name)) then
      allocate(character(len=10) :: g_component_name)
      g_component_name = 'SFINCS BMI'
    end if

    if (.not. allocated(g_time_units)) then
      allocate(character(len=1) :: g_time_units)
      g_time_units = 's'
    end if

    this%component_name => g_component_name
    this%is_initialized = .true.
    status = BMI_SUCCESS
  end function sfincs_bmi_initialize

  function sfincs_bmi_update(this) result(status)
    class(sfincs_bmi), intent(inout) :: this
    integer :: status
    integer :: ierr
    double precision :: dtrange

    if (.not. this%is_initialized) then
      status = BMI_FAILURE
      return
    end if

    dtrange = this%dt
    ierr = sfincs_update(dtrange)
    if (ierr /= 0) then
      status = BMI_FAILURE
      return
    end if

    this%t = t
    status = BMI_SUCCESS
  end function sfincs_bmi_update

  function sfincs_bmi_update_until(this, time) result(status)
    class(sfincs_bmi), intent(inout) :: this
    double precision,  intent(in)    :: time
    integer :: status
    integer :: ierr
    double precision :: dtrange

    if (.not. this%is_initialized) then
      status = BMI_FAILURE
      return
    end if

    if (time <= this%t + 1.0d-12) then
      status = BMI_SUCCESS
      return
    end if

    dtrange = time - this%t
    ierr = sfincs_update(dtrange)
    if (ierr /= 0) then
      status = BMI_FAILURE
      return
    end if

    this%t = t
    status = BMI_SUCCESS
  end function sfincs_bmi_update_until

  function sfincs_bmi_finalize(this) result(status)
    class(sfincs_bmi), intent(inout) :: this
    integer :: status
    integer :: ierr

    if (.not. this%is_initialized) then
      status = BMI_SUCCESS
      return
    end if

    ierr = sfincs_finalize()
    nullify(this%component_name)
    this%is_initialized = .false.

    if (ierr /= 0) then
      status = BMI_FAILURE
    else
      status = BMI_SUCCESS
    end if
  end function sfincs_bmi_finalize

  function sfincs_bmi_get_component_name(this, name) result(status)
    class(sfincs_bmi),         intent(in)  :: this
    character(len=:), pointer, intent(out) :: name
    integer :: status

    if (associated(this%component_name)) then
      name => this%component_name
      status = BMI_SUCCESS
    else
      nullify(name)
      status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_component_name

  function sfincs_bmi_get_input_var_names(this, names) result(status)
    class(sfincs_bmi),         intent(in)  :: this
    character(len=:), pointer, intent(out) :: names(:)
    integer :: status

    names => g_input_names
    status = BMI_SUCCESS
  end function sfincs_bmi_get_input_var_names

  function sfincs_bmi_get_output_var_names(this, names) result(status)
    class(sfincs_bmi),         intent(in)  :: this
    character(len=:), pointer, intent(out) :: names(:)
    integer :: status

    names => g_output_names
    status = BMI_SUCCESS
  end function sfincs_bmi_get_output_var_names

  function sfincs_bmi_get_input_item_count(this, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer, intent(out) :: count
    integer :: status

    count = size(g_input_names)
    status = BMI_SUCCESS
  end function sfincs_bmi_get_input_item_count

  function sfincs_bmi_get_output_item_count(this, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer, intent(out) :: count
    integer :: status

    count = size(g_output_names)
    status = BMI_SUCCESS
  end function sfincs_bmi_get_output_item_count

  function sfincs_bmi_get_start_time(this, time) result(status)
    class(sfincs_bmi), intent(in)  :: this
    double precision,  intent(out) :: time
    integer :: status
    time = this%t_start
    status = BMI_SUCCESS
  end function sfincs_bmi_get_start_time

  function sfincs_bmi_get_end_time(this, time) result(status)
    class(sfincs_bmi), intent(in)  :: this
    double precision,  intent(out) :: time
    integer :: status
    time = this%t_end
    status = BMI_SUCCESS
  end function sfincs_bmi_get_end_time

  function sfincs_bmi_get_current_time(this, time) result(status)
    class(sfincs_bmi), intent(in)  :: this
    double precision,  intent(out) :: time
    integer :: status
    time = this%t
    status = BMI_SUCCESS
  end function sfincs_bmi_get_current_time

  function sfincs_bmi_get_time_step(this, time_step) result(status)
    class(sfincs_bmi), intent(in)  :: this
    double precision,  intent(out) :: time_step
    integer :: status
    time_step = this%dt
    status = BMI_SUCCESS
  end function sfincs_bmi_get_time_step

  function sfincs_bmi_get_time_units(this, units) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(out) :: units
    integer :: status
    units = 's'
    status = BMI_SUCCESS
  end function sfincs_bmi_get_time_units

  function sfincs_bmi_get_var_grid(this, name, grid) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    integer,           intent(out) :: grid
    integer :: status
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)

    select case (trim(cname))
    case (VAR_ZS, VAR_ETA2, VAR_TROUTE_ETA2, VAR_ZSMAX, VAR_ZVOL, VAR_QEXT, VAR_PRCP, VAR_WINDU, VAR_WINDV, &
          VAR_PATM, VAR_UORB, VAR_Z_XZ, VAR_Z_YZ, VAR_KCS, VAR_Z_IREF)
      grid = GRID_CELL
      status = BMI_SUCCESS
    case (VAR_Q, VAR_UV, VAR_UV_NM, VAR_UV_NMU, VAR_UV_DIR, VAR_UV_TYPE)
      grid = GRID_EDGE
      status = BMI_SUCCESS
    case default
      grid = -1
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_var_grid

  function sfincs_bmi_get_var_type(this, name, type) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    character(len=*),  intent(out) :: type
    integer :: status
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)

    select case (trim(cname))
    case (VAR_ZS, VAR_ETA2, VAR_TROUTE_ETA2, VAR_ZVOL)
      type = 'double precision'
      status = BMI_SUCCESS
    case (VAR_Q, VAR_UV, VAR_ZSMAX, VAR_QEXT, VAR_PRCP, VAR_WINDU, VAR_WINDV, VAR_PATM, VAR_UORB, &
          VAR_Z_XZ, VAR_Z_YZ)
      type = 'real'
      status = BMI_SUCCESS
    case (VAR_KCS, VAR_Z_IREF, VAR_UV_NM, VAR_UV_NMU, VAR_UV_DIR, VAR_UV_TYPE)
      type = 'integer'
      status = BMI_SUCCESS
    case default
      type = ''
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_var_type

  function sfincs_bmi_get_var_units(this, name, units) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    character(len=*),  intent(out) :: units
    integer :: status
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)

    select case (trim(cname))
    case (VAR_ZS, VAR_ETA2, VAR_TROUTE_ETA2, VAR_ZSMAX, VAR_ZVOL)
      units = 'm'
      status = BMI_SUCCESS
    case (VAR_Q, VAR_QEXT)
      units = 'm3 s-1'
      status = BMI_SUCCESS
    case (VAR_UV, VAR_WINDU, VAR_WINDV, VAR_UORB)
      units = 'm s-1'
      status = BMI_SUCCESS
    case (VAR_PRCP)
      units = 'm s-1'
      status = BMI_SUCCESS
    case (VAR_PATM)
      units = 'Pa'
      status = BMI_SUCCESS
    case (VAR_Z_XZ, VAR_Z_YZ)
      units = 'm'
      status = BMI_SUCCESS
    case (VAR_KCS, VAR_Z_IREF, VAR_UV_NM, VAR_UV_NMU, VAR_UV_DIR, VAR_UV_TYPE)
      units = '1'
      status = BMI_SUCCESS
    case default
      units = ''
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_var_units

  function sfincs_bmi_get_var_itemsize(this, name, size) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    integer,           intent(out) :: size
    integer :: status
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)

    select case (trim(cname))
    case (VAR_ZS, VAR_ETA2, VAR_TROUTE_ETA2, VAR_ZVOL)
      size = 8
      status = BMI_SUCCESS
    case (VAR_Q, VAR_UV, VAR_ZSMAX, VAR_QEXT, VAR_PRCP, VAR_WINDU, VAR_WINDV, VAR_PATM, VAR_UORB, &
          VAR_Z_XZ, VAR_Z_YZ)
      size = 4
      status = BMI_SUCCESS
    case (VAR_KCS, VAR_Z_IREF, VAR_UV_NM, VAR_UV_NMU, VAR_UV_DIR, VAR_UV_TYPE)
      size = 4
      status = BMI_SUCCESS
    case default
      size = 0
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_var_itemsize

  function sfincs_bmi_get_var_nbytes(this, name, nbytes) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    integer,           intent(out) :: nbytes
    integer :: status, itemsize, n
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)
    status = this%get_var_itemsize(cname, itemsize)
    if (status /= BMI_SUCCESS) then
      nbytes = 0
      return
    end if

    select case (trim(cname))
    case (VAR_ZS, VAR_ETA2, VAR_TROUTE_ETA2, VAR_ZSMAX, VAR_ZVOL, VAR_QEXT, VAR_PRCP, VAR_WINDU, VAR_WINDV, &
          VAR_PATM, VAR_UORB, VAR_Z_XZ, VAR_Z_YZ, VAR_KCS, VAR_Z_IREF)
      n = np
    case (VAR_Q, VAR_UV, VAR_UV_NM, VAR_UV_NMU, VAR_UV_DIR, VAR_UV_TYPE)
      n = npuv
    case default
      nbytes = 0
      status = BMI_FAILURE
      return
    end select

    nbytes = itemsize * n
    status = BMI_SUCCESS
  end function sfincs_bmi_get_var_nbytes

  function sfincs_bmi_get_var_location(this, name, location) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    character(len=*),  intent(out) :: location
    integer :: status
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)

    select case (trim(cname))
    case (VAR_ZS, VAR_ETA2, VAR_TROUTE_ETA2, VAR_ZSMAX, VAR_ZVOL, VAR_QEXT, VAR_PRCP, VAR_WINDU, VAR_WINDV, &
          VAR_PATM, VAR_UORB, VAR_Z_XZ, VAR_Z_YZ, VAR_KCS, VAR_Z_IREF)
      location = 'node'
      status = BMI_SUCCESS
    case (VAR_Q, VAR_UV, VAR_UV_NM, VAR_UV_NMU, VAR_UV_DIR, VAR_UV_TYPE)
      location = 'edge'
      status = BMI_SUCCESS
    case default
      location = ''
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_var_location

  function sfincs_bmi_get_grid_rank(this, grid, rank) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: rank
    integer :: status
    if (grid == GRID_CELL .or. grid == GRID_EDGE) then
      rank = 2
      status = BMI_SUCCESS
    else
      rank = 0
      status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_grid_rank

  function sfincs_bmi_get_grid_size(this, grid, size) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: size
    integer :: status
    select case (grid)
    case (GRID_CELL)
      size = np
      status = BMI_SUCCESS
    case (GRID_EDGE)
      size = npuv
      status = BMI_SUCCESS
    case default
      size = 0
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_grid_size

  function sfincs_bmi_get_grid_type(this, grid, type) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    character(len=*),  intent(out) :: type
    integer :: status
    if (grid == GRID_CELL .or. grid == GRID_EDGE) then
      type = 'unstructured'
      status = BMI_SUCCESS
    else
      type = ''
      status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_grid_type

  function sfincs_bmi_get_grid_shape(this, grid, shape) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: shape(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_shape

  function sfincs_bmi_get_grid_spacing(this, grid, spacing) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: spacing(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_spacing

  function sfincs_bmi_get_grid_origin(this, grid, origin) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: origin(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_origin

  function sfincs_bmi_get_grid_x(this, grid, x) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: x(:)
    integer :: status
    if (grid /= GRID_CELL) then
      status = BMI_FAILURE
      return
    end if
    if (size(x) < np) then
      status = BMI_FAILURE
      return
    end if
    x(1:np) = real(z_xz(1:np), kind=real64)
    status = BMI_SUCCESS
  end function sfincs_bmi_get_grid_x

  function sfincs_bmi_get_grid_y(this, grid, y) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: y(:)
    integer :: status
    if (grid /= GRID_CELL) then
      status = BMI_FAILURE
      return
    end if
    if (size(y) < np) then
      status = BMI_FAILURE
      return
    end if
    y(1:np) = real(z_yz(1:np), kind=real64)
    status = BMI_SUCCESS
  end function sfincs_bmi_get_grid_y

  function sfincs_bmi_get_grid_z(this, grid, z) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: z(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_z

  function sfincs_bmi_get_grid_node_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: count
    integer :: status
    if (grid == GRID_CELL) then
      count = np
      status = BMI_SUCCESS
    else
      count = 0
      status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_grid_node_count

  function sfincs_bmi_get_grid_edge_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: count
    integer :: status
    if (grid == GRID_EDGE) then
      count = npuv
      status = BMI_SUCCESS
    else
      count = 0
      status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_grid_edge_count

  function sfincs_bmi_get_grid_face_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: count
    integer :: status
    if (grid == GRID_CELL) then
      count = np
      status = BMI_SUCCESS
    else
      count = 0
      status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_grid_face_count

  function sfincs_bmi_get_grid_edge_nodes(this, grid, edge_nodes) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: edge_nodes(:)
    integer :: status
    integer :: ip

    if (grid /= GRID_EDGE) then
      status = BMI_FAILURE
      return
    end if
    if (size(edge_nodes) < 2*npuv) then
      status = BMI_FAILURE
      return
    end if

    do ip = 1, npuv
      edge_nodes(2*ip - 1) = uv_index_z_nm(ip)  - 1
      edge_nodes(2*ip    ) = uv_index_z_nmu(ip) - 1
    end do

    status = BMI_SUCCESS
  end function sfincs_bmi_get_grid_edge_nodes

  function sfincs_bmi_get_grid_face_nodes(this, grid, face_nodes) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: face_nodes(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_face_nodes

  function sfincs_bmi_get_grid_nodes_per_face(this, grid, nodes_per_face) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: nodes_per_face(:)
    integer :: status
    if (grid /= GRID_CELL) then
      status = BMI_FAILURE
      return
    end if
    if (size(nodes_per_face) < np) then
      status = BMI_FAILURE
      return
    end if
    nodes_per_face(1:np) = 4
    status = BMI_SUCCESS
  end function sfincs_bmi_get_grid_nodes_per_face

  function sfincs_bmi_get_grid_face_edges(this, grid, face_edges) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: face_edges(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_face_edges

  function sfincs_bmi_get_value_float(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(real32),      intent(inout) :: dest(:)
    integer :: status
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)

    select case (trim(cname))
    case (VAR_ZS)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = real(zs(1:np), kind=real32)
      status = BMI_SUCCESS
    case (VAR_Q)
      if (size(dest) < npuv) then
        status = BMI_FAILURE; return
      end if
      dest(1:npuv) = q(1:npuv)
      status = BMI_SUCCESS
    case (VAR_UV)
      if (size(dest) < npuv) then
        status = BMI_FAILURE; return
      end if
      dest(1:npuv) = uv(1:npuv)
      status = BMI_SUCCESS
    case (VAR_ZSMAX)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = zsmax(1:np)
      status = BMI_SUCCESS
    case (VAR_QEXT)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = qext(1:np)
      status = BMI_SUCCESS
    case (VAR_PRCP)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = prcp(1:np)
      status = BMI_SUCCESS
    case (VAR_WINDU)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = windu(1:np)
      status = BMI_SUCCESS
    case (VAR_WINDV)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = windv(1:np)
      status = BMI_SUCCESS
    case (VAR_PATM)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = patm(1:np)
      status = BMI_SUCCESS
    case (VAR_UORB)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = uorb(1:np)
      status = BMI_SUCCESS
    case (VAR_Z_XZ)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = z_xz(1:np)
      status = BMI_SUCCESS
    case (VAR_Z_YZ)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = z_yz(1:np)
      status = BMI_SUCCESS
    case default
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_value_float

  function sfincs_bmi_set_value_float(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    real(real32),      intent(in)    :: src(:)
    integer :: status
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)

    select case (trim(cname))
    case (VAR_QEXT)
      if (size(src) < np) then
        status = BMI_FAILURE; return
      end if
      qext(1:np) = src(1:np)
      status = BMI_SUCCESS
    case (VAR_PRCP)
      if (size(src) < np) then
        status = BMI_FAILURE; return
      end if
      prcp(1:np) = src(1:np)
      status = BMI_SUCCESS
    case (VAR_WINDU)
      if (size(src) < np) then
        status = BMI_FAILURE; return
      end if
      windu(1:np) = src(1:np)
      status = BMI_SUCCESS
    case (VAR_WINDV)
      if (size(src) < np) then
        status = BMI_FAILURE; return
      end if
      windv(1:np) = src(1:np)
      status = BMI_SUCCESS
    case (VAR_PATM)
      if (size(src) < np) then
        status = BMI_FAILURE; return
      end if
      patm(1:np) = src(1:np)
      status = BMI_SUCCESS
    case (VAR_UORB)
      if (size(src) < np) then
        status = BMI_FAILURE; return
      end if
      uorb(1:np) = src(1:np)
      status = BMI_SUCCESS
    case default
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_set_value_float

  function sfincs_bmi_get_value_ptr_float(this, name, dest_ptr) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    real(real32),      pointer, intent(inout) :: dest_ptr(:)
    integer :: status
    nullify(dest_ptr)
    status = BMI_FAILURE
  end function sfincs_bmi_get_value_ptr_float

  function sfincs_bmi_get_value_at_indices_float(this, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(real32),      intent(inout) :: dest(:)
    integer,           intent(in)    :: inds(:)
    integer :: status, n, k
    real(real32), allocatable :: tmp(:)
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)
    n = size(inds)

    if (size(dest) < n) then
      status = BMI_FAILURE
      return
    end if

    select case (trim(cname))
    case (VAR_Q, VAR_UV)
      allocate(tmp(npuv))
    case default
      allocate(tmp(np))
    end select

    status = this%get_value_float(cname, tmp)
    if (status /= BMI_SUCCESS) then
      deallocate(tmp)
      return
    end if

    do k = 1, n
      dest(k) = tmp(inds(k))
    end do
    deallocate(tmp)
    status = BMI_SUCCESS
  end function sfincs_bmi_get_value_at_indices_float

  function sfincs_bmi_set_value_at_indices_float(this, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: inds(:)
    real(real32),      intent(in)    :: src(:)
    integer :: status, k, n
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)
    n = size(inds)

    if (size(src) < n) then
      status = BMI_FAILURE
      return
    end if

    select case (trim(cname))
    case (VAR_QEXT)
      do k = 1, n
        qext(inds(k)) = src(k)
      end do
      status = BMI_SUCCESS
    case (VAR_PRCP)
      do k = 1, n
        prcp(inds(k)) = src(k)
      end do
      status = BMI_SUCCESS
    case (VAR_WINDU)
      do k = 1, n
        windu(inds(k)) = src(k)
      end do
      status = BMI_SUCCESS
    case (VAR_WINDV)
      do k = 1, n
        windv(inds(k)) = src(k)
      end do
      status = BMI_SUCCESS
    case (VAR_PATM)
      do k = 1, n
        patm(inds(k)) = src(k)
      end do
      status = BMI_SUCCESS
    case (VAR_UORB)
      do k = 1, n
        uorb(inds(k)) = src(k)
      end do
      status = BMI_SUCCESS
    case default
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_set_value_at_indices_float

  function sfincs_bmi_get_value_double(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(real64),      intent(inout) :: dest(:)
    integer :: status
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)

    select case (trim(cname))
    case (VAR_ZS)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = zs(1:np)
      status = BMI_SUCCESS
    case (VAR_ZVOL)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = z_volume(1:np)
      status = BMI_SUCCESS
    case default
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_value_double

  function sfincs_bmi_set_value_double(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    real(real64),      intent(in)    :: src(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_set_value_double

  function sfincs_bmi_get_value_ptr_double(this, name, dest_ptr) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    real(real64),      pointer, intent(inout) :: dest_ptr(:)
    integer :: status
    nullify(dest_ptr)
    status = BMI_FAILURE
  end function sfincs_bmi_get_value_ptr_double

  function sfincs_bmi_get_value_at_indices_double(this, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(real64),      intent(inout) :: dest(:)
    integer,           intent(in)    :: inds(:)
    integer :: status, n, k
    real(real64), allocatable :: tmp(:)
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)
    n = size(inds)

    if (size(dest) < n) then
      status = BMI_FAILURE
      return
    end if

    allocate(tmp(np))
    status = this%get_value_double(cname, tmp)
    if (status /= BMI_SUCCESS) then
      deallocate(tmp)
      return
    end if

    do k = 1, n
      dest(k) = tmp(inds(k))
    end do
    deallocate(tmp)
    status = BMI_SUCCESS
  end function sfincs_bmi_get_value_at_indices_double

  function sfincs_bmi_set_value_at_indices_double(this, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: inds(:)
    real(real64),      intent(in)    :: src(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_set_value_at_indices_double

  function sfincs_bmi_get_value_int(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(inout) :: dest(:)
    integer :: status
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)

    select case (trim(cname))
    case (VAR_KCS)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = int(kcs(1:np))
      status = BMI_SUCCESS
    case (VAR_Z_IREF)
      if (size(dest) < np) then
        status = BMI_FAILURE; return
      end if
      dest(1:np) = int(z_flags_iref(1:np))
      status = BMI_SUCCESS
    case (VAR_UV_NM)
      if (size(dest) < npuv) then
        status = BMI_FAILURE; return
      end if
      dest(1:npuv) = uv_index_z_nm(1:npuv)
      status = BMI_SUCCESS
    case (VAR_UV_NMU)
      if (size(dest) < npuv) then
        status = BMI_FAILURE; return
      end if
      dest(1:npuv) = uv_index_z_nmu(1:npuv)
      status = BMI_SUCCESS
    case (VAR_UV_DIR)
      if (size(dest) < npuv) then
        status = BMI_FAILURE; return
      end if
      dest(1:npuv) = int(uv_flags_dir(1:npuv))
      status = BMI_SUCCESS
    case (VAR_UV_TYPE)
      if (size(dest) < npuv) then
        status = BMI_FAILURE; return
      end if
      dest(1:npuv) = int(uv_flags_type(1:npuv))
      status = BMI_SUCCESS
    case default
      status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_value_int

  function sfincs_bmi_set_value_int(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: src(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_set_value_int

  function sfincs_bmi_get_value_ptr_int(this, name, dest_ptr) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    integer,           pointer, intent(inout) :: dest_ptr(:)
    integer :: status
    nullify(dest_ptr)
    status = BMI_FAILURE
  end function sfincs_bmi_get_value_ptr_int

  function sfincs_bmi_get_value_at_indices_int(this, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(inout) :: dest(:)
    integer,           intent(in)    :: inds(:)
    integer :: status, n, k
    integer, allocatable :: tmp(:)
    character(len=:), allocatable :: cname

    cname = canon_var_name(name)
    n = size(inds)

    if (size(dest) < n) then
      status = BMI_FAILURE
      return
    end if

    select case (trim(cname))
    case (VAR_UV_NM, VAR_UV_NMU, VAR_UV_DIR, VAR_UV_TYPE)
      allocate(tmp(npuv))
    case default
      allocate(tmp(np))
    end select

    status = this%get_value_int(cname, tmp)
    if (status /= BMI_SUCCESS) then
      deallocate(tmp)
      return
    end if

    do k = 1, n
      dest(k) = tmp(inds(k))
    end do
    deallocate(tmp)
    status = BMI_SUCCESS
  end function sfincs_bmi_get_value_at_indices_int

  function sfincs_bmi_set_value_at_indices_int(this, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: inds(:)
    integer,           intent(in)    :: src(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_set_value_at_indices_int

  function lower_str(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
    integer :: i

    out = str
    do i = 1, len(str)
      if (iachar(str(i:i)) >= iachar('A') .and. iachar(str(i:i)) <= iachar('Z')) then
        out(i:i) = achar(iachar(str(i:i)) + 32)
      end if
    end do
  end function lower_str

  function canon_var_name(name) result(canon)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: canon
    character(len=:), allocatable :: cname

    cname = lower_str(trim(name))

    select case (trim(cname))
    case ('zs', 'eta2', 'troute_eta2', 'troute-eta2', 'trouteeta2')
      canon = VAR_ZS
    case ('q')
      canon = VAR_Q
    case ('uv')
      canon = VAR_UV
    case ('zsmax')
      canon = VAR_ZSMAX
    case ('z_volume')
      canon = VAR_ZVOL
    case ('qext')
      canon = VAR_QEXT
    case ('prcp')
      canon = VAR_PRCP
    case ('windu')
      canon = VAR_WINDU
    case ('windv')
      canon = VAR_WINDV
    case ('patm')
      canon = VAR_PATM
    case ('uorb')
      canon = VAR_UORB
    case ('z_xz')
      canon = VAR_Z_XZ
    case ('z_yz')
      canon = VAR_Z_YZ
    case ('kcs')
      canon = VAR_KCS
    case ('uv_index_z_nm')
      canon = VAR_UV_NM
    case ('uv_index_z_nmu')
      canon = VAR_UV_NMU
    case ('uv_flags_dir')
      canon = VAR_UV_DIR
    case ('uv_flags_type')
      canon = VAR_UV_TYPE
    case ('z_flags_iref')
      canon = VAR_Z_IREF
    case default
      canon = cname
    end select
  end function canon_var_name

end module sfincs_bmi2
