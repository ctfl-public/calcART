! Minimal forpy example: call Python calcART.calc_dq_ED from Fortran
!
! This example uses forpy (Fortran-Python interface) to import the Python
! module `calcART` directly, call `calc_dq_ED` with an array of depths (t),
! then convert the returned NumPy ndarray to a Python list and copy it into
! a Fortran array for printing.
!
! NOTE: forpy's API names can differ slightly by version. The calls used here
! (forpy_initialize, import_module, getattr, tuple_create, tuple_setitem, call,
! to_double, forpy_finalize) are common; adjust names if your forpy version uses
! e.g. create_tuple/tuple_setitem or call_object.

program forpy_calc_dq_example
   use iso_c_binding
   use forpy_mod
   implicit none

   integer :: ierror
   type(module_py) :: calcART
   type(object) :: result
   type(ndarray) :: dq_nd    ! numpy ndarray wrapper for returned result
   type(tuple) :: args
   type(dict) :: kwargs

   real(8) :: q, beta, omega, g1
   real(8), pointer :: dq_ptr(:)  ! pointer view into NumPy data (no copy)
   real(8), allocatable :: dq(:)  ! optional Fortran copy of results
   real(8), allocatable :: t(:)
   integer :: n_t, i
   type(ndarray) :: py_t
   character(len=:), allocatable :: SF

   ! Minimal example: import Python module `calcART` and call `calc_dq_ED`
   ierror = forpy_initialize()
   ierror = import_py(calcART, "calcART")

   ! Prepare arguments for calc_dq_ED(q, t, beta, omega, SF, g1)
   q = 500.0e4
   beta = 10000.0d0
   omega = 0.5d0
   g1 = 0.0d0
   SF = "HG"

   ! build a small t array and wrap it as a NumPy ndarray
   n_t = 10
   allocate(t(n_t))
   do i = 1, n_t
      t(i) = (i-1) * 0.0001d0  ! depths in meters (example)
   end do

   ierror = ndarray_create(py_t, t)

   ierror = tuple_create(args, 6)
   ierror = args%setitem(0, q)
   ierror = args%setitem(1, py_t)
   ierror = args%setitem(2, beta)
   ierror = args%setitem(3, omega)
   ierror = args%setitem(4, SF)
   ierror = args%setitem(5, g1)

   ierror = dict_create(kwargs) ! empty kwargs

   ! Call calc_dq_ED from the calcART module
   ierror = call_py(result, calcART, "calc_dq_ED", args, kwargs)
   ! ! Print the returned Python object (e.g. NumPy array)
   ! ierror = print_py(result)

   ! Cast returned Python object to numpy ndarray and obtain Fortran pointer
   ierror = cast(dq_nd, result)
   if (ierror /= 0) then
      call err_print()
   else
      ierror = dq_nd%get_data(dq_ptr)  ! default requires Fortran-order, will transpose if C-order
      if (ierror /= 0) then
         call err_print()
      else
         ! Make an allocatable Fortran array copy (optional)
         allocate(dq(size(dq_ptr)))
         dq = dq_ptr
         write(*,'(A,I0)') 'dq length = ', size(dq)
         write(*,'(A)') 'dq values:'
         do i = 1, size(dq)
            write(*,'(I4,1X,ES15.7)') i, dq(i)
         end do
      end if
   end if



   call args%destroy
   call py_t%destroy
   call kwargs%destroy
   call dq_nd%destroy
   call result%destroy
   call calcART%destroy

   deallocate(t)

   call forpy_finalize

end program forpy_calc_dq_example
