subroutine create_hdf5_file(file_name, NX, Xgrid, m, dx, file_id)
    use hdf5
    use iso_fortran_env, only: RK => real64

    implicit none
  
    integer,                 intent(in) :: NX
    real(RK),                intent(in) :: m, dx
    real(RK), dimension(NX), intent(in) :: Xgrid
  
    character(len=*), intent(in)  :: file_name
    integer(hid_t)  , intent(out) :: file_id
    
    integer(hid_t) :: dataspace_id, group_id, dataset_id
    integer(C_INT) :: rank
    integer(8) :: dims(1)
    integer :: status
    character(len=100) :: dataset_name
  
    ! Open the HDF5 library
    call h5open_f(status)
  
    ! Create a new HDF5 file
    call h5fcreate_f(trim(file_name), H5F_ACC_TRUNC_F, file_id, status)
  
    ! Check for errors
    if (status /= 0) then
      write(*,*) 'Error creating HDF5 file:', trim(file_name)
      return
    end if
  
    ! SAVE THE GRID

    dims = shape(Xgrid)
    rank = size(dims)
  
    call h5screate_simple_f(rank, dims, dataspace_id, status)
    write(dataset_name, "(A)") "Xgrid"
    call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, Xgrid, dims, status) 

    ! SAVE SOME PARAMETERS

    dims = 1
    rank = 0

    call h5screate_simple_f(rank, dims, dataspace_id, status)
    write(dataset_name, "(A)") "Mass"
    call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, m, dims, status)

    write(dataset_name, "(A)") "dx"
    call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dx, dims, status)
  
    ! Close the datasets
    call h5dclose_f(dataset_id, status)
    ! Close the dataspace
    call h5sclose_f(dataspace_id, status)
  
    ! Check for errors
    if (status /= 0) then
        write(*,*) 'Error saving to HDF5 file:', trim(file_name)
        return
    end if
end subroutine create_hdf5_file


subroutine close_hdf5_file(file_id)
    use hdf5
    implicit none
    
    integer(hid_t), intent(in) :: file_id
    integer :: status
    
    ! Close the HDF5 file
    call h5fclose_f(file_id, status)
    
    ! Check for errors
    if (status /= 0) then
      write(*,*) 'Error closing HDF5 file'
    end if
    
    ! Close the HDF5 library
    call h5close_f(status)
    
end subroutine close_hdf5_file


subroutine write_arrays_to_hdf5(file_id, NX, B, E, rho, t, dt, iteration)

    use iso_fortran_env, only: RK => real64
    use hdf5
    implicit none
    
    integer,                 intent(in) :: iteration, NX
    integer(hid_t),          intent(in) :: file_id
    real(RK), dimension(NX), intent(in) :: B, E, rho
    real(RK),                intent(in) :: t, dt

    integer(hid_t) :: dataspace_id, group_id, dataset_id
    integer(C_INT) :: rank, status
    integer(8) :: dims(1)
    character(len=100) :: group_name, dataset_name
    

    write(group_name, "(I0)") iteration
    call h5gcreate_f(file_id, trim(group_name), group_id, status)

    ! SAVE ARRAYS

    dims = shape(B)
    rank = size(dims)
    
    call h5screate_simple_f(rank, dims, dataspace_id, status)

    write(dataset_name, "(A)") "B"
    call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, B, dims, status) 

    write(dataset_name, "(A)") "e^b"    
    call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, E, dims, status)

    write(dataset_name, "(A)") "rho"    
    call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, rho, dims, status)

    ! SAVE SCALARS

    dims = 1
    rank = 0

    call h5screate_simple_f(rank, dims, dataspace_id, status)
    write(dataset_name, "(A)") "t"
    call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, t, dims, status)

    write(dataset_name, "(A)") "dt"
    call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dt, dims, status)


    ! Close the datasets
    call h5dclose_f(dataset_id, status)
    ! Close the dataspace
    call h5sclose_f(dataspace_id, status)
    ! Close the group
    call h5gclose_f(group_id, status)

    if (status /= 0) then
        write(*,*) 'Error writing HDF5 file'
    end if

  end subroutine write_arrays_to_hdf5
  
