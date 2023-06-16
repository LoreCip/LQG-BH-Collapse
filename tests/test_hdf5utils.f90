subroutine test_hdf5utilsAll(T)

    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK), intent(out) :: T
    real(RK) :: T1, T2

    call cpu_time(T1)

    call test_generichdf5()

    call cpu_time(T2)
    T = T2 - T1

end subroutine test_hdf5utilsAll

subroutine test_generichdf5()

    use hdf5
    use iso_fortran_env, only: RK => real64
    use assert_statments

    implicit none

    integer :: status
    integer(8) :: dims(1)
    integer(C_INT) :: rank
    integer(hid_t) :: file_id, dataset_id, dataspace_id, plist_id
    integer(hsize_t) :: cdims(1)

    real(RK) :: ddata(5), rdata(5)
    
    
    real(RK), parameter :: tol = 1E-10

    ! Initialize data
    ddata = (/ 1, 2, 3, 4, 5 /)
    dims(1) = 1
    rank = size(dims)

    ! Create HDF5 file
    call h5open_f(status)
    call h5fcreate_f("test.h5", H5F_ACC_TRUNC_F, file_id, status)

    ! Create dataspace
    dims = 5
    rank = size(dims)
    cdims(1) = int(5./2.)

    call h5screate_simple_f(1, dims, dataspace_id, status)

    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, status)
    call h5pset_chunk_f(plist_id, rank, cdims, status)
    call h5pset_deflate_f(plist_id, 9, status)

    ! Create dataset
    call h5dcreate_f(file_id, "data", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)

    ! Write data to dataset
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, ddata, dims, status)

    ! Close dataset, dataspace, and file
    call h5dclose_f(dataset_id, status)
    call h5sclose_f(dataspace_id, status)
    call h5fclose_f(file_id, status)

    ! Read data from HDF5 file
    call h5fopen_f("test.h5", H5F_ACC_RDONLY_F, file_id, status)
    call h5dopen_f(file_id, "data", dataset_id, status)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, rdata, dims, status)

    ! Verify data
    call assertArrayEqual(ddata, rdata, tol, "General test of HDF5")

    ! Close dataset and file
    call h5dclose_f(dataset_id, status)
    call h5fclose_f(file_id, status)
    call h5close_f(status)

    ! Delete the HDF5 file
    call system("rm -f test.h5")
    
end subroutine test_generichdf5