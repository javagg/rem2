//! Raw FFI bindings to system MPI library (Open MPI / MPICH).
//!
//! Enabled via `native-mpi` feature. These are thin wrappers around the
//! C MPI API, using the same pattern as `jsmpi` for WASM targets.

/// C MPI function signatures (linked at build time via `-lmpi`).
extern "C" {
    fn MPI_Init(argc: *mut i32, argv: *mut *mut *mut u8) -> i32;
    fn MPI_Finalize() -> i32;
    fn MPI_Finalized(flag: *mut i32) -> i32;
    fn MPI_Comm_rank(comm: u32, rank: *mut i32) -> i32;
    fn MPI_Comm_size(comm: u32, size: *mut i32) -> i32;
    fn MPI_Barrier(comm: u32) -> i32;
    fn MPI_Bcast(
        buffer: *mut std::ffi::c_void,
        count: i32,
        datatype: u32,
        root: i32,
        comm: u32,
    ) -> i32;
    fn MPI_Allreduce(
        sendbuf: *const std::ffi::c_void,
        recvbuf: *mut std::ffi::c_void,
        count: i32,
        datatype: u32,
        op: u32,
        comm: u32,
    ) -> i32;
}

/// MPI type constants (matching Open MPI / MPICH).
pub const MPI_COMM_WORLD: u32 = 0x44000000;
pub const MPI_BYTE: u32       = 0x4c00010d;
pub const MPI_DOUBLE: u32     = 0x4c00080b;
pub const MPI_SUM: u32        = 0x58000003;

/// User-friendly native MPI initialisation.
pub(crate) fn mpi_init() {
    unsafe {
        let mut flag: i32 = 0;
        MPI_Finalized(&mut flag);
        if flag == 0 {
            MPI_Init(std::ptr::null_mut(), std::ptr::null_mut());
        }
    }
}

/// MPI finalisation.
pub(crate) fn mpi_finalize() {
    unsafe {
        let mut flag: i32 = 0;
        MPI_Finalized(&mut flag);
        if flag == 0 {
            MPI_Finalize();
        }
    }
}

/// Return the rank of the calling process in `comm`.
pub(crate) fn mpi_comm_rank(comm: u32, rank: &mut i32) {
    unsafe { MPI_Comm_rank(comm, rank); }
}

/// Return the number of processes in `comm`.
pub(crate) fn mpi_comm_size(comm: u32, size: &mut i32) {
    unsafe { MPI_Comm_size(comm, size); }
}

/// Barrier synchronisation.
pub(crate) fn mpi_barrier(comm: u32) {
    unsafe { MPI_Barrier(comm); }
}

/// Broadcast `data` from `root` to all ranks in `comm`.
pub(crate) fn mpi_bcast(data: *mut std::ffi::c_void, count: i32, datatype: u32, root: i32, comm: u32) {
    unsafe { MPI_Bcast(data, count, datatype, root, comm); }
}

/// Element-wise sum across all ranks in `comm`.
pub(crate) fn mpi_allreduce(
    send: *const std::ffi::c_void,
    recv: *mut std::ffi::c_void,
    count: i32,
    datatype: u32,
    op: u32,
    comm: u32,
) {
    unsafe { MPI_Allreduce(send, recv, count, datatype, op, comm); }
}
