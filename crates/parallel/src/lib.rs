pub mod jsmpi;
#[cfg(feature = "native-mpi")]
mod native_mpi;

use crate::jsmpi::*;

/// MPI-style communicator trait for REM parallel solvers.
pub trait Comm {
    fn rank(&self) -> i32;
    fn size(&self) -> i32;
    fn barrier(&self);
    fn bcast_u8(&self, data: &mut [u8], root: i32);
    fn allreduce_f64(&self, val: f64) -> f64;
    fn allreduce_f64_vec(&self, data: &mut [f64]);
}

/// World communicator implementation using jsmpi (WASM target).
pub struct WorldComm;

impl WorldComm {
    pub fn new() -> Self {
        mpi_init();
        Self
    }
}

impl Drop for WorldComm {
    fn drop(&mut self) {
        mpi_finalize();
    }
}

impl Comm for WorldComm {
    fn rank(&self) -> i32 {
        let mut rank = 0;
        mpi_comm_rank(MPI_COMM_WORLD, &mut rank);
        rank
    }

    fn size(&self) -> i32 {
        let mut size = 0;
        mpi_comm_size(MPI_COMM_WORLD, &mut size);
        size
    }

    fn barrier(&self) {
        mpi_barrier(MPI_COMM_WORLD);
    }

    fn bcast_u8(&self, data: &mut [u8], root: i32) {
        mpi_bcast(
            data.as_mut_ptr(),
            data.len() as i32,
            MPI_BYTE,
            root,
            MPI_COMM_WORLD,
        );
    }

    fn allreduce_f64(&self, val: f64) -> f64 {
        let mut res = 0.0f64;
        mpi_allreduce(
            &val as *const f64 as *const u8,
            &mut res as *mut f64 as *mut u8,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            MPI_COMM_WORLD,
        );
        res
    }

    fn allreduce_f64_vec(&self, data: &mut [f64]) {
        let mut res = data.to_vec();
        mpi_allreduce(
            data.as_ptr() as *const u8,
            res.as_mut_ptr() as *mut u8,
            data.len() as i32,
            MPI_DOUBLE,
            MPI_SUM,
            MPI_COMM_WORLD,
        );
        data.copy_from_slice(&res);
    }
}

/// Native MPI world communicator (feature-gated behind `native-mpi`).
///
/// Uses raw MPI FFI to the system MPI library (Open MPI / MPICH)
/// for true multi-process parallelism on HPC clusters.
#[cfg(feature = "native-mpi")]
pub struct MpiWorldComm;

#[cfg(feature = "native-mpi")]
impl MpiWorldComm {
    pub fn new() -> Self {
        native_mpi::mpi_init();
        Self
    }
}

#[cfg(feature = "native-mpi")]
impl Drop for MpiWorldComm {
    fn drop(&mut self) {
        native_mpi::mpi_finalize();
    }
}

#[cfg(feature = "native-mpi")]
impl Comm for MpiWorldComm {
    fn rank(&self) -> i32 {
        let mut rank: i32 = 0;
        native_mpi::mpi_comm_rank(native_mpi::MPI_COMM_WORLD, &mut rank);
        rank
    }

    fn size(&self) -> i32 {
        let mut size: i32 = 0;
        native_mpi::mpi_comm_size(native_mpi::MPI_COMM_WORLD, &mut size);
        size
    }

    fn barrier(&self) {
        native_mpi::mpi_barrier(native_mpi::MPI_COMM_WORLD);
    }

    fn bcast_u8(&self, data: &mut [u8], root: i32) {
        native_mpi::mpi_bcast(
            data.as_mut_ptr() as *mut std::ffi::c_void,
            data.len() as i32,
            native_mpi::MPI_BYTE,
            root,
            native_mpi::MPI_COMM_WORLD,
        );
    }

    fn allreduce_f64(&self, val: f64) -> f64 {
        let mut res: f64 = 0.0;
        native_mpi::mpi_allreduce(
            &val as *const f64 as *const std::ffi::c_void,
            &mut res as *mut f64 as *mut std::ffi::c_void,
            1,
            native_mpi::MPI_DOUBLE,
            native_mpi::MPI_SUM,
            native_mpi::MPI_COMM_WORLD,
        );
        res
    }

    fn allreduce_f64_vec(&self, data: &mut [f64]) {
        let n = data.len();
        let mut recv = vec![0.0f64; n];
        native_mpi::mpi_allreduce(
            data.as_ptr() as *const std::ffi::c_void,
            recv.as_mut_ptr() as *mut std::ffi::c_void,
            n as i32,
            native_mpi::MPI_DOUBLE,
            native_mpi::MPI_SUM,
            native_mpi::MPI_COMM_WORLD,
        );
        data.copy_from_slice(&recv);
    }
}

/// Dummy communicator for single-threaded/non-MPI targets.
pub struct NoComm;

impl Comm for NoComm {
    fn rank(&self) -> i32 { 0 }
    fn size(&self) -> i32 { 1 }
    fn barrier(&self) {}
    fn bcast_u8(&self, _data: &mut [u8], _root: i32) {}
    fn allreduce_f64(&self, val: f64) -> f64 { val }
    fn allreduce_f64_vec(&self, _data: &mut [f64]) {}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_comm_is_single_rank() {
        let comm = NoComm;
        assert_eq!(comm.rank(), 0);
        assert_eq!(comm.size(), 1);
    }

    #[test]
    fn no_comm_allreduce_identity() {
        let comm = NoComm;
        assert_eq!(comm.allreduce_f64(42.0), 42.0);
        let mut v = vec![1.0, 2.0, 3.0];
        comm.allreduce_f64_vec(&mut v);
        assert_eq!(v, vec![1.0, 2.0, 3.0]);
    }
}
