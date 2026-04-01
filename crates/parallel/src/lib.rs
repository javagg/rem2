pub mod jsmpi;

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

/// World communicator implementation using jsmpi.
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
