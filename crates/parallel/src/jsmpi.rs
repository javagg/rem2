#[cfg(target_arch = "wasm32")]
use js_sys::{global, Reflect};

#[cfg(target_arch = "wasm32")]
fn global_i32(name: &str, default: i32) -> i32 {
    Reflect::get(&global(), &wasm_bindgen::JsValue::from_str(name))
        .ok()
        .and_then(|v| v.as_f64())
        .map(|v| v as i32)
        .unwrap_or(default)
}

pub fn mpi_init() -> i32 {
    0
}

pub fn mpi_finalize() -> i32 {
    0
}

pub fn mpi_comm_size(_comm: i32, size: *mut i32) -> i32 {
    #[cfg(target_arch = "wasm32")]
    let detected = global_i32("__jsmpi_size", 1);
    #[cfg(not(target_arch = "wasm32"))]
    let detected = 1;

    unsafe {
        if !size.is_null() {
            *size = detected;
        }
    }
    0
}

pub fn mpi_comm_rank(_comm: i32, rank: *mut i32) -> i32 {
    #[cfg(target_arch = "wasm32")]
    let detected = global_i32("__jsmpi_rank", 0);
    #[cfg(not(target_arch = "wasm32"))]
    let detected = 0;

    unsafe {
        if !rank.is_null() {
            *rank = detected;
        }
    }
    0
}

pub fn mpi_send(_buf: *const u8, _count: i32, _datatype: i32, _dest: i32, _tag: i32, _comm: i32) -> i32 {
    0
}

pub fn mpi_recv(_buf: *mut u8, _count: i32, _datatype: i32, _source: i32, _tag: i32, _comm: i32, _status: *mut i32) -> i32 {
    0
}

pub fn mpi_bcast(_buffer: *mut u8, _count: i32, _datatype: i32, _root: i32, _comm: i32) -> i32 {
    0
}

pub fn mpi_barrier(_comm: i32) -> i32 {
    0
}

pub fn mpi_allreduce(sendbuf: *const u8, recvbuf: *mut u8, count: i32, datatype: i32, _op: i32, _comm: i32) -> i32 {
    // Single-rank fallback: recv = send.
    if sendbuf.is_null() || recvbuf.is_null() || count <= 0 {
        return 0;
    }

    let elem_size = match datatype {
        MPI_BYTE => 1usize,
        MPI_INT => core::mem::size_of::<i32>(),
        MPI_DOUBLE => core::mem::size_of::<f64>(),
        _ => 1usize,
    };

    let byte_len = (count as usize).saturating_mul(elem_size);
    unsafe {
        core::ptr::copy_nonoverlapping(sendbuf, recvbuf, byte_len);
    }
    0
}

// MPI Ops
pub const MPI_SUM: i32 = 0;

// MPI Constants defined in jsmpi
pub const MPI_COMM_WORLD: i32 = 0;
pub const MPI_BYTE: i32 = 1;
pub const MPI_INT: i32 = 2;
pub const MPI_DOUBLE: i32 = 3;
pub const MPI_ANY_SOURCE: i32 = -1;
pub const MPI_ANY_TAG: i32 = -1;
