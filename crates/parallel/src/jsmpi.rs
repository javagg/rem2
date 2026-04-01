use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = jsmpi, js_name = Init)]
    pub fn mpi_init() -> i32;

    #[wasm_bindgen(js_namespace = jsmpi, js_name = Finalize)]
    pub fn mpi_finalize() -> i32;

    #[wasm_bindgen(js_namespace = jsmpi, js_name = Comm_size)]
    pub fn mpi_comm_size(comm: i32, size: *mut i32) -> i32;

    #[wasm_bindgen(js_namespace = jsmpi, js_name = Comm_rank)]
    pub fn mpi_comm_rank(comm: i32, rank: *mut i32) -> i32;

    #[wasm_bindgen(js_namespace = jsmpi, js_name = Send)]
    pub fn mpi_send(buf: *const u8, count: i32, datatype: i32, dest: i32, tag: i32, comm: i32) -> i32;

    #[wasm_bindgen(js_namespace = jsmpi, js_name = Recv)]
    pub fn mpi_recv(buf: *mut u8, count: i32, datatype: i32, source: i32, tag: i32, comm: i32, status: *mut i32) -> i32;

    #[wasm_bindgen(js_namespace = jsmpi, js_name = Bcast)]
    pub fn mpi_bcast(buffer: *mut u8, count: i32, datatype: i32, root: i32, comm: i32) -> i32;

    #[wasm_bindgen(js_namespace = jsmpi, js_name = Barrier)]
    pub fn mpi_barrier(comm: i32) -> i32;

    #[wasm_bindgen(js_namespace = jsmpi, js_name = Allreduce)]
    pub fn mpi_allreduce(sendbuf: *const u8, recvbuf: *mut u8, count: i32, datatype: i32, op: i32, comm: i32) -> i32;
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
