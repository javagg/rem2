use rem_parallel::{Comm, NoComm};

#[test]
fn test_nocomm() {
    let comm = NoComm;
    assert_eq!(comm.rank(), 0);
    assert_eq!(comm.size(), 1);

    let mut data = [1, 2, 3];
    comm.bcast_u8(&mut data, 0);
    assert_eq!(data, [1, 2, 3]);
}

// jsmpi tests require a WASM environment and cannot run as native unit tests.
// They will be verified during the browser-based integration phase.
