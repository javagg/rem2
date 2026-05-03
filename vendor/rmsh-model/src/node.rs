use serde::{Deserialize, Serialize};

use crate::Point3;

/// A mesh node with a unique ID and 3D coordinates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node {
    pub id: u64,
    pub position: Point3,
}

impl Node {
    pub fn new(id: u64, x: f64, y: f64, z: f64) -> Self {
        Self {
            id,
            position: Point3::new(x, y, z),
        }
    }
}
