use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::element::Element;
use crate::node::Node;

/// A finite element mesh containing nodes and elements.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Mesh {
    /// All nodes indexed by their ID.
    pub nodes: HashMap<u64, Node>,
    /// All elements.
    pub elements: Vec<Element>,
    /// Physical group names: (dimension, tag) -> name.
    pub physical_names: HashMap<(i32, i32), String>,
}

impl Mesh {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_node(&mut self, node: Node) {
        self.nodes.insert(node.id, node);
    }

    pub fn add_element(&mut self, element: Element) {
        self.elements.push(element);
    }

    /// Get the bounding box of all nodes: (min, max).
    pub fn bounding_box(&self) -> Option<(crate::Point3, crate::Point3)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut min = crate::Point3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut max = crate::Point3::new(f64::MIN, f64::MIN, f64::MIN);
        for node in self.nodes.values() {
            min.x = min.x.min(node.position.x);
            min.y = min.y.min(node.position.y);
            min.z = min.z.min(node.position.z);
            max.x = max.x.max(node.position.x);
            max.y = max.y.max(node.position.y);
            max.z = max.z.max(node.position.z);
        }
        Some((min, max))
    }

    /// Get the center of the bounding box.
    pub fn center(&self) -> crate::Point3 {
        match self.bounding_box() {
            Some((min, max)) => nalgebra::center(&min, &max),
            None => crate::Point3::origin(),
        }
    }

    /// Get the diagonal length of the bounding box.
    pub fn diagonal_length(&self) -> f64 {
        match self.bounding_box() {
            Some((min, max)) => nalgebra::distance(&min, &max),
            None => 1.0,
        }
    }

    /// Filter elements by topological dimension.
    pub fn elements_by_dimension(&self, dim: u8) -> Vec<&Element> {
        self.elements
            .iter()
            .filter(|e| e.dimension() == dim)
            .collect()
    }

    /// Get node count.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Get element count.
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }
}
