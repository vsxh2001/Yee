//! ML surrogate models for parameterized electromagnetic simulation
//! outputs.
//!
//! ## Phase 3 walking-skeleton scope
//!
//! - `Sample`: a single (parameters, output) pair recorded from a full solver run.
//! - `Dataset`: a growing collection of samples.
//! - `Surrogate` trait: `train(dataset)` + `predict(params)`.
//! - `NearestNeighbor`: trivial baseline surrogate. Returns the output of the
//!   closest sample by Euclidean distance in parameter space.
//!
//! Phase 3.1+ will add Gaussian-process, MLP, and Fourier neural operator
//! backends behind the same `Surrogate` trait.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use num_complex::Complex64;

/// Single labeled training sample: parameters -> output.
#[derive(Debug, Clone)]
pub struct Sample {
    /// Continuous design parameters (e.g. patch length, width, dielectric thickness).
    pub params: Vec<f64>,
    /// Output complex S-parameter vector across all frequencies (length F).
    pub output: Vec<Complex64>,
}

/// Collection of labeled samples.
#[derive(Debug, Clone, Default)]
pub struct Dataset {
    /// Underlying sample storage.
    pub samples: Vec<Sample>,
}

impl Dataset {
    /// Construct an empty dataset.
    pub fn new() -> Self {
        Self::default()
    }
    /// Append a sample to the dataset.
    pub fn push(&mut self, sample: Sample) {
        self.samples.push(sample);
    }
    /// Number of samples in the dataset.
    pub fn len(&self) -> usize {
        self.samples.len()
    }
    /// Whether the dataset has zero samples.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }
}

/// Surrogate-layer errors.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Training with an empty dataset.
    #[error("training failed: dataset is empty")]
    EmptyDataset,
    /// Prediction parameter vector length does not match training samples.
    #[error("parameter dimensionality mismatch: training={training}, query={query}")]
    DimMismatch {
        /// Parameter dimensionality the surrogate was trained on.
        training: usize,
        /// Parameter dimensionality of the query.
        query: usize,
    },
    /// Output vectors in the dataset have inconsistent lengths.
    #[error("inconsistent output length: expected {expected}, got {got}")]
    InconsistentOutput {
        /// Output length implied by the first sample.
        expected: usize,
        /// Output length of the offending sample.
        got: usize,
    },
}

/// Convenience `Result` alias for surrogate operations.
pub type Result<T> = core::result::Result<T, Error>;

/// Surrogate model trait: train on a `Dataset`, predict for new parameter
/// vectors.
pub trait Surrogate {
    /// Fit the surrogate to the dataset.
    fn train(&mut self, dataset: &Dataset) -> Result<()>;
    /// Predict the output vector for the given parameter vector.
    fn predict(&self, params: &[f64]) -> Result<Vec<Complex64>>;
}
