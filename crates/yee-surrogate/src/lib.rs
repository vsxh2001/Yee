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

/// Trivial baseline: returns the output of the training sample whose
/// parameter vector is closest in Euclidean L2 distance.
#[derive(Debug, Clone, Default)]
pub struct NearestNeighbor {
    samples: Vec<Sample>,
    param_dim: usize,
    output_dim: usize,
}

impl NearestNeighbor {
    /// Construct an untrained nearest-neighbor surrogate.
    pub fn new() -> Self {
        Self::default()
    }
}

impl Surrogate for NearestNeighbor {
    fn train(&mut self, dataset: &Dataset) -> Result<()> {
        if dataset.is_empty() {
            return Err(Error::EmptyDataset);
        }
        let param_dim = dataset.samples[0].params.len();
        let output_dim = dataset.samples[0].output.len();
        for s in &dataset.samples {
            if s.output.len() != output_dim {
                return Err(Error::InconsistentOutput {
                    expected: output_dim,
                    got: s.output.len(),
                });
            }
        }
        self.samples = dataset.samples.clone();
        self.param_dim = param_dim;
        self.output_dim = output_dim;
        Ok(())
    }

    fn predict(&self, params: &[f64]) -> Result<Vec<Complex64>> {
        if params.len() != self.param_dim {
            return Err(Error::DimMismatch {
                training: self.param_dim,
                query: params.len(),
            });
        }
        let nearest = self
            .samples
            .iter()
            .min_by(|a, b| {
                let da = a
                    .params
                    .iter()
                    .zip(params)
                    .map(|(x, y)| (x - y).powi(2))
                    .sum::<f64>();
                let db = b
                    .params
                    .iter()
                    .zip(params)
                    .map(|(x, y)| (x - y).powi(2))
                    .sum::<f64>();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("at least one sample after train()");
        Ok(nearest.output.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(params: Vec<f64>, output: Vec<Complex64>) -> Sample {
        Sample { params, output }
    }

    #[test]
    fn empty_dataset_fails_to_train() {
        let mut surr = NearestNeighbor::new();
        let ds = Dataset::new();
        assert!(matches!(surr.train(&ds), Err(Error::EmptyDataset)));
    }

    #[test]
    fn dim_mismatch_fails_predict() {
        let mut surr = NearestNeighbor::new();
        let mut ds = Dataset::new();
        ds.push(sample(vec![1.0, 2.0], vec![Complex64::new(0.5, 0.0)]));
        surr.train(&ds).unwrap();
        assert!(matches!(surr.predict(&[1.0]), Err(Error::DimMismatch { .. })));
    }

    #[test]
    fn inconsistent_output_fails() {
        let mut surr = NearestNeighbor::new();
        let mut ds = Dataset::new();
        ds.push(sample(vec![1.0], vec![Complex64::new(0.1, 0.0); 5]));
        ds.push(sample(vec![2.0], vec![Complex64::new(0.2, 0.0); 7]));
        assert!(matches!(
            surr.train(&ds),
            Err(Error::InconsistentOutput { .. })
        ));
    }

    #[test]
    fn nearest_neighbor_returns_closest_sample() {
        let mut surr = NearestNeighbor::new();
        let mut ds = Dataset::new();
        ds.push(sample(vec![0.0, 0.0], vec![Complex64::new(0.1, 0.0)]));
        ds.push(sample(vec![1.0, 1.0], vec![Complex64::new(0.9, 0.0)]));
        ds.push(sample(vec![10.0, 10.0], vec![Complex64::new(0.5, 0.5)]));
        surr.train(&ds).unwrap();
        let pred = surr.predict(&[0.1, 0.1]).unwrap();
        assert_eq!(pred[0], Complex64::new(0.1, 0.0));
    }
}
