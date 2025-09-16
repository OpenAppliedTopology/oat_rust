//! Vec-of-Vec sparse matrices (row-major)
//! 
//! A "vec-of-vec" is a vector of vectors representing the nonzero entries of a sparse matrix.
//! For example, if `A` is the following matrix
//! 
//! ```text
//! 5 6
//! 7 0
//! ```
//! 
//! Then the row-major vec-of-vec representation of `A` would be 
//! 
//! ```
//! vec![
//!     vec![ (0, 5), (1, 6) ],
//!     vec![ (0, 7) ],
//! ];
//! ```
//! 
//! Compare this with the column-major vec-of-vec representation of `A`, which would be
//! 
//! ```
//! vec![
//!     vec![ (0, 5), (1, 7) ],
//!     vec![ (0, 6) ],
//! ];
//! ```
//! 
//! This module contains several variants of the vec-of-vec data structure for use
//! with the [matrix oracle traits](crate::algebra::matrices::matrix_oracle_traits).
//! 
//! # Design elements
//! 
//! #### Why `VecOfVec` doesn't implement [MatrixOracle](crate::algebra::matrices::query::MatrixOracle) but `&'a VecOfVec` does.
//! 
//! In brief, because Rust's lifetime checker does not allow `VecOfVec` to implement [MatrixOracle](crate::algebra::matrices::query::MatrixOracle).
//! The developers have found no way around this hurdle without resorting to methods that are clearly memory-inefficient.
//! 
//! #### Sorting
//! 
//! All sparse matrix data structure store lists of (strictly) sorted lists. Sorting is enforced in
//! order to compliance with the sorting requirements in the [MatrixOracle](crate::algebra::matrices::query::MatrixOracle) trait.
//! 
//! #### Modifying entries
//! 
//! Due to the sorting requirements, the user is highly limited in their ability to modify that data stored internally by a `VecOfVec`.
//!         The most general route is to retrieve the internally stored `Vec<Vec<EntryType>>`, modify it, then 
//!         wrap it in a new `VecOfVec` struct.  Note, however, that this will incur the cost of (i) re-sorting
//!         each internal vector, or (ii) verifying that each internal vector is sorted, if it is already sorted.
//! 
//! 
//! 

pub mod sorted;
// pub mod sorted_ref;
pub mod sorted_custom;