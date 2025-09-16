//! Umatch factorization
//!
//! U-match is a form of matrix factorization described in the paper [U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co)homology](https://arxiv.org/abs/2108.08831), by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.
//! It is a matrix equation
//! 
//! ```text
//! TM = DS
//! ```
//! where `D` is an `m` x `n` matrix, `S` and `T` are upper triangular matrices with 1's on the diagonal, and `M` is a generalized matching matrix.  We call
//! 
//! - `D` the factored matrix
//! - `T` the target COMB
//! - `S` the source COMB
//! - `M` the matching matrix
//! 
//! # Solving kernels, images, and systems of equations
//! 
//! U-match factorization applies to several core problems in linear algebra:
//! 
//! - solving [`Dx = b`](crate::algebra::matrices::operations::umatch::row_major::Umatch::solve_dx_equals_b) and [`xD = b`](crate::algebra::matrices::operations::umatch::row_major::Umatch::solve_xd_equals_b) for `x`
//! - computing bases for the [image](crate::algebra::matrices::operations::umatch::row_major::Umatch::image) and [kernel](crate::algebra::matrices::operations::umatch::row_major::Umatch::kernel) of `D`
//! 
//! for example, `{T[:,i] : M[i,j] != 0 }` is a basis for column space of `D`, and `{ S[:,j] : M[i,j] != 0 }` is a basis for the kernel of `D`.
//! 
//! OAT provides functionality to compute a U-match factorization, via the function  `Umatch::new`.
//! This produces an `Umatch` struct that provides the user with access to `M`, `S`, `S^{-1}`, `T`, and `T^{-1}`.
//! 
//! # Memory-efficiency
//! 
//! When OAT computes a U-match factorization, it stores only three items in memory: the factored matrix `D`, the generalized matching array `M`, and the square submatrix of `T^{-1}` indexed by row-pivot indices.
//! It can be [shown](https://arxiv.org/abs/2108.08831) that this information suffices to reconstruct any row or column of `M`, `S`, `S^{-1}`, `T`, and `T^{-1}` efficiently (row access is generally much faster than column access).
//! This makes U-match factorization highly memory efficient, in practice.
//! 

use crate::algebra::matrices::{query::MatrixOracle, types::matching::{GeneralizedMatchingMatrix, GeneralizedMatchingMatrixWithSequentialOrder}};

use std::hash::Hash;


// pub mod row_major_only;
pub mod gimbled;
pub mod row_major;
pub mod differential;
// pub mod developer;
// pub mod diagnostics;






















// pub trait UmatchOracle< MatrixToFactor > 
//     where
//         MatrixToFactor:     MatrixOracle<
//                                 ColumnIndex: Hash, // required for the generalized matching matrix
//                                 RowIndex:    Hash, // required for the generalized matching matrix
//                             >,
// {
//     type SourceComb;
//     type TargetComb;
//     type SourceCombInverse;
//     type TargetCombInverse;

//     fn source_comb(&self) -> Self::SourceComb;
//     fn target_comb(&self) -> Self::TargetComb;
//     fn source_comb_inverse(&self) -> Self::SourceCombInverse;
//     fn target_comb_inverse(&self) -> Self::TargetCombInverse;
//     fn generalized_matching_matrix<'a>(&'a self) -> &'a GeneralizedMatchingMatrixWithSequentialOrder<
//         MatrixToFactor::ColumnIndex,
//         MatrixToFactor::RowIndex,
//         MatrixToFactor::Coefficient
//     >;
// }