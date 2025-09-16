//! Chain complexes.
//! 
//! This module provides a pair of traits which define a standardized set of
//! functions/methods for working with (filtered) chain complexes. It also houses
//! the [barcode] module, which provides tools to work with the barcode of a persistent homology module.
//! 
//! Most of the methods available to analyze chain complexes are currently provided through
//! the [differential Umatch](crate::algebra::matrices::operations::umatch::differential)
//! module; additional tools may be available in the future.
//! 
//! 
//! 

use std::{hash::Hash, iter::{Chain, Flatten}};

use crate::{algebra::{matrices::query::MatrixAlgebra, vectors::entries::KeyValGet}, utilities::sequences_and_ordinals::BijectiveSequence};



pub mod barcode;









/// A trait representing a (based) chain complex.
/// 
/// A based chain complex is a chain complex 
/// 
/// ```text
/// ... <-- C_{n-1} <-- C_n <-- C_{n+1} <-- ... 
/// ```
/// equiped with a basis for each chain group `C_n`. The **differntial operator**
/// of the chain complex is the linear operator on the direct sum of the chain groups
/// which maps each vector in `C_n` to its boundary in `C_{n-1}`. Because bases are
/// given for each chain group, this operator can be represented by a matrix,
/// called the **differential matrix**.
/// 
/// 
/// To implement this trait, a type `T` must provide lookup methods to access both the
/// differential matrix and the (indices of the) basis vectors for each chain group.
/// Concretely, we require `T` to provide
/// 
/// - lookup methods to access the rows/columns/entries of the differential matrix.
///   This requirement is formalized by requiring `T` to implement the [MatrixOracle] and
///   [MatrixAlgebra] traits.
/// - a lookup method to access the complete list of (indices of the) basis vectors for 
///   each chain group.
/// - a lookup method to determine the homological dimension of each basis vector (
///   i.e. the chain group to which it belongs).
pub trait ChainComplex 
    
    where 
        Self: MatrixAlgebra
{

    type BasisVectorIndicesIterable: IntoIterator<Item = Self::RowIndex>;

    /// Returns an iterable that runs over the indices of the basis vectors.
    /// 
    /// Indices should be returned sorted in the same order they are asigned in the rows/columns of the differential matrix.
    /// That is, strictly increasing order with respect to `self.order_operator_for_row_indices()`.
    fn basis_vector_indices_for_dimension( &self, dimension: isize ) -> Self::BasisVectorIndicesIterable;


    /// Returns an iterable that runs over the indices of the basis vectors
    /// in the given dimensions.
    /// 
    /// **Dimensions are automatically sorted and deduplicated.** Internally, this function
    /// 
    /// - sorts the values `[d_1, d_2, ... ]` in `dimensions`
    /// - removes duplicates
    /// - calls `self.basis_vector_indices_for_dimension(d_i)` for each `d_i`
    /// - chains the resulting iterables together
    /// - returns the resulting iterable
    fn basis_vector_indices_for_dimensions< I >( &self, dimensions: I ) -> Flatten< std::vec::IntoIter< Self::BasisVectorIndicesIterable > >
    where
        I: IntoIterator<Item = isize>,
    {
        let mut dimensions: Vec<isize> = dimensions.into_iter().collect();
        dimensions.sort();
        dimensions.dedup();        
        dimensions.into_iter()
            .map( |d| self.basis_vector_indices_for_dimension(d) )
            .collect::<Vec<_>>()
            .into_iter()
            .flatten()
    }


    /// Returns a bijective sequence that runs over the indices of the basis vectors
    /// in the given dimensions.
    /// 
    /// **Dimensions are automatically sorted and deduplicated.** Internally, this function
    /// 
    /// - sorts the values `[d_1, d_2, ... ]` in `dimensions`
    /// - removes duplicates
    /// - calls `self.basis_vector_indices_for_dimension(d_i)` for each `d_i`
    /// - chains the resulting iterables together
    /// 
    /// The resulting sequence of indices `[index_1, index_2, ... ]` is then placed in a [BijectiveSequence],
    /// which allows for easy mapping between indices and their positions in the sequence.
    /// 
    /// # Errors
    /// 
    /// If the resulting sequence of indices contains duplicates, then this function returns `Err(duplicate_index)`.
    /// Otherwise, it returns `Ok(bijective_sequence)`.
    fn basis_vector_index_bimap_for_dimensions< I >( &self, dimensions: I ) 
        -> 
        Result< 
            BijectiveSequence< Self::RowIndex >,
            Self::RowIndex,
        >
    where
        I: IntoIterator<Item = isize>,
        Self::RowIndex: Hash,
    {
        let indices: Vec< Self::RowIndex > = self.basis_vector_indices_for_dimensions(dimensions).into_iter().collect();
        BijectiveSequence::from_iter(indices)
    }


    /// Returns the dimension of the basis vector with the given index.
    fn dimension_for_basis_vector_with_index( &self, index: & Self::RowIndex ) -> Result<isize, Self::RowIndex>;



}    


/// A trait representing a (based) filtered chain complex.
/// 
/// This is a special "subtrait" of [ChainComplex] that adds the ability to
/// retrieve the filtration value for each basis vector in the chain complex.
pub trait FilteredChainComplex: ChainComplex {

    type FiltrationValue;

    /// Returns the filtration value for the basis vector with the given index.
    fn filtration_value_for_basis_vector_with_index( 
        &self, 
        index: & Self::RowIndex 
    ) ->    
        Result<
            Self::FiltrationValue, 
            Self::RowIndex 
        >;

    /// Returns the filtration value for a vector
    /// 
    /// Concertely, if the vector is represented by a collection of index-coefficient pairs
    /// `[(index_1, coefficient_1), (index_2, coefficient_2), ...]`,
    /// then the filtration value is the maximum of the filtration values of the indices.
    /// In this case the function returns `Ok(max_value)`.
    /// 
    /// # Errors
    /// 
    /// - If the list `index_1, index_2, .. ` contains an invalide index, then this function returns `Err(Some(index))`, where
    ///   `index` is the first invalid index in the list.
    /// - If the list is empty, then this function returns `Err(None)`.
    /// 
    fn filtration_value_for_vector< V >( 
        &self, 
        vector: V 
    ) -> 
        Result<
            Self::FiltrationValue, 
            Option< Self::RowIndex > 
        >   
    where
        V:                          IntoIterator< Item: KeyValGet< Key = Self::RowIndex > >,
        Self::FiltrationValue:      Ord,
    {

        let mut filtration_value = None;

        for entry in vector.into_iter() {
            match self.filtration_value_for_basis_vector_with_index(&entry.key()) {
                Ok(value) => {
                    filtration_value = filtration_value.max(Some(value));
                },
                Err(invalid_index) => return Err(Some(invalid_index)),
            }
        }

        match filtration_value {
            Some(value) => Ok(value),
            None => Err(None), // empty vector
        }
    }
}