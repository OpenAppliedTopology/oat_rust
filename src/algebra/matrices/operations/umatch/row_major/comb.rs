//! Columar ordered matching bases
//! 
//! The notion of a columnar ordered matching basis (COMB) was introduced in [U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co)homology](https://arxiv.org/abs/2108.08831), by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.
//! 
//! This module defines structs for the domain and target COMBs, as well as their inverses.
//! 
//! # Design notes
//! 
//! We also define wrapper structs for rows and columns of the COMBs.  We use structs and not type
//! aliases because structs allow one to take advantage of type bounds on type parameters
//! (specifically, bounds that ensure existence of associated types), whereas parameters on type
//! aliases can't currently be constrained with `where` clauses.


// use debugit::DebugIt as D;
// use debugit::debugit;

//  =========================================================================================================
//  COMB'S (COLUMNAR ORDERED MATCHING BASES -- UNCOMPRESSED)
//  =========================================================================================================

use super::Umatch;
use super::construction::*;

use crate::algebra::matrices::operations::solve::echelon::{RowEchelonSolverReindexed};
use crate::algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix, multiply_column_vector_with_matrix_and_return_reversed};
use crate::algebra::matrices::operations::transform_vector_wise::{OnlyColumnIndicesInsideCollection};
use crate::algebra::matrices::operations::solve::triangle::{TriangularSolveForColumnVectorReverse};





use crate::algebra::matrices::operations::combine_rows_and_columns::LinearCombinationOfColumnsReverse;
use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::matrices::query::{MatrixAlgebra, MatrixOracle};
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet, KeyValNew, KeyValPair, ReindexEntry};
use crate::algebra::rings::traits::{SemiringOperations, RingOperations, DivisionRingOperations};


use crate::algebra::vectors::operations::RemapEntryTuple;
use crate::utilities::functions::evaluate::{ IdentityFunction, EvaluateFunctionFnMutWrapper };
use crate::utilities::iterators::general::{TwoTypeIterator, IterWrappedVec, OncePeekable, MapByTransform};
use crate::utilities::iterators::merge::hit::{IteratorsMergedInSortedOrder, hit_merge_by_predicate};
use crate::utilities::iterators::merge::two_type::MergeTwoIteratorsByOrderOperator;
use crate::utilities::order::{JudgeOrder, JudgePartialOrder, OrderOperatorByKey, ReverseOrder};


use crate::algebra::vectors::operations::{Scale, VectorOperations, Simplify, OnlyIndicesInsideCollection, OnlyIndicesOutsideCollection, ChangeIndexSimple, ChangeEntryType, LinearCombinationSimplified};

use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;
use std::iter::Rev;
use std::iter::{Cloned, Peekable, Once, Chain};
use std::marker::PhantomData;
use std::slice::{Iter};
use std::vec::IntoIter;

use derive_getters::Dissolve;
use derive_new::new;
use itertools::Itertools;





// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the target COMB in a lazy fashion.
#[derive(Copy, Clone, new, Dissolve)]
pub struct TargetComb< 'a, MatrixToFactor, > 
    where   
        MatrixToFactor:                        MatrixAlgebra,
        MatrixToFactor::ColumnIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:              Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct 
{
    pub umatch:     &'a Umatch< MatrixToFactor > ,
}

// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the inverse of the target COMB in a lazy fashion.
#[derive(Copy, Clone, new, Dissolve)]
pub struct TargetCombInverse< 'a, MatrixToFactor, > 
    where   
        MatrixToFactor:                        MatrixAlgebra,
        MatrixToFactor::ColumnIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:              Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct 
{
    pub umatch:     &'a Umatch< MatrixToFactor > ,
}

// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the source COMB in a lazy fashion.
#[derive(Copy, Clone, new, Dissolve)]
pub struct SourceComb< 'a, MatrixToFactor, > 
    where   
        MatrixToFactor:                        MatrixAlgebra,
        MatrixToFactor::ColumnIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:              Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct 
{
    pub umatch:     &'a Umatch< MatrixToFactor > ,
}

// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the inverse of the source COMB in a lazy fashion.
#[derive(Copy, Clone, new, Dissolve)]
pub struct SourceCombInverse< 'a, MatrixToFactor, > 
    where   
        MatrixToFactor:                        MatrixAlgebra,
        MatrixToFactor::ColumnIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:              Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct 
{
    pub umatch:     &'a Umatch< MatrixToFactor > ,
}




//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: SOURCE
//  ---------------------------------------------------------------------------------------------------------




//  MATRIX ORACLE
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, MatrixToFactor >
    
    MatrixOracle for 
    
    SourceComb
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:                         MatrixAlgebra<
                                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                                    RingOperator:           DivisionRingOperations,
                                                    RowEntry:               KeyValPair,
                                                    ColumnEntry:            KeyValPair,        
                                                >,

{   
    type Coefficient            =   MatrixToFactor::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex               =   MatrixToFactor::ColumnIndex;       // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   MatrixToFactor::ColumnIndex;       // The type of column indices

    type RowEntry               =   MatrixToFactor::RowEntry   ;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   MatrixToFactor::RowEntry   ;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                    =   SourceCombRow< 'a, MatrixToFactor >;  // What you get when you ask for a row.
    type RowReverse             =   IntoIter< MatrixToFactor::RowEntry >;        // What you get when you ask for a row with the order of entries reversed
    type Column                 =   IntoIter< MatrixToFactor::RowEntry >;            // What you get when you ask for a column   
    type ColumnReverse          =   SourceCombColumnReverse< 'a, MatrixToFactor >;     // What you get when you ask for a column with the order of entries reversed                                

    /// The source COMB (and its inverse) has a column for index `j` iff the factored matrix has a column for index `j`    
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.umatch.matrix_to_factor.has_column_for_index(index)
    }

    /// The source COMB (and its inverse) has a row for index `j` iff the factored matrix has a column for index `j`
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        self.umatch.matrix_to_factor.has_column_for_index(index)
    }

    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > 
        {
            let row = self.row( row );
            let order_operator_for_column_indices = self.umatch.order_operator_for_column_indices();
            for entry in row {
                match order_operator_for_column_indices.judge_cmp( & entry.key(), column  ) {
                    Ordering::Equal => { return Some( entry.val() )  },
                    Ordering::Greater => { return None },
                    Ordering::Less => { continue }
                }
            }
            return None
        }  
    
    fn row(                     & self,  lookup_index: & MatrixToFactor::ColumnIndex    )       -> Self::Row
    {
        match self.umatch.matching.has_a_match_for_column_index( &lookup_index ) { 
            false => { 

                SourceCombRow{
                                iter_unwrapped:         TwoTypeIterator::Version1(
                                                                std::iter::once(  MatrixToFactor::RowEntry::new( lookup_index.clone(), MatrixToFactor::RingOperator::one() )  ) 
                                                            ),
                                phantom_arraymapping:   PhantomData,
                            }
            }
            true => { 

                let ring_operator = self.umatch.ring_operator();

                // The matrix A from Hang et al., "Umatch factorization ...", with rows indexed by row index ordinals
                // Note: this struct only contains a reference
                let seed = self.umatch.matched_blocks_of_target_comb_inverse_times_matrix_to_factor_oc();


                // Struct that encodes the matching from column indices to row index ordinals
                let generalized_matching_matrix_ref = self.umatch.generalized_matching_matrix_ref();
                let column_index_to_row_ordinal_wrapped = generalized_matching_matrix_ref.bijection_column_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal();

                // obtain a row of A^{-1}
                let seed_inv_row_vec = 
                    RowEchelonSolverReindexed::solve(
                            std::iter::once( MatrixToFactor::RowEntry::new( lookup_index.clone(), MatrixToFactor::RingOperator::one() ) ), // the standard unit vector supported on `column_index`
                            seed, // matrix A
                            column_index_to_row_ordinal_wrapped.clone(),
                            ring_operator.clone(),
                            self.umatch.order_operator_for_row_entries(),
                        )
                        .solution().unwrap();

                // obtain a copy of R_{\rho \rho}
                let comb_target_inv_matched_block = self.umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row();

                // obtain a copy of the column submatrix of the factored matrix indexed by unmatched column indices
                let mapping_npcols = self.umatch.matrix_to_factor_matchless_columns_only();

                // y = - seed_inv_row_vec * R_{\rho \rho}^{-1} * D_{ \rho \bar \kappa}
                let x 
                    =   multiply_row_vector_with_matrix(
                                seed_inv_row_vec
                                    .iter()
                                    .map(   |x| 
                                            (   
                                                generalized_matching_matrix_ref.ordinal_for_column_index(& x.key() ).unwrap(),  // replace column indices with row index indices
                                                ring_operator.negate( x.val() ) // multiply by -1 (recally that we are computing *MINUS* seed_inv_row_vec * R_{\rho \rho}^{-1} * D_{ \rho \bar \kappa})
                                            ) 
                                        ),
                                comb_target_inv_matched_block,
                                ring_operator.clone(),
                                OrderOperatorByKey::new(),
                            );

                let y
                    =   multiply_row_vector_with_matrix(
                                x.map( 
                                        | ( row_ordinal, snzval ) |
                                        ( 
                                                generalized_matching_matrix_ref.row_index_for_ordinal( row_ordinal ),
                                                snzval
                                            )
                                    ),
                                mapping_npcols,
                                ring_operator.clone(),
                                self.umatch.order_operator_for_row_entries(),                                
                            );

                // rescale entries of seed_inv_row_vec, transforming it into row `column_index` of `A^{-1} M_{\rho \kappa}`
                // NB: this has to occur AFTER we've used `seed_inv_row_vec` to construct `y`
                let mut seed_inv_row_vec_times_matching = seed_inv_row_vec;
                for entry in seed_inv_row_vec_times_matching.iter_mut() {
                    entry.set_val( 
                            ring_operator.multiply( 
                                    entry.val(),  
                                    self.umatch.generalized_matching_matrix_ref().coefficient_for_column_index( &entry.key() )
                                )                         
                        )
                }

                // merge seed_inv_row_vec
                let z = MergeTwoIteratorsByOrderOperator::new( 
                                                            IterWrappedVec::new( seed_inv_row_vec_times_matching ).peekable(), 
                                                            y.peekable(), 
                                                            self.umatch.order_operator_for_row_entries()
                                                        );

          
                SourceCombRow{
                                iter_unwrapped:         TwoTypeIterator::Version2( z ),
                                phantom_arraymapping:   PhantomData,
                            }

            }
        }
    }
    fn row_reverse(             &self,  index: &Self::RowIndex    )   -> Self::RowReverse           { 
        let mut vec = self.row(&index).collect_vec();
        (&mut vec).reverse();
        vec.into_iter() 
    }
    
    fn column(                  &self,  index: &Self::ColumnIndex )   -> Self::Column               { 
        let mut vec = self.column_reverse(&index).collect_vec();
        (&mut vec).reverse();
        vec.into_iter() 
    }
    fn column_reverse(          &self,  index: &Self::ColumnIndex )   -> Self::ColumnReverse
    {
        match self.umatch.matching.bijection_column_indices_to_ordinals_and_inverse().ordinal_for_element( & index ) {

            Some( row_ordinal ) => {                     
                let matched_snzval  =   self.umatch.matching.coefficient_for_ordinal(row_ordinal);
                // let matched_row_index = self.umatch.matching.row_row_index_for_ordinal(row_ordinal);
                let unit_vec = OncePeekable::new( MatrixToFactor::RowEntry::new( index.clone(), matched_snzval ) );
                let col = 
                    TriangularSolveForColumnVectorReverse::solve(
                        unit_vec,
                        self.umatch.target_comb_inverse_times_matrix_to_factor_matched_block_with_rows_indexed_by_matched_column_index(), // matrix A
                    ).ok().unwrap();
                // let col = ChangeEntryType::new( col ); // change entries from tuples to ArrayMappingRowEntry
                SourceCombColumnReverse{ 
                                iter_unwrapped:     TwoTypeIterator::Version1( col ) 
                            }
            }
            None => {
                // a column c of D[ matched_rows, : ] 
                let col_matched_rows_only 
                    =   self.umatch.matrix_to_factor_matched_rows_only().column_reverse( & index )
                            .negate( self.umatch.ring_operator() ); // this multiplies the column by -1
                            // .map( | entry | (self.umatch.matching.bimap_row.ordinal_for_element( entry.key() ), entry.val()) ); // reindex so that entries are indexed by ordinals of matched row indices
                // get an object that represents the REVERSED total order on column indices
                let order_operator_for_row_entries_reverse_total = self.umatch.order_operator_for_row_entries_reverse();
                // a linear combination v of columns of R_{\rho \rho}^{-1}; we reindex by column index, then sort by column index
                let mut lin_comb_r 
                    =   multiply_column_vector_with_matrix_and_return_reversed(
                                col_matched_rows_only, 
                                self.umatch.matched_block_of_target_comb_inverse(), 
                                self.umatch.ring_operator(), 
                                self.umatch.order_operator_for_column_entries()
                            )
                            .map(   |x| //|(key,val)| 
                                    MatrixToFactor::RowEntry::new(
                                            self.umatch.matching.column_index_for_row_index(&x.key() ).unwrap(),
                                            x.val(),
                                        )
                                )
                            .collect_vec(); // collect into a vector
                lin_comb_r.sort_by( |a,b| order_operator_for_row_entries_reverse_total.judge_partial_cmp( a, b ).unwrap() );  // sort the vector according to column index
                // let lin_comb_r = lin_comb_r.into_iter();
                // A^{-1} v  <-- this is equal to the matched part of the column we want to construct
                let matched_part = 
                    TriangularSolveForColumnVectorReverse::solve(
                        lin_comb_r, // the vector b in "solve Ax = b"
                        self.umatch.target_comb_inverse_times_matrix_to_factor_matched_block_with_rows_indexed_by_matched_column_index(), // matrix A
                    ).ok().unwrap();
                // change the entry type of this iterator from (RowIndex, Coefficient) to MatrixToFactor::ColumnEntry
                let matched_part =   ChangeEntryType::new( matched_part, ); // rust infers the new entry type that we want, automatically
                let col =   MergeTwoIteratorsByOrderOperator::new(
                                    matched_part.peekable(),
                                    OncePeekable::new( MatrixToFactor::RowEntry::new( index.clone(), MatrixToFactor::RingOperator::one() ) ),
                                    self.umatch.order_operator_for_row_entries_reverse(), // recall that we want entries returned in *descending* order
                                );                             
                SourceCombColumnReverse{ 
                                    iter_unwrapped:     TwoTypeIterator::Version2( col )
                                }                                
            }
        }
    }
} 




impl < 'a, MatrixToFactor >
    
    MatrixAlgebra for 
    
    SourceComb
        < 'a, MatrixToFactor >


    where   
        MatrixToFactor:                         MatrixAlgebra<
                                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                                    RingOperator:           DivisionRingOperations,
                                                    RowEntry:               KeyValPair,
                                                    ColumnEntry:            KeyValPair,        
                                                >,
{
    type RingOperator                           =   MatrixToFactor::RingOperator;

    type OrderOperatorForRowEntries             =   MatrixToFactor::OrderOperatorForRowEntries;

    type OrderOperatorForRowIndices             =   MatrixToFactor::OrderOperatorForColumnIndices;

    type OrderOperatorForColumnEntries          =   MatrixToFactor::OrderOperatorForRowEntries;

    type OrderOperatorForColumnIndices          =   MatrixToFactor::OrderOperatorForColumnIndices;

    /// The ring operator for the source COMB is the same as for the factored matrix
    fn ring_operator( &self ) -> Self::RingOperator {
        self.umatch.matrix_to_factor.ring_operator()
    }

    /// The order operators for row and column entries for the source COMB are the same as the order operator for the row entries of the factored matrix.
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.umatch.matrix_to_factor.order_operator_for_row_entries()
    }

    /// The order operators for row and column indices for the source COMB are the same as the order operator for the column indices of the factored matrix.    
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.umatch.matrix_to_factor.order_operator_for_column_indices()
    }

    /// The order operators for row and column entries for the source COMB are the same as the order operator for the row entries of the factored matrix.    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.umatch.matrix_to_factor.order_operator_for_row_entries()
    }

    /// The order operators for row and column indices for the source COMB are the same as the order operator for the column indices of the factored matrix.        
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.umatch.matrix_to_factor.order_operator_for_column_indices()
    }
}       






impl < 'a, MatrixToFactor >
    
    MatrixOracleOperations for 
    
    SourceComb
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:                     MatrixAlgebra,
        MatrixToFactor::ColumnIndex:        Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct        
{}        




//  DEFINE ROW
//  ---------------------------------------------------------------------------------------------------------

pub struct SourceCombRow< 'a, MatrixToFactor, >


    where   
        MatrixToFactor:                         MatrixAlgebra<
                                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                                    RingOperator:           DivisionRingOperations,
                                                    RowEntry:               KeyValPair,
                                                    ColumnEntry:            KeyValPair,        
                                                >,

{
    iter_unwrapped:     TwoTypeIterator< 
                                Once< MatrixToFactor::RowEntry >,
                                MergeTwoIteratorsByOrderOperator<
                                        Peekable< IterWrappedVec< MatrixToFactor::RowEntry > >, 
                                        Peekable< 
                                                LinearCombinationSimplified<
                                                        OnlyIndicesOutsideCollection< MatrixToFactor::Row, &'a HashMap<MatrixToFactor::ColumnIndex, usize> >, 
                                                        MatrixToFactor::RingOperator, MatrixToFactor::OrderOperatorForRowEntries
                                                    >, 
                                            >,
                                        MatrixToFactor::OrderOperatorForRowEntries
                                    >,
                            >,
    phantom_arraymapping:   PhantomData< MatrixToFactor >
}  


impl < 'a, MatrixToFactor, >

    Iterator for

    SourceCombRow< 
            'a, MatrixToFactor, 
        >   

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,        
{
    type Item = MatrixToFactor::RowEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}         


//  DEFINE (REVERSE) COLUMN
//  ---------------------------------------------------------------------------------------------------------

pub struct SourceCombColumnReverse< 
                    'a, MatrixToFactor, 
                >
    where   
        MatrixToFactor:                         MatrixAlgebra<
                                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                                    RingOperator:           DivisionRingOperations,
                                                    RowEntry:               KeyValPair,
                                                    ColumnEntry:            KeyValPair,        
                                                >,                
{
    iter_unwrapped:    TwoTypeIterator<     
                                TriangularSolveForColumnVectorReverse<
                                        Vec< MatrixToFactor::RowEntry >, 
                                        TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex
                                            <'a, MatrixToFactor >,                                                                                                                             
                                    >,                                                  
                                MergeTwoIteratorsByOrderOperator<
                                        Peekable< 
                                                ChangeEntryType<
                                                        TriangularSolveForColumnVectorReverse<
                                                                Vec< MatrixToFactor::RowEntry >,
                                                                TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex
                                                                    <'a, MatrixToFactor >,                                                                                                                           
                                                            >, 
                                                        MatrixToFactor::RowEntry,                                   
                                                    >,
                                            >,
                                        OncePeekable<MatrixToFactor::RowEntry>, 
                                        ReverseOrder< MatrixToFactor::OrderOperatorForRowEntries >
                                    >,        
                            >,
}

impl < 'a, MatrixToFactor, >

    Iterator for 

    SourceCombColumnReverse< 
            'a, MatrixToFactor, 
        >        
    where   
        MatrixToFactor:                         MatrixAlgebra<
                                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                                    RingOperator:           DivisionRingOperations,
                                                    RowEntry:               KeyValPair,
                                                    ColumnEntry:            KeyValPair,        
                                                >,                         
{
    type Item = MatrixToFactor::RowEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}            










//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: SOURCE INVERSE
//  ---------------------------------------------------------------------------------------------------------


//  IMPLEMENT MATRIX ORACLE
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, MatrixToFactor >
    
    MatrixOracle for 
    
    SourceCombInverse
        < 'a, MatrixToFactor >


    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,

{   
    type Coefficient            =   MatrixToFactor::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex               =   MatrixToFactor::ColumnIndex;       // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   MatrixToFactor::ColumnIndex;       // The type of column indices

    type RowEntry               =   MatrixToFactor::RowEntry   ;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   MatrixToFactor::RowEntry   ;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                    =   SourceCombInverseRow< MatrixToFactor >;  // What you get when you ask for a row.
    type RowReverse             =   IntoIter< MatrixToFactor::RowEntry >;        // What you get when you ask for a row with the order of entries reversed
    type Column                 =   IntoIter< MatrixToFactor::RowEntry >;            // What you get when you ask for a column   
    type ColumnReverse          =   IntoIter< MatrixToFactor::RowEntry >;     // What you get when you ask for a column with the order of entries reversed

    /// The source COMB (and its inverse) has a column for index `j` iff the factored matrix has a column for index `j`    
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.umatch.matrix_to_factor.has_column_for_index(index)
    }
    /// The source COMB (and its inverse) has a column for index `j` iff the factored matrix has a column for index `j`        
    fn has_column_for_index_result(  &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnIndex, Self::ColumnIndex > {
        self.umatch.matrix_to_factor.has_column_for_index_result(index)        
    }
    /// The source COMB (and its inverse) has a row for index `j` iff the factored matrix has a column for index `j`
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        self.umatch.matrix_to_factor.has_column_for_index(index)
    }
    /// The source COMB (and its inverse) has a row for index `j` iff the factored matrix has a column for index `j`    
    fn has_row_for_index_result(     &   self, index: & Self::RowIndex   )   -> Result< Self::RowIndex   , Self::RowIndex    > {
        self.umatch.matrix_to_factor.has_column_for_index_result(index)        
    }

    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { 
            let row = self.row( row );
            let order_operator_for_column_indices = self.umatch.order_operator_for_column_indices();
            for entry in row {
                match order_operator_for_column_indices.judge_cmp( & entry.key(), column  ) {
                    Ordering::Equal => { return Some( entry.val() )  },
                    Ordering::Greater => { return None },
                    Ordering::Less => { continue }
                }
            }
            return None
        }
    
    fn row(                     & self,  column_index: & MatrixToFactor::ColumnIndex    )       -> Self::Row
    {
        match self.umatch.matching.ordinal_for_column_index( &column_index ) { // this function call looks strange b/c the row indices of the source COMB have type MatrixToFactor::ColumnIndex!  
            None => { 
                // this is a non-pivot row, so return a unit vector
                SourceCombInverseRow{
                                iter_unwrapped:     TwoTypeIterator::Version1(
                                                            std::iter::once(  MatrixToFactor::RowEntry::new( column_index.clone(), MatrixToFactor::RingOperator::one() )  ) 
                                                        )
                    }
            }
            Some( row_ordinal ) => {  
                let scalar  =   self.umatch.matching.coefficient_for_ordinal( row_ordinal );
                let scalar_inv      =   self.umatch.ring_operator().invert( scalar );
                // let key_to_target_comb    =   self.umatch.matching.row_row_index_for_ordinal(row_ordinal);

                // this equals row `row_index` of M_{rk}^{-1} * R^{-1}_{\rho \rho}
                let lefthand_factor_vec     =   
                        self.umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row()                                                    
                            .row( &row_ordinal )
                            .scale_by( scalar_inv, self.umatch.ring_operator() )
                            .map(   |(x,y)|   // re-index the entries, so that indices align with the row indices of the factored matrix
                                    (self.umatch.matching.row_index_for_ordinal(x), y) 
                                );

                let iter = multiply_row_vector_with_matrix ( 
                            lefthand_factor_vec,
                            & self.umatch.matrix_to_factor,
                            self.umatch.ring_operator(),
                            self.umatch.order_operator_for_row_entries(),
                        );

                SourceCombInverseRow{
                                iter_unwrapped:     TwoTypeIterator::Version2(  iter  )
                            }                        
            }
        }
    }   

    fn row_reverse(             &self,  index: &Self::RowIndex    )   -> Self::RowReverse               { 
        let mut vec = self.row(&index).collect_vec();
        (&mut vec).reverse();
        vec.into_iter()         
    }
    
    fn column(                  &self,  index: &Self::ColumnIndex )   -> Self::Column                   {
        let mut vec = self.column_reverse(&index).collect_vec();
        (&mut vec).reverse();
        vec.into_iter()                 
    }

    fn column_reverse(          &self,  index: &Self::ColumnIndex )   -> Self::ColumnReverse
    {
        let ring_operator = self.umatch.ring_operator();
        let col     =   self.umatch.matrix_to_factor_matched_rows_only().column_reverse( index );
        let mut solution    =   
            multiply_column_vector_with_matrix_and_return_reversed(
                    col, 
                    self.umatch.matched_block_of_target_comb_inverse(), 
                       self.umatch.ring_operator(), 
                    self.umatch.order_operator_for_column_entries(),
                )
                .map(
                    |x|
                    {
                        let row_ordinal = self.umatch.matching.ordinal_for_row_index( &x.key() ).unwrap();
                        let matched_column_index = self.umatch.matching.column_index_for_ordinal( row_ordinal );
                        let matched_snzval = self.umatch.matching.coefficient_for_ordinal( row_ordinal );
                        let coefficient = ring_operator.divide( x.val(), matched_snzval );
                        MatrixToFactor::RowEntry::new( matched_column_index, coefficient )
                    }
                )
                .collect_vec();
        // if `index` is not a matched column index, then append an entry of form (index, 1)
        if ! self.umatch.matching.has_a_match_for_column_index( & index ) {
            solution.push( MatrixToFactor::RowEntry::new( index.clone(), MatrixToFactor::RingOperator::one() ) );
        }

        // sort the entries in the solution in descending order
        let order_operator_for_row_entries_reverse_extended = self.umatch.order_operator_for_row_entries_reverse(); //_extended_functionality();
        solution.sort_by( |a,b| order_operator_for_row_entries_reverse_extended.judge_partial_cmp(a, b).unwrap() ); // sort the entries of the vector
        
        solution.into_iter()

    }
} 








impl < 'a, MatrixToFactor >
    
    MatrixAlgebra for 
    
    SourceCombInverse
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,
{
    type RingOperator                           =   MatrixToFactor::RingOperator;

    type OrderOperatorForRowEntries             =   MatrixToFactor::OrderOperatorForRowEntries;

    type OrderOperatorForRowIndices             =   MatrixToFactor::OrderOperatorForColumnIndices;

    type OrderOperatorForColumnEntries          =   MatrixToFactor::OrderOperatorForRowEntries;

    type OrderOperatorForColumnIndices          =   MatrixToFactor::OrderOperatorForColumnIndices;

    /// The ring operator for the source COMB is the same as for the factored matrix
    fn ring_operator( &self ) -> Self::RingOperator {
        self.umatch.matrix_to_factor.ring_operator()
    }

    /// The order operators for row and column entries for the source COMB (and its inverse) are the same as the order operator for the row entries of the factored matrix.
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.umatch.matrix_to_factor.order_operator_for_row_entries()
    }

    /// The order operators for row and column indices for the source COMB (and its inverse) are the same as the order operator for the column indices of the factored matrix.    
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.umatch.matrix_to_factor.order_operator_for_column_indices()
    }

    /// The order operators for row and column entries for the source COMB (and its inverse) are the same as the order operator for the row entries of the factored matrix.    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.umatch.matrix_to_factor.order_operator_for_row_entries()
    }

    /// The order operators for row and column indices for the source COMB (and its inverse) are the same as the order operator for the column indices of the factored matrix.        
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.umatch.matrix_to_factor.order_operator_for_column_indices()
    }
}







impl < 'a, MatrixToFactor >
    
    MatrixOracleOperations for 
    
    SourceCombInverse
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:                     MatrixAlgebra,
        MatrixToFactor::ColumnIndex:        Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct        
{}        









//  DEFINE ROW
//  ---------------------------------------------------------------------------------------------------------

#[derive(Debug, Clone, )]
pub struct SourceCombInverseRow
                < MatrixToFactor >
    where   
        MatrixToFactor:                                 MatrixAlgebra,        
        MatrixToFactor::ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
        MatrixToFactor::RowIndex:                       Hash, // required for the hashing performed by the generalized matching array
        MatrixToFactor::RowEntry:                       KeyValPair, 
        MatrixToFactor::RingOperator:                   DivisionRingOperations,   
{
    iter_unwrapped:     TwoTypeIterator< 
                                Once
                                    < MatrixToFactor::RowEntry >,
                                LinearCombinationSimplified
                                    < MatrixToFactor::Row, MatrixToFactor::RingOperator, MatrixToFactor::OrderOperatorForRowEntries >,
                            >,
}

impl < MatrixToFactor >

    Iterator for

    SourceCombInverseRow
        < MatrixToFactor >
        
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,             
        
{
    type Item = MatrixToFactor::RowEntry;
    
    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}



























//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: TARGET
//  ---------------------------------------------------------------------------------------------------------







//  IMPLEMENT MATRIX ORACLE
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, MatrixToFactor >
    
    MatrixOracle for 
    
    TargetComb
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,

{   
    type Coefficient            =   MatrixToFactor::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex               =   MatrixToFactor::RowIndex;       // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   MatrixToFactor::RowIndex;       // The type of column indices

    type RowEntry               =   MatrixToFactor::ColumnEntry   ;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   MatrixToFactor::ColumnEntry   ;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                    =   TargetCombRow< MatrixToFactor::ColumnEntry >;  // What you get when you ask for a row.
    type RowReverse             =   IntoIter< MatrixToFactor::ColumnEntry >;        // What you get when you ask for a row with the order of entries reversed
    type Column                 =   IntoIter< MatrixToFactor::ColumnEntry >;            // What you get when you ask for a column   
    type ColumnReverse          =   TargetCombColumnReverse< MatrixToFactor >;     // What you get when you ask for a column with the order of entries reversed                                

    
    
    /// The target COMB (and its inverse) has a column for index `j` iff the factored matrix has a row for index `j`    
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.umatch.matrix_to_factor.has_row_for_index(index)
    }
    /// The target COMB (and its inverse) has a column for index `j` iff the factored matrix has a row for index `j`        
    fn has_column_for_index_result(  &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnIndex, Self::ColumnIndex > {
        self.umatch.matrix_to_factor.has_row_for_index_result(index)        
    }
    /// The target COMB (and its inverse) has a row for index `j` iff the factored matrix has a row for index `j`
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        self.umatch.matrix_to_factor.has_row_for_index(index)
    }
    /// The target COMB (and its inverse) has a row for index `j` iff the factored matrix has a row for index `j`    
    fn has_row_for_index_result(     &   self, index: & Self::RowIndex   )   -> Result< Self::RowIndex   , Self::RowIndex    > {
        self.umatch.matrix_to_factor.has_row_for_index_result(index)        
    }    
    
    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        {
            let row = self.row( row );
            let order_operator_for_column_indices = self.umatch.order_operator_for_row_indices();
            for entry in row {
                match order_operator_for_column_indices.judge_cmp( & entry.key(), column  ) {
                    Ordering::Equal => { return Some( entry.val() )  },
                    Ordering::Greater => { return None },
                    Ordering::Less => { continue }
                }
            }
            return None
        }
    
    fn row(                     & self,  index: & MatrixToFactor::RowIndex    )       -> Self::Row
    {

        // define the matrix A that will fit into an equation xA = b
        let seed_with_integer_indexed_rows = self.umatch.matched_blocks_of_target_comb_inverse_times_matrix_to_factor_oc();

        // define the problem vector b that will fit into an equation xA = b
        let matrix_to_factor_matched_columns_only = self.umatch.matrix_to_factor_matched_columns_only();
        let problem_vector = matrix_to_factor_matched_columns_only.row( index );

        // Struct that encodes the matching from column indices of A to row index ordinals of A
        let generalized_matching_matrix_ref = self.umatch.generalized_matching_matrix_ref();
        let column_index_to_row_ordinal = | column_index: MatrixToFactor::ColumnIndex | -> Option< usize > { generalized_matching_matrix_ref.ordinal_for_column_index( &column_index ) };
        let column_index_to_row_ordinal_wrapped = EvaluateFunctionFnMutWrapper::new( column_index_to_row_ordinal );        

        // Solve xA = b for x
        let mut solution_vec = 
        RowEchelonSolverReindexed::solve( // solution to xA = b
                problem_vector, // mabrix b
                seed_with_integer_indexed_rows, // matrix A
                column_index_to_row_ordinal_wrapped,
                self.umatch.ring_operator(),
                self.umatch.order_operator_for_row_entries(),
            )
            .solution().unwrap();

        // Sort solution vector x
        // println!("!!! REDO THIS SECTION OF CODE AFTER REFACTORING <CompareOrder>");
        let order_operator_for_row_entries_clone = self.umatch.order_operator_for_row_entries(); // we have to clone the order comparator in order to compare order, since the `lt` method takes `& mut self` as an argument -- but to get a mutable reference in the present context we would have to borrow the umatch struct as mutable
        let compare_order_for_vector_sort 
            = | x: &MatrixToFactor::RowEntry, y: &MatrixToFactor::RowEntry |-> Ordering { // we have "repackage" order_operator_for_column_entries so that it returns an std::cmp::Ordering rather than a bool
            match order_operator_for_row_entries_clone.judge_lt( x, y ) {
                true => { Ordering::Less },
                false => { 
                    match order_operator_for_row_entries_clone.judge_lt( y, x ) {
                        true => { Ordering::Greater },
                        false => { Ordering:: Equal }
                    }
                }
            }
        };
        solution_vec.sort_by( compare_order_for_vector_sort );

        // Reindex solution vector x
        let mut solution_vec_reindexed = 
            solution_vec.into_iter()
                .map(   | item |
                        MatrixToFactor::ColumnEntry::new( 
                            self.umatch.matching.row_index_for_column_index( &item.key() ).unwrap(),
                            item.val(),
                        )
                    )
                .collect_vec();

        // a function that sends (a,b) to `Less` if a < b and to `Greater` otherwise
        let order_operator_for_column_entries = self.umatch.order_operator_for_column_entries(); // we have to make a clone because, as currently written, `self` is behind a mutable reference and `order_operator_for_column_entries` may mutate itself when it compares two objects
        solution_vec_reindexed.sort_by( 
                |a,b| {  
                    if order_operator_for_column_entries.judge_lt( a, b ) { Ordering::Less }
                    else { Ordering::Greater }
                    }
            ); 
        let solution_vec_reindexed_and_sorted = IterWrappedVec::new( solution_vec_reindexed );                            

        // Possibly append an entry with coefficient 1
        match self.umatch.generalized_matching_matrix_ref().has_a_match_for_row_index( & index ){
            false =>    { // if row_index is UNMATCHED then merge in a copy of ( row_index, 1 )
                TargetCombRow{
                                iter_unwrapped:     TwoTypeIterator::Version1(
                                                            std::iter::once( MatrixToFactor::ColumnEntry::new(index.clone(), MatrixToFactor::RingOperator::one() ) )
                                                                .chain( solution_vec_reindexed_and_sorted ) 
                                                        )
                            }

            }
            true =>     { // otherwise just return the reindexed, sorted solution vector
                TargetCombRow{
                                iter_unwrapped:     TwoTypeIterator::Version2(
                                                            solution_vec_reindexed_and_sorted
                                                        )
                            }
            }
        }
    }   
    fn row_reverse(             &self,  index: &Self::RowIndex    )   -> Self::RowReverse               {
        let mut vec = self.row(index).collect_vec();
        (&mut vec).reverse();
        vec.into_iter() 
    }    
    fn column(                  &self,  index: &Self::ColumnIndex )   -> Self::Column                   {
        let mut vec = self.column_reverse(index).collect_vec();
        (&mut vec).reverse();
        vec.into_iter() 
    }
    fn column_reverse(          &self,  index: &Self::ColumnIndex )   -> Self::ColumnReverse
    {
        match self.umatch.generalized_matching_matrix_ref().column_index_for_row_index( & index ) {

            // follow this branch if row_index is matched to a column index, column_index
            Some( column_index ) => {

                // 
                let unit_vector_column_index = std::iter::once( MatrixToFactor::RowEntry::new( column_index, MatrixToFactor::RingOperator::one() ) );   // the standard unit vector supported on `row_index`                                

                // the matched block of R^{-1} D, indexed by column indices along rows and columns
                let seed = self.umatch.target_comb_inverse_times_matrix_to_factor_matched_block_with_rows_indexed_by_matched_column_index();

                // DEPRECATED (OK TO DELETE)
                // the matching relation from rows of A to columns of A
                // let column_index_for_row_index = self.umatch.generalized_matching_matrix_ref().bijection_row_index_to_column_index();

            
                // a column c of A^{-1}
                let column_of_a_inv = TriangularSolveForColumnVectorReverse::solve(
                        unit_vector_column_index,
                        seed, // matrix A
                    ).ok().unwrap();

                // the product vector D_{m k} * c  (refer to the matching identities in the Umatch facotrization paper for explanation)
                let column_of_comb  =   multiply_column_vector_with_matrix_and_return_reversed(
                                    column_of_a_inv,
                                                self.umatch.matrix_to_factor_ref(),
                                                self.umatch.ring_operator(),
                                                self.umatch.order_operator_for_column_entries(), // the constructor handles reversing the order for us, so we don't have to
                                            );                                          
                
                // wrap the product vector in an enum
                TargetCombColumnReverse{
                                iter_unwrapped:     TwoTypeIterator::Version1( column_of_comb )
                            }

            }

            // follow this branch if row_index is not matched to a column index
            None => {

                let unit_vector_row_index = std::iter::once( MatrixToFactor::ColumnEntry::new( index.clone(), MatrixToFactor::RingOperator::one() ) );   // the standard unit vector supported on `row_index`                
                TargetCombColumnReverse{
                                iter_unwrapped:     TwoTypeIterator::Version2( unit_vector_row_index )
                            }                                
            }
        }
    }
} 










impl < 'a, MatrixToFactor >
    
    MatrixAlgebra for 
    
    TargetComb
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,

{
    type RingOperator                           =   MatrixToFactor::RingOperator;

    type OrderOperatorForRowEntries             =   MatrixToFactor::OrderOperatorForColumnEntries;

    type OrderOperatorForRowIndices             =   MatrixToFactor::OrderOperatorForRowIndices;

    type OrderOperatorForColumnEntries          =   MatrixToFactor::OrderOperatorForColumnEntries;

    type OrderOperatorForColumnIndices          =   MatrixToFactor::OrderOperatorForRowIndices;

    /// The ring operator for the target COMB is the same as for the factored matrix
    fn ring_operator( &self ) -> Self::RingOperator {
        self.umatch.matrix_to_factor.ring_operator()
    }

    /// The order operators for row and column entries for the target COMB are the same as the order operator for the column entries of the factored matrix.
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.umatch.matrix_to_factor.order_operator_for_column_entries()
    }

    /// The order operators for row and column indices for the target COMB are the same as the order operator for the row indices of the factored matrix.    
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.umatch.matrix_to_factor.order_operator_for_row_indices()
    }

    /// The order operators for row and column entries for the target COMB are the same as the order operator for the column entries of the factored matrix.    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.umatch.matrix_to_factor.order_operator_for_column_entries()
    }

    /// The order operators for row and column indices for the target COMB are the same as the order operator for the row indices of the factored matrix.        
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.umatch.matrix_to_factor.order_operator_for_row_indices()
    }
}      






impl < 'a, MatrixToFactor >
    
    MatrixOracleOperations for 
    
    TargetComb
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:                     MatrixAlgebra,
        MatrixToFactor::ColumnIndex:        Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct        
{}        









//  DEFINE ROW
//  ---------------------------------------------------------------------------------------------------------




pub struct TargetCombRow< T: Clone >
{
    iter_unwrapped:     TwoTypeIterator< 
                                Chain<
                                        Once
                                            < T >,
                                        IterWrappedVec
                                            < T >,
                                    >,
                                IterWrappedVec
                                    < T >,
                            >
}

impl < T: Clone >

    Iterator for

    TargetCombRow< T >
{
    type Item = T;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}




//  DEFINE COLUMN
//  ---------------------------------------------------------------------------------------------------------



pub struct TargetCombColumnReverse
                < MatrixToFactor >

    where   
        MatrixToFactor:                                 MatrixAlgebra,        
        MatrixToFactor::ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
        MatrixToFactor::RowIndex:                       Hash, // required for the hashing performed by the generalized matching array

        MatrixToFactor::ColumnEntry:                    KeyValPair, 
        MatrixToFactor::RingOperator:                   DivisionRingOperations,                         
{
    iter_unwrapped:     TwoTypeIterator<
                                LinearCombinationOfColumnsReverse< MatrixToFactor >,
                                std::iter::Once< MatrixToFactor::ColumnEntry >,
                            >
}             

impl    < MatrixToFactor >

        Iterator for 

        TargetCombColumnReverse
                < MatrixToFactor >

    where   
        MatrixToFactor:                                 MatrixAlgebra,        
        MatrixToFactor::ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
        MatrixToFactor::RowIndex:                       Hash, // required for the hashing performed by the generalized matching array
        MatrixToFactor::ColumnEntry:                    KeyValPair, 
        MatrixToFactor::RingOperator:                   DivisionRingOperations,                      
{
    type Item = MatrixToFactor::ColumnEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}
































//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: TARGET INV
//  ---------------------------------------------------------------------------------------------------------




//  IMPLEMENT MATRIX ORACLE
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, MatrixToFactor >
    
    MatrixOracle for 
    
    TargetCombInverse
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >, 

{   
    type Coefficient            =   MatrixToFactor::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex               =   MatrixToFactor::RowIndex;       // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   MatrixToFactor::RowIndex;       // The type of column indices

    type RowEntry               =   MatrixToFactor::ColumnEntry   ;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   MatrixToFactor::ColumnEntry   ;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                    =   TargetCombInverseRow< 'a, MatrixToFactor >;  // What you get when you ask for a row.
    type RowReverse             =   IntoIter< MatrixToFactor::ColumnEntry >;        // What you get when you ask for a row with the order of entries reversed
    type Column                 =   IntoIter< MatrixToFactor::ColumnEntry >;            // What you get when you ask for a column   
    type ColumnReverse          =   TargetCombInverseColumnReverse< 'a, MatrixToFactor >;     // What you get when you ask for a column with the order of entries reversed                                


    /// The target COMB (and its inverse) has a column for index `j` iff the factored matrix has a row for index `j`    
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.umatch.matrix_to_factor.has_row_for_index(index)
    }
    /// The target COMB (and its inverse) has a column for index `j` iff the factored matrix has a row for index `j`        
    fn has_column_for_index_result(  &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnIndex, Self::ColumnIndex > {
        self.umatch.matrix_to_factor.has_row_for_index_result(index)        
    }
    /// The target COMB (and its inverse) has a row for index `j` iff the factored matrix has a row for index `j`
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        self.umatch.matrix_to_factor.has_row_for_index(index)
    }
    /// The target COMB (and its inverse) has a row for index `j` iff the factored matrix has a row for index `j`    
    fn has_row_for_index_result(     &   self, index: & Self::RowIndex   )   -> Result< Self::RowIndex   , Self::RowIndex    > {
        self.umatch.matrix_to_factor.has_row_for_index_result(index)        
    }   

    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        {
            let row = self.row( row );
            let order_operator_for_column_indices = self.umatch.order_operator_for_row_indices();
            for entry in row {
                match order_operator_for_column_indices.judge_cmp( & entry.key(), column  ) {
                    Ordering::Equal => { return Some( entry.val() )  },
                    Ordering::Greater => { return None },
                    Ordering::Less => { continue }
                }
            }
            return None
        }
    
    fn row(                     & self,  index: & MatrixToFactor::RowIndex    )       -> Self::Row
    {
        // a struct that realizes the function sending (row_ordinal: usize, coefficient) --> t: MatrixToFactor::ColumnEntry,
        // where t has 
        //   - key equal to the row_ordinal'th row index in the matching matrix
        //   - val equal to coefficient
        let entry_changer = RemapEntryTuple::new(
                                    self.umatch.generalized_matching_matrix_ref().bijection_row_indices_to_ordinals_and_inverse().vec_elements_in_order(),
                                    IdentityFunction,
                                );

        match self.umatch.generalized_matching_matrix_ref().ordinal_for_row_index( & index ) {
            // unmatched row index
            None => {
                // define the matrix A that will fit into an equation xA = b
                let seed_with_integer_indexed_rows = self.umatch.matched_blocks_of_target_comb_inverse_times_matrix_to_factor_oc();
                // let seed_with_integer_indexed_rows_ref = & seed_with_integer_indexed_rows;

                // define the problem vector b that will fit into an equation xA = b
                let matrix_to_factor_matched_columns_only = self.umatch.matrix_to_factor_matched_columns_only();
                let problem_vector 
                    = matrix_to_factor_matched_columns_only
                        .row( index )
                        .negate( self.umatch.ring_operator() );  // recall that there is a MINUS SIGN in the formula in Theorem 6 of Hang et al., "Umatch factorization: ..."

                // Struct that encodes the matching from column indices of A to row index ordinals of A
                let generalized_matching_matrix_ref = self.umatch.generalized_matching_matrix_ref();
                let column_index_to_row_ordinal = | column_index: MatrixToFactor::ColumnIndex | -> Option< usize > { generalized_matching_matrix_ref.ordinal_for_column_index( &column_index ) };
                let column_index_to_row_ordinal_wrapped = EvaluateFunctionFnMutWrapper::new( column_index_to_row_ordinal );        

                // Solve xA = b for x
                let solution_vec = 
                RowEchelonSolverReindexed::solve( // solution to xA = b
                        problem_vector, // mabrix b
                        seed_with_integer_indexed_rows, // matrix A
                        column_index_to_row_ordinal_wrapped,
                        self.umatch.ring_operator(),
                        self.umatch.order_operator_for_row_entries(),
                    )
                    .quotient()
                    .collect::<Vec<_>>()   // NB: we dump the entries returned by this iterator into a Rust Vec becuase the ChangeIndexSimple::new(..) function keeps asking for an exact type specification for `solution_vec` ... which we can't otherwise give, because the solver contains a closure (which has no specifiable type)
                    .into_iter();
                let solution_vec_integer_indexed:   ChangeIndexSimple<
                                                        IntoIter<_>,                                              // reindexed vector
                                                        &HashMap<MatrixToFactor::ColumnIndex, usize>,   // re-indexing function
                                                        usize,                                          // new index type
                                                    >
                    = ChangeIndexSimple::new( solution_vec, self.umatch.generalized_matching_matrix_ref().bijection_column_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal() );
                    
                // Multiply the solution x by matrix R_{\rho \rho}^{-1}
                let comb_target_inv_matched_block = self.umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row();

                // Compute the portion of desired matrix row indexed by matched row indices -- however, this vector will actually have usize indices which we must change next
                let vec_with_matched_indices 
                    = multiply_row_vector_with_matrix( 
                            solution_vec_integer_indexed, 
                            comb_target_inv_matched_block, 
                            self.umatch.ring_operator(),
                            OrderOperatorByKey::new(),
                        );

                // reindex the preceding vector so that it is indeed indexed by row indices
                let vec_with_matched_indices_reindexed 
                    = MapByTransform::new( 
                            vec_with_matched_indices,
                            entry_changer.clone(),             
                        );
                
                // Compute the portion of desired matrix row indexed by unmatched row indices -- this portion of the vector has exactly one nonzero entry, and it is the first entry in the vector, because it's the diagonal element of a row of a triangular matrix
                let vec_with_unmatched_index = OncePeekable::new(  MatrixToFactor::ColumnEntry::new(index.clone(), MatrixToFactor::RingOperator::one() )  ); 

                // Merge the two portions of the vector together, to form a whole
                let merged  
                    =   vec_with_unmatched_index.chain( vec_with_matched_indices_reindexed );                  
                
                return  TargetCombInverseRow::new( TwoTypeIterator::Version1( merged ) ) //TargetCombInverseRow{ iter_unwrapped: ChangeEntryType::new( TwoTypeIterator::Version1( merged ) ) }
            }
            // matched row index
            Some( row_ordinal ) => {
                let reindexed_iter = 
                    MapByTransform::new( 
                            self.umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row().row( & row_ordinal ),
                            entry_changer,
                        );
                return  TargetCombInverseRow::new( TwoTypeIterator::Version2( reindexed_iter ) ) // TargetCombInverseRow{ iter_unwrapped: ChangeEntryType::new( TwoTypeIterator::Version2( reindexed_iter ) ) }
            }
        }
    }         
    fn row_reverse(             &self,  index: &Self::RowIndex    )   -> Self::RowReverse               { 
        let mut vec = self.row(index).collect_vec();
        (&mut vec).reverse();
        vec.into_iter() 
    }
    
    fn column(                  &self,  index: &Self::ColumnIndex )   -> Self::Column                   {
        let mut vec = self.column_reverse(index).collect_vec();
        (&mut vec).reverse();
        vec.into_iter()     
    }
    fn column_reverse(          &self,  index: &Self::ColumnIndex )   -> Self::ColumnReverse
    {
        match self.umatch.matching.has_a_match_for_row_index( &index ) {
            true => {
                // get the ordinal of the matched row index
                let row_ordinal = self.umatch.matching.ordinal_for_row_index( &index ).unwrap();
                // a column, c, of (R_{\rho \rho}^{-1})  [THIS NOTATION COMES FROM INNER IDENTITIES, C.F. UMATCH FACTORIZATION]
                // we have to reindex this column with column indices, then sort it according to the order on column indices, because our next step will be to run a triangular solve, and the triangular matrix is only triangular with respect to the order on column indices                                
                let mut comb_target_inv_matched_block_col     = 
                    self.umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row()
                        .column_reverse( & row_ordinal )
                        .map(   |(key,val)|   // here we convert indices to column indices
                                MatrixToFactor::RowEntry::new(
                                        self.umatch.matching.column_index_for_ordinal(key), 
                                        val
                                    )  
                            )
                        .collect_vec(); // collect entries into a vector to simplify sorting
                // sort the entries of the column
                // let cmp_style_comparator = InferTotalOrderFromJudgePartialOrder::new( // we have to construct an order comparator that spits out values of type Ord rather than type bool
                //         self.umatch.order_operator_for_row_entries() 
                //     );    
                comb_target_inv_matched_block_col.sort_by( |x,y| self.umatch.order_operator_for_row_entries_reverse().judge_partial_cmp( x, y ).unwrap() );
                let A = self.umatch.target_comb_inverse_times_matrix_to_factor_matched_block_with_rows_indexed_by_matched_column_index();
                // A^{-1} * c
                let echelon_solution    =   TriangularSolveForColumnVectorReverse::solve(
                                                    comb_target_inv_matched_block_col.into_iter(),
                                                    A,
                                                ).ok().unwrap();
                let echelon_solution_minus = echelon_solution.negate( self.umatch.ring_operator() );
                let unmatched_part_of_solution = multiply_column_vector_with_matrix_and_return_reversed(
                        echelon_solution_minus,
                        self.umatch.matrix_to_factor_matchless_rows_only(),
                        self.umatch.ring_operator(),
                        self.umatch.order_operator_for_column_entries(),
                    );

                // this equals the variable comb_target_inv_matched_block_col defined above; we need a second copy, and calling the constructor a second time avoids the need to impose "Clone" on the struct
                let matched_part_of_solution     = 
                    self.umatch.matched_block_of_target_comb_inverse()
                        .column_reverse( index );
                      
                // let   matched_part_of_solution = TwoTypeIterator::Version1( matched_part_of_solution );
                // let unmatched_part_of_solution = TwoTypeIterator::Version2( unmatched_part_of_solution );

                let solution = MergeTwoIteratorsByOrderOperator::new(
                            matched_part_of_solution.peekable(),
                            unmatched_part_of_solution.peekable(),
                            self.umatch.order_operator_for_column_entries_reverse(),
                );
                let solution = ChangeEntryType::new( TwoTypeIterator::Version1( solution ) );

                // wrap the output in a struct that has a simpler type signature
                TargetCombInverseColumnReverse{ iter_unwrapped: solution }
                
            }
            false => {
                let solution = TwoTypeIterator::Version2( 
                        std::iter::once( 
                                Self::ColumnEntry::new( index.clone(), MatrixToFactor::RingOperator::one()  ) 
                            ),
                    );
                let solution = ChangeEntryType::new( solution );
                // wrap the output in a struct that has a simpler type signature                    
                TargetCombInverseColumnReverse{ iter_unwrapped: solution }
            }
        }
    }
} 






impl < 'a, MatrixToFactor >
    
    MatrixAlgebra for 
    
    TargetCombInverse
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,
{
    type RingOperator                           =   MatrixToFactor::RingOperator;

    type OrderOperatorForRowEntries             =   MatrixToFactor::OrderOperatorForColumnEntries;

    type OrderOperatorForRowIndices             =   MatrixToFactor::OrderOperatorForRowIndices;

    type OrderOperatorForColumnEntries          =   MatrixToFactor::OrderOperatorForColumnEntries;

    type OrderOperatorForColumnIndices          =   MatrixToFactor::OrderOperatorForRowIndices;

    /// The ring operator for the target COMB is the same as for the factored matrix
    fn ring_operator( &self ) -> Self::RingOperator {
        self.umatch.matrix_to_factor.ring_operator()
    }

    /// The order operators for row and column entries for the target COMB are the same as the order operator for the column entries of the factored matrix.
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.umatch.matrix_to_factor.order_operator_for_column_entries()
    }

    /// The order operators for row and column indices for the target COMB are the same as the order operator for the row indices of the factored matrix.    
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.umatch.matrix_to_factor.order_operator_for_row_indices()
    }

    /// The order operators for row and column entries for the target COMB are the same as the order operator for the column entries of the factored matrix.    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.umatch.matrix_to_factor.order_operator_for_column_entries()
    }

    /// The order operators for row and column indices for the target COMB are the same as the order operator for the row indices of the factored matrix.        
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.umatch.matrix_to_factor.order_operator_for_row_indices()
    }
}      







impl < 'a, MatrixToFactor >
    
    MatrixOracleOperations for 
    
    TargetCombInverse
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:                     MatrixAlgebra,
        MatrixToFactor::ColumnIndex:        Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct        
{}        








//  DEFINE ROW
//  ---------------------------------------------------------------------------------------------------------




#[derive(Clone, Dissolve, new)]
pub struct TargetCombInverseRow
            < 'a, MatrixToFactor,  >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,

{
    iter_unwrapped:
    TwoTypeIterator<
            std::iter::Chain<
                    OncePeekable< MatrixToFactor::ColumnEntry >, 
                    MapByTransform<
                            Simplify<
                                    IteratorsMergedInSortedOrder<
                                            Scale<
                                                    std::iter::Chain<
                                                            std::iter::Once<(usize, MatrixToFactor::Coefficient)>, 
                                                            Cloned<std::slice::Iter<'a, (usize, MatrixToFactor::Coefficient)>>
                                                        >, 
                                                    MatrixToFactor::RingOperator, 
                                                >, 
                                                OrderOperatorByKey,
                                        >, 
                                    MatrixToFactor::RingOperator, 
                                >, 
                            MatrixToFactor::ColumnEntry, 
                            RemapEntryTuple<
                                    usize, 
                                    MatrixToFactor::RowIndex, 
                                    &'a Vec<MatrixToFactor::RowIndex>, 
                                    MatrixToFactor::Coefficient, 
                                    MatrixToFactor::Coefficient, 
                                    IdentityFunction, 
                                    MatrixToFactor::ColumnEntry,
                                >
                        >
                >, 
            MapByTransform<
                    Chain<
                            Once<(usize, MatrixToFactor::Coefficient)>, 
                            Cloned<Iter<'a, (usize, MatrixToFactor::Coefficient)>>
                        >,
                    MatrixToFactor::ColumnEntry,
                    RemapEntryTuple<
                            usize, 
                            MatrixToFactor::RowIndex, 
                            &'a Vec<MatrixToFactor::RowIndex>, 
                            MatrixToFactor::Coefficient, 
                            MatrixToFactor::Coefficient, 
                            IdentityFunction, 
                            MatrixToFactor::ColumnEntry,
                        >                      
                >          
    >
}    



impl < 'a, MatrixToFactor,  >

    Iterator for

    TargetCombInverseRow
            < 'a, MatrixToFactor,  >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,
{
    type Item = MatrixToFactor::ColumnEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}       





//  DEFINE COLUMN
//  ---------------------------------------------------------------------------------------------------------


pub struct TargetCombInverseColumnReverse
                < 'a, MatrixToFactor >
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,             
{
    iter_unwrapped: 
                    // <   SparseVector, EntryNew, Index, RingElement, >
                    ChangeEntryType<
                            TwoTypeIterator<
                                    MergeTwoIteratorsByOrderOperator<
                                            Peekable<
                                                    MapByTransform<
                                                            Chain<
                                                                    Once<(usize, MatrixToFactor::Coefficient)>, 
                                                                    // VecOfVecMatrixColumnReverse<'a, usize, MatrixToFactor::Coefficient>
                                                                    Cloned<Rev<std::slice::Iter<'a, (usize, MatrixToFactor::Coefficient)>>>,
                                                                >, 
                                                            MatrixToFactor::ColumnEntry, 
                                                            ReindexEntry< &'a Vec<MatrixToFactor::RowIndex> >
                                                        >,                                                                
                                                >,
                                            Peekable< 
                                                LinearCombinationSimplified<
                                                        OnlyIndicesOutsideCollection<
                                                                MatrixToFactor::ColumnReverse, 
                                                                &'a HashMap< MatrixToFactor::RowIndex, usize>, 
                                                            >,
                                                            MatrixToFactor::RingOperator, 
                                                            ReverseOrder<MatrixToFactor::OrderOperatorForColumnEntries>
                                                    >,                                                                                                                             
                                                >,                                                        
                                            ReverseOrder< 
                                                    MatrixToFactor::OrderOperatorForColumnEntries
                                                >,                                                                
                                        >,
                                    std::iter::Once< MatrixToFactor::ColumnEntry >,
                                >,
                            MatrixToFactor::ColumnEntry,
                        >,
}

impl < 'a, MatrixToFactor >

    TargetCombInverseColumnReverse
                < 'a, MatrixToFactor >
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,            
{
    pub fn unwrap_iter( self ) -> 
                    ChangeEntryType<
                            TwoTypeIterator<
                                    MergeTwoIteratorsByOrderOperator<
                                            Peekable<
                                                    MapByTransform<
                                                            Chain<
                                                                    Once<(usize, MatrixToFactor::Coefficient)>, 
                                                                    // VecOfVecMatrixColumnReverse<'a, usize, MatrixToFactor::Coefficient>
                                                                    Cloned<Rev<std::slice::Iter<'a, (usize, MatrixToFactor::Coefficient)>>>,
                                                                >, 
                                                            MatrixToFactor::ColumnEntry, 
                                                            ReindexEntry< &'a Vec<MatrixToFactor::RowIndex> >
                                                        >,                                                                
                                                >,
                                            Peekable< 
                                                LinearCombinationSimplified<
                                                        OnlyIndicesOutsideCollection<
                                                                MatrixToFactor::ColumnReverse, 
                                                                &'a HashMap< MatrixToFactor::RowIndex, usize>, 
                                                            >, 
                                                            MatrixToFactor::RingOperator, 
                                                            ReverseOrder<MatrixToFactor::OrderOperatorForColumnEntries>
                                                    >,                                                                                                                             
                                                >,                                                        
                                            ReverseOrder< 
                                                    MatrixToFactor::OrderOperatorForColumnEntries
                                                >,                                                                
                                        >,
                                    std::iter::Once< MatrixToFactor::ColumnEntry >,
                                >,
                            MatrixToFactor::ColumnEntry,
                        >
    {
        self.iter_unwrapped
    }    
}



impl < 'a, MatrixToFactor >

    Iterator for

    TargetCombInverseColumnReverse
                < 'a, MatrixToFactor >
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,            
{
    type Item = MatrixToFactor::ColumnEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}

































//  =========================================================================================================
//  MATCHED BLOCK OF (TARGET COMB INV * MAPPING ARRAY)
//  =========================================================================================================





/// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// 
/// Concretely, this matrix equals (the pivot block of the inverse of the target COMB) * (the pivot block of the matching array).
/// It also equals the matched block of (SM^-), where S is the source COMB and M^- is the generalized inverse of the generalized
/// matching array M obtained by inverting nonzero entries and transposing.
#[derive(Copy, Clone, Debug, PartialEq, Dissolve)]
pub struct TargetCombInverseTimesMatrixToFactorMatchedBlock< 'a, MatrixToFactor, > 
    where   // these are the bare minimum requirements for a Umatch object
        MatrixToFactor:                     MatrixAlgebra,
        MatrixToFactor::ColumnIndex:        Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
{
    pub umatch:     & 'a Umatch< MatrixToFactor > ,  
}

impl < 'a, MatrixToFactor, > 
    
    TargetCombInverseTimesMatrixToFactorMatchedBlock
        < 'a, MatrixToFactor, > 
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,       
{
    // Make a new [`SourceCombInverseMatchedBlock`].
    pub fn new( umatch: &'a Umatch< MatrixToFactor >  ) -> Self {
        TargetCombInverseTimesMatrixToFactorMatchedBlock{ umatch }
    }
}



// Implement `Eq` if matrix coefficients implement `Eq`.
// (in fact, Eq can be implemented for this struct if and only if matrix coefficients implement Eq)
impl < 'a, MatrixToFactor, > 

    Eq for
    
    TargetCombInverseTimesMatrixToFactorMatchedBlock
        < 'a, MatrixToFactor, > 
    where   
        MatrixToFactor:                 MatrixAlgebra<
                                            ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                            RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                            RingOperator:           DivisionRingOperations,
                                            RowEntry:               KeyValPair,
                                            ColumnEntry:            KeyValPair,        
                                        >,
        MatrixToFactor:                 PartialEq,
        MatrixToFactor::Coefficient:    Eq, 
{}






//  IMPLEMENT MATRIX ORACLE
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, MatrixToFactor >
    
    MatrixOracle for 
    
    TargetCombInverseTimesMatrixToFactorMatchedBlock
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,

{   
    type Coefficient            =   MatrixToFactor::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex               =   MatrixToFactor::RowIndex;       // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   MatrixToFactor::ColumnIndex;       // The type of column indices

    type RowEntry               =   MatrixToFactor::RowEntry   ;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   MatrixToFactor::ColumnEntry   ;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                    =   LinearCombinationSimplified< 
                                        OnlyIndicesInsideCollection< MatrixToFactor::Row, &'a HashMap<MatrixToFactor::ColumnIndex, usize> >, 
                                        MatrixToFactor::RingOperator, MatrixToFactor::OrderOperatorForRowEntries 
                                    >;  // What you get when you ask for a row.
    type RowReverse             =   IntoIter< MatrixToFactor::RowEntry >;        // What you get when you ask for a row with the order of entries reversed
    type Column                 =   IntoIter< MatrixToFactor::ColumnEntry >;            // What you get when you ask for a column   
    type ColumnReverse          =   SourceCombInverseMatchedBlockColumnReverse< 'a, MatrixToFactor >;     // What you get when you ask for a column with the order of entries reversed                                

    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        {
            let row = self.row( row );
            let order_operator_for_column_indices = self.umatch.order_operator_for_column_indices();
            for entry in row {
                match order_operator_for_column_indices.judge_cmp( & entry.key(), column  ) {
                    Ordering::Equal => { return Some( entry.val() )  },
                    Ordering::Greater => { return None },
                    Ordering::Less => { continue }
                }
            }
            return None
        }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool 
        { self.umatch.generalized_matching_matrix_ref().has_a_match_for_column_index( index ) }
    fn has_row_for_index(  &   self, index: & Self::RowIndex)   -> bool 
        { self.umatch.generalized_matching_matrix_ref().has_a_match_for_row_index( index ) }        
    
    fn row(                     & self,  index: & MatrixToFactor::RowIndex    )       -> Self::Row
    {

        // get a vector representing the row of T^{-1} indexed by `index`
        // --------------------------------------------------------------

        // first we have to look up the right index, in order to find the proper row of T^{-1}, where T is the target COMB
        // NB: we have to translate the index `row_index` which has type `MatrixToFactor::RowIndex` into the index `row_ordinal` which has type `usize`, because `self.umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal` is indexed by unsigned integers 
        let row_ordinal  =   self.umatch.matching.ordinal_for_row_index( &index ).unwrap();

        // here we get the row; denote this vector `r`
        let combining_coefficients 
            = self.umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row().row( & row_ordinal );     

        // get a submatrix of the factored matrix that automatically filters out / removes entries that appear in non-pivot columns; denote this matrix A
        // ------------------------------------------------------------------------------------------------------------------------------------------

        // the matched columns of the factored matrix
        let matched_cols_of_matrix_to_factor : OnlyColumnIndicesInsideCollection< &'a MatrixToFactor, &'a HashMap<MatrixToFactor::ColumnIndex, usize>, >
            =   self.umatch.matrix_to_factor_matched_columns_only();
        // let matched_cols_of_matrix_to_factor // : OnlyColumnIndicesInsideCollection<'a, 'b, MatrixToFactor, MatrixToFactor::Row, HashMap<MatrixToFactor::ColumnIndex, usize>, MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient>
        //     =   self.umatch.matrix_to_factor_ref();            

        // compute the product rA
        // ----------------------

        // a collection of terms; their sum equals the product of the vector with the matrix
        let iter_over_scaled_rows = 
            combining_coefficients
                    // .iter()
                    .map(   
                            |( row_ordinal, snzval)|  
                            matched_cols_of_matrix_to_factor.row( 
                                    & self.umatch.matching.row_index_for_ordinal( row_ordinal ) 
                                )
                                .scale_by( snzval, self.umatch.ring_operator() )                                
                        );                 
                                                         

        // sum the terms
        hit_merge_by_predicate( iter_over_scaled_rows, self.umatch.order_operator_for_row_entries() )
                    .simplify( self.umatch.ring_operator() )
    }
    fn row_reverse(             &self,  index: &Self::RowIndex    )   -> Self::RowReverse               { 
        let mut vec = self.row(index).collect_vec();
        vec.reverse();
        vec.into_iter() 
    }    
    fn column(                  &self,  index: &Self::ColumnIndex )   -> Self::Column                   {
        let mut vec = self.column_reverse(index).collect_vec();
        vec.reverse();
        vec.into_iter() 
    }
    fn column_reverse(          &self,  index: &Self::ColumnIndex )   -> Self::ColumnReverse
    {

        // get a column of the factored matrix
        let matched_row_submatrix_of_matrix_to_factor = self.umatch.matrix_to_factor_matched_rows_only();
        let column      =   matched_row_submatrix_of_matrix_to_factor.column_reverse( index );

        // get an oracle for the inverse of the target COMB, with non-pivot rows removed
        let comb_target_inv_matched_block = self.umatch.matched_block_of_target_comb_inverse();
        
        // compute the product
        let column_new  = multiply_column_vector_with_matrix_and_return_reversed(
                    column,
                    comb_target_inv_matched_block,
                    self.umatch.ring_operator(),
                    self.umatch.order_operator_for_column_entries(),
                );

        // return the product in a wrapper
        SourceCombInverseMatchedBlockColumnReverse{ 
                column_reverse_of_target_comb_inverse_times_matrix_to_factor_matched_block: column_new 
            }
        
    }
} 









impl < 'a, MatrixToFactor >
    
    MatrixOracleOperations for 
    
    TargetCombInverseTimesMatrixToFactorMatchedBlock
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:                     MatrixAlgebra,
        MatrixToFactor::ColumnIndex:        Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct        
{}    






//  DEFINE COLUMN
//  --------------------------------------------------------------------------------------------------------------


/// A wrapper struct for descending columns of `SourceCombInverseMatchedBlock` 
pub struct  SourceCombInverseMatchedBlockColumnReverse
                < 'a, MatrixToFactor >
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,
    { 
        column_reverse_of_target_comb_inverse_times_matrix_to_factor_matched_block:
            LinearCombinationSimplified<
                    MapByTransform<
                            Chain<
                                    Once<(usize, MatrixToFactor::Coefficient)>, 
                                    Cloned<
                                            Rev<
                                                    std::slice::Iter<
                                                            'a, 
                                                            (usize, MatrixToFactor::Coefficient)
                                                        >
                                                >
                                        >,                                            
                                >, 
                            MatrixToFactor::ColumnEntry, 
                            ReindexEntry<
                                    &'a Vec<MatrixToFactor::RowIndex>
                                >
                        >, 
                    MatrixToFactor::RingOperator, 
                    ReverseOrder<MatrixToFactor::OrderOperatorForColumnEntries>,
                >,           
    }

impl < 'a, MatrixToFactor >

    Iterator for 

    SourceCombInverseMatchedBlockColumnReverse
        < 'a, MatrixToFactor >
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                    RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                    RingOperator:           DivisionRingOperations,
                                    RowEntry:               KeyValPair,
                                    ColumnEntry:            KeyValPair,        
                                >,
{
    type Item = MatrixToFactor::ColumnEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.column_reverse_of_target_comb_inverse_times_matrix_to_factor_matched_block.next()
    }
}

