

//! U-match factorization of differential matrices, with functionality for (co)homology, cycles, boundaries, etc.
//! 
//! A differential U-match decomposition is a matrix factorization technique specialized for
//! chain complexes. It can be used to compute 
//! 
//! - barcodes in persistent (co)homology
//! - homology groups
//! - bounding chains
//! - (co)cycle representatives
//! 
//! and more!  
//! 
//! # Definitions
//! 
//! A **differential U-match decomposition** is a special type of U-match decomposition which consists of an equation 
//! 
//! ```text
//! JM = DJ
//! ```
//! 
//! where 
//! 
//! - `D` is the **differential matrix** of a chain complex. For example, if the chain complex comes from a simplicial
//!   complex, then `D` has a row (respectively, column) for every simplex (including all simplices of all dimensions), and the column for simplex `s`
//!   encodes the boundary of `s`. In this sense, `D` essentially contains every boundary matrix of the simplicial complex
//!   as a submatrix.
//! 
//!   We also require the matrix `J` to be *homogeneous*, in the sense that each nonzero column represents a linear combination of
//!   simplices (or cubes, etc.) which have the same dimension. 
//! 
//! - `M` is a generalized matching matrix
//! - `J` is an upper triangular matrix with 1's on the diagonal.
//! 
//! We call `J` the **differential COMB**.
//! 
//! 
//! # Examples
//! 
//! For examples on how to use a differential U-match decomposition, see the modules for
//! 
//! - [Vietoris-Rips complexes](crate::topology::simplicial::from::graph_weighted::VietorisRipsComplex)
//! - [Dowker complexes](crate::topology::simplicial::from::relation)

use itertools::Itertools;
use ordered_float::OrderedFloat;

use crate::algebra::chain_complexes::barcode::{get_barcode, Barcode};
use crate::algebra::chain_complexes::{ChainComplex, FilteredChainComplex};
use crate::algebra::matrices::operations::iterate_rows_and_columns::SequenceOfRows;
use crate::algebra::matrices::operations::umatch::gimbled::GimbledUmatch;
use crate::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder;
use crate::algebra::matrices::types::transpose::OrderAntiTranspose;
use crate::utilities::iterators::general::{TwoTypeIterator, IterWrappedVec, IterWrappedVecReverse};
use crate::algebra::matrices::{operations::umatch::row_major::Umatch};
use crate::algebra::matrices::operations::umatch::row_major::comb::{SourceComb, SourceCombInverse, TargetComb, TargetCombInverse};
use crate::algebra::matrices::operations::{MatrixOracleOperations, SequenceOfColumns};
use crate::algebra::matrices::query::{MatrixAlgebra, MatrixOracle};
use crate::algebra::rings::traits::{SemiringOperations, DivisionRingOperations};
use crate::utilities::order::{JudgeOrder, JudgePartialOrder};
use crate::algebra::vectors::entries::{KeyValPair};

use std::fmt::Debug;
use std::cmp::Ordering;
use std::hash::Hash;
use std::collections::{HashMap, HashSet};
use std::iter::Cloned;
use std::vec::IntoIter;

use derive_new::new;
use derive_getters::Dissolve;


//  ---------------------------------------------------------------------------------------
//  COMB
//  ---------------------------------------------------------------------------------------



/// The COMB of a differential U-match decomposition
/// 
/// This is the matrix `J` in a differential U-match decomposition `JM = DJ`. For details
/// on this decomposition, see the documentation for [differential U-match module](crate::algebra::matrices::operations::umatch::differential).
/// 
/// # Construction
/// 
/// We obtain the differential COMB `J` from a secondary U-match decomposition `TM = DS`.
/// Concretely, the construction proceeds as follows
/// 
/// - First let `J = S`, then 
/// - For each nonzero entry `M[i,j]` in the generalized matching matrix `M`, replace column  `i` of `J` with column `i` of `T`
/// 
/// It can then be shown, mathematically, that `JM = DJ`, so `J` is a valid differential COMB.
/// 
/// # Examples
/// 
/// For examples on how to use a differential U-match decomposition, see the modules for
/// 
/// - [Vietoris-Rips complexes](crate::topology::simplicial::from::graph_weighted::VietorisRipsComplex)
/// - [Dowker complexes](crate::topology::simplicial::from::relation)
#[derive(new,Clone,Copy,Debug,PartialEq)]
pub struct DifferentialComb< 
                    'a, 
                    MatrixToFactor,                 
                >      
    where
        MatrixToFactor:                    MatrixAlgebra,
        MatrixToFactor::ColumnIndex:       Hash + Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:          Hash + Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
{
    umatch: &'a GimbledUmatch<  MatrixToFactor >,                 
}




//  ---------------------------------------------------------------------------------------
//  IMPLEMENT MATRIX ORACLE
//  ---------------------------------------------------------------------------------------


impl < 
        'a,
        BoundaryMatrix, 
        OrderOperatorForRowAndColumnEntries, 
        OrderOperatorForRowAndColumnIndices, 
        IndexForRowsAndColumns,
        EntryForRowsAndColumns,
    > 

    MatrixOracle for
    
    DifferentialComb
        < 'a, BoundaryMatrix, >   

    where
        BoundaryMatrix:                         MatrixAlgebra< 
                                                        RowEntry                        =   EntryForRowsAndColumns, 
                                                        RowIndex                        =   IndexForRowsAndColumns, 
                                                        OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                        ColumnEntry                     =   EntryForRowsAndColumns, 
                                                        ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                        OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                    >,
        IndexForRowsAndColumns:                     Clone + Debug + Eq + Hash, // required for the hashing performed by the generalized matching array   
        EntryForRowsAndColumns:                     Clone + Debug + PartialEq + KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,
        OrderOperatorForRowAndColumnEntries:        JudgePartialOrder< EntryForRowsAndColumns >,
        BoundaryMatrix::RingOperator:           DivisionRingOperations,   
{
    type Coefficient            =   BoundaryMatrix::Coefficient;    // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =   IndexForRowsAndColumns;    // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   IndexForRowsAndColumns;    // The type of column indices    
    type RowEntry               =   EntryForRowsAndColumns;  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   EntryForRowsAndColumns;  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`                 =   
    type Row                    =   IterWrappedVec< EntryForRowsAndColumns >;  // What you get when you ask for a row
    type RowReverse             =   IterWrappedVecReverse< EntryForRowsAndColumns >;  // What you get when you ask for a row with the order of entries reversed
    type Column                 =   IterWrappedVec< EntryForRowsAndColumns >;
    type ColumnReverse          =   IterWrappedVecReverse< EntryForRowsAndColumns >; 


    fn structural_nonzero_entry(                   &   self, row: &Self::RowIndex, column: &Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        
        // take the corresponding nonzero entry from the target comb if `column` is a matched row index, and from the source comb otherwise
        match self.umatch.generalized_matching_matrix_ref().has_a_match_for_row_index( column ) {
            false   =>  { self.umatch.source_comb().structural_nonzero_entry( row, column ) },            
            true    =>  { self.umatch.target_comb().structural_nonzero_entry( row, column ) },
        }
    }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.umatch.matrix_to_factor_ref().has_column_for_index( index ) 
    }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool {
        self.umatch.matrix_to_factor_ref().has_row_for_index( index )
    }
    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row                 
    { 
        // here we assume a proper Umatch decomposition T * M = D * S
        // we are tryig to build a row of the differntial COMB which is by definition
        // 
        // J = T * Delta_t + S * (I - Delta_t)
        //
        // see the section on Differential COMBs in the Umatch paper for details
        //
        // if `index` is not a matched column index then because the decomposition is
        // proper row `index` of `S` is a standard unit vector. when we swap the 
        // set of all non-unit-vector columns of `T` into `S`, we therefore transform
        // this unit vector into the row vector of `T` at index `index`.

        let matching_matrix     =   self.umatch.generalized_matching_matrix_ref();
        if matching_matrix.lacks_a_match_for_column_index( index ) { 
            return  IterWrappedVec::new( self.umatch.target_comb().row( index ).collect() );
        }

        // otherwise let's take the parts we need from each row of S and T and combine them
        
        // because the decomposition is proper all the off-diagonal nonzero elements of T occur in
        // columns whose indices are matched-row-indices. thus we don't throw away any elements from
        // the target row
        let target_row          =   self.umatch.target_comb().row( index );
        
        // we do have to throw out some elements of the source row
        let mut source_row      =   self.umatch.source_comb().row( index );
        source_row.next(); // throw out the leading entry (we get that for free from the target row)
        let source_row          =   source_row.filter( // throw out entries indexed by matched-row-indices
                                        |entry| 
                                        matching_matrix.lacks_a_match_for_row_index( & entry.key() ) 
                                    );

        // merge the two rows
        let order_operator      =   self.umatch.matrix_to_factor_ref().order_operator_for_row_entries();
        let merged              =   target_row.merge_by( 
                                        source_row,  
                                        | a,b | {  order_operator.judge_lt( a, b )  }
                                    )
                                    .collect::<Vec<_>>();
        return  IterWrappedVec::new( merged );
        
    }
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse          {
        self.row( index ).reset_and_run_backward()
    }     
    fn column(                  &   self, index: &Self::ColumnIndex)   -> Self::Column              
    {
        match self.umatch.generalized_matching_matrix_ref().has_a_match_for_row_index( index ) {
            false   =>  IterWrappedVec::new( self.umatch.source_comb().column( index ).collect() ),
            true    =>  IterWrappedVec::new( self.umatch.target_comb().column( index ).collect() ),
        }
    }       
    fn column_reverse(          &   self, index: &Self::ColumnIndex)   -> Self::ColumnReverse       
    {  
        self.column( index ).reset_and_run_backward()
    }
}




impl < 
        'a,
        BoundaryMatrix, 
        OrderOperatorForRowAndColumnEntries, 
        OrderOperatorForRowAndColumnIndices, 
        IndexForRowsAndColumns,
        EntryForRowsAndColumns,
    > 

    MatrixAlgebra for
    
    DifferentialComb
        < 'a, BoundaryMatrix, >   

    where
        BoundaryMatrix:                         MatrixAlgebra< 
                                                        RowEntry                        =   EntryForRowsAndColumns, 
                                                        RowIndex                        =   IndexForRowsAndColumns, 
                                                        OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                        ColumnEntry                     =   EntryForRowsAndColumns, 
                                                        ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                        OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                    >,
        IndexForRowsAndColumns:                     Clone + Debug + Eq + Hash, // required for the hashing performed by the generalized matching array   
        EntryForRowsAndColumns:                     Clone + Debug + PartialEq + KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,
        OrderOperatorForRowAndColumnEntries:        Clone + JudgeOrder< EntryForRowsAndColumns >,
        OrderOperatorForRowAndColumnIndices:        Clone + JudgeOrder< IndexForRowsAndColumns >,        
        BoundaryMatrix::RingOperator:           DivisionRingOperations,   
{
    type RingOperator                       =   BoundaryMatrix::RingOperator; // The type of ring used to store coefficients in the matrix

    type OrderOperatorForRowEntries         =   OrderOperatorForRowAndColumnEntries;

    type OrderOperatorForRowIndices         =   OrderOperatorForRowAndColumnIndices;

    type OrderOperatorForColumnEntries      =   OrderOperatorForRowAndColumnEntries;

    type OrderOperatorForColumnIndices      =   OrderOperatorForRowAndColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator {
        self.umatch.ring_operator()
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.umatch.order_operator_for_row_entries()
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.umatch.order_operator_for_row_indices()
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.umatch.order_operator_for_column_entries()
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.umatch.order_operator_for_column_indices()
    }
}













impl < 
        'a,
        BoundaryMatrix, 
        OrderOperatorForRowAndColumnEntries, 
        OrderOperatorForRowAndColumnIndices, 
        IndexForRowsAndColumns,
        EntryForRowsAndColumns,
    > 

    MatrixOracleOperations for
    
    DifferentialComb
        < 'a, BoundaryMatrix, >   

    where
        BoundaryMatrix:                         MatrixAlgebra< 
                                                        RowEntry                        =   EntryForRowsAndColumns, 
                                                        RowIndex                        =   IndexForRowsAndColumns, 
                                                        OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                        ColumnEntry                     =   EntryForRowsAndColumns, 
                                                        ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                        OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                    >,
        IndexForRowsAndColumns:                     Clone + Debug + Eq + Hash, // required for the hashing performed by the generalized matching array   
        EntryForRowsAndColumns:                     Clone + Debug + PartialEq + KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,
        OrderOperatorForRowAndColumnEntries:        JudgePartialOrder< EntryForRowsAndColumns >,
        BoundaryMatrix::RingOperator:           DivisionRingOperations,   
{}        

















//  ==============================================================================
//  COMB INVERSE
//  ==============================================================================



/// The inverse of the COMB of a differential U-match decomposition
/// 
/// The differential COMB is the matrix `J` in a differential U-match decomposition `JM = DJ`. For details
/// on this decomposition, see the documentation for [differential U-match module](crate::algebra::matrices::operations::umatch::differential).
/// 
/// # Construction
/// 
/// We obtain the differential COMB `J` from a secondary U-match decomposition `TM = DS`.
/// Concretely, the construction proceeds as follows
/// 
/// - First let `J = S`, then 
/// - For each nonzero entry `M[i,j]` in the generalized matching matrix `M`, replace column  `i` of `J` with column `i` of `T`
/// 
/// It can then be shown, mathematically, that `JM = DJ`, so `J` is a valid differential COMB.
/// 
/// We compute the inverse as follows:
/// 
/// - First let `J^{-1} = T^{-1}`, then 
/// - For each nonzero entry `M[i,j]` in the generalized matching matrix `M`, replace column  `j` of `J^{-1}` with column `j` of `S^{-1}`
/// 
/// It can then be shown, mathematically, that `J^{-1}` is the inverse of `J`.
/// 
/// # Examples
/// 
/// For examples on how to use a differential U-match decomposition, see the modules for
/// 
/// - [Vietoris-Rips complexes](crate::topology::simplicial::from::graph_weighted::VietorisRipsComplex)
/// - [Dowker complexes](crate::topology::simplicial::from::relation)
#[derive(new,Clone,Copy,Debug,PartialEq)]
pub struct DifferentialCombInverse< 
                    'a, 
                    MatrixToFactor,                 
                >      
    where
        MatrixToFactor:                    MatrixAlgebra,
        MatrixToFactor::ColumnIndex:       Hash + Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:          Hash + Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
{
    umatch: &'a GimbledUmatch<  MatrixToFactor >,                 
}




//  ---------------------------------------------------------------------------------------
//  IMPLEMENT MATRIX ORACLE
//  ---------------------------------------------------------------------------------------


impl < 
        'a,
        BoundaryMatrix, 
        OrderOperatorForRowAndColumnEntries, 
        OrderOperatorForRowAndColumnIndices, 
        IndexForRowsAndColumns,
        EntryForRowsAndColumns,
    > 

    MatrixOracle for
    
    DifferentialCombInverse
        < 'a, BoundaryMatrix, >   

    where
        BoundaryMatrix:                         MatrixAlgebra< 
                                                        RowEntry                        =   EntryForRowsAndColumns, 
                                                        RowIndex                        =   IndexForRowsAndColumns, 
                                                        OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                        ColumnEntry                     =   EntryForRowsAndColumns, 
                                                        ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                        OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                    >,
        IndexForRowsAndColumns:                     Clone + Debug + Eq + Hash + 'a, // required for the hashing performed by the generalized matching array   
        EntryForRowsAndColumns:                     Clone + Debug + PartialEq + KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,
        OrderOperatorForRowAndColumnEntries:        JudgePartialOrder< EntryForRowsAndColumns >,
        BoundaryMatrix::RingOperator:           DivisionRingOperations,   
{
    type Coefficient            =   BoundaryMatrix::Coefficient;    // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =   IndexForRowsAndColumns;    // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   IndexForRowsAndColumns;    // The type of column indices    
    type RowEntry               =   EntryForRowsAndColumns;  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   EntryForRowsAndColumns;  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`                 =   
    // type ColumnReverse          =   TwoTypeIterator< 
    //                                     IterWrappedVec< EntryForRowsAndColumns >,
    //                                     <TargetCombInverse<'a, BoundaryMatrix> as MatrixOracle>::Column,
    //                                 >;      
    // type Column                 =   IterWrappedVecReverse< EntryForRowsAndColumns >;  // What you get when you ask for a row with the order of entries reversed
    // type Row                    =   TwoTypeIterator< 
    //                                     <SourceCombInverse<'a, BoundaryMatrix> as MatrixOracle>::Row, 
    //                                     <TargetCombInverse<'a, BoundaryMatrix> as MatrixOracle>::Row,
    //                                 >;  
    // type RowReverse             =   TwoTypeIterator< 
    //                                     <SourceCombInverse<'a, BoundaryMatrix> as MatrixOracle>::RowReverse, 
    //                                     <TargetCombInverse<'a, BoundaryMatrix> as MatrixOracle>::RowReverse,
    //                                 >;  
    type Row                        =   IterWrappedVec< EntryForRowsAndColumns >;
    type RowReverse                 =   IterWrappedVec< EntryForRowsAndColumns >;    
    type Column                     =   IterWrappedVecReverse< EntryForRowsAndColumns >;
    type ColumnReverse              =   IterWrappedVec< EntryForRowsAndColumns >;        
    


    fn structural_nonzero_entry(                   &   self, row: &Self::RowIndex, column: &Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        
        // take the corresponding nonzero entry from the inverse source comb if `row` is a matched column index (yes, column index), and from the inverse target comb otherwise
        match self.umatch.generalized_matching_matrix_ref().has_a_match_for_column_index( row ) {
            false   =>  { self.umatch.target_comb_inverse().structural_nonzero_entry( row, column ) },            
            true    =>  { self.umatch.source_comb_inverse().structural_nonzero_entry( row, column ) },
        }
    }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.umatch.matrix_to_factor_ref().has_column_for_index( index ) 
    }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool {
        self.umatch.matrix_to_factor_ref().has_row_for_index( index )
    }
    fn column_reverse(                &   self, index: &Self::RowIndex   )   -> Self::ColumnReverse               
    { 
        // here we assume a proper Umatch decomposition T * M = D * S
        // we are tryig to build a row of the differntial COMB which is by definition
        // 
        // J = T * Delta_t + S * (I - Delta_t)
        //
        // see the section on Differential COMBs in the Umatch paper for details
        //
        // if `index` is not a matched row index then because the decomposition is
        // proper column `index` of `T` is a standard unit vector. it can therefore
        // be argued that column `index` of `T^{-1}` is a standard unit vector (consider
        // the process of computing the inverse by doing column clearing operations
        // working low-to-high, right-to-left). when we swap the 
        // set of all non-unit-vector rows of `S^{-1}` into `T^{-1}`, we therefore transform
        // this unit vector into the row vector of `S^{-1}` at index `index`.

        let matching_matrix     =   self.umatch.generalized_matching_matrix_ref();
        if matching_matrix.lacks_a_match_for_row_index( index ) { 
            return      
                        // TwoTypeIterator::Version2( 
                        //     self.umatch.source_comb_inverse().column_reverse( index ) 
                        // )
                        IterWrappedVec::new(
                            self.umatch.source_comb_inverse().column_reverse( index ).collect()
                        )
        }

        // otherwise let's take the parts we need from each row of T^{-1} and S^{-1} and combine them
        
        // because the decomposition is proper all the off-diagonal nonzero elements of S^{-1} occur in
        // rows whose indices are matched-column-indices. thus we don't throw away any elements from
        // the S^{-1] column
        let source_inverse_column_rev          =   self.umatch.source_comb_inverse().column_reverse( index );
        
        // we do have to throw out some elements of the target column, corresponding to rows that we swap out for rows of S^{-1}
        let mut target_inverse_column_rev      =   self.umatch.target_comb_inverse().column_reverse( index );
        target_inverse_column_rev.next(); // throw out the leading entry (we get that for free from the source column)
        let target_inverse_column_rev          =   target_inverse_column_rev.filter( // throw out entries indexed by matched-column-indices
                                        |entry| 
                                        matching_matrix.lacks_a_match_for_column_index( & entry.key() ) 
                                    );

        // merge the two columns
        let order_operator      =   self.umatch.matrix_to_factor_ref()
                                        .order_operator_for_column_entries_reverse();
        let merged              =   source_inverse_column_rev.merge_by( 
                                        target_inverse_column_rev,  
                                        | a,b | {  order_operator.judge_lt( a, b )  }
                                    )
                                    .collect::<Vec<_>>();
        return  
                // TwoTypeIterator::Version1(
                    IterWrappedVec::new( merged )
                // );
        
    }
    fn column(             &   self, index: &Self::RowIndex   )   -> Self::Column          {

        return self.column_reverse( index ).reset_and_run_backward();
        // match self.column_reverse( index ) {
        //     TwoTypeIterator::Version1( iter_wrapped_vec ) => {
        //         return iter_wrapped_vec.reset_and_run_backward()
        //     } TwoTypeIterator::Version2( iterator ) => {
        //         let iter_wrapped_vecc = IterWrappedVec::new( iterator.collect() );
        //         return iter_wrapped_vecc.reset_and_run_backward()
        //     }
        // } 
    }     
    fn row(                  &   self, index: &Self::RowIndex)   -> Self::Row              
    {
        match self.umatch.generalized_matching_matrix_ref().has_a_match_for_column_index( index ) {
            false   =>  IterWrappedVec::new( self.umatch.target_comb_inverse().row( index ).collect() ),
            true    =>  IterWrappedVec::new( self.umatch.source_comb_inverse().row( index ).collect() ),
        }
    }       
    fn row_reverse(          &   self, index: &Self::ColumnIndex)   -> Self::RowReverse       
    {  
        match self.umatch.generalized_matching_matrix_ref().has_a_match_for_column_index( index ) {
            false   =>  IterWrappedVec::new( self.umatch.target_comb_inverse().row_reverse( index ).collect() ),
            true    =>  IterWrappedVec::new( self.umatch.source_comb_inverse().row_reverse( index ).collect() ),
        }
    }
}




impl < 
        'a,
        BoundaryMatrix, 
        OrderOperatorForRowAndColumnEntries, 
        OrderOperatorForRowAndColumnIndices, 
        IndexForRowsAndColumns,
        EntryForRowsAndColumns,
    > 

    MatrixAlgebra for
    
    DifferentialCombInverse
        < 'a, BoundaryMatrix, >   

    where
        BoundaryMatrix:                         MatrixAlgebra< 
                                                        RowEntry                        =   EntryForRowsAndColumns, 
                                                        RowIndex                        =   IndexForRowsAndColumns, 
                                                        OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                        ColumnEntry                     =   EntryForRowsAndColumns, 
                                                        ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                        OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                    >,
        IndexForRowsAndColumns:                     Clone + Debug + Eq + Hash + 'a, // required for the hashing performed by the generalized matching array   
        EntryForRowsAndColumns:                     Clone + Debug + PartialEq + KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,
        OrderOperatorForRowAndColumnEntries:        Clone + JudgeOrder< EntryForRowsAndColumns >,
        OrderOperatorForRowAndColumnIndices:        Clone + JudgeOrder< IndexForRowsAndColumns >,        
        BoundaryMatrix::RingOperator:           DivisionRingOperations,   
{
    type RingOperator                       =   BoundaryMatrix::RingOperator; // The type of ring used to store coefficients in the matrix

    type OrderOperatorForRowEntries         =   OrderOperatorForRowAndColumnEntries;

    type OrderOperatorForRowIndices         =   OrderOperatorForRowAndColumnIndices;

    type OrderOperatorForColumnEntries      =   OrderOperatorForRowAndColumnEntries;

    type OrderOperatorForColumnIndices      =   OrderOperatorForRowAndColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator {
        self.umatch.ring_operator()
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.umatch.order_operator_for_row_entries()
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.umatch.order_operator_for_row_indices()
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.umatch.order_operator_for_column_entries()
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.umatch.order_operator_for_column_indices()
    }
}













impl < 
        'a,
        BoundaryMatrix, 
        OrderOperatorForRowAndColumnEntries, 
        OrderOperatorForRowAndColumnIndices, 
        IndexForRowsAndColumns,
        EntryForRowsAndColumns,
    > 

    MatrixOracleOperations for
    
    DifferentialCombInverse
        < 'a, BoundaryMatrix, >   

    where
        BoundaryMatrix:                         MatrixAlgebra< 
                                                        RowEntry                        =   EntryForRowsAndColumns, 
                                                        RowIndex                        =   IndexForRowsAndColumns, 
                                                        OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                        ColumnEntry                     =   EntryForRowsAndColumns, 
                                                        ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                        OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                    >,
        IndexForRowsAndColumns:                     Clone + Debug + Eq + Hash + 'a, // required for the hashing performed by the generalized matching array   
        EntryForRowsAndColumns:                     Clone + Debug + PartialEq + KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,
        OrderOperatorForRowAndColumnEntries:        JudgePartialOrder< EntryForRowsAndColumns >,
        BoundaryMatrix::RingOperator:           DivisionRingOperations, 
{}        





















//  ==============================================================================
//  DIFFERENTIAL U-MATCH STRUCT
//  ==============================================================================


/// A differential U-match decomposition
/// 
/// For definitions see [differential U-match decomposition](crate::algebra::matrices::operations::umatch::differential).
/// 
/// This struct can be used to compute bases for homology, cycle space, boundary spaces,
/// and more.  See the methods on this struct for details.
/// 
/// 
/// 
/// # Data stored internally
/// 
/// This object contains
/// 
/// - A pair of integers `self.min_homology_dimension` and `self.max_homology_dimension` that specify the range of dimensions for which 
///   the decomposition is valid.
/// 
/// - A [Umatch] struct for the submatrix of the differential matrix `D` indexed by row and column indices of 
///   dimension `self.min_homology_dimension - 1 <= d <= self.max_homology_dimension + 1`.
/// 
/// 
/// This is not itself a *differential* U-match decomposition; it contains the same essential information, but formatted in a different form.
/// 
/// 
/// # Valid only for dimensions 0 ..= N
/// 
/// When constructing this object, the [constructor function](DifferentialUmatch::new) takes the differential matrix `D` and the iterable `row_indices` as arguments. 
/// If the user attempts to look up information about the differential U-match decomposition (e.g. a row or column of the differential COMB) using an index
/// of dimension strictly greater than N+1, then the result is not guaranteed to be accurate; it may result in an error or return an invalid result.
///  
/// 
/// 
/// 
/// 
/// # Examples
/// 
/// For examples on how to use a differential U-match decomposition, see the modules for
/// 
/// - [Vietoris-Rips complexes](crate::topology::simplicial::from::graph_weighted::VietorisRipsComplex)
/// - [Dowker complexes](crate::topology::simplicial::from::relation)
#[derive(Debug,Clone,Dissolve)]
pub struct DifferentialUmatch< 
                BoundaryMatrix,                                   
            > 
    where 
        BoundaryMatrix:                     MatrixAlgebra,
        BoundaryMatrix::ColumnIndex:        Hash,
        BoundaryMatrix::RowIndex:           Hash,                     
{
    umatch:                                     GimbledUmatch<  BoundaryMatrix >, 
    min_homology_dimension:                     isize,
    max_homology_dimension:                     isize,
}


impl < 
        BoundaryMatrix, 
        OrderOperatorForRowAndColumnEntries, 
        OrderOperatorForRowAndColumnIndices, 
        IndexForRowsAndColumns,
        EntryForRowsAndColumns,
    > 

    DifferentialUmatch< 
                    BoundaryMatrix,              
                > 
    where
        BoundaryMatrix:                         MatrixAlgebra< 
                                                        RowEntry                        =   EntryForRowsAndColumns, 
                                                        RowIndex                        =   IndexForRowsAndColumns, 
                                                        OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                        ColumnEntry                     =   EntryForRowsAndColumns, 
                                                        ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                        OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                    >
                                                    + MatrixOracleOperations
                                                    + ChainComplex,
        OrderOperatorForRowAndColumnEntries:        JudgePartialOrder< EntryForRowsAndColumns >,
        OrderOperatorForRowAndColumnIndices:        JudgeOrder< IndexForRowsAndColumns >,
        IndexForRowsAndColumns:                     Clone + Debug + Hash + Eq, // required for the hashing performed by the generalized matching array   
        BoundaryMatrix::Coefficient:                Debug,       
        BoundaryMatrix::RowEntry:                   KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >, 
        BoundaryMatrix::ColumnEntry:                KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,  
        BoundaryMatrix::RingOperator:               DivisionRingOperations,   
        EntryForRowsAndColumns:                     Clone + Debug + PartialEq,     
{    



    /// Factors a differential matrix, returning a differential U-match decomposition
    /// 
    /// # Arguments
    /// 
    /// - `boundary_matrix`: the [differential matrix](crate::algebra::matrices::operations::umatch::differential#definition)
    /// - `row_indices`: an iterator that runs over row indices in
    ///    (first) ascending order of dimension, and (second) descending order,
    ///    as determined by the `boundary_matrix.order_operator_for_row_indices()` struct. For example, when working with a Vietoris-Rips
    ///    complex, it should iterate over
    ///    - 0-simplices in descending order of filtration, then
    ///    - 1-simplices in descending order of filtration, then
    ///    - 2-simplices in descending order of filtration,
    ///    et cetera.
    /// 
    /// # Max homology dimension
    /// 
    /// It is not necessary for the `row_indices` iterator to run over all row indices of the differential matrix.
    /// The user can pick a dimension `d`, and then iterate over all row indices of dimension `d` and below. 
    /// We call `d` the **maximum homology dimension** in this use case.
    /// This is useful for
    /// computing persistent homology, since we typically only compute homology in low dimensions (and the number of higher
    /// dimensional simplices can be very high).  If 
    /// - `row_indices` iterates over all row indices of dimension `d` and below, and
    /// - `u = DifferentialUmatch::new( boundary_matrix, row_indices )` and 
    /// - `JM = DJ` is the differential U-match decomposition of the entire differential matrix `D`,
    /// 
    /// then the struct `u` will provide correct information about
    /// - any row `J[i,:]` provided that `i` has dimension `<= d`
    /// - any column `J[:,j]` provided that `j` has dimension `<= d+1`
    /// - (co)homology in dimensions 0 .. d
    /// - (co)cycles and (co)boundaries in dimensions 0 .. d
    /// 
    /// **Nota bene** 
    /// 
    /// - It is also acceptable for the `row_indices` argument to iterate
    /// over row indices in descending order of filtration value, ignoring dimension (so that
    /// 0-simplices, 1-simplices, 2-simplices, etc. can appear in interleaved order). This should
    /// produce an identical output, but we expect time/memory use to be much higher. This is
    /// because grouping/ordering simplices by dimension allows us to take advantage of the 
    /// compress optimization for persistent homology. 
    /// - Whatever order you choose to iterate over row indices, it should respect the 
    ///   [order operator](crate::utilities::order) `boundary_matrix.order_operator_for_row_indices()`. In particular,
    ///   the sequence of row indices shoudl be strictly decreasing, with respect to this
    ///   order operator.
    /// 
    /// # Returns
    /// 
    /// Returns a [DifferentialUmatch], which provides access to bases for
    /// homology, cycles, boundaries, and more.
    pub fn new(
                boundary_matrix:                        BoundaryMatrix,
                min_homology_dimension:                     isize,
                max_homology_dimension:                     isize,
            ) 
            ->  DifferentialUmatch
                    < BoundaryMatrix > 
    {

        // get the row indices in the correct order
        let mut buffer = Vec::new();
        let mut row_indices_in_decreasing_order = Vec::new();

        for degree in min_homology_dimension-1 .. max_homology_dimension+1 {
            // get the row indices for this degree
            buffer.extend(  
                boundary_matrix
                    .basis_vector_indices_for_dimension( degree )
            );
            // reverse the order
            buffer.reverse(); 
            // move from the buffer into the row indices vector
            row_indices_in_decreasing_order.append( &mut buffer );
        }

        // compute the inner U-match decomposition
        let umatch = Umatch::new_with_compression(
                // the matrix we want to factor
                boundary_matrix, 
                // the relevant rows of the matrix, sorted first by dimension and then in reverse lexicographic order
                // computing homology of dimension d, we only need to iterate over simplices of dimension d and below.
                row_indices_in_decreasing_order, 
            );       

        // return the differential U-match decomposition
        DifferentialUmatch { 
            umatch:                     GimbledUmatch::RowMajor(umatch), 
            min_homology_dimension,
            max_homology_dimension,
        }
    }


    /// Returns a column-major version of the differential U-match decomposition, or `None` if the decomposition is already in column-major form.
    /// 
    /// The column-major version may have different source and target COMB's, but will have the same generalized matching matrix.
    pub fn column_major( & self ) -> Option< Self >

        where
            BoundaryMatrix: Clone,
    {
        match self.umatch.column_major() {
            Some(umatch) => {
                Some( DifferentialUmatch {
                    umatch,
                    min_homology_dimension:     self.min_homology_dimension,
                    max_homology_dimension:     self.max_homology_dimension,
                })
            },
            None => None, // already in column-major form
        }
    }

    /// Returns `true` if the differential U-match decomposition is in column-major form
    /// 
    /// Concretely, row-major means that internally this object stores a row-major [Umatch]
    /// decomposition of the differential matrix `D`, and column-major means that it stores
    /// row-major [Umatch] decomposition of the anti-transpose of `D`.
    pub fn is_column_major( & self ) -> bool {
        self.umatch.is_column_major()
    }


    /// Returns a reference to the differential matrix `D` of the differential U-match decomposition
    pub fn boundary_matrix( &self ) -> & BoundaryMatrix 
    { 
        // the boundary matrix is the matrix D in the equation JM = DJ
        self.umatch.matrix_to_factor_ref()
    }

    /// Returns a reference to the generalized matching matrix `M` of the differential U-match decomposition `JM = DJ`
    pub fn generalized_matching_matrix( &self ) -> & GeneralizedMatchingMatrixWithSequentialOrder< 
                                                        IndexForRowsAndColumns, 
                                                        IndexForRowsAndColumns,                                                         
                                                        BoundaryMatrix::Coefficient,
                                                    > 
    { 
        self.asymmetric_umatch().generalized_matching_matrix_ref()
    }


    /// Returns the minimum dimension `d` for which homology is computed
    /// 
    /// In concrete operational terms, this means that the differential COMB `J` in the differential U-match decomposition `JM = DJ`
    /// has rows and columns for indices of dimension `d-1` and above.
    pub fn min_homology_dimension( &self ) -> isize 
    { 
        self.min_homology_dimension 
    }

    /// Returns the maximum dimension `d` for which homology is computed
    /// 
    /// In concrete operational terms, this means that the differential COMB `J` in the differential U-match decomposition `JM = DJ`
    /// has rows and columns for indices of dimension `d+1` and below.
    pub fn max_homology_dimension( &self ) -> isize 
    { 
        self.max_homology_dimension 
    }    

    /// Returns the range of dimensions where homology can be computed using the differential U-match decomposition
    /// 
    /// Concretely, this is `self.min_homology_dimension .. self.max_homology_dimension + 1`.
    pub fn dimensions_where_homology_is_computed( &self ) -> std::ops::Range<isize> {
        self.min_homology_dimension .. self.max_homology_dimension + 1
    }


    // Differential COMB
    // ---------------

    /// The differential COMB
    /// 
    /// Concretely, this is the matrix `J` in the equation `JM = DJ`.  The indices of this matrix are the same as for the differential matrix `D`.
    /// 
    /// For details, see the documenation for [DifferentialComb](super::differential_comb::DifferentialComb)
    pub fn differential_comb( &self ) -> DifferentialComb< '_, BoundaryMatrix, >  
    { DifferentialComb::new( & self.umatch ) }


    /// The inverse of the differential COMB
    /// 
    /// Concretely, this is the inverse of the matrix `J` in the equation `JM = DJ`.  The indices of this matrix are the same as for the differential matrix `D`.
    /// 
    /// For details, see the documenation for [DifferentialComb](super::differential_comb::DifferentialComb)
    pub fn differential_comb_inverse( &self ) -> DifferentialCombInverse< '_, BoundaryMatrix, >  
    { DifferentialCombInverse::new( & self.umatch ) }    


    /// Reference to the internally stored U-match factorization
    /// 
    /// This internally stored U-match is **not** the differential Umatch; it is the `raw data` from
    /// which the differential U-match is constructed.
    /// 
    /// In general, this internal U-match is assymetric, in the sense that it has form `(T,M,D,S)` where `T != S`.
    pub fn asymmetric_umatch(&self) -> & GimbledUmatch<  BoundaryMatrix >
    { & self.umatch }



    /// Returns the ring operator used to perform arithmetic on the coefficients of the differential matrix
    pub fn ring_operator( &self ) -> BoundaryMatrix::RingOperator {
        self.umatch.ring_operator()
    }



    // Index iterators
    // ---------------



    /// Returns the indices of the differential matrix `D` in dimensions where homology has been computed
    /// 
    /// Conceretely, this means each dimension `d` such that `self.min_homology_dimension <= d <= self.max_homology_dimension`.
    /// 
    /// Indices are returned in order of dimension. Indices of equal dimension are returned in order.
    pub fn indices_in_homologically_valid_dimensions( & self ) -> Vec< IndexForRowsAndColumns > {
                
        let mut indices = Vec::new();
        for degree in self.min_homology_dimension .. self.max_homology_dimension + 1 {
            indices.extend(  
                self.boundary_matrix()
                    .basis_vector_indices_for_dimension( degree )
            );
        }
        indices
    }


    /// Runs over the subset of row indices of the boundary matrix visited during row-reduction, in order
    /// 
    /// This includes all indices of dimension `self.min_homology_dimension -1 <= d <= self.max_homology_dimension`
    /// 
    /// Indices are listed in order of dimension. Indices of equal dimension are returned in reverse order.
    pub fn row_reduction_indices( &self ) -> Vec< IndexForRowsAndColumns > { 

        // get the row indices in the correct order
        let mut buffer = Vec::new();
        let mut row_indices_in_decreasing_order = Vec::new();

        for degree in self.min_homology_dimension - 1 .. self.max_homology_dimension + 1 {
            // get the row indices for this degree
            buffer.extend(  
                self.boundary_matrix()
                    .basis_vector_indices_for_dimension( degree )
            );
            // reverse the order
            buffer.reverse(); 
            // move from the buffer into the row indices vector
            row_indices_in_decreasing_order.append( &mut buffer );
        }        
        row_indices_in_decreasing_order
    }  


    /// A hashset containing (a) the [DifferentialUmatch::row_reduction_indices], and (b) all column indices that are matched to row indices in the U-match decomposition
    /// 
    /// In practice we usually force `row_reduction_indices` to run over all simplices (or other types of indices)
    /// in dimension 0 .. d, where d is the maximum homology dimension of interest. Therefore this vector contains
    /// all indices of dimension 0 .. d, plus all **negative simplices** of dimension d+1.
    pub fn row_reduction_indices_plus_negatives( &self ) -> HashSet< IndexForRowsAndColumns > {

        let mut indices = HashSet::new();
        for dimension in self.min_homology_dimension - 1 .. self.max_homology_dimension + 1 {
            // get the row indices for this degree
            indices.extend(  
                self.boundary_matrix()
                    .basis_vector_indices_for_dimension( dimension )
            );
        }
        // add all negative simplices of dimension d+1
        let matching_matrix = self.umatch.generalized_matching_matrix_ref();
        for column_index in matching_matrix.matched_column_indices_in_sequence() {
            indices.insert( column_index.clone() );
        }
        indices
    }        

    /// Indices for a basis of homology
    /// 
    /// Returns a sequence of indices i0 .. ik such that the columns `J[:,i0] .. J[:,ik]`
    /// form a basis for the homology of the chain complex 
    /// - in each dimension `d` such that `self.min_homology_dimension() <= d <= self.max_homology_dimension()`
    /// - here `J` is the [differential COMB](crate::algebra::matrices::operations::umatch::differential)
    /// 
    /// # How this is calculated
    /// 
    /// The indices i0 .. ik are the elements of `self.indices_in_homologically_valid_dimensions()` such that both row `M[ip,:]` and column `M[:,ip]`
    /// are zero for all `p`. That is, i0 .. ik are the indices which select zero-rows and zero-columns of the generalized matching
    /// matrix `M` in the differential U-match decomposition `JM = DJ`.
    /// 
    /// # Identical to cohomology
    /// 
    /// This function returns the same set of indices as `self.cohomology_indices()`.
    pub fn homology_indices( &self ) -> 
        Vec< IndexForRowsAndColumns >
        // SelectedIndices    < '_, BoundaryMatrix,  DecreasingRowIndexIterator, >  
    { 
        let v = self.indices_in_homologically_valid_dimensions();
        let a: SelectedIndices<'_, BoundaryMatrix, Vec<IndexForRowsAndColumns>> = SelectedIndices{ 
            umatch: &self.umatch, 
            row_index_iterator: v.into_iter(), 
            criteria: IndexSelectionCriterion::for_homology() 
        };
        a.collect()
    }   


    /// Indices for a basis of cohomology
    /// 
    /// Returns a sequence of indices i0 .. ik such that the rows `K[i0,:] .. K[ik,:]`
    /// form a basis for the cohomology of the chain complex 
    /// - in each dimension `d` such that `self.min_homology_dimension() <= d <= self.max_homology_dimension()`
    /// - here here K is the *inverse* of the [differential COMB](crate::algebra::matrices::operations::umatch::differential) J in the equation `JM = DJ`
    /// 
    /// # How this is calculated
    /// 
    /// The indices i0 .. ik are the elements of `self.indices_in_homologically_valid_dimensions()` such that row `M[ip,:]` and column `M[:,ip]`
    /// are zero for all `p`. That is, i0 .. ik are the indices which select zero-rows and zero-columns of the generalized matching
    /// matrix `M` in the differential U-match decomposition `JM = DJ`.
    /// 
    /// # Identical to cohomology
    /// 
    /// This function returns the same set of indices as `self.homology_indices()`.
    pub fn cohomology_indices( &self ) -> 
        // SelectedIndices    < '_, BoundaryMatrix,  DecreasingRowIndexIterator, >  
        Vec< IndexForRowsAndColumns >        
    { 
        self.homology_indices()
    }    
    
    /// Cycle indices of the Differential COMB
    /// 
    /// Returns a sequence of indices i0 .. ik such that the columns `J[:,i0] .. J[:,ik]`
    /// form a basis for the cycle space of the chain complex
    /// - in each dimension `d` such that `self.min_homology_dimension() <= d <= self.max_homology_dimension()`
    /// - here `J` is the [differential COMB](crate::algebra::matrices::operations::umatch::differential)
    pub fn cycle_space_indices( &self ) -> 
        Vec< IndexForRowsAndColumns >
        // SelectedIndices    < '_, BoundaryMatrix,  DecreasingRowIndexIterator, >  
    { 
        let v = self.indices_in_homologically_valid_dimensions();
        let a: SelectedIndices<'_, BoundaryMatrix, Vec<IndexForRowsAndColumns>> = SelectedIndices{ 
            umatch: &self.umatch, 
            row_index_iterator: v.into_iter(), 
            criteria: IndexSelectionCriterion::for_cycle_space() 
        };
        a.collect()
    }
    
    /// Indices for a basis of cocycle space
    /// 
    /// Returns a sequence of indices i0 .. ik such that the rows `K[i0,:] .. K[ik,:]`
    /// form a basis for the cocycle space of the chain complex 
    /// - in each dimension `d` such that `self.min_homology_dimension() <= d <= self.max_homology_dimension()`
    /// - here K is the *inverse* of the [differential COMB](crate::algebra::matrices::operations::umatch::differential)
    /// 
    /// # How this is calculated
    /// 
    /// This is the sequence of row indices i0 .. ik such that
    /// - each ip has dimension `<=d`
    /// - each row `M[ip,:]` is zero
    /// 
    /// where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`
    pub fn cocycle_space_indices( &self ) -> 
        // SelectedIndices    < '_, BoundaryMatrix,  DecreasingRowIndexIterator, >  
        Vec<IndexForRowsAndColumns>
    { 
        // SelectedIndices{ umatch: &self.umatch, row_index_iterator: self.row_reduction_indices().into_iter(), criteria: IndexSelectionCriterion::for_cocycle_space() } 
        let v = self.indices_in_homologically_valid_dimensions();
        let a: SelectedIndices<'_, BoundaryMatrix, Vec<IndexForRowsAndColumns>> = SelectedIndices{ 
            umatch: &self.umatch, 
            row_index_iterator: v.into_iter(), 
            criteria: IndexSelectionCriterion::for_cocycle_space()
        };
        a.collect()        
    }
    

    /// Indices for a basis of boundary space
    /// 
    /// Returns a sequence of indices i0 .. ik such that the columns `J[:,i0] .. J[:,ik]`
    /// form a basis for the boundary space of the chain complex 
    /// - in each dimension `d` such that `self.min_homology_dimension() <= d <= self.max_homology_dimension()`
    /// - here J represents the [differential COMB](crate::algebra::matrices::operations::umatch::differential)
    /// 
    /// # How this is calculated
    /// 
    /// This is the sequence of column indices i0 .. ik such that
    /// - each ip has dimension `<=d`
    /// - each column `M[:,ip]` is nonzero
    /// 
    /// where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`
    pub fn boundary_space_indices( &self ) -> 
        // SelectedIndices    < '_, BoundaryMatrix,  DecreasingRowIndexIterator, >  
        Vec< IndexForRowsAndColumns >
    { 
        // SelectedIndices{ umatch: &self.umatch, row_index_iterator: self.row_reduction_indices().into_iter(), criteria: IndexSelectionCriterion::for_cocycle_space() } 
        let v = self.indices_in_homologically_valid_dimensions();
        let a: SelectedIndices<'_, BoundaryMatrix, Vec<IndexForRowsAndColumns>> = SelectedIndices{ 
            umatch: &self.umatch, 
            row_index_iterator: v.into_iter(), 
            criteria: IndexSelectionCriterion::for_boundary_space(),
        };
        a.collect()        
    }        
    
    /// Indices for a basis of coboundary space
    /// 
    /// Returns a sequence of indices i0 .. ik such that the rows `K[i0,:] .. K[ik,:]`
    /// form a basis for the coboundary space of the chain complex 
    /// - in each dimension `d` such that `self.min_homology_dimension() <= d <= self.max_homology_dimension()`
    /// - here K is the *inverse* of the [differential COMB](crate::algebra::matrices::operations::umatch::differential)
    /// 
    /// # How this is calculated
    /// 
    /// This is the sequence of row indices i0 .. ik such that
    /// - each ip has dimension `<=d`
    /// - each row `M[ip,:]` is nonzero
    /// 
    /// where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`
    pub fn coboundary_space_indices( &self ) -> 
        Vec< IndexForRowsAndColumns >
    { 
        // self.umatch
        //     .generalized_matching_matrix_ref()
        //     . bijection_column_indices_to_ordinals_and_inverse()
        //     .vec_elements_in_order() 
        let v = self.indices_in_homologically_valid_dimensions();
        let a: SelectedIndices<'_, BoundaryMatrix, Vec<IndexForRowsAndColumns>> = SelectedIndices{ 
            umatch: &self.umatch, 
            row_index_iterator: v.into_iter(), 
            criteria: IndexSelectionCriterion::for_coboundary_space(),
        };
        a.collect()          
    }

    /// Returns all elements of `self.indices_in_homologically_valid_dimensions()` not contained in `self.homology_indices()`
    /// 
    /// # How this is calculated
    /// 
    /// This is the sequence of column indices i0 .. ik in `self.indices_in_homologically_valid_dimensions()` such that
    /// either column `M[:,ip]` is nonzero or row `M[ip,:]` is nonzero, 
    /// where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`.
    pub fn non_homology_indices( &self ) -> 
        // SelectedIndices    < '_, BoundaryMatrix,  DecreasingRowIndexIterator, >  
        Vec< IndexForRowsAndColumns >
    // { SelectedIndices{ umatch: &self.umatch, row_index_iterator: self.row_reduction_indices().into_iter(), criteria: IndexSelectionCriterion::for_non_homology() } }   
    { 
        let v = self.indices_in_homologically_valid_dimensions();
        let a: SelectedIndices<'_, BoundaryMatrix, Vec<IndexForRowsAndColumns>> = SelectedIndices{ 
            umatch: &self.umatch, 
            row_index_iterator: v.into_iter(), 
            criteria: IndexSelectionCriterion::for_non_homology(),
        };
        a.collect()          
    }  

    /// Returns a list of indices corresponding to a basis for persistent homology
    /// 
    /// Concretely, returns a sequence of indices i0 .. ik such that the columns `J[:,i0] .. J[:,ik]`
    /// form a basis for the persistent homology of the filtered chain complex 
    /// - in dimensions `self.min_homology_dimension() .. self.max_homology_dimension()`, including `self.max_homology_dimension()`
    /// - here `J` represents the differential COMB in the differential U-match decomposition `JM = DJ`
    pub fn persistent_homology_indices( &self ) -> 
        Vec< IndexForRowsAndColumns >

        where
            BoundaryMatrix:   FilteredChainComplex< FiltrationValue: Ord >
    { 
        let mut ph_indices = Vec::new();
        let generalized_matching_matrix = self.generalized_matching_matrix();
        let boundary_matrix = self.boundary_matrix();

        for index in self.indices_in_homologically_valid_dimensions() {
            if generalized_matching_matrix.has_a_match_for_column_index( & index ) {
                continue // in this case the index doesn't correspond to a cycle
            }
            if let Some( column_index ) = generalized_matching_matrix.column_index_for_row_index(&index) {
                let birth_filtration = boundary_matrix.filtration_value_for_basis_vector_with_index( &index ).ok().unwrap();
                let death_filtration = boundary_matrix.filtration_value_for_basis_vector_with_index( &column_index ).ok().unwrap(); 
                if birth_filtration < death_filtration {
                    ph_indices.push( index );
                }               
            }
        }
        ph_indices        
    }  




    /// Returns `Some(bounding_index)` if there exists a `bounding_index` such that `D * J[:,bounding_index] = J[:,index]` up to nonzero scaling
    /// 
    /// Such a `bounding_index` exists iff `M[index, bounding_index] != 0`,
    /// where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`. If this condition holds,
    /// then
    /// ```text
    /// J[:,index] * M[index, bounding_index] = D * J[:,bounding_index]
    /// ```
    /// as a special case of the U-match equation `JM = DJ`.
    pub fn bounding_index_for( &self, index: &IndexForRowsAndColumns ) -> Option< IndexForRowsAndColumns > 
        where 
            BoundaryMatrix:   ChainComplex,
    {
        self.generalized_matching_matrix()
            .column_index_for_row_index( index )
    }


    /// Returns `Some(bounded_index)` if there exists a `bounded_index` such that `D * J[:,index] = J[:,bounded_index]` up to nonzero scaling
    /// 
    /// Such a `bounded_index` exists iff `M[bounded_index,index] != 0`,
    /// where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`. If this condition holds,
    /// then
    /// ```text
    /// J[:,bounded_index] * M[bounded_index, index] = D * J[:,index]
    /// ```
    /// as a special case of the U-match equation `JM = DJ`.
    pub fn bounded_index_for( &self, index: &IndexForRowsAndColumns ) -> Option< IndexForRowsAndColumns > 
        where 
            BoundaryMatrix:   ChainComplex,
    {
        self.generalized_matching_matrix()
            .row_index_for_column_index( index )
    }    



    /// Returns `Some(cobounding_index)` if there exists a `cobounding_index` such that `K[cobounding_index,:] * D = K[index,:]` up to nonzero scaling
    /// 
    /// Such a `cobounding_index` exists iff `M[cobounding_index,index] != 0`,
    /// where 
    /// - `K` is the inverse of the differential COMB `J` in the differential U-match decomposition `JM = DJ`. Therefore `MK = KD`
    /// - `M` is the generalized matching matrix 
    /// If this condition holds, then
    /// ```text
    /// K[cobounding_index,:] * D = M[cobounding_index,index] * K[index,:]
    /// ```
    /// as a special case of the equation `MK = KD`
    pub fn cobounding_index_for( &self, index: &IndexForRowsAndColumns ) -> Option< IndexForRowsAndColumns > 
        where 
            BoundaryMatrix:   ChainComplex,
    {
        self.generalized_matching_matrix()
            .row_index_for_column_index( index )
    }


    /// Returns `Some(cobounded_index)` if there exists a `cobounded_index` such that `K[index,:] * D = K[cobounded_index,:]` up to nonzero scaling
    /// 
    /// Such a `cobounded_index` exists iff `M[index,cobounded_index] != 0`,
    /// where 
    /// - `K` is the inverse of the differential COMB `J` in the differential U-match decomposition `JM = DJ`. Therefore `MK = KD`
    /// - `M` is the generalized matching matrix 
    /// If this condition holds, then
    /// ```text
    /// K[index,:] * D = M[index,cobounded_index] * K[cobounded_index,:]
    /// ```
    /// as a special case of the equation `MK = KD`
    pub fn cobounded_index_for( &self, index: &IndexForRowsAndColumns ) -> Option< IndexForRowsAndColumns > 
        where 
            BoundaryMatrix:   ChainComplex,
    {
        self.generalized_matching_matrix()
            .column_index_for_row_index( index )
    }       







    // /// Checks if the index is a cycle index
    // /// 
    // /// Concretely, returns `true` if 
    // /// - column `J[:,index]` is zero is a cycle, where `J` is the differential COMB in the differential U-match decomposition `JM = DJ`
    // /// - equivalently, if `M[:,index]` is zero, where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`
    // /// 
    // /// # Returns
    // /// 
    // /// - `Ok(true/false)` if the index is valid
    // /// - `Err(())` if the index is not valid, either because it is not in the range of dimensions where homology is computed,
    // ///    or because it is not a valid index for the boundary matrix
    // pub fn is_cycle_index( &self, index: &IndexForRowsAndColumns ) -> Result< bool, () > 
    //     where 
    //         BoundaryMatrix:   ChainComplex,
    // {

    //     let boundary_matrix = self.boundary_matrix();

    //     // get the dimension of the index
    //     let dim = match boundary_matrix.dimension_for_basis_vector_with_index(index) {
    //         Ok(dim) => dim,
    //         Err(_) => return Err(()), // return an error if the index is not valid
    //     };

    //     // ensure the index is in a valid range
    //     let valid =     ( dim >= self.min_homology_dimension ) 
    //                           && 
    //                           ( dim <= self.max_homology_dimension );
    //     if ! valid {
    //         return Err( () ); // return an error if the index is not in-bounds
    //     }        

    //     // check if the column is zero
    //     Ok( 
    //         self.generalized_matching_matrix()
    //             .lacks_a_match_for_column_index( index )            
    //     )
    // }




    // /// Checks if the index is a boundary index
    // /// 
    // /// Concretely, returns `true` if 
    // /// - column `J[:,index]` is zero is a boundary, where `J` is the differential COMB in the differential U-match decomposition `JM = DJ`
    // /// - equivalently, if `M[:,index]` is nonzero, where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`
    // /// 
    // /// # Returns
    // /// 
    // /// - `Ok(true/false)` if the index is valid
    // /// - `Err(())` if the index is not valid, either because it is not in the range of dimensions where homology is computed,
    // ///    or because it is not a valid index for the boundary matrix
    // pub fn is_boundary_index( &self, index: &IndexForRowsAndColumns ) -> Result< bool, () > 
    //     where 
    //         BoundaryMatrix:   ChainComplex,
    // {

    //     let boundary_matrix = self.boundary_matrix();

    //     // get the dimension of the index
    //     let dim = match boundary_matrix.dimension_for_basis_vector_with_index(index) {
    //         Ok(dim) => dim,
    //         Err(_) => return Err(()), // return an error if the index is not valid
    //     };

    //     // ensure the index is in a valid range
    //     let valid =     ( dim >= self.min_homology_dimension ) 
    //                           && 
    //                           ( dim <= self.max_homology_dimension );
    //     if ! valid {
    //         return Err( () ); // return an error if the index is not in-bounds
    //     }        

    //     // check if the column is zero
    //     Ok( 
    //         self.generalized_matching_matrix()
    //             .lacks_a_match_for_column_index( index )            
    //     )
    // }




    // /// Checks if the index is a cocycle index
    // /// 
    // /// Concretely, returns `true` if 
    // /// - row `K[index]` is zero is a cocycle, where `J` is the differential COMB in the differential U-match decomposition `JM = DJ`
    // /// - equivalently, if `M[:,index]` is zero, where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`
    // /// 
    // /// # Returns
    // /// 
    // /// - `Ok(true/false)` if the index is valid
    // /// - `Err(())` if the index is not valid, either because it is not in the range of dimensions where homology is computed,
    // ///    or because it is not a valid index for the boundary matrix
    // pub fn is_cocycle_index( &self, index: &IndexForRowsAndColumns ) -> Result< bool, () > 
    //     where 
    //         BoundaryMatrix:   ChainComplex,
    // {

    //     let boundary_matrix = self.boundary_matrix();

    //     // get the dimension of the index
    //     let dim = match boundary_matrix.dimension_for_basis_vector_with_index(index) {
    //         Ok(dim) => dim,
    //         Err(_) => return Err(()), // return an error if the index is not valid
    //     };

    //     // ensure the index is in a valid range
    //     let valid =     ( dim >= self.min_homology_dimension ) 
    //                           && 
    //                           ( dim <= self.max_homology_dimension );
    //     if ! valid {
    //         return Err( () ); // return an error if the index is not in-bounds
    //     }        

    //     // check if the column is zero
    //     Ok( 
    //         self.generalized_matching_matrix()
    //             .lacks_a_match_for_column_index( index )            
    //     )
    // }    














    // Chain iterators
    // ---------------

    /// A cycle basis for homology 
    /// 
    /// The output is an iterator that runs over a collection of cycle representatives 
    /// - for the homology groups in dimensions `self.min_homology_dimension .. self.max_homology_dimension`, including `self.max_homology_dimension`
    /// 
    /// # How this is calculated
    /// 
    /// The iterator returns the columns `J[:,i0] .. J[:,ik]` of the differential COMB `J` in the differential U-match decomposition `JM = DJ`,
    /// where `i0 .. ik` are the indices in `self.homology_indices()`.
    pub fn homology_basis( &self ) -> 
        // SequenceOfDifferentialCombColumns
        //     < '_, BoundaryMatrix,  IntoIter<IndexForRowsAndColumns>, >  
        SequenceOfColumns
            < 
                DifferentialComb<'_, BoundaryMatrix>,  
                IntoIter<IndexForRowsAndColumns>, 
            >                  
        { 
            SequenceOfColumns::new(
                self.differential_comb(),
                self.homology_indices().into_iter(),
            )
        }


    /// A cocycle basis for cohomology 
    /// 
    /// The output is an iterator that runs over a collection of cocycle representatives for the cohomology groups
    /// in dimensions `self.min_homology_dimension .. self.max_homology_dimension`, including `self.max_homology_dimension`.
    /// 
    /// # How this is calculated
    /// 
    /// The iterator returns the columns `J[:,i0] .. J[:,ik]` of the differential COMB `J` in the differential U-match decomposition `JM = DJ`,
    /// where `i0 .. ik` are the indices in `self.homology_indices()`.
    pub fn cohomology_basis( &self ) -> 
        SequenceOfRows<
            DifferentialCombInverse< '_, BoundaryMatrix, >, 
            IntoIter<IndexForRowsAndColumns>
        >
        { 
            SequenceOfRows::new(
                self.differential_comb_inverse(),
                self.homology_indices().into_iter(),
            )
        }        
    
    /// Indices for a basis of cycle space
    /// 
    /// Returns a sequence of indices i0 .. ik such that the columns `J[:,i0] .. J[:,ik]`
    /// form a basis for the cycle space of the chain complex in dimensions 0 .. d, where
    /// - here J represents the [differential COMB](crate::algebra::matrices::operations::umatch::differential)
    /// - d is the [maximum homology dimension](DifferentialUmatch::new)
    /// 
    /// # How this is calculated
    /// 
    /// This is the sequence of column indices i0 .. ik such that
    /// - each ip has dimension `<=d`
    /// - each column `M[:,ip]` is zero
    /// 
    /// where `M` is the generalized matching matrix in the differential U-match decomposition `JM = DJ`
    pub fn cycle_space_basis( &self ) -> 
        // SequenceOfDifferentialCombColumns    < '_, BoundaryMatrix,  IntoIter<IndexForRowsAndColumns>, >  
        SequenceOfColumns
            < 
                DifferentialComb<'_, BoundaryMatrix>,  
                IntoIter<IndexForRowsAndColumns>, 
            >           
        { 
            SequenceOfColumns::new(
                self.differential_comb(),
                self.cycle_space_indices().into_iter(),
            )
        }
    
    /// A basis for the space of cocycles in dimensions
    /// 
    /// Iterates over `self.cocycle_indices()`, converting each index to the corresponding row of the **inverse** Differential COMB.
    /// 
    /// This forms a basis for the space of cocycles in dimensions `self.min_homology_dimension .. self.max_homology_dimension` (including `self.max_homology_dimension`)
    pub fn cocycle_space_basis( &self ) -> 
        SequenceOfRows<
            DifferentialCombInverse< '_, BoundaryMatrix, >, 
            IntoIter<IndexForRowsAndColumns>
        >
        { 
            SequenceOfRows::new(
                self.differential_comb_inverse(),
                self.cocycle_space_indices().into_iter(),
            )
        } 
    
    /// A basis for the space of boundaries
    /// 
    /// Iterates over `self.boundary_indices()`, converting each index to the corresponding column of the differrential COMB.
    /// 
    /// This forms a basis for the space of boundaries in dimensions `self.min_homology_dimension .. self.max_homology_dimension` (including `self.max_homology_dimension`)
    pub fn boundary_space_basis( &self ) -> 
        SequenceOfColumns
            < 
                DifferentialComb<'_, BoundaryMatrix>,  
                IntoIter<IndexForRowsAndColumns>, 
            >   
        { 
            SequenceOfColumns::new(
                self.differential_comb(),
                self.boundary_space_indices().into_iter(),
            )
        }
    
    /// A basis for the space of coboundaries
    /// 
    /// Iterates over `self.cocycle_indices()`, converting each index to the corresponding row of the **inverse** Differential COMB.
    /// 
    /// This forms a basis for the space of cocycles in dimensions `self.min_homology_dimension .. self.max_homology_dimension` (including `self.max_homology_dimension`)
    pub fn coboundary_space_basis( &self ) -> 
        SequenceOfRows<
            DifferentialCombInverse< '_, BoundaryMatrix, >, 
            IntoIter<IndexForRowsAndColumns>
        >
        { 
            SequenceOfRows::new(
                self.differential_comb_inverse(),
                self.coboundary_space_indices().into_iter(),
            )
        }    

    /// Non-homology indices of the Differential COMB
    /// 
    /// Iterates over `self.non_homology_indices()`, converting each index to the corresponding column of the differrential COMB.
    pub fn non_homology_basis( &self ) -> 
        SequenceOfColumns
            < 
                DifferentialComb<'_, BoundaryMatrix>,  
                IntoIter<IndexForRowsAndColumns>, 
            >   
        { 
            SequenceOfColumns::new(
                self.differential_comb(),
                self.non_homology_indices().into_iter(),
            )
        }


    /// Prints the homology basis vectors to the console
    pub fn print_homology_basis( & self ) {
        for ( counter, basis_vector ) in self.homology_basis().into_iter().enumerate() {
            println!("    Homology basis vector {}: {:?}", counter, basis_vector);
        }
    }


    /// Prints the cohomology basis vectors to the console
    pub fn print_cohomology_basis( & self ) {
        for ( counter, basis_vector ) in self.cohomology_basis().into_iter().enumerate() {
            println!("    Coomology basis vector {}: {:?}", counter, basis_vector);
        }
    }    


    // Barcode
    // -------

    /// Returns the homological barcode of the filtered chain complex
    /// 
    /// This method can only be applied if the `BoundaryMatrix` implements the `FilteredChainComplex` trait.
    /// 
    /// It returns a [Barcode](crate::algebra::matrices::operations::umatch::barcode) for persistent homology
    /// in dimensions `self.min_homology_dimension .. self.max_homology_dimension` (including `self.max_homology_dimension`).
    pub fn barcode
    ( 
        &self,
        return_cycle_representatives:   bool,
        return_bounding_chains:         bool,
    ) 
    -> Barcode<
                BoundaryMatrix::RowIndex, 
                BoundaryMatrix::ColumnEntry 
            > 
    where

        BoundaryMatrix:                         FilteredChainComplex< FiltrationValue = OrderedFloat< f64 > >,
        OrderOperatorForRowAndColumnEntries:        Clone,
        OrderOperatorForRowAndColumnIndices:        Clone,
    {
        get_barcode(
            &self, 
            return_cycle_representatives, 
            return_bounding_chains
        )
    }



    // Betti numbers
    // -------------
    /// The betti numbers of the chain complex
    /// 
    /// - The output is a hashmap `H` such that `H[i]` equals the dimension of homology in dimension `i`.
    /// - The hashmap `H` has keys `self.min_homology_dimension .. self.max_homology_dimension`, including `self.max_homology_dimension`
    pub fn betti_numbers( &self ) -> HashMap< isize, isize > 
        where
            BoundaryMatrix: ChainComplex
    {
        let boundary_matrix = self.boundary_matrix();
        let mut betti_numbers = HashMap::new();
        for dimension in self.dimensions_where_homology_is_computed() {
            betti_numbers.insert(dimension, 0);
        }
        for key in self.homology_indices() {
            let dim = boundary_matrix
                .dimension_for_basis_vector_with_index( &key ).ok().unwrap();
            * betti_numbers.entry( dim ).or_insert(0) += 1;
        }
        betti_numbers
    }    

    // Cycle numbers
    // -------------
    /// The dimensions of the cycle spaces of the chain complex
    /// 
    /// - The output is a hashmap `H` such that `H[i]` equals the dimension of the cycle space in dimension `i`.
    /// - The hashmap `H` has keys `self.min_homology_dimension .. self.max_homology_dimension`, including `self.max_homology_dimension`
    pub fn cycle_space_dimensions( &self ) -> HashMap< isize, isize > 
        where
            BoundaryMatrix: ChainComplex
    {    
        let boundary_matrix = self.boundary_matrix();
        let mut subspace_dimensions = HashMap::new();
        for dimension in self.dimensions_where_homology_is_computed() {
            subspace_dimensions.insert(dimension, 0);
        }
        for key in self.cycle_space_indices() {
            let dim = boundary_matrix
                .dimension_for_basis_vector_with_index( &key ).ok().unwrap();
            * subspace_dimensions.entry( dim ).or_insert(0) += 1;
        }
        subspace_dimensions
    }  


    // Cocycle numbers
    // -------------
    /// The dimensions of the cocycle spaces of the chain complex
    ///
    /// - The output is a hashmap `H` such that `H[i]` equals the dimension of the cocycle space in dimension `i`.
    /// - The hashmap `H` has keys `self.min_homology_dimension .. self.max_homology_dimension`, including `self.max_homology_dimension`
    pub fn cocycle_space_dimensions( &self ) -> HashMap< isize, isize > 
        where
            BoundaryMatrix: ChainComplex
    {    
        let boundary_matrix = self.boundary_matrix();
        let mut subspace_dimensions = HashMap::new();
        for dimension in self.dimensions_where_homology_is_computed() {
            subspace_dimensions.insert(dimension, 0);
        }
        for key in self.cocycle_space_indices() {
            let dim = boundary_matrix
                .dimension_for_basis_vector_with_index( &key ).ok().unwrap();
            * subspace_dimensions.entry( dim ).or_insert(0) += 1;
        }
        subspace_dimensions
    }        

    // Boundary numbers
    // ----------------
    /// The dimensions of the boundary spaces of the chain complex
    /// 
    /// - The output is a hashmap `H` such that `H[i]` equals the dimension of the boundary space in dimension `i`.
    /// - The hashmap `H` has keys `self.min_homology_dimension .. self.max_homology_dimension`, including `self.max_homology_dimension`
    pub fn boundary_space_dimensions( &self ) -> HashMap< isize, isize > 
        where
            BoundaryMatrix: ChainComplex
    {    
        let boundary_matrix = self.boundary_matrix();
        let mut subspace_dimensions = HashMap::new();
        for dimension in self.dimensions_where_homology_is_computed() {
            subspace_dimensions.insert(dimension, 0);
        }
        for key in self.boundary_space_indices() {
            let dim = boundary_matrix
                .dimension_for_basis_vector_with_index( &key ).ok().unwrap();
            * subspace_dimensions.entry( dim ).or_insert(0) += 1;
        }
        subspace_dimensions
    }  


    // Coboundary numbers
    // ----------------
    /// The dimensions of the coboundary spaces of the chain complex
    /// 
    /// - The output is a hashmap `H` such that `H[i]` equals the dimension of the coboundary space in dimension `i`.
    /// - The hashmap `H` has keys `self.min_homology_dimension .. self.max_homology_dimension`, including `self.max_homology_dimension`
    pub fn coboundary_space_dimensions( &self ) -> HashMap< isize, isize > 
        where
            BoundaryMatrix: ChainComplex
    {    
        let boundary_matrix = self.boundary_matrix();
        let mut subspace_dimensions = HashMap::new();
        for dimension in self.dimensions_where_homology_is_computed() {
            subspace_dimensions.insert(dimension, 0);
        }
        for key in self.coboundary_space_indices() {
            let dim = boundary_matrix
                .dimension_for_basis_vector_with_index( &key ).ok().unwrap();
            * subspace_dimensions.entry( dim ).or_insert(0) += 1;
        }
        subspace_dimensions
    }      







    /// Columns indices for the cycle minimization constraint matrix in Escolar Hiraoka 2015 (with extended functionality).
    /// 
    /// Let `birth_column` be the index of a column in the boundary matrix.  This function returns
    /// every column index `i` such that 
    /// - `i` has the same dimension in homology as `birth_column`
    /// - `i` strictly precedes `birth_column` in the linear order on column indices
    /// - column `i` of the Differential COMB is a cycle
    /// - if column `i` ever becomes a boundary, then it does so no later than `birth_column` (if `birth_column` indexes an essential class then this condition is always satisfied)
    /// 
    /// # Arguments
    /// 
    /// - `birth_column`: The column we wish to optimize.
    /// - `get_dimension`: A function which returns the homology dimension of each column index of the boundary matrix.
    ///   This argument is necessary because indices can have multiple different types (simplices, cubes, etc.), so the user
    ///   must specify how dimension is calculated.
    /// 
    /// # Implications for computation
    /// 
    /// This function returns a collection of column indices `i0 .. im` for the [differential COMB](crate::algebra::matrices::operations::umatch::differential), `J`.
    /// If `u = J[:,birth_column]` is the column of the differential COMB indexed `birth_column`, then adding any linear combination of the columns `J[:,i0] .. J[:,im]`
    /// to u  will yield
    /// a cycle representative `v` for a persistent homology class with the same birth and death time as the cycle `u`.  Moreover,
    /// swapping `u` for `v` will not create any new linear dependencies in the basis represented by `J`; that is, it will produce a
    /// basis for cycle space which includes a basis for persistent homology. We can perform this swap on every column
    /// of `J`, simultaneously, and the result will be a valid basis for persistent homology.

    pub fn escolar_hiraoka_indices< DimensionEvaluator >( &self, birth_column: & IndexForRowsAndColumns, mut get_dimension: DimensionEvaluator, ) 
        -> Vec< IndexForRowsAndColumns > 
        where
            DimensionEvaluator:  FnMut( & IndexForRowsAndColumns)-> isize,       
    {        
        // define some parameters
        let order_operator = self.asymmetric_umatch().order_operator_for_row_entries();
        let generalized_matching_matrix_ref      =   self.umatch.generalized_matching_matrix_ref();
        
        // a convenience function; we want to use our order operator, but it technically only operates
        // on matrix entries, not on indices; so we embed indices in entries
        let zero = BoundaryMatrix::RingOperator::zero();
        let to_entry = |i: IndexForRowsAndColumns| EntryForRowsAndColumns::new( i, zero.clone() );        
        
        // calculate some baselines for comparison
        let objective_dimension     =   get_dimension( birth_column );
        let death_column_opt         =   generalized_matching_matrix_ref                              //  death filtration value
                                            .column_index_for_row_index( birth_column );
        
        // now we can enumerate the desired column indices
        let mut columns = Vec::new();
        for row_index in self.row_reduction_indices() {

            // must have equal dimension
            if get_dimension( & row_index ) != objective_dimension { continue }
            // must be born strictly before birth_column
            if order_operator.judge_ge( &to_entry(row_index.clone()), &to_entry(birth_column.clone()) ) { continue }
            // must be a cycle
            if generalized_matching_matrix_ref.has_a_match_for_column_index( & row_index ) { continue }
            // cannot die strictly later than birth_column
            let comp_death_opt          =   generalized_matching_matrix_ref                                      
                                            .column_index_for_row_index( birth_column );
            if let Some( death_column ) = death_column_opt.as_ref() {
                if let Some( comp_death ) = comp_death_opt {
                    if order_operator.judge_gt( &to_entry(comp_death.clone()), &to_entry(death_column.clone()) ) {
                        continue
                    }
                } else {
                    continue
                }
            }
            
            columns.push(row_index);
        }        
        columns
    }


    /// Columns of cycle minimization constraint matrix in Escolar Hiraoka 2015 (**relaxed** and with extended functionality).
    /// 
    /// Similar to [DifferentialUmatch::escolar_hiraoka_indices], but we relax the constraint on which indices to include,
    /// by allowing some indices with equal filtration value.
    /// 
    /// Let `birth_column` be the index of a column in the boundary matrix.  This function returns
    /// every column index `i` such that 
    /// - column `i` of the differential COMB is a cycle with the same dimension in homology as `birth_column`
    /// - birth time of the PH class represented by `i` ≤ birth time of the PH class represented by `birth_column`
    /// - death time of the PH class represented by `i` ≤ death time of the PH class represented by `birth_column`
    /// 
    ///
    /// # Implications for computation
    /// 
    /// Adding any linear combination of these columns to a the column, `u`, of the differential COMB `J` indexed by `birth_column` will yield
    /// a cycle representative `v` for a persistent homology class with the same birth and death time as the cycle `u`.  Moreover,
    /// swapping `u` for `v` will not create any new linear dependencies in the basis; that is, it will produce a
    /// basis for cycle space which includes a basis for persistent homology. **HOWEVER**, this gaurantee only holds
    /// for a single substitution; were we to optimize two different columns and swap them both into the Differential COMB,
    /// we might create linear dependence.
    /// 
    /// # Arguments
    /// 
    /// - `birth_column`: the column we wish to compare
    /// - `get_dimension`: returns the homology dimension of each column index of the boundary matrix
    /// - `filtration_order`: a closure operator that compares two indices on the basis of their **filtration alone**
    pub fn escolar_hiraoka_indices_relaxed< DimensionEvaluator, FiltrationOrder >( &self, birth_column: & IndexForRowsAndColumns, mut get_dimension: DimensionEvaluator, mut filtration_order: FiltrationOrder ) 
        -> Vec< IndexForRowsAndColumns > 
        where
            DimensionEvaluator:  FnMut( & IndexForRowsAndColumns)-> isize,
            FiltrationOrder:  FnMut( & IndexForRowsAndColumns, & IndexForRowsAndColumns )-> Ordering,
    {        
        // define some parameters
        let generalized_matching_matrix_ref      =   self.umatch.generalized_matching_matrix_ref();      
        
        // calculate some baselines for comparison
        let objective_dimension     =   get_dimension( birth_column );
        let death_column_opt         =   generalized_matching_matrix_ref                              //  death filtration value
                                            .column_index_for_row_index( birth_column );
        
        // now we can enumerate the desired column indices
        let mut columns = Vec::new();
        for row_index in self.row_reduction_indices() {

            if row_index == *birth_column { continue }

            // must have equal dimension
            if get_dimension( & row_index ) != objective_dimension { continue }
            // must be born strictly before birth_column
            if filtration_order( &row_index, birth_column) == Ordering::Greater { continue }
            // must be a cycle
            if generalized_matching_matrix_ref.has_a_match_for_column_index( & row_index ) { continue }
            // cannot die strictly later than birth_column
            let comp_death_opt          =   generalized_matching_matrix_ref                                      
                                            .column_index_for_row_index( birth_column );
            if let Some( death_column ) = death_column_opt.as_ref() {
                if let Some( comp_death ) = comp_death_opt {
                    if filtration_order( &comp_death, death_column) == Ordering::Greater { continue }
                } else {
                    continue
                }
            }
            
            columns.push(row_index);
        }        
        columns
    }    

}


//  ----------------------------------------------------------------------------------------
//  ITERATOR OF INDICES
//  ----------------------------------------------------------------------------------------


/// Flag used to indicate a desired type of index
/// 
/// We use this object as an input token to functions that return sets of indices. Different values of the token
/// indicate different sets of indices. 
/// 
/// # Example
/// 
/// To illustrate, suppose we have a differential matrix `D` and we want to compute a basis
/// for the space of boundaries, i.e. the image of `D`.
/// 
/// One way to obtain this basis is to compute a differential Umatch decomposition of `D`, 
/// i.e. a matrix equation
/// 
/// ```text
/// JM = DJ
/// ```
/// 
/// Because `J` is invertible, the image of `D` equals the image of `DJ = JM`. Therefore the image of
/// `D` equals the span of the nonzero columns of `JM`. Because `M` is a generalized matching matrix,
/// this is the same as asking for the span of some subset of the columns of `J`.
/// 
/// In OAT, we can compute a differential Umatch decomposition and store it in a [DifferentialUmatch] object;
/// let's call it `decomp`. We can access the matrix `J` by calling `decomp.differential_comb()`.
/// And we can get a list of the correct column indices to choose by calling `decomp.boundary_indices()`.
/// (These steps are combined by the function call `decomp.boundary_basis()`, which first calls `decomp.boundary_indices()`
/// to get a list of indices, then returns an iterator that runs over the corresponding columns of `J`).
/// 
/// When we call `decomp.boundary_indices()`, it constructs an iterator with the following command:
/// 
/// ```text
/// SelectedIndices{ 
///     umatch:                 decomp.umatch, 
///     row_index_iterator:     decomp.row_reduction_indices().into_iter(), 
///     criteria:               IndexSelectionCriterion::for_boundary_space() 
/// }
/// ```
/// 
/// The [SelectedIndices] object is essentially a wrapper around the `row_index_iterator`, which filters out
/// certain row indices based on a user-provided criterion. In this case, the criterion is `IndexSelectionCriterion::for_boundary_space() `.
#[derive(Debug,Clone)]
pub struct IndexSelectionCriterion{ boundary: bool, bounding: bool, harmonic: bool}

impl IndexSelectionCriterion{
    /// Homology
    fn for_homology() -> Self { IndexSelectionCriterion { boundary: false, bounding: false, harmonic: true } }    
    /// Cycle space
    fn for_cycle_space() -> Self { IndexSelectionCriterion { boundary: true, bounding: false, harmonic: true } }
    /// Cocycle space
    fn for_cocycle_space() -> Self { IndexSelectionCriterion { boundary: false, bounding: true, harmonic: true } }    
    /// Boundary space
    fn for_boundary_space() -> Self { IndexSelectionCriterion { boundary: true, bounding: false, harmonic: false } }
    /// Coboundary space
    fn for_coboundary_space() -> Self { IndexSelectionCriterion { boundary: false, bounding: true, harmonic: false } }
    /// Spans the sum of two subspaces: (1) boundaries, (2) a complement of the space of cycles within the space of all chains.
    fn for_non_homology() -> Self { IndexSelectionCriterion { boundary: true, bounding: true, harmonic: false } }
}

/// Indices of (non)zero rows and columns of `M`
/// 
/// This object is an iterator; it is primarily used as a convenient package for the outputs of functions that return lists
/// of indices related to a differential U-match factorization `JM=DJ`.  For example, [DifferentialUmatch::boundaryIndices]
/// returns a [SelectedIndices] object.
/// 
/// 
/// Depending on how parameters are chosen, this iterator will run over the zero (respectively, nonzero) rows (respectively, column) of `M`.
/// It can also run over any union of this sets. 
#[derive(Debug,Clone)]
pub struct SelectedIndices< 
                    'a,
                    BoundaryMatrix,                  
                    DecreasingRowIndexIterator, 
                > 
    where
        BoundaryMatrix:                     MatrixAlgebra,
        BoundaryMatrix::RowIndex:           Hash,
        BoundaryMatrix::ColumnIndex:        Hash,        
        DecreasingRowIndexIterator:         Clone + IntoIterator,
{
    umatch:                 &'a GimbledUmatch<  BoundaryMatrix >,
    row_index_iterator:     DecreasingRowIndexIterator::IntoIter,
    criteria:               IndexSelectionCriterion,
}      




impl    < 
            'a,
            BoundaryMatrix, 
            OrderOperatorForRowAndColumnEntries, 
            OrderOperatorForRowAndColumnIndices, 
            IndexForRowsAndColumns,
            EntryForRowsAndColumns,
            DecreasingRowIndexIterator, 
        > 

    Iterator for     

    SelectedIndices< 
                    'a,
                    BoundaryMatrix,              
                    DecreasingRowIndexIterator, 
                > 
    where
        DecreasingRowIndexIterator:                 Clone + IntoIterator< Item = BoundaryMatrix::RowIndex >, 
        BoundaryMatrix:                             MatrixAlgebra< 
                                                        RowEntry                        =   EntryForRowsAndColumns, 
                                                        RowIndex                        =   IndexForRowsAndColumns, 
                                                        OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                        ColumnEntry                     =   EntryForRowsAndColumns, 
                                                        ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                        OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                    >,
        IndexForRowsAndColumns:                         Clone + Debug + Hash + Eq, // required for the hashing performed by the generalized matching array   
        BoundaryMatrix::RowEntry:                   KeyValPair< Key = BoundaryMatrix::RowIndex, Val = BoundaryMatrix::Coefficient >, 
        BoundaryMatrix::ColumnEntry:                KeyValPair< Key = BoundaryMatrix::RowIndex, Val = BoundaryMatrix::Coefficient >,  
        BoundaryMatrix::RingOperator:               DivisionRingOperations,        

{
    type Item = IndexForRowsAndColumns;

    fn next(&mut self) -> Option<Self::Item> {
        let matching            =   self.umatch.generalized_matching_matrix_ref();
        let criterion           =   &self.criteria;
        self.row_index_iterator
            .find(  |x|  
                        (   criterion.bounding  &&  matching.has_a_match_for_column_index(x) )
                    ||
                        (   criterion.boundary  &&  matching.has_a_match_for_row_index(x) )
                    ||
                        (   criterion.harmonic  &&  matching.lacks_match_for_index(x)       )            
                )                
    }
}


/// [Type alias](https://doc.rust-lang.org/rust-by-example/types/alias.html) for a sequence of columns of the 
type SequenceOfDifferentialCombColumns< 'a, BoundaryMatrix, IndexIterator >  = 
    SequenceOfColumns<
        DifferentialComb< 'a, BoundaryMatrix >,
        SelectedIndices<  'a, BoundaryMatrix, IndexIterator >
    >;














































mod tests {
    

    
    
    

    
    
    
    
    
    
    
    

    

    // Note this useful idiom: importing names from outer (for mod tests) scope.

    use itertools::Itertools;

    use crate::algebra::chain_complexes::ChainComplex;
    use crate::algebra::matrices::debug::{matrix_oracle_is_internally_consistent, matrix_order_operators_are_internally_consistent, product_is_identity_matrix};
    use crate::algebra::matrices::operations::umatch::differential::DifferentialUmatch;
    use crate::algebra::matrices::operations::MatrixOracleOperations;
    use crate::algebra::matrices::query::{MatrixAlgebra, MatrixOracle};
    use crate::algebra::matrices::types::product::ProductMatrix;
    use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    use crate::algebra::rings::traits::DivisionRingOperations;
    use crate::algebra::vectors::entries::KeyValPair;
    use crate::topology::simplicial::from::graph_weighted::DiagonalEntryIterator;
    use crate::utilities::order::{JudgeOrder, JudgePartialOrder};


    use std::fmt::Debug;
    use std::hash::Hash;



    


    #[test]
    fn test_differential_umatch_random_symmetric_matrix() {

        use crate::algebra::vectors::entries::KeyValGet;
        use crate::algebra::matrices::types::third_party::IntoCSR;
        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::vectors::operations::VectorOperations;        
        use crate::algebra::matrices::operations::umatch::differential::DifferentialUmatch;        
        use crate::topology::simplicial::from::graph_weighted::VietorisRipsComplex;
        use crate::utilities::iterators::general::minmax;        
        use ordered_float::OrderedFloat;  
        use std::sync::Arc;   
        use itertools::Itertools;                   
        
        let number_of_points = 5;
        let min_homology_dimension = 0; 
        let max_homology_dimension = 1; 

        // initialize a random dissimilarity matrix
        let dissimilarity_matrix_vecvec = VecOfVec::random_symmetric_zero_diagonal_with_enclosing_radius_threshold(number_of_points);
        let dissimilarity_matrix =  Arc::new( dissimilarity_matrix_vecvec.into_csr(number_of_points, number_of_points) );                

     
        // initialize the Vietoris-Rips complex
        let ring_operator = crate::algebra::rings::types::native::FieldRationalSize::new();
        let boundary_matrix_data = VietorisRipsComplex::new( dissimilarity_matrix.clone(), number_of_points, ring_operator ).ok().unwrap();
        let boundary_matrix = Arc::new(boundary_matrix_data);  


        // get a comprehensive list of indices; runs over all simplices of dimension <= d in SORTED ORDER
        let mut indices_to_check = boundary_matrix.cliques_in_row_reduction_order(2);
        let order_operator = boundary_matrix.order_operator_for_row_indices();
        indices_to_check.sort_unstable_by( |a,b| order_operator.judge_cmp(a,b) );

            
        // get the differential umatch decomposition
        let differential_umatch = DifferentialUmatch::new(
            boundary_matrix, 
            min_homology_dimension,
            max_homology_dimension,
        );


        // check the differential umatch
        validate_differential_umatch(
            & differential_umatch,
            & indices_to_check,
        );     


        // check the anti-transpose umatch
        validate_differential_umatch(
            & differential_umatch.column_major().unwrap(),
            & indices_to_check,
        ); 


             


    }    













    #[allow(dead_code)]
    fn validate_differential_umatch
        <
            BoundaryMatrix, 
            OrderOperatorForRowAndColumnEntries, 
            OrderOperatorForRowAndColumnIndices, 
            IndexForRowsAndColumns,
            EntryForRowsAndColumns,
        > 
    (
        differential_umatch:    & DifferentialUmatch< BoundaryMatrix >,
        indices_to_check:       & Vec< BoundaryMatrix::RowIndex >
    )        
        where
            BoundaryMatrix:                         MatrixAlgebra< 
                                                            RowEntry                        =   EntryForRowsAndColumns, 
                                                            RowIndex                        =   IndexForRowsAndColumns, 
                                                            OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                            OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                            ColumnEntry                     =   EntryForRowsAndColumns, 
                                                            ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                            OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                            OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                        >
                                                        + MatrixOracleOperations
                                                        + ChainComplex
                                                        + Clone,

            BoundaryMatrix::Coefficient:                Debug,       
            BoundaryMatrix::RowEntry:                   KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >, 
            BoundaryMatrix::ColumnEntry:                KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,  
            BoundaryMatrix::RingOperator:               DivisionRingOperations,      


            IndexForRowsAndColumns:                     Clone + Debug + Eq + Hash, // required for the hashing performed by the generalized matching array   
            EntryForRowsAndColumns:                     Clone + Debug + PartialEq + KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,
            OrderOperatorForRowAndColumnEntries:        Clone + JudgeOrder< EntryForRowsAndColumns >,
            OrderOperatorForRowAndColumnIndices:        Clone + JudgeOrder< IndexForRowsAndColumns >,                    
    {


        let matching_packet         =   differential_umatch.asymmetric_umatch().generalized_matching_matrix_ref_packet();
        let matching_matrix         =   differential_umatch.generalized_matching_matrix();
        let diff_comb               =   differential_umatch.differential_comb();
        let diff_comb_inverse       =   differential_umatch.differential_comb_inverse();        
        let boundary_matrix         =   differential_umatch.boundary_matrix();

        let dj  =   ProductMatrix::new( boundary_matrix, diff_comb.clone() ); // product of differential with comb
        let jm  =   ProductMatrix::new( diff_comb.clone(), matching_packet.clone() ); // product of comb with matching



        // check that the differential comb and its inverse have valid order operators
        assert!(
            matrix_order_operators_are_internally_consistent(
                & diff_comb, 
                indices_to_check.iter().cloned(), 
                indices_to_check.iter().cloned(),
            ).is_ok()
        );

        matrix_order_operators_are_internally_consistent(
            &diff_comb_inverse,
            indices_to_check.iter().cloned(),
            indices_to_check.iter().cloned(),
        ).unwrap_or_else(|e| panic!("{}", e));




        // check that the differential comb and its inverse are internally valid oracles
        assert!(
            matrix_oracle_is_internally_consistent(
                differential_umatch.asymmetric_umatch().source_comb_inverse(), 
                indices_to_check.iter().cloned(), 
                indices_to_check.iter().cloned(),
            )
        );
        assert!(
            matrix_oracle_is_internally_consistent(
                differential_umatch.asymmetric_umatch().target_comb_inverse(), 
                indices_to_check.iter().cloned(), 
                indices_to_check.iter().cloned(),
            )
        );
        assert!(
            matrix_oracle_is_internally_consistent(
                & diff_comb, 
                indices_to_check.iter().cloned(), 
                indices_to_check.iter().cloned(),
            )
        );
        assert!(
            matrix_oracle_is_internally_consistent(
                & diff_comb_inverse, 
                indices_to_check.iter().cloned(), 
                indices_to_check.iter().cloned(),
            )
        );    



        // check that the differential comb and its inverse are upper triangular
        for index in indices_to_check.iter() {
            assert_eq!(
                index,
                & diff_comb
                    .row(index)
                    .next() // get the first nonzero entry of the row, which should be the diagonal
                    .unwrap()
                    .key()
            );
        }

        for index in indices_to_check.iter() {
            assert_eq!(
                index,
                & diff_comb_inverse
                    .row(index)
                    .next() // get the first nonzero entry of the row, which should be the diagonal
                    .unwrap()
                    .key()
            );
        }      



        // check that differential comb inverse comforms to the formula
        for index in indices_to_check.iter() {
            let row = diff_comb_inverse.row(index);
            match matching_matrix.has_a_match_for_column_index(index) {
                true => {
                    assert!( row.eq( 
                        differential_umatch.umatch.source_comb_inverse().row(index) 
                    ) );
                },
                false => {
                    // if the index is not matched, then the differential comb inverse column is nonzero
                    assert!( row.eq( 
                        differential_umatch.umatch.target_comb_inverse().row(index) 
                    ) );
                }
            }
        }    




        // check that the differential comb comforms to the formula
        for index in indices_to_check.iter() {
            let column = diff_comb.column(index);
            match matching_matrix.has_a_match_for_row_index(index) {
                true => {
                    assert!( column.eq( 
                        differential_umatch.umatch.target_comb().column(index) 
                    ) );
                },
                false => {
                    // if the index is not matched, then the differential comb inverse column is nonzero
                    assert!( column.eq( 
                        differential_umatch.umatch.source_comb().column(index) 
                    ) );
                }
            }
        }         




        // check that all unmatched columns of the differential comb are cycles
        for row_index in differential_umatch.row_reduction_indices() {
            if ! matching_matrix.has_match_for_index(&row_index) { 
                assert!( dj.column(&row_index).next().is_none() );
            }
        }



        // check that JM = DJ, where J is the differential comb and M is the generalized matching matrix
        for column_index in matching_matrix.matched_column_indices_in_sequence() {   
            if !                 dj.column(column_index).eq( 
                    jm.column(column_index) 
                ) {
                println!("Column index: {:?}", column_index);
                println!("DJ column: {:#?}", dj.column(column_index).collect_vec());
                println!("JM column: {:#?}", jm.column(column_index).collect_vec());
            }

            assert!( 
                dj.column(column_index).eq( 
                    jm.column(column_index) 
                )  
            );
        }          




        // check that the differential comb and its inverse are indeed inverses        
        assert!(
            product_is_identity_matrix(
                & diff_comb,
                & diff_comb_inverse,
                indices_to_check.iter().cloned()
            )
        );          




   




    

    }













}