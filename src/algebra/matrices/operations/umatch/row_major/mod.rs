//! U-match factorization (uses row-major storage)
//! 
//! # Background
//! 
//! To learn more about U-match factorization, check out the [umatch](crate::algebra::matrices::operations::umatch) module, or see the paper [U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co)homology](https://arxiv.org/abs/2108.08831), by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.
//! All COMB's are calculated as described in this paper.
//! 
//! # Quick start
//!
//! To factor a matrix, call [Umatch::new].  This will produce a [Umatch], which you can use to obtain copies of the generalized matching matrix, the COMB's, and more.
//! Here's an example:
//! 
//! ```
//! use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
//! use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
//! use oat_rust::algebra::matrices::operations::umatch::row_major::Umatch;
//! use oat_rust::algebra::matrices::types::product::ProductMatrix;
//! use oat_rust::algebra::matrices::debug::product_is_identity_matrix;
//! use oat_rust::algebra::matrices::query::MatrixOracle;
//! use oat_rust::algebra::matrices::display::print_indexed_columns;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::order::{OrderOperatorAuto, ReverseOrder};
//! use itertools::Itertools;
//! 
//! // DEFINE INPUTS
//! // ===============================
//! 
//! // define the ring operator and order operator
//! let modulus             =   5;
//! let ring_operator       =   PrimeOrderField::new( modulus );        
//! let order_operator      =   OrderOperatorAuto;
//! 
//! // define the matrix we wish to factor
//! let matrix_to_factor       =   & VecOfVec::new( 
//!                                 vec![   
//!                                     vec![(0,1), (1,1), (2,1)],  // row 0
//!                                     vec![                   ],  // row 1
//!                                     vec![              (2,1)],  // row 2
//!                                 ] 
//!                             ).ok().unwrap();
//! let matrix_to_factor       =   MatrixAlgebraPacket::with_default_order( &matrix_to_factor, ring_operator.clone());
//!                                 
//! // COMPUTE U-MATCH
//! // ===============================
//!                                 
//! let umatch
//!     =   Umatch::new(
//!             matrix_to_factor,  // the matrix we wish to factor
//!             (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
//!         );
//!     
//!     
//! // INSPECT FACTORIZATION
//! // ===============================
//!     
//! // extract T, T^{-1}, S, S^{-1}, and M
//! let t           =   umatch.target_comb();        // the target COMB
//! let tinv        =   umatch.target_comb_inverse();    // inverse of the the target COMB
//! let s           =   umatch.source_comb();          // the source COMB
//! let sinv        =   umatch.source_comb_inverse();      // inverse of the source COMB
//! let m           =   umatch.generalized_matching_matrix_ref();         // the generalized matching matrix
//!     
//!     
//! println!("\nColumns of the target COMB");     print_indexed_columns( &t, 0..3 ); 
//! println!("\nColumns of the source COMB");     print_indexed_columns( &s, 0..3 ); 
//! println!("\nColumns of the generalized matching matrix"); print_indexed_columns( &m, 0..3 ); 
//!     
//! // this will print the following:
//! //
//! // Columns of the target COMB
//! // column 0: [(0, 1)]
//! // column 1: [(1, 1)]
//! // 
//! // Columns of the   source COMB
//! // column 0: [(0, 1)]
//! // column 1: [(1, 1), (0, 3)]
//! // column 2: [(2, 1), (0, 2)]
//! // 
//! // Columns of the generalized matching matrix
//! // column 0: [(0, 1)]
//! // column 1: []
//! // column 2: [(1, 1)]
//! 
//! // SOLVE Ax = b FOR x
//! // ===============================
//! 
//! let b   =   [ (0,1), (2,1) ];
//! let x   =   umatch.solve_dx_equals_b( b.clone() ).unwrap();
//! let dx  =   umatch.multiply_dx(x);
//! assert!( dx.eq( b ) );
//!     
//!     
//! // VERIFY THE CALCULATION
//! // ===============================
//!     
//! // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
//! product_is_identity_matrix( &s, &sinv, 0..3 );
//!     
//! // check that the product of the target COMB with its inverse is identity: T * T^{-1} = I
//! product_is_identity_matrix( &t, &tinv, 0..3 );
//!     
//! // check the factorization: T^{-1} * D * S = M
//! let rinv_d   = ProductMatrix::new( &tinv,   &matrix_to_factor   );      
//! let rinv_d_c = ProductMatrix::new( &rinv_d, &s                  );                
//! for row_index in 0 .. 3 { 
//!     assert_eq!(
//!         rinv_d_c.row( &row_index ).collect_vec(),
//!         m.row( &row_index ).collect_vec()
//!     ) 
//! }   
//! ```
//! 
//! 

pub mod comb;
pub mod construction;


//  ==================================



use itertools::Itertools;



use comb::*;
use construction::*;
use num::rational::Ratio;
use ordered_float::OrderedFloat;
use sprs::linalg::ordering::order;

use crate::algebra::matrices::types::product::ProductMatrix;
use crate::algebra::matrices::operations::solve::triangle::{TriangularSolveForColumnVectorReverse, TriangularSolveForRowVector};
use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::matrices::types::bimajor::MatrixBimajorData;
use crate::algebra::matrices::types::packet::MatrixAlgebraPacket;
use crate::algebra::vectors::entries::KeyValNew;
use crate::algebra::matrices::operations::transform_entry_wise::{ReindexMatrixColumns, ReindexSquareMatrix};
use crate::algebra::matrices::types::matching::{GeneralizedMatchingMatrixWithSequentialOrder};
use crate::algebra::matrices::types::scalar_diagonal_triangle::SumOfScalarAndStrictlyUpperTriangularMatrices;

use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
use crate::algebra::matrices::query::{ MatrixAlgebra, MatrixOracle, };
use crate::algebra::matrices::operations::iterate_rows_and_columns::SequenceOfReverseColumns;
use crate::algebra::matrices::operations::combine_rows_and_columns::{LinearCombinationOfColumns, LinearCombinationOfColumnsReverse, LinearCombinationOfRows };
use crate::topology::simplicial::simplices::weighted::WeightedSimplex;
use crate::utilities::functions::evaluate::{ EvaluateFunctionFnMutWrapper };

use crate::utilities::iterators::general::{FilterOutMembers};
use crate::utilities::iterators::merge::hit::{IteratorsMergedInSortedOrder};
use crate::algebra::vectors::entries::{KeyValGet, KeyValPair};
use crate::algebra::rings::traits::{SemiringOperations, DivisionRingOperations};

use crate::utilities::order::{JudgeOrder, OrderOperatorAuto, OrderOperatorByKey, OrderOperatorByKeyCustom, ReverseOrder};


use crate::algebra::vectors::operations::{Scale, Simplify, VectorOperations};

use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;
use std::iter::Cloned;



use crate::algebra::matrices::operations::solve::echelon::{RowEchelonSolver};

use crate::algebra::matrices::operations::transform_vector_wise::{OnlyColumnIndicesInsideCollection, OnlyColumnIndicesOutsideCollection, OnlyRowIndicesInsideCollection, OnlyRowIndicesOutsideCollection};





use derive_getters::Dissolve;






//  =========================================================================================================
//  U-MATCH OBJECT
//  =========================================================================================================





/// A [U-match decomposition](https://arxiv.org/abs/2108.08831)
/// 
/// This object represents a U-match decomposition `TM=DS` (equivalently `RM = DC`) of a matrix`D`. Internally, it stores three pieces of information
/// 
/// - The matrix to be factored, `D`. (This must implement the [MatrixOracle] and [MatrixAlgebra] traits. These traits
///   are automatically inherited by references, so you can also pass in a reference `&D`).
/// - A copy of the sparse matrix `R_{\rho \rho}` defined in [Hang et al. 2021](https://arxiv.org/abs/2108.08831). 
///   This is an upper triangular matrix with 1's on the diagonal. The off-diagonal entries are stored in a [VecOfVec](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec)
///   data structure.
/// - The matrix `M`, stored as a [GeneralizedMatchingMatrixWithSequentialOrder].
///   This data structure essentially encodes a list of tuples `(r0,c0,x0), .., (rN,cN,xN)` such that `M[ri,ci] = xi` for all `i`.
///   We call `i` the *ordinal* of `ri` (respectively, of `ci`). The data is organized such that `r0 < .. < rN`, where order is
///   determined by the order operator `D.order_operator_for_row_indices()`.  The column indices `c0, .., cN` are NOT sorted in
///   ascending order. Because `M` is a generalized matching matrix, there are no repeats in the sequence of row indices; nor are
///   ther repeats in the sequence of column indices.
/// 
///   As explained by the Inner Identities in [Hang et al. 2021](https://arxiv.org/abs/2108.08831), this data is enough to rapidly calculate
///   any row or column of `T, S`, or their inverses.
#[derive(Clone,Debug,Dissolve,Eq,PartialEq)] // can't automatically derive PartialOrd or Ord for this struct because can't derive PartialOrd or Ord for GeneralizedMatchingMatrixWithSequentialOrder (because GeneralizedMatchingMatrixWithSequentialOrder contains hashmaps, which don't impelment PartialOrd or Ord)
pub struct Umatch< MatrixToFactor > 
    where   
        MatrixToFactor:                     MatrixAlgebra,
        MatrixToFactor::ColumnIndex:        Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
{
    matrix_to_factor:                       MatrixToFactor,
    matching:                               GeneralizedMatchingMatrixWithSequentialOrder< 
                                                MatrixToFactor::ColumnIndex, 
                                                MatrixToFactor::RowIndex, 
                                                MatrixToFactor::Coefficient 
                                            >,
    matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal:   
                                            MatrixBimajorData<
                                                VecOfVec< usize, MatrixToFactor::Coefficient >,
                                                VecOfVec< usize, MatrixToFactor::Coefficient >,
                                            >,  
}


//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- CONSTRUCTORS
//  ---------------------------------------------------------------------------------------------------------



impl < MatrixToFactor >  

    Umatch 
        < MatrixToFactor >  
    
    where   
        MatrixToFactor:                        MatrixAlgebra,
        MatrixToFactor::ColumnIndex:           Hash + std::cmp::Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:              Hash + std::cmp::Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct 
{

    /// Generate a new U-match factorization
    /// 
    /// # Arguments
    /// 
    /// - `matrix_to_factor`: matrix you want to factor
    /// - `row_indices_in_reverse_order`: an iterator that runs over the row indices of the matrix, in *strictly descending order*
    pub fn new
            < RowIndicesInReverseOrder > ( 
                matrix_to_factor:                       MatrixToFactor, 
                row_indices_in_reverse_order:           RowIndicesInReverseOrder,        
            ) 
        -> 
        Self

    where   
            MatrixToFactor::RingOperator:       DivisionRingOperations< Element = MatrixToFactor::Coefficient >,           
            MatrixToFactor::RowEntry:           KeyValPair,
            MatrixToFactor::ColumnEntry:        KeyValPair,            
            RowIndicesInReverseOrder:           Iterator< Item = MatrixToFactor::RowIndex >,
            
    {
        
        let ( comb_target_inv_off_diag_pivot_block, matching ) : ( VecOfVec<usize, MatrixToFactor::Coefficient>, GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient > )
            = get_pivot_block_of_target_comb_inverse_with_deleted_diagonal( 
                    & matrix_to_factor, 
                    row_indices_in_reverse_order, 
                );

        let comb_target_inv_off_diag_pivot_block
            =   MatrixBimajorData { 
                    matrix_columns_data:    comb_target_inv_off_diag_pivot_block
                                                .transpose_deep( matching.number_of_structural_nonzeros() ) // the number of rows we specify is gauranteed to be correct; it makes the matrix square
                                                .unwrap(),                     
                    matrix_rows_data:       comb_target_inv_off_diag_pivot_block, 
                };
        
        Umatch{ 
                matrix_to_factor, 
                matching, 
                matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal:     comb_target_inv_off_diag_pivot_block,   
            }
        
    }



    /// Same as [Umatch::new], but applies the compress optimization.
    /// 
    /// Concretely, this means that when iterating over row indices, the solver skips indices that have
    /// already been identified as nonzero columns in the matching array.
    /// 
    /// This method **only works correctly** if the set of matched row indices is disjoint from the set
    /// of matched column indices in the actual U-match decomposition of the matrix to be factored.
    /// It can be shown mathematically that this condition is always satisfied when the matrix to
    /// be factored is the differential matrix of a chain complex.
    pub fn new_with_compression
            < RowIndicesInReverseOrder, IndexForRowsAndColumns, EntryForRowsAndColumns, Coefficient > ( 
                matrix_to_factor:              MatrixToFactor, 
                row_indices_in_reverse_order:                RowIndicesInReverseOrder,
            ) 
        -> 
        Umatch< MatrixToFactor, > 

    where   
        IndexForRowsAndColumns:     Clone + Debug + Eq + Hash, // hash is required for the hashing performed by the generalized matching array
        EntryForRowsAndColumns:     PartialEq + KeyValPair< Key=IndexForRowsAndColumns, Val=Coefficient >,
        MatrixToFactor:             MatrixAlgebra<
                                        ColumnIndex=                    IndexForRowsAndColumns,  // for the pareto short circuit to work, rows and columns must have the same index type
                                        RowIndex=                       IndexForRowsAndColumns,  // for the pareto short circuit to work, rows and columns must have the same index type
                                        RowEntry=                       EntryForRowsAndColumns,
                                        ColumnEntry=                    EntryForRowsAndColumns,                        
                                        RingOperator:                   DivisionRingOperations< Element =  Coefficient >, // the ring operator for the coefficient ring
                                        Coefficient=                    Coefficient,  // the coefficient type        
                                    >
                                    + MatrixOracleOperations,  
        RowIndicesInReverseOrder:               IntoIterator< Item = MatrixToFactor::RowIndex >,                    
        Coefficient:                Clone + Debug + PartialEq,                
    {        
        let ( comb_target_inv_off_diag_pivot_block, matching ) : 
            ( 
                VecOfVec<usize, 
                MatrixToFactor::Coefficient>, 
                GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, 
                MatrixToFactor::RowIndex, 
                MatrixToFactor::Coefficient > 
            )
            =   target_comb_inv_off_diag_pivot_block_skipmatched( 
                    & matrix_to_factor, 
                    row_indices_in_reverse_order, 
                );

        // for computational efficiency, store a copy of the pivot block in AND a copy of its transpose
        let comb_target_inv_off_diag_pivot_block
            =   MatrixBimajorData { 
                    matrix_columns_data: comb_target_inv_off_diag_pivot_block.transpose_deep( matching.number_of_structural_nonzeros() ).unwrap(), // the number of rows we specify is gauranteed to be correct; it makes the matrix square                
                    matrix_rows_data: comb_target_inv_off_diag_pivot_block, 
                };                
        
        Umatch{ 
                matrix_to_factor, 
                matching, 
                matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal:     comb_target_inv_off_diag_pivot_block,   
            }
        
    }


    /// Returns the [ring operator](crate::algebra::rings) for the coefficient ring used in the factorization.
    pub fn ring_operator( &self ) -> MatrixToFactor::RingOperator {
        self.matrix_to_factor.ring_operator()
    }
}





//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- GENERAL IMPLEMENTATIONS (THERE ARE SPECIFIC IMPLEMENTATIONS FOR ColumnIndex=RowIndex BELOW)
//  ---------------------------------------------------------------------------------------------------------



impl < MatrixToFactor >  

    Umatch 
    < MatrixToFactor >  
    
    where   
        MatrixToFactor:                             MatrixAlgebra,    
        MatrixToFactor::RingOperator:               DivisionRingOperations,        
        MatrixToFactor::ColumnIndex:                Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:                   Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct 
        MatrixToFactor::RowEntry:                   KeyValPair,
        MatrixToFactor::ColumnEntry:                KeyValPair,

{


    //  =========================================================================================================
    //  U-MATCH REF OBJECT
    //  =========================================================================================================


    /// Rank of the factored matrix
    /// 
    /// Equivalently, 
    /// - the dimension of the image of the linear map represented by the matrix
    /// - the number of nonzero entries in the generalized matching matrix of the U-match factorization
    pub fn rank( &self ) -> usize        
    {
        self.matching.number_of_structural_nonzeros()
    }  



    /// Returns a copy of the order comparator for (column-index, coefficient) pairs
    pub fn order_operator_for_row_entries( &self ) -> MatrixToFactor::OrderOperatorForRowEntries 
    { self.matrix_to_factor.order_operator_for_row_entries() } 
    
    /// Returns a copy of the **inverted** order comparator for (column-index, coefficient) pairs
    pub fn order_operator_for_row_entries_reverse( &self ) -> ReverseOrder< MatrixToFactor::OrderOperatorForRowEntries >
    { ReverseOrder::new(self.matrix_to_factor.order_operator_for_row_entries()) }    

    /// Returns a copy of the order comparator for row indices
    pub fn order_operator_for_row_indices( &self ) -> MatrixToFactor::OrderOperatorForRowIndices
    { self.matrix_to_factor.order_operator_for_row_indices() } 
    
    /// Returns a copy of the **inverted** order comparator for row indices
    pub fn order_operator_for_row_indices_reverse( &self ) -> ReverseOrder< MatrixToFactor::OrderOperatorForRowIndices >
    { ReverseOrder::new(self.matrix_to_factor.order_operator_for_row_indices()) }                

    /// Returns a copy of the order comparator for (row-index, coefficient) pairs
    pub fn order_operator_for_column_entries( &self ) -> MatrixToFactor::OrderOperatorForColumnEntries 
    { self.matrix_to_factor.order_operator_for_column_entries() }  
    
    /// Returns a copy of the **inverted** order comparator for (row-index, coefficient) pairs
    pub fn order_operator_for_column_entries_reverse( &self ) -> ReverseOrder< MatrixToFactor::OrderOperatorForColumnEntries >
    { ReverseOrder::new( self.matrix_to_factor.order_operator_for_column_entries() ) } 

    /// Returns a copy of the order comparator for (row-index, coefficient) pairs
    pub fn order_operator_for_column_indices( &self ) -> MatrixToFactor::OrderOperatorForColumnIndices 
    { self.matrix_to_factor.order_operator_for_column_indices() }  
    
    /// Returns a copy of the **inverted** order comparator for (row-index, coefficient) pairs
    pub fn order_operator_for_column_indices_reverse( &self ) -> ReverseOrder< MatrixToFactor::OrderOperatorForColumnIndices >
    { ReverseOrder::new( self.matrix_to_factor.order_operator_for_column_indices() ) }       

    
    /// The sequence of matched row indices in *ascending order*
    /// 
    /// Concretely, this is the sequence of matched row indices `r_0 < .. < r_k`, where
    /// order is deteremined by the order operator for row indices associated with the factored matrix.
    pub fn matched_row_indices_in_ascending_order( &self ) -> &Vec< MatrixToFactor::RowIndex > {
        self.matching.matched_row_indices_in_sequence()
    }

    /// The sequence of matched column indices, ordered according to the associated row indices
    /// 
    /// Concretely, this is the sequence of matched column indices `c_0, .., c_k`, obtained by 
    /// ordering the sequence of matched row-column index pairs `(r0,c0), .., (rk,ck)` 
    /// such that `r_0 < .. < r_k`.
    /// 
    /// **In particular, there is no guarantee that `c_0 < .. < c_k`**.
    pub fn matched_column_indices_in_matched_row_order( &self ) -> &Vec< MatrixToFactor::ColumnIndex > {
        self.matching.matched_column_indices_in_sequence()
    }  

    /// The sequence of matched column indices in *ascending order*
    /// 
    /// Concretely, this is the sequence of matched column indices `c_0 < .. < c_k`, where
    /// order is deteremined by the order operator for column indices associated with the factored matrix.
    /// 
    /// # Performance
    /// 
    /// The U-match data structure stores matched column indices in a different order. Thus to obtain this
    /// sequence, we must copy the stored data, and sort the column indices according to the order operator.
    /// If all you need is the sequence of matched column indices, use [Umatch::matched_column_indices_in_matched_row_order] instead.
    pub fn matched_column_indices_in_ascending_order( &self ) -> Vec< MatrixToFactor::ColumnIndex > {
        let mut indices = self.matching.matched_column_indices_in_sequence().clone();
        let order_operator = self.matrix_to_factor.order_operator_for_column_indices();
        indices.sort_by( |a,b| order_operator.judge_cmp( a, b ) );
        indices
    }      
    
   

    /// Returns the  target COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the factored matrix.
    pub fn target_comb( &self ) -> TargetComb< '_, MatrixToFactor >  {
        TargetComb{ umatch: self }
    }

    /// Returns the  inverse of the target COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the factored matrix.
    pub fn target_comb_inverse( &self ) -> TargetCombInverse< '_, MatrixToFactor >  {
        TargetCombInverse{ umatch: self }
    }  
    
    /// Returns the  source COMB, indexed by `MatrixToFactor::ColumnIndex`
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the factored matrix.
    pub fn source_comb( &self ) -> SourceComb< '_, MatrixToFactor >  {
        SourceComb{ umatch: self }
    }

    /// Returns the  inverse of the source COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the factored matrix.
    pub fn source_comb_inverse( &self ) -> SourceCombInverse< '_, MatrixToFactor >  {
        SourceCombInverse{ umatch: self }
    }      

    /// Returns a reference to the matching array of the internally stored  U-match factorization.
    pub fn generalized_matching_matrix_ref( &self ) -> & GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient > { & self.matching }

    /// Returns a reference to the factored matrix of the internally stored U-match factorization.
    pub fn matrix_to_factor_ref( &self ) -> & MatrixToFactor { & self.matrix_to_factor }    


    // pub fn packet_comb_target< 'a >( &'a self ) -> 
    //     MatrixAlgebraPacket< TargetComb< 'a, MatrixToFactor >, MatrixToFactor::RingOperator, MatrixToFactor::OrderOperatorForRowEntries, MatrixToFactor::OrderOperatorForColumnEntries >
    //     where
    //         MatrixToFactor::RingOperator:               Clone,
    //         MatrixToFactor::OrderOperatorForRowEntries:    Clone,
    //         MatrixToFactor::OrderOperatorForColumnEntries:    Clone,            
    // {
    //     MatrixAlgebraPacket{ matrix: self.target_comb(), ring: self.ring_operator(), row_entry_order: self.order_operator_for_row_entries(), col_entry_order: self.order_operator_for_column_entries() }
    // }

    // /// Returns a convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket) for the  target COMB
    // pub fn comb_target_packet( &self ) -> 
    //     MatrixAlgebraPacket< TargetComb< '_, MatrixToFactor, MatrixToFactor::RingOperator, MatrixToFactor::OrderOperatorForRowEntries, OrderOperatorForRowIndices, MatrixToFactor::OrderOperatorForColumnEntries >, MatrixToFactor::RingOperator, MatrixToFactor::OrderOperatorForColumnEntries, MatrixToFactor::OrderOperatorForColumnEntries >
    //     where
    //         MatrixToFactor::RingOperator:               Clone,
    //         MatrixToFactor::OrderOperatorForRowEntries:    Clone,
    //         MatrixToFactor::OrderOperatorForColumnEntries:    Clone,            
    // {
    //     MatrixAlgebraPacket{ matrix: self.target_comb(), ring: self.ring_operator(), order_operator_for_row_entries: self.order_operator_for_column_entries(), order_operator_for_column_entries: self.order_operator_for_column_entries() }
    // }

    // /// Returns a convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket) for the  inverse of the source COMB
    // pub fn comb_target_inv( &self ) -> 
    //     MatrixAlgebraPacket< TargetCombInverse< '_, MatrixToFactor >, MatrixToFactor::RingOperator, MatrixToFactor::OrderOperatorForColumnEntries, MatrixToFactor::OrderOperatorForColumnEntries >
    //     where
    //         MatrixToFactor::RingOperator:               Clone,
    //         MatrixToFactor::OrderOperatorForRowEntries:    Clone,
    //         MatrixToFactor::OrderOperatorForColumnEntries:    Clone,            
    // {
    //     MatrixAlgebraPacket{ matrix: self.target_comb_inverse(), ring: self.ring_operator(), order_operator_for_row_entries: self.order_operator_for_column_entries(), order_operator_for_column_entries: self.order_operator_for_column_entries() }
    // }

    // /// Returns a convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket) for the  source COMB
    // pub fn comb_source_packet( &self ) -> 
    //     MatrixAlgebraPacket< 
    //             SourceComb< '_, MatrixToFactor >, 
    //             MatrixToFactor::RingOperator, 
    //             MatrixToFactor::OrderOperatorForRowEntries, 
    //             MatrixToFactor::OrderOperatorForRowEntries
    //         >
    //     where
    //         MatrixToFactor::RingOperator:               Clone,
    //         MatrixToFactor::OrderOperatorForRowEntries:    Clone,         
    // {
    //     MatrixAlgebraPacket{ matrix: self.source_comb(), ring: self.ring_operator(), order_operator_for_row_entries: self.order_operator_for_row_entries(), order_operator_for_column_entries: self.order_operator_for_row_entries() }
    // }

    // /// Returns a convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket) for the inverse of the  source COMB
    // pub fn comb_source_inv_packet( &self ) -> 
    //     MatrixAlgebraPacket< 
    //             SourceCombInverse< '_, MatrixToFactor >, 
    //             MatrixToFactor::RingOperator, 
    //             MatrixToFactor::OrderOperatorForRowEntries, 
    //             MatrixToFactor::OrderOperatorForRowEntries,
    //         >
    //     where
    //         MatrixToFactor::RingOperator:               Clone,
    //         MatrixToFactor::OrderOperatorForRowEntries:    Clone,          
    // {
    //     MatrixAlgebraPacket{ matrix: self.source_comb_inverse(), ring: self.ring_operator(), order_operator_for_row_entries: self.order_operator_for_row_entries(), order_operator_for_column_entries: self.order_operator_for_row_entries() }
    // }       

    /// Returns a reference to the matching array of the internally stored  U-match factorization, wrapped in a convenient convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket)
    pub fn generalized_matching_matrix_ref_packet( &self ) 
        -> MatrixAlgebraPacket< 
            & GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient >,
            MatrixToFactor::RingOperator, 
            OrderOperatorByKeyCustom < MatrixToFactor::OrderOperatorForColumnIndices >, // order operator for row entries
            MatrixToFactor::OrderOperatorForRowIndices, // order operator for column indices            
            OrderOperatorByKeyCustom< MatrixToFactor::OrderOperatorForRowIndices >, // order operator for column entries            
            MatrixToFactor::OrderOperatorForColumnIndices, // order operator for column indices
        >
    {
        MatrixAlgebraPacket{ 
            matrix: self.generalized_matching_matrix_ref(), 
            ring_operator: self.ring_operator(), 
            order_operator_for_row_entries:     OrderOperatorByKeyCustom::< MatrixToFactor::OrderOperatorForColumnIndices >::new(  // note: we have to use this instead of `matrix_to_factor_ref().order_operator_for_row_entries()` because the order operator for row entries is specific to the type of row entries in the matrix
                                                    self.matrix_to_factor.order_operator_for_column_indices() 
                                                ),
            order_operator_for_row_indices:     self.matrix_to_factor.order_operator_for_row_indices(),
            order_operator_for_column_entries:  OrderOperatorByKeyCustom::< MatrixToFactor::OrderOperatorForRowIndices >::new( 
                                                    self.matrix_to_factor.order_operator_for_row_indices()  // note: we have to use this instead of `matrix_to_factor_ref().order_operator_for_column_entries()` because the order operator for column entries is specific to the type of column entries in the matrix
                                                ),            
            order_operator_for_column_indices:  self.matrix_to_factor.order_operator_for_column_indices(),
        }
    }  

             

    /// The column submatrix of the factored matrix indexed by matched column indices.
    pub fn matrix_to_factor_matched_columns_only( &self ) -> OnlyColumnIndicesInsideCollection< &MatrixToFactor, &HashMap< MatrixToFactor::ColumnIndex, usize >, >
    {
        OnlyColumnIndicesInsideCollection::new( & self.matrix_to_factor, self.matching.bijection_column_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal() )
    }     

    /// The column submatrix of the factored matrix indexed by unmatched column indices.
    pub fn matrix_to_factor_matchless_columns_only( &self ) -> OnlyColumnIndicesOutsideCollection< &MatrixToFactor, &HashMap< MatrixToFactor::ColumnIndex, usize >, >
    {
        OnlyColumnIndicesOutsideCollection::new( & self.matrix_to_factor, self.matching.bijection_column_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal() )
    }    

    /// The row submatrix of the factored matrix indexed by matched row indices.
    pub fn matrix_to_factor_matched_rows_only( &self ) -> OnlyRowIndicesInsideCollection< &MatrixToFactor, &HashMap< MatrixToFactor::RowIndex, usize >, >
        // where 'a: 'b,
    {
        OnlyRowIndicesInsideCollection::new( & self.matrix_to_factor, self.matching.bijection_row_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal() )
    }    

    /// The row submatrix of the factored matrix indexed by unmatched row indices.
    pub fn matrix_to_factor_matchless_rows_only( &self ) -> OnlyRowIndicesOutsideCollection< &MatrixToFactor, &HashMap< MatrixToFactor::RowIndex, usize >, >
    {
        OnlyRowIndicesOutsideCollection::new( & self.matrix_to_factor, self.matching.bijection_row_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal() )
    }  

 


    /// The square submatrix of the factored matrix indexed by matched rows and columns.
    /// 
    /// This matrix is indexed by the same type of row and column indices as the factored matrix. We simply exclude
    /// entries indexed by row and column indices that are unmatched.
    pub fn matched_block_of_matrix_to_factor( &self ) -> 
        OnlyRowIndicesInsideCollection< 
            OnlyColumnIndicesInsideCollection<
                &MatrixToFactor,
                &HashMap< MatrixToFactor::ColumnIndex, usize >, 
            >, 
            &HashMap< MatrixToFactor::RowIndex, usize >, 
        >
    {
        let matched_row_collection      =   self.matching.bijection_row_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal();
        let matched_column_collection   =   self.matching.bijection_column_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal();        
        OnlyRowIndicesInsideCollection::new(
            OnlyColumnIndicesInsideCollection::new(
                & self.matrix_to_factor,
                matched_column_collection,
            ),
            matched_row_collection,
        )
    }




    /// Returns a reference to the internally stored compressed representation of the inverse of the target COMB;
    /// this representation consists of a `VecOfVec` which encodes the off-diagonal entries of the square submatrix of the 
    /// inverse of the target COMB which is indexed by matched (i.e. pivot) indices.
    /// 
    /// Concretely, the represented matrix is obtained from the inverse of the target COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the factored matrix, (ii) deleting the diagonal elements
    /// (each of which is equal to 1), and (iii)
    /// replacing the index `r_i` with `i` for all `i`, where `r_0 < .. < r_k` is the sequence of matched row indices
    /// in sorted order, according to the user-provided [order operator](crate::utilities::order).
    pub fn matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal_ref( &self ) 
        -> 
        & MatrixBimajorData<
                VecOfVec< usize, MatrixToFactor::Coefficient >,
                VecOfVec< usize, MatrixToFactor::Coefficient >,                
            >
        { & self.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal }   

        
    /// Returns a nested double reference to the internally stored compressed representation of the inverse of the target COMB;
    /// this representation consists of a `VecOfVec` which encodes the off-diagonal entries of the square submatrix of the 
    /// inverse of the target COMB which is indexed by matched (i.e. pivot) indices.
    /// 
    /// Concretely, the represented matrix is obtained from the inverse of the target COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the factored matrix, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).    
    /// 
    /// # Design note
    /// 
    /// This function exists because many parts of the `OAT` library use *references* to objects that implement
    /// matrix oracle traits.  A `VecOfVec` simple does not implement oracle traits in general, but a reference
    /// `& VecOfVec` does.  Therefore we often need to work with objects of form `&'b &'b VecOfVec`.  In
    /// practice, we find that Rust is prone to inferring the wrong lifetime if we simply write 
    /// `& self.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal_ref()` (for example, one finds errors alluding to
    /// dropped temprorary values).  This function has succeeded in sidestepping such errors in the past; please 
    /// let us know if it fails to do so successfully in future examples.
    


    /// The square, block submatrix of the inverse of the target COMB which is indexed by matched row indices.
    /// 
    /// Concretely, this object is obtained from the inverse of the target COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the factored matrix, and (ii) 
    /// replacing the index `r_i` with `i` for all `i`, where `r_0 < .. < r_k` is the sequence of matched row indices
    /// in sorted order, according to the user-provided [order operator](crate::utilities::order).
    pub fn matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row< 'a >( &'a self ) 
        -> 
        SumOfScalarAndStrictlyUpperTriangularMatrices<
                &'a  MatrixBimajorData<
                            VecOfVec< usize, MatrixToFactor::Coefficient >,
                            VecOfVec< usize, MatrixToFactor::Coefficient >,                
                        >
            >
        {   

            let prepended //: SumOfScalarAndStrictlyUpperTriangularMatrices< &'a VecOfVec<usize, MatrixToFactor::Coefficient> >
            = SumOfScalarAndStrictlyUpperTriangularMatrices::new( 
                        self.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal_ref(), 
                        MatrixToFactor::RingOperator::one() 
                    );  
            
            prepended
        }  
        
    /// The square, block submatrix of the inverse of the target COMB which is indexed by matched row indices.
    /// 
    /// Concretely, this object is obtained from the inverse of the target COMB by restricting to the square
    /// submatrix indexed by the pivot-row indices of the factored matrix.
    /// 
    /// This submatrix `S` is upper triangular in the sense that `S[i,j] != 0` implies `i <= j`, where order
    /// is determined by [order operator](crate::utilities::order) on row indices which is returned by the
    /// factored matrix using the `self.order_operator_for_row_indices()` method in the [MatrixAlgebra] trait.
    pub fn matched_block_of_target_comb_inverse( &self ) 
        -> 
        MatrixAlgebraPacket<
            ReindexSquareMatrix< 
                SumOfScalarAndStrictlyUpperTriangularMatrices<
                        & MatrixBimajorData<
                                VecOfVec< usize, MatrixToFactor::Coefficient >,
                                VecOfVec< usize, MatrixToFactor::Coefficient >,                
                            >
                    >,                
                &Vec< MatrixToFactor::RowIndex >,
                &HashMap< MatrixToFactor::RowIndex, usize >,
                usize,
                MatrixToFactor::RowIndex,
                MatrixToFactor::ColumnEntry,
            >,
            MatrixToFactor::RingOperator,
            MatrixToFactor::OrderOperatorForRowEntries,
            MatrixToFactor::OrderOperatorForRowIndices,
            MatrixToFactor::OrderOperatorForRowEntries,
            MatrixToFactor::OrderOperatorForRowIndices,                  
        >
    {   

        let matrix_integer_indexed = self.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row();
        let matrix =    ReindexSquareMatrix::new(
                            matrix_integer_indexed,
                            self.matching.bijection_row_indices_to_ordinals_and_inverse().vec_elements_in_order(),
                            self.matching.bijection_row_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal(),                    
                        );
        MatrixAlgebraPacket { 
            matrix, 
            ring_operator:                      self.ring_operator(),
            order_operator_for_row_entries:     self.order_operator_for_row_entries(),
            order_operator_for_row_indices:     self.order_operator_for_row_indices(),
            order_operator_for_column_entries:  self.order_operator_for_row_entries(),
            order_operator_for_column_indices:  self.order_operator_for_row_indices(),
        }
    }      


    /// The matched part of the inverse target COMB, with rows indxed by rank order and columns indexed by native row indices.
    /// 
    /// If the matched row indices are `rho_0 < .. < rho_k` (arranged in sorted order according to `MatrixToFactor::OrderOpeartorForRowIndices`),
    /// then this matrix, `M` has rows indexed by `0 .. k` and columns indexed by `rho_0 .. rho_k`. Moreover, `M[i,\rho_j] = T[\rho_i, \rho_j]`,
    /// where `T` is the target COMB.
    /// 
    /// # Comment on notation
    /// 
    /// The suffix `_or` at the end of this method name refers to the fact that rows are indexed by (o)rdinal, while
    /// columns are indxed by (r)ow indices of the matrix to be factored..
    pub fn matched_block_of_target_comb_inverse_or( &self ) ->
        MatrixAlgebraPacket<
            ReindexMatrixColumns< 
                SumOfScalarAndStrictlyUpperTriangularMatrices<
                    & MatrixBimajorData<
                        VecOfVec< usize, MatrixToFactor::Coefficient >,
                        VecOfVec< usize, MatrixToFactor::Coefficient >,                
                    >
                >,                
                &Vec< MatrixToFactor::RowIndex >,
                &HashMap< MatrixToFactor::RowIndex, usize >,
                usize,
                MatrixToFactor::RowIndex,
                MatrixToFactor::ColumnEntry,
            >,
            MatrixToFactor::RingOperator,
            MatrixToFactor::OrderOperatorForColumnEntries,      // order operator for row entries (which look like column entries of the factored matrix)
            OrderOperatorAuto,                                  // order opeartor for row indices (which are usize)
            OrderOperatorByKey,                                 // order operator for column entries
            MatrixToFactor::OrderOperatorForRowIndices,         // order operator for column indices (which are row indices of the factored matrix)   
        >
    {
        let matrix_integer_indexed = self.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row();
        let matrix =    ReindexMatrixColumns::new(
                            matrix_integer_indexed,
                            self.matching.bijection_row_indices_to_ordinals_and_inverse().vec_elements_in_order(),
                            self.matching.bijection_row_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal(),                    
                        );
        MatrixAlgebraPacket { 
            matrix, 
            ring_operator:                      self.ring_operator(),
            order_operator_for_row_entries:     self.order_operator_for_column_entries(),
            order_operator_for_row_indices:     OrderOperatorAuto,                          
            order_operator_for_column_entries:  OrderOperatorByKey,
            order_operator_for_column_indices:  self.order_operator_for_row_indices(),      
        }                        
    }


    /// Product of the matched blocks of the target COMB and the factored matrix.
    /// 
    /// Let `M` be the matrix returned by this function, and let `A = Tinv D` be the product of the inverse of the target COMB with the factored matrix `D`. If
    /// 
    /// - the matched row indices are `r0 < .. < rN` (arranged in sorted order according to `MatrixToFactor::OrderOpeartorForRowIndices`)
    /// 
    /// - the matched column indices are `k0 < .. < kN` (arranged in sorted order according to `MatrixToFactor::OrderOpeartorForColumnIndices`;
    /// note that this ordering *does not imply that `rho_i` matches with `kappa_i`*)
    /// 
    /// then `M[ i, kj] = A[ ri, kj]`
    /// 
    /// # Naming convention
    /// 
    /// The suffix `oc` at the end of this function name refers to the fact that rows are indxed by (o)rdinal, while columns
    /// are indexed by the (c)olumn indices of the matrix to be factored.
    pub fn matched_blocks_of_target_comb_inverse_times_matrix_to_factor_oc( & self ) ->
        ProductMatrix<
            MatrixAlgebraPacket< // this packet contains the matched block of the target comb:
                ReindexMatrixColumns< 
                    SumOfScalarAndStrictlyUpperTriangularMatrices<
                        & MatrixBimajorData<
                            VecOfVec< usize, MatrixToFactor::Coefficient >,
                            VecOfVec< usize, MatrixToFactor::Coefficient >,                
                        >
                    >,                
                    &Vec< MatrixToFactor::RowIndex >,                   // mapping from old index to new index
                    &HashMap< MatrixToFactor::RowIndex, usize >,        // mapping from new index to old index
                    usize,                                              // old index
                    MatrixToFactor::RowIndex,                           // new index
                    MatrixToFactor::ColumnEntry,                        // new entry
                >,
                MatrixToFactor::RingOperator,
                MatrixToFactor::OrderOperatorForColumnEntries,          // order operator for row entries 
                OrderOperatorAuto,                                      // order opeartor for row indices (which are usize)
                OrderOperatorByKey,                                     // order operator for column entries
                MatrixToFactor::OrderOperatorForRowIndices,             // order operator for column indices (which are row indices of the factored matrix)   
            >,            
            OnlyRowIndicesInsideCollection<                             // the submatrix is the matched part of the factored matrix
                OnlyColumnIndicesInsideCollection<
                    &MatrixToFactor,
                    &HashMap< MatrixToFactor::ColumnIndex, usize >, 
                >, 
                &HashMap< MatrixToFactor::RowIndex, usize >, 
            >            
        >
    {
        let target_block        =   self.matched_block_of_target_comb_inverse_or();
        let factored_block      =   self.matched_block_of_matrix_to_factor();
        ProductMatrix::new( target_block, factored_block )
    }    


    /// Returns a matrix with rows indexed by `MatrixToFactor::RowIndex` and columns indexed by `MatrixToFactor::ColumnIndex`
    pub fn target_comb_inverse_times_matrix_to_factor_matched_block( &self ) 
        -> 
        TargetCombInverseTimesMatrixToFactorMatchedBlock< '_, MatrixToFactor > 
    {
        TargetCombInverseTimesMatrixToFactorMatchedBlock{ umatch: self } // SourceCombInverseMatchedBlock::new(self )
    }

 

    pub fn target_comb_inverse_times_matrix_to_factor_matched_block_with_rows_indexed_by_matched_column_index( &self ) 
        -> 
        TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex< '_, MatrixToFactor > 
    {
        TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex{ umatch: self } // SourceCombInverseMatchedBlock::new(self )
    }

    /// Solve `Tx = b`, where `T` is the target COMB.
    /// 
    /// Solution returns entries in strictly descending order.
    /// 
    /// `b` must iterate over entries in strictly descending order.
    /// 
    /// # Example
    /// 
    /// In this case the target COMB is a `2x2` identity matrix.
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::matrices::operations::{umatch::row_major::Umatch, MatrixOracleOperations};
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// 
    /// // define inputs
    /// let b               =   vec![ (1, 6.), (0, 6.), ];  // note: entries bust appear in descending order of index
    /// let data            =   vec![   
    ///                             vec![ (0, 1.),   (1, 2.)  ],
    ///                             vec![            (1, 1.)  ],
    ///                         ];
    /// let data            =   VecOfVec::new( data ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &data );
    /// let umatch          =   Umatch::new( 
    ///                             &matrix, 
    ///                             (0..2).rev() 
    ///                         );               
    /// 
    /// // solve
    /// let x               =   umatch
    ///                             .solve_tx_equals_b( b.clone() )             // solve the equation `Tx = b`
    ///                             .ok().unwrap()                              // unwrap the solution from its enclosing Result struct
    ///                             .collect::< Vec<(usize, f64)> >();          // collect the sparse vector entries into a Rust Vec struct
    /// 
    /// // calculat Tx
    /// let tx              =   umatch
    ///                             .target_comb()
    ///                             .multiply_with_column_vector_reverse( x )   // adding _reverse ensures that entries of Tx are returned in reverse order of index
    ///                             .collect::< Vec<(usize, f64)> >();          // collect the sparse vector entries into a Rust Vec struct
    ///
    /// // verify Tx = b
    /// assert!( tx.eq(  &b  ) );
    /// ```
    pub fn solve_tx_equals_b< I >( &self, b: I ) 
        -> 
        // Vec< MatrixToFactor::ColumnEntry >
        Result< 
            TriangularSolveForColumnVectorReverse<
                Vec< MatrixToFactor::ColumnEntry >, 
                TargetComb< MatrixToFactor >,
            >,
            ()
        >         
    where
        I:      IntoIterator<Item=MatrixToFactor::ColumnEntry>,
    {
        TriangularSolveForColumnVectorReverse::solve( 
            b, 
            self.target_comb(), 
        )
    }




    /// Returns `x = b * Sinv`, where `Sinv` is the inverse of the source COMB.
    /// 
    /// Solution returns entries in strictly ascending order.
    /// 
    /// The entries of `b` do not have to appear in ascending order.
    /// 
    /// This product is computed by simplying multiplying `b * Sinv`; in the past other methods were used, see below.
    /// 
    /// # Design notes
    /// 
    /// In other situations it might have made sense to compute this vector by solving `x * S = b` via
    /// back substitution. However, the rows of `Sinv` are easy to look up, compared to `S`
    /// (because up to some basic transformations `Sinv` can be recovered efficiently from `Tinv *D`, where
    /// `Tinv` is the inverse of the target COMB and `D` is the matrix we want to factor.). So we
    /// posited that simply multiplying `x * Sinv` would be more efficient. **Future users may wish to compare these two techniques experimentally**.
    pub fn solve_x_equals_b_times_source_comb_inverse< I >( &self, b: I ) 
        -> 
        // Vec< MatrixToFactor::ColumnEntry >
        // TriangularSolveForRowVector<
        //         I, 
        //         SourceComb< MatrixToFactor >, 
        //     >
        LinearCombinationOfRows< SourceCombInverse< MatrixToFactor > >
    where
        I:      IntoIterator<Item=MatrixToFactor::RowEntry>,  // NB: we actually want the type to be *ROW* entries, because row entries are the onces with *COLUMN INDICES*
    {
        // TriangularSolveForRowVector::solve( 
        //     b, 
        //     self.source_comb(), 
        // )
        b.multiply_self_as_a_row_vector_with_matrix( self.source_comb_inverse() )
    }



    /// Return `x = bS` by solving `x * Sinv = b` where `S` is the source COMB and `Sinv` is its inverse.
    /// 
    /// Solution returns entries in strictly ascending order.
    /// 
    /// `b` must iterate over entries in strictly ascending order.
    /// 
    /// **Note** It is probably more efficient than computing `bS` directly, because
    /// the rows of `Sinv` are easier to compute than the rows of `S`.
    /// 
    /// # Errors
    /// 
    /// Returns `Err(())` if the entries of `b` are not in strictly ascending order.
    pub fn solve_x_equals_b_times_source_comb< I >( &self, b: I ) 
        -> 
        // Vec< MatrixToFactor::ColumnEntry >
        Result<
            TriangularSolveForRowVector<
                Vec< MatrixToFactor::RowEntry >, 
                SourceCombInverse< MatrixToFactor >, 
            >,
            ()
        >
    where
        I:      IntoIterator<Item=MatrixToFactor::RowEntry>,                                                            
    {
        TriangularSolveForRowVector::solve( 
            b, 
            self.source_comb_inverse(), 
        )
    }




    /// Solve `Dx = b`, where `D` is the factored matrix.
    /// 
    /// Returns `None` if there is no solution.
    /// 
    /// # Arguments
    /// 
    /// The iterable `b` can iterate over entries in any order. 
    /// If `b` contains multiple entries with the same index, these entries will be summed.
    /// 
    /// # Calculation
    /// 
    /// The U-match equation implies `TM = DS` implies that `TMS^{-1} = D`. It can therefore be shown that
    /// `SM^{^-1}T^{-1}` is a generalized inverse of `D`, where `M^{-1}` is the generalized
    /// inverse of `M` obtained by transposing and inverting nonzero entries. Therefore the solution `x`
    /// can be computed as `SM^{^-1}T^{-1}b`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::operations::umatch::row_major::Umatch;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;   
    /// use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket; 
    /// use itertools::Itertools;    
    /// 
    /// // DEFINE THE MATRIX
    /// // ===============================
    /// let matrix          =   VecOfVec::new( 
    ///                                         vec![   
    ///                                                     vec![(0,true), (1,true), (2,true)],
    ///                                                     vec![                            ], 
    ///                                                     vec![                    (2,true)], 
    ///                                         ] 
    ///                                     ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( &matrix );
    ///                                 
    /// // COMPUTE U-MATCH
    /// // ===============================
    ///                                 
    /// let umatch
    ///     =   Umatch::new(
    ///             & matrix,  // the matrix we wish to factor
    ///             (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
    ///         );        
    /// 
    /// // SOLVE Dx = b FOR x
    /// // ===============================
    /// 
    /// let b   =   [ (0,true), (2,true) ]; 
    /// let x   =   umatch.solve_dx_equals_b( b.clone() ).unwrap();
    /// let dx  =   umatch.multiply_dx(x);
    /// assert!( dx.eq( b ) );  
    /// 
    /// // SOLVE Dx = b FOR x (WHEN NO SOLUTION EXISTS)
    /// // ===============================
    /// 
    /// let b   =   [ (1,true) ]; 
    /// assert!( umatch.solve_dx_equals_b( b ).is_none() ); // no solution exists
    /// ```
    pub fn solve_dx_equals_b< Vector >( &self, b: Vector ) 
        -> 
        Option< Vec< MatrixToFactor::RowEntry > >
    where
        Vector:     IntoIterator<Item=MatrixToFactor::ColumnEntry>,                            
    {
        let matching_inverse = self.generalized_matching_matrix_ref().generalized_inverse(self.ring_operator());
        let matching_inverse = MatrixAlgebraPacket{
            matrix: matching_inverse,
            ring_operator: self.ring_operator(),
            order_operator_for_row_entries:  OrderOperatorByKeyCustom::new(self.order_operator_for_row_indices()), // self.order_operator_for_column_entries(),//
            order_operator_for_row_indices: self.order_operator_for_column_indices(),
            order_operator_for_column_entries:  OrderOperatorByKeyCustom::new(self.order_operator_for_column_indices()), // self.order_operator_for_row_entries(),//
            order_operator_for_column_indices: self.order_operator_for_row_indices(),
        };
        let comb_target_inverse = self.target_comb_inverse();
        let comb_source = self.source_comb();

        // multiply with T^{-1}
        let tinv_b = comb_target_inverse.multiply_with_column_vector( b ).collect_vec();
        
        // check that all nonzeros occur in matched indices
        for entry in tinv_b.iter() {
            if matching_inverse.matrix_ref().lacks_a_match_for_column_index( & entry.key() ) {
                // if the matching matrix does not have a match for this column index, then we cannot solve
                // the equation `Dx = b` for this `b`
                return None;
            }
        }
        // multiply with M^{-1}
        let minv_tinv_b = matching_inverse.multiply_with_column_vector(tinv_b).collect_vec();

        // multiply with S^{-1}
        let x = comb_source.multiply_with_column_vector( minv_tinv_b ).collect_vec();
        
        return Some( x )
    }




    /// Solve `xD = b`, where `D` is the factored matrix.
    /// 
    /// # Arguments
    /// 
    /// `b` must iterate over entries in strictly ascending order.
    /// 
    /// # Returns
    /// 
    ///  `None` if there is no solution, otherwise `Some( iter )`, where `iter` is an iterator that runs over the entries of a solution in strictly ascending order
    /// 
    /// # Solution strategry
    /// 
    /// 
    /// We have the U-match factorization `TM = DS`, where `T` is the target COMB, `M` is the factored matrix, and `S` is the source COMB.
    /// This equation implies that `MS^{-1} = T^{-1}D`.
    /// Because `S^{-1}` is upper triangular, the matrix `M Sinv` is in a variation of row echelon form. 
    /// So too is `T^{-1}D = MS^{-1}`. Thus we can solve `yT^{-1}D = b` efficiently using back-substitution (or verify that no solution exists).
    /// We then set `x = yT^{-1}`.
    /// 
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket;
    /// use oat_rust::algebra::matrices::operations::umatch::row_major::Umatch;
    /// 
    /// // Define the matrix
    /// // -----------------
    /// let d = VecOfVec::new(
    ///             vec![
    ///                 vec![  (0,true), (1,true),           ],
    ///                 vec![            (1,true), (2,true), ],
    ///             ]
    ///         ).ok().unwrap();
    /// let d = MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( &d );
    /// 
    /// // Obtain a u-match factorization
    /// // ------------------------------
    /// let umatch  =   Umatch::new( 
    ///     &d, 
    ///     (0..2).rev(),             
    /// );
    /// 
    /// // Solve xD = b for x
    /// // ------------------
    /// 
    /// // Case 1: a solution exists; in this case we are gauaranteed to find one
    /// let x = umatch.solve_xd_equals_b( vec![ (0,true), (2,true), ] );        
    /// assert!( x.is_some() );
    /// assert!( x.unwrap().eq( & vec![ (0,true), (1,true), ] ) );
    /// 
    /// // Case 2: no solution exists; in this case we get a certificate that no solution exists
    /// let x = umatch.solve_xd_equals_b( vec![ (0,true), (1,true), (2,true) ] );        
    /// assert!( x.is_none() );
    /// ```
    pub fn solve_xd_equals_b< Vector >( &self, b: Vector ) 
        -> 
        Option< Vec< MatrixToFactor::ColumnEntry > >
        where 
            Vector:                     IntoIterator<Item=MatrixToFactor::RowEntry>,                             
    {

        let prod        =   self.target_comb_inverse()
                            .multiply_on_the_left_of( self.matrix_to_factor_ref() );
        let matching = self.generalized_matching_matrix_ref();
        let key_matching = |x| { matching.row_index_for_column_index(&x) };
        // let key_matching    =   self.generalized_matching_matrix_ref().bijection_column_indices_to_row_indices();
        RowEchelonSolver::solve(
                    b.into_iter(),
                    prod,
                    EvaluateFunctionFnMutWrapper::new( key_matching ),
                    self.ring_operator(),
                    self.order_operator_for_row_entries(),
                )
                .solution() // returns Some(solution) if a solution exists, otherwise None
                .map(|x| x.multiply_self_as_a_row_vector_with_matrix( self.target_comb_inverse() ).collect_vec() )
    }




    
    
    /// A basis for the kernel of the factored matrix ("factored matrix" is a term from U-match factorization)
    /// 
    /// # Inputs
    /// 
    /// `column_indices` is an iterator that runs over the columns of the factored matrix, without repeats
    /// 
    /// # Returns
    /// 
    /// An iterator that wraps around `column_indices`.  The wrapper filters out every element
    /// of `column_indices` that indexes a nonzero column of the generalized matching matrix; every other
    /// element is mapped to the corresponding column of the source COMB.
    /// 
    /// # Why this works
    /// 
    /// Every vector returned by this iterator lies in the kernel of the factored matrix, by
    /// the identity `TM = DS`.  The set of iterators is linearly independent because they form
    /// a subset of the columns of an invertible matrix.  The number of vectors equals the nullity
    /// of the matrix, again by the identity `TM = DS`.  Therefore the iterator runs over a basis for the kernel.
    pub fn kernel< ColumnIndices >( &self, column_indices: ColumnIndices ) 
        -> 
        SequenceOfReverseColumns<
                SourceComb
                    < MatrixToFactor >,
                FilterOutMembers
                    < ColumnIndices::IntoIter, & HashMap< MatrixToFactor::ColumnIndex, usize > >,
            >
    where   
        ColumnIndices:              IntoIterator< Item = MatrixToFactor::ColumnIndex >,          

    {
        SequenceOfReverseColumns::new(
            self.source_comb(),  
            self.matching.filter_out_matched_column_indices( column_indices )
        )
    }  



    /// A basis for the kernel of the factored matrix
    /// 
    /// 
    /// # Returns
    /// 
    /// An iterator that returns a subset of the columns of the target COMB, which form a basis for the image.
    /// 
    /// Specifically, the iterator returns `T[:,i]` for all `i` such that `M[i,:]` is nonzero, where `T`
    /// is the target COMB and `M` is the generalized matching matrix.
    /// 
    /// # Why this works
    /// 
    /// Every vector returned by this iterator lies in the image of the factored matrix, by
    /// the identity `TM = DS`.  The set of iterators in linearly independent because they form
    /// a subset of the columns of an invertible matrix.  The number of vectors equals the rank
    /// of the matrix, again by the identity `TM = DS`.  Therefore the iterator runs over a basis for the image.
    pub fn image( &self ) 
        -> 
        SequenceOfReverseColumns<
                TargetComb
                    < MatrixToFactor >,
                Cloned< std::slice::Iter< MatrixToFactor::RowIndex > >,
            >          
    {
        SequenceOfReverseColumns::new(
            self.target_comb(),  
            self.matching.bijection_row_indices_to_ordinals_and_inverse().vec_elements_in_order().iter().cloned()
        )
    }      



    //  MATRIX-VECTOR MULTIPLICATION
    //  ---------------------------------------------------------------------------------------------------------


    /// Calculate the product `vD`, where `D` is the factored matrix and `v` is a vector.
    pub fn multiply_xd< Vector, VectorEntry >( & self, v: Vector ) 
        -> 
        LinearCombinationOfRows< MatrixToFactor >
        where 
            Vector:                         IntoIterator<Item=VectorEntry>,
            VectorEntry:                    KeyValGet < Key = MatrixToFactor::RowIndex, Val = MatrixToFactor::Coefficient >,        
                             
    {
        let matrix = |i| self.matrix_to_factor.row( & i );
        v.multiply_matrix_fnmut( matrix, self.ring_operator(), self.order_operator_for_row_entries() )
    }    

    /// Calculate the product `Dv`, where `D` is the factored matrix and `v` is a vector.
    pub fn multiply_dx< Vector, VectorEntry >( &self, v: Vector ) 
        -> 
        LinearCombinationOfColumns< MatrixToFactor >
        where 
            Vector:                         IntoIterator<Item=VectorEntry>,
            VectorEntry:                    KeyValGet < Key = MatrixToFactor::ColumnIndex, Val = MatrixToFactor::Coefficient >,                             
    {
        let matrix = |i| self.matrix_to_factor.column( & i );
        v.multiply_matrix_fnmut( matrix, self.ring_operator(), self.order_operator_for_column_entries() )
    }



    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //
    //  RETURN WHEN POSSIBLE -- THE PROBLEM IS LINE 2254
    //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //
    // /// Calculate the product `xC`, where `S` is the source COMB
    // pub fn multiply_vs< Vector, VectorEntry >( & self, v: Vector ) 
    //     -> 
    //     Simplify<
    //             IteratorsMergedInSortedOrder<
    //                     Scale< 
    //                             SourceCombRow<
    //                                     MatrixToFactor, 
    //                                     MatrixToFactor::RingOperator, 
    //                                     MatrixToFactor::OrderOperatorForRowEntries, 
    //                                 >,
    //                             MatrixToFactor::ColumnIndex, 
    //                             MatrixToFactor::RingOperator, 
    //                             MatrixToFactor::Coefficient, 
    //                         >,
    //                     MatrixToFactor::OrderOperatorForRowEntries,
    //                 >,
    //             MatrixToFactor::ColumnIndex,
    //             MatrixToFactor::RingOperator,
    //             MatrixToFactor::Coefficient,
    //         >
    //     where 
    //         Vector:                         IntoIterator<Item=VectorEntry>,
    //         VectorEntry:                    KeyValGet < Key = MatrixToFactor::RowIndex, Val = MatrixToFactor::Coefficient >,                            
    // {
    //     let comb_source = self.source_comb();        
    //     let matrix = |i| comb_source.row(i);
    //     return v.multiply_matrix( matrix, self.ring_operator(), self.order_operator_for_row_entries() )
    // }    



    // /// Calculate the product `Sv`, where `S` is the source COMB and `v` is a vector.
    // pub fn multiply_sv< 'a, Vector, VectorEntry >( &'a self, v: Vector ) 
    //     -> 
    //     LinearCombinationOfColumnsReverse< SourceComb<'a, _> >
 
    //     where 
    //         Vector:         IntoIterator<Item=VectorEntry>,
    //         VectorEntry:    KeyValGet < Key = MatrixToFactor::ColumnIndex, Val = MatrixToFactor::Coefficient >,                          
    // {
    //     let comb_source = self.source_comb();
    //     let matrix = |i| comb_source.column_reverse(i);
    //     v.multiply_matrix_fnmut( matrix, self.ring_operator(), self.order_operator_for_row_entries_reverse() )
    // }


    // /// Calculate the product `vS`, where `S` is the source COMB and `v` is a vector.
    // pub fn multiply_vs< 'a, Vector, VectorEntry >( &'a self, v: Vector ) 
    //     -> 
    //     LinearCombinationOfRows< SourceComb< 'a, _> >
    //     where 
    //         Vector:         IntoIterator<Item=VectorEntry>,
    //         VectorEntry:    KeyValGet < Key = MatrixToFactor::ColumnIndex, Val = MatrixToFactor::Coefficient >,                          
    // {
    //     let comb_source = self.source_comb();
    //     let matrix = |i| comb_source.row(i);
    //     v.multiply_matrix_fnmut( matrix, self.ring_operator(), self.order_operator_for_row_entries_reverse() )
    // }



}


//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- IMPLEMENTATIONS FOR ColumnIndex == RowIndex
//  ---------------------------------------------------------------------------------------------------------


impl < MatrixToFactor, EntryForRowsAndColumns, IndexForRowsAndColumns >  

    Umatch 
    < MatrixToFactor >  
    
    where   
        MatrixToFactor:                         MatrixAlgebra + MatrixOracle< RowEntry=EntryForRowsAndColumns, ColumnEntry=EntryForRowsAndColumns, RowIndex=IndexForRowsAndColumns >,    
        MatrixToFactor::ColumnIndex:            Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:               Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct 
        MatrixToFactor::RowEntry:               KeyValGet < Key = MatrixToFactor::ColumnIndex, Val = MatrixToFactor::Coefficient >,
        // OrderOperatorByKey<usize, MatrixToFactor::Coefficient, (usize, MatrixToFactor::Coefficient)>: JudgePartialOrder< (usize, MatrixToFactor::Coefficient)>
{
    
     
}































//  =========================================================================================================
//  TESTS
//  =========================================================================================================


//  ---------------------------------------------------------------------
//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    use crate::algebra::{matrices::types::packet::MatrixAlgebraPacket, rings::types::field_prime_order::BooleanField};



    #[test]
    fn doc_test_umatchrowmajor_comprehensive_tiny() {

        // import packages
        use crate::algebra::matrices::operations::umatch::row_major::doc_test_drafts::BooleanField;
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
        
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::types::product::ProductMatrix;
        use crate::algebra::matrices::query::MatrixOracle;
        
        use itertools::Itertools;

        // define the coefficient ring
        let ring_operator                   =   BooleanField::new();

        // define the matrix we wish to factor
        let num_indices_row           =   1;
        let num_indices_col           =   1;
        let matrix_to_factor_data              =   VecOfVec::new( 
            vec![   
                                vec![(0,true),], 
                        ] ).ok().unwrap();
        let matrix_to_factor = MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( & matrix_to_factor_data );

        // compute the U-match factorization
        let umatch
            =   Umatch::new(
                    matrix_to_factor, 
                    (0..num_indices_row).rev(), 
                );
        
        // extract T, T^{-1}, S, S^{-1}, and M
        let matching = umatch.generalized_matching_matrix_ref();
        let comb_target = umatch.target_comb();
        let comb_target_inv = umatch.target_comb_inverse();        
        let comb_source = umatch.source_comb();        
        let comb_source_inv = umatch.source_comb_inverse(); 

        // get references to T, T^{-1}, S, S^{-1}, and M        
        let comb_target_ref         =   & comb_target;
        let comb_target_inv_ref         =   & comb_target_inv;
        let comb_source_ref         =   & comb_source;
        let comb_source_inv_ref         =   & comb_source_inv;   
        
        // compute some products
        let product_source = ProductMatrix::new( comb_source_ref, comb_source_inv_ref );
        let product_target = ProductMatrix::new( comb_target_ref, comb_target_inv_ref );        
        let product_target_comb_inv_times_matrix_to_factor = ProductMatrix::new( comb_target_inv_ref, matrix_to_factor );      
        let product_target_comb_inv_times_matrix_to_factor_times_source_comb = ProductMatrix::new( product_target_comb_inv_times_matrix_to_factor, comb_source_ref );                


        // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
        for column_index in 0 .. num_indices_col { 
            assert_eq!(
                product_source.row( & column_index ).collect_vec(),
                vec![ (column_index, true) ]
            ) 
        }

        // check that the product of the target COMB with its inverse is identity T * T^{-1} = I
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target.row( &row_index ).collect_vec(),
                vec![ (row_index, true) ]
            ) 
        }    
        
        // check the factorization T^{-1} * D * S = M
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target_comb_inv_times_matrix_to_factor_times_source_comb.row( &row_index ).collect_vec(),
                matching.row( &row_index ).collect_vec()
            ) 
        }     

    }



    #[test]
    fn doc_test_umatchrowmajor_comprehensive_tiny_waist() {

        // import packages
        use crate::algebra::matrices::operations::umatch::row_major::doc_test_drafts::BooleanField;
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
        
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::types::product::ProductMatrix;
        use crate::algebra::matrices::query::MatrixOracle;
        
        use itertools::Itertools;

        // define the coefficient ring
        let ring_operator                   =   BooleanField::new();

        // define the matrix we wish to factor
        let num_indices_row           =   2;
        let num_indices_col           =   1;
        let matrix_to_factor_data              =   VecOfVec::new( 
            vec![   
                                vec![(0,true),], 
                                vec![(0,true),],                                 
                            ] ).ok().unwrap();
        let matrix_to_factor = MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( & matrix_to_factor_data );

        // compute the U-match factorization
        let umatch
            =   Umatch::new(
                    matrix_to_factor, 
                    (0..num_indices_row).rev(), 
                );
        
        // extract T, T^{-1}, S, S^{-1}, and M
        let matching = umatch.generalized_matching_matrix_ref();
        let comb_target = umatch.target_comb();
        let comb_target_inv = umatch.target_comb_inverse();        
        let comb_source = umatch.source_comb();        
        let comb_source_inv = umatch.source_comb_inverse(); 

        // get references to T, T^{-1}, S, S^{-1}, and M        
        let comb_target_ref         =   & comb_target;
        let comb_target_inv_ref         =   & comb_target_inv;
        let comb_source_ref         =   & comb_source;
        let comb_source_inv_ref         =   & comb_source_inv;   
        
        // compute some products
        let product_source = ProductMatrix::new( comb_source_ref, comb_source_inv_ref );
        let product_target = ProductMatrix::new( comb_target_ref, comb_target_inv_ref );        
        let product_target_comb_inv_times_matrix_to_factor = ProductMatrix::new( comb_target_inv_ref, matrix_to_factor );      
        let product_target_comb_inv_times_matrix_to_factor_times_source_comb = ProductMatrix::new( product_target_comb_inv_times_matrix_to_factor, comb_source_ref );                


        // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
        for column_index in 0 .. num_indices_col { 
            assert_eq!(
                product_source.row( & column_index ).collect_vec(),
                vec![ (column_index, true) ]
            ) 
        }

        // check that the product of the target COMB with its inverse is identity T * T^{-1} = I
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target.row( & row_index ).collect_vec(),
                vec![ (row_index, true) ]
            ) 
        }    
        
        // check the factorization T^{-1} * D * S = M
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target_comb_inv_times_matrix_to_factor_times_source_comb.row( &row_index ).collect_vec(),
                matching.row( & row_index ).collect_vec()
            ) 
        }     

    }    





    #[test]
    fn doc_test_umatchrowmajor_comprehensive_tiny_height() {

        // import packages
        use crate::algebra::matrices::operations::umatch::row_major::doc_test_drafts::BooleanField;
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
        
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::types::product::ProductMatrix;
        use crate::algebra::matrices::query::MatrixOracle;
        
        use itertools::Itertools;

        // define the coefficient ring
        let ring_operator                   =   BooleanField::new();

        // define the matrix we wish to factor
        let num_indices_row           =   1;
        let num_indices_col           =   2;
        let matrix_to_factor_data              =   VecOfVec::new( 
            vec![   
                                vec![(0,true), (1,true)], 
                            ] ).ok().unwrap();
        let matrix_to_factor                 =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( & matrix_to_factor_data );

        // compute the U-match factorization
        let umatch
            =   Umatch::new(
                    matrix_to_factor, 
                    (0..num_indices_row).rev(), 
                );
        
        // extract T, T^{-1}, S, S^{-1}, and M
        let matching = umatch.generalized_matching_matrix_ref();
        let comb_target = umatch.target_comb();
        let comb_target_inv = umatch.target_comb_inverse();        
        let comb_source = umatch.source_comb();        
        let comb_source_inv = umatch.source_comb_inverse(); 

        // get references to T, T^{-1}, S, S^{-1}, and M        
        let comb_target_ref         =   & comb_target;
        let comb_target_inv_ref         =   & comb_target_inv;
        let comb_source_ref         =   & comb_source;
        let comb_source_inv_ref         =   & comb_source_inv;   
        
        // compute some products
        let product_source = ProductMatrix::new( comb_source_ref, comb_source_inv_ref );
        let product_target = ProductMatrix::new( comb_target_ref, comb_target_inv_ref );        
        let product_target_comb_inv_times_matrix_to_factor = ProductMatrix::new( comb_target_inv_ref, matrix_to_factor );      
        let product_target_comb_inv_times_matrix_to_factor_times_source_comb = ProductMatrix::new( product_target_comb_inv_times_matrix_to_factor, comb_source_ref );                


        // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
        for column_index in 0 .. num_indices_col { 
            assert_eq!(
                product_source.row( & column_index ).collect_vec(),
                vec![ (column_index, true) ]
            ) 
        }

        // check that the product of the target COMB with its inverse is identity T * T^{-1} = I
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target.row( &row_index ).collect_vec(),
                vec![ (row_index, true) ]
            ) 
        }    
        
        // check the factorization T^{-1} * D * S = M
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target_comb_inv_times_matrix_to_factor_times_source_comb.row( &row_index ).collect_vec(),
                matching.row( &row_index ).collect_vec()
            ) 
        }     

    }





    #[test]
    fn doc_test_umatchrowmajor_comprehensive_small() {

        // import packages
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
        
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::types::product::ProductMatrix;
        use crate::algebra::matrices::query::MatrixOracle;
        
        use itertools::Itertools;

        // define the coefficient ring
        let modulus                     =   5;
        let ring_operator                   =   PrimeOrderField::new( modulus );        

        // define the matrix we wish to factor
        let num_indices_row           =   2;
        let num_indices_col           =   3;
        let matrix_to_factor_data              =   VecOfVec::new( 
            vec![   
                                vec![(0,1), (1,2), (2,3)], 
                                vec![              (2,1)]  
                            ] ).ok().unwrap();
        let matrix_to_factor = MatrixAlgebraPacket::with_default_order( & matrix_to_factor_data, ring_operator.clone() );

        // compute the U-match factorization
        let umatch
            =   Umatch::new(
                    matrix_to_factor, 
                    (0..num_indices_row).rev(), 
                );
        
        // extract T, T^{-1}, S, S^{-1}, and M
        let matching = umatch.generalized_matching_matrix_ref();
        let comb_target = umatch.target_comb();
        let comb_target_inv = umatch.target_comb_inverse();        
        let comb_source = umatch.source_comb();        
        let comb_source_inv = umatch.source_comb_inverse(); 

        // get references to T, T^{-1}, S, S^{-1}, and M        
        let comb_target_ref         =   & comb_target;
        let comb_target_inv_ref         =   & comb_target_inv;
        let comb_source_ref         =   & comb_source;
        let comb_source_inv_ref         =   & comb_source_inv;   
        
        // compute some products
        let product_source = ProductMatrix::new( comb_source_ref, comb_source_inv_ref );
        let product_target = ProductMatrix::new( comb_target_ref, comb_target_inv_ref );        
        let product_target_comb_inv_times_matrix_to_factor = ProductMatrix::new( comb_target_inv_ref, matrix_to_factor );      
        let product_target_comb_inv_times_matrix_to_factor_times_source_comb = ProductMatrix::new( product_target_comb_inv_times_matrix_to_factor, comb_source_ref );                


        // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
        for column_index in 0 .. num_indices_col { 
            assert_eq!(
                product_source.row( & column_index ).collect_vec(),
                vec![ (column_index, 1) ]
            ) 
        }

        // check that the product of the target COMB with its inverse is identity T * T^{-1} = I
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target.row( &row_index ).collect_vec(),
                vec![ (row_index, 1) ]
            ) 
        }    
        
        // check the factorization T^{-1} * D * S = M
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target_comb_inv_times_matrix_to_factor_times_source_comb.row( &row_index ).collect_vec(),
                matching.row( &row_index ).collect_vec()
            ) 
        }     

    }


}    


//  ---------------------------------------------------------------------
//  Unit tests
//  ---------------------------------------------------------------------

#[cfg(test)]
mod unit_tests {
    

    use itertools::{assert_equal, Itertools};

    use crate::algebra::matrices::{debug::{matrix_oracle_is_internally_consistent, matrix_order_operators_are_internally_consistent}, operations::{multiply::multiply_row_vector_with_matrix, MatrixOracleOperations} };
    use crate::algebra::matrices::debug::verify_rows_compatible_with_columns;
    
    use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
     
    use crate::algebra::rings::types::field_prime_order::PrimeOrderField;



//     //  ===================================================================================
//     //  CONSTRUCTION OF PIVOT BLOCK OF INVERSE OF THE CODOMAIN COMB
//     //  ===================================================================================    

    

    /// This test targets the initial computation of the pivot block of the inverse of the target COMB.
    #[test]
    fn test_initial_decomposition() {
        use crate::algebra::matrices::operations::umatch::row_major::get_pivot_block_of_target_comb_inverse_with_deleted_diagonal;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;

        let matrix_to_factor       =   VecOfVec::new(
                                            vec![
                                                vec![ (0, 1), (1, 1) ],
                                                vec![ (0, 1), (1, 1) ],
                                                vec![ (0, 1), (1, 1) ],                                                
                                            ]
                                        ).ok().unwrap();
        let ring_operator       =   PrimeOrderField::new(5);
        let matrix_to_factor_packet     =   MatrixAlgebraPacket::with_default_order( & matrix_to_factor, ring_operator);

        let row_indices_in_reverse_order         = (0..2).rev();
        
        let _target_comb_inv_pivot_block = get_pivot_block_of_target_comb_inverse_with_deleted_diagonal(
                                            & matrix_to_factor_packet,
                                            row_indices_in_reverse_order,
                                        );       
                                        
        println!("{:#?}", & _target_comb_inv_pivot_block.0);
        println!("{:#?}", & _target_comb_inv_pivot_block.1);        
    }


    /// This test targets the initial computation of the pivot block of the inverse of the target COMB.    
    #[test]
    fn test_initial_decomposition_another_example() {
        use itertools::Itertools;

        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::matrices::operations::umatch::row_major::get_pivot_block_of_target_comb_inverse_with_deleted_diagonal;
                
        use crate::algebra::matrices::operations::invert::InverseUpperTriangularMatrix;
        
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::utilities::order::{OrderOperatorByKey};


        // set some parameters
        let matrix_size = 5;
        let max_column_index = matrix_size - 1;
        let modulus = 7;
        let ring_operator       =   PrimeOrderField::new( modulus );        

        // define the matrix
        let matrix_to_factor = VecOfVec::new(
                                    vec![
                                        vec![(0, 1), (1, 2),                 (4, 0)],
                                        vec![        (1, 1), (2, 0),         (4, 1)],
                                        vec![                (2, 1), (3, 0), (4, 0)],
                                        vec![                        (3, 1), (4, 0)],
                                        vec![                                (4, 1)],
                                    ]
                                ).ok().unwrap();
        let matrix_to_factor_packet  =   MatrixAlgebraPacket::with_default_order( &matrix_to_factor, ring_operator.clone() );

        // make a copy of the matrix with same set of columns but in reversed order
        let mut flipped_vertically = matrix_to_factor.clone();
        flipped_vertically.reverse_the_sequence_of_columns_in_place(max_column_index);
        let flipped_vertically_packet  =   MatrixAlgebraPacket::with_default_order( &flipped_vertically, ring_operator.clone() );


        let row_indices_in_reverse_order         = (0 .. matrix_size).rev();

        // compute the inverse
        let inverse =   InverseUpperTriangularMatrix::new( & matrix_to_factor_packet );    
        

        // compute the target COMB of `matrix_to_factor_transformed`
        //
        // NOTE: the target COMB of `matrix_to_factor_transformed` is the inverse of `matrix_to_factor`, where -- by contrast -- 
        //       the target COMB of `matrix_to_factor` is the identity matrix (which is less interesting)
        let (array_comb, matching) = get_pivot_block_of_target_comb_inverse_with_deleted_diagonal(
            & flipped_vertically_packet, // the matrix matrix_to_factor_transformed ought to implement Copy
            row_indices_in_reverse_order,
        );          


        for row_index in 0 .. matrix_size {

            // obtain a row of the inverse
            let inv_row     =    inverse.row(&row_index).collect_vec();

            // obtain a row of the target COMB
            let ordinal    =   matching.ordinal_for_row_index( &row_index ).unwrap();  

            let comb_off_diag_view =    (& array_comb)
                                        .row( &ordinal )
                                        // .iter()
                                        .map(   | (x,y)| 
                                                ( matching.row_index_for_ordinal( x ), y ) // reindex the row from ordinals to row indices
                                            );
            let mut comb_view   =   comb_off_diag_view.collect_vec();
            comb_view.push( (row_index, 1) );
            comb_view.sort();       
            
            // compute products
            let product_inv=   multiply_row_vector_with_matrix(
                                                                        inv_row.clone(),
                                                                        & flipped_vertically_packet,
                                                                        ring_operator,
                                                                        OrderOperatorByKey::new()
                                                                    );
            let product_umatch=   multiply_row_vector_with_matrix(
                                                                        comb_view.clone(),
                                                                        & flipped_vertically_packet,
                                                                        ring_operator,
                                                                        OrderOperatorByKey::new()
                                                                    );                                                                    
                                                                

            // compare
            // println!("MATRIX:");
            for k in 0 .. matrix_size { println!("{:?}", matrix_to_factor_packet.row(&k).collect_vec()  ) }
            // println!("MATRIX TRANSFORMED:");
            for k in 0 .. matrix_size { println!("{:?}", flipped_vertically_packet.row(&k)  ) }            
            // println!("product_inv: {:?}", product_inv.clone().collect_vec() );
            // println!("product_umatch: {:?}", product_umatch.clone().collect_vec());            
            // println!("COMB row len: {:?}, inv row len{:?}", comb_view.len(), inv_row.len());
            assert_eq!(comb_view, inv_row);
                            
        } 
        
        
    }    


    /// This test targets the initial computation of the pivot block of the inverse of the target COMB.
    ///
    /// The key idea of this test is the following fact: let M be a square upper-unitriangular matrix,
    /// and let N be the matrix obtained by reversing the order of columns of M.  Then the standard
    /// "cohomology algorithm," applied to N, produces a target COMB equal to M^{-1}.
    /// 
    /// This test applies the standard cohomology algorithm to compute a target COMB of N.  We 
    /// check to ensure that this target COMB equals M^{-1}.
    #[test]
    fn test_initial_decomposition_larger() {
        use itertools::Itertools;

        use crate::algebra::matrices::operations::umatch::row_major::get_pivot_block_of_target_comb_inverse_with_deleted_diagonal;
                
        use crate::algebra::matrices::operations::invert::InverseUpperTriangularMatrix;
        use crate::algebra::matrices::query::MatrixOracle;        
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::utilities::order::OrderOperatorByKey;

        // set parameters
        let matrix_size     =   10;
        let max_column_index    =   matrix_size - 1;
        let modulus             =   7;
        let ring_operator       =   PrimeOrderField::new( modulus );        
        
        // randomly generate a matrix
        let matrix_to_factor = VecOfVec::random_mod_p_upper_unitriangular( matrix_size, modulus );
        let matrix_to_factor_packet  =   MatrixAlgebraPacket::with_default_order( &matrix_to_factor, ring_operator.clone() );

        // make a copy of the matrix with same set of columns but in reversed order
        let mut flipped_vertically = matrix_to_factor.clone();
        flipped_vertically.reverse_the_sequence_of_columns_in_place(max_column_index);
        let flipped_vertically_packet  =   MatrixAlgebraPacket::with_default_order( &flipped_vertically, ring_operator.clone() );


        // compute the inverse
        let inverse =   InverseUpperTriangularMatrix::new( & matrix_to_factor_packet );    


        // compute the target COMB of `matrix_to_factor_transformed`
        //
        // NOTE: the target COMB of `matrix_to_factor_transformed` is the inverse of `matrix_to_factor`, where -- by contrast -- 
        //       the target COMB of `matrix_to_factor` is the identity matrix (which is less interesting)
        let row_indices_in_reverse_order         = (0 .. matrix_size).rev();                
        let (array_comb, matching) = get_pivot_block_of_target_comb_inverse_with_deleted_diagonal(
            & flipped_vertically_packet,
            row_indices_in_reverse_order,
        );          


        for row_index in 0 .. matrix_size {

            // obtain a row of the inverse
            let inv_row     =    inverse.row(&row_index).collect_vec();

            // obtain a row of the target COMB
            let ordinal    =   matching.ordinal_for_row_index( &row_index ).unwrap();

            let comb_off_diag_view =    (& array_comb)
                                        .row( & ordinal )
                                        // .iter()
                                        .map(   | (x,y)| 
                                                ( matching.row_index_for_ordinal( x ), y ) // reindex the row from ordinals to row indices
                                            );
            let mut comb_view   =   comb_off_diag_view.collect_vec();
            comb_view.push( (row_index, 1) );
            comb_view.sort();       
            
            // compute products
            let product_inv=   multiply_row_vector_with_matrix(
                                                                        inv_row.clone(),
                                                                        & matrix_to_factor_packet,
                                                                        ring_operator,
                                                                        OrderOperatorByKey::new()
                                                                    );
            let product_umatch=   multiply_row_vector_with_matrix(
                                                                        comb_view.clone(),
                                                                        & matrix_to_factor_packet,
                                                                        ring_operator,
                                                                        OrderOperatorByKey::new()
                                                                    );                                                                    
                                                                

            // compare
            // println!("MATRIX:");
            for k in 0 .. matrix_size { println!("{:?}", matrix_to_factor_packet.row(&k).collect_vec()  ) }
            // println!("MATRIX TRANSFORMED:");
            for k in 0 .. matrix_size { println!("{:?}", flipped_vertically_packet.row(&k)  ) }            
            // println!("product_inv: {:?}", product_inv.clone().collect_vec() );
            // println!("product_umatch: {:?}", product_umatch.clone().collect_vec());            
            // println!("COMB row len: {:?}, inv row len{:?}", comb_view.len(), inv_row.len());
            assert_eq!(comb_view, inv_row);
                            
        } 
        
        
    }    



    //  ===================================================================================
    //  RECOVERY OF ROWS OF COMBS -- COMPREHENSIVE
    //  ===================================================================================    


    /// Checks that Umatch decomposition is correct (using a small example matrix, D) in the following sense:
    /// T^{-1} * T = I
    /// S^{-1} * S = I
    /// T^{-1} * D * S = M
    /// And the rows of T, T^{-1}, S, and S^{-1} appear in strictly ascending order
    #[test]
    fn test_umatchrowmajor_comprehensive_small() {
        
        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
        use crate::algebra::matrices::types::product::ProductMatrix;
        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::matrices::display::print_indexed_rows;           

        let num_indices_row             =   2;
        let num_indices_col             =   3;
        let modulus                     =   3;

        let ring_operator               =   PrimeOrderField::new( modulus );
        let matrix_to_factor_data        =   VecOfVec::new( 
            vec![   
                                vec![(0,1), (1,2), (2,0)], 
                                vec![              (2,0)]  
                            ] ).ok().unwrap();
        let matrix_to_factor = MatrixAlgebraPacket::with_default_order( & matrix_to_factor_data, ring_operator );

        
        let umatch  =   Umatch::new(
                            matrix_to_factor, 
                            (0..num_indices_row).rev(),
                        );


        


        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch );
        let matching = umatch.generalized_matching_matrix_ref();

        let comb_target = umatch.target_comb();
        let comb_target_inv = umatch.target_comb_inverse();        
        let comb_source = umatch.source_comb();        
        let comb_source_inv = umatch.source_comb_inverse(); 

        let comb_target_ref             =   & comb_target;
        let comb_target_inv_ref         =   & comb_target_inv;
        let comb_source_ref             =   & comb_source;
        let comb_source_inv_ref         =   & comb_source_inv;                        
        
        
        let product_source = ProductMatrix::new( comb_source_ref, comb_source_inv_ref );
        let product_target = ProductMatrix::new( comb_target_ref, comb_target_inv_ref );        
        let product_target_comb_inv_times_matrix_to_factor = ProductMatrix::new( comb_target_inv_ref, matrix_to_factor );      
        let product_target_comb_inv_times_matrix_to_factor_times_source_comb = ProductMatrix::new( product_target_comb_inv_times_matrix_to_factor, comb_source_ref );                


        // println!("matrix_to_factor:");
        // print_indexed_rows( & matrix_to_factor, 0 .. num_indices_row );
        // println!("matching:");
        // print_indexed_rows( & matching, 0 .. num_indices_row );        
        // println!("comb_source:");
        // print_indexed_rows( & comb_source, 0 .. num_indices_col );        
        // println!("comb_source_inv:");
        // print_indexed_rows( & comb_source_inv, 0 .. num_indices_col );     
        // println!("comb_target:");
        // print_indexed_rows( & comb_target, 0 .. num_indices_row );        
        // println!("comb_target_inv:");
        // print_indexed_rows( & comb_target_inv, 0 .. num_indices_row );    
        // println!("comb_target_inv * matrix_to_factor * comb_source:");
        // print_indexed_rows( & product_target_comb_inv_times_matrix_to_factor_times_source_comb, 0 .. num_indices_row );                                
        for column_index in 0 .. num_indices_col {
            // println!("{:?}", product_source.row( column_index ).collect_vec() );
            itertools::assert_equal( product_source.row( & column_index ), std::iter::once( (column_index, 1) )   );
        }


        // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
        for column_index in 0 .. num_indices_col { 
            assert_eq!(
                product_source.row( & column_index ).collect_vec(),
                vec![ (column_index, 1) ]
            ) 
        }

        // check that the product of the target COMB with its inverse is identity T * T^{-1} = I
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target.row( &row_index ).collect_vec(),
                vec![ (row_index, 1) ]
            ) 
        }    
        
        // check the factorization T^{-1} * D * S = M
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target_comb_inv_times_matrix_to_factor_times_source_comb.row( &row_index ).collect_vec(),
                matching.row( &row_index ).collect_vec()
            ) 
        }     
        
        // check that rows are sorted in strictly ascending order
        for row_index in 0 .. num_indices_row { 
            assert!(    matrix_to_factor.row( & row_index             ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_target.row( & row_index       ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_target_inv.row( & row_index   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_source.row( & row_index         ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_source_inv.row( & row_index     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
        }        
        
        // check that the major and columns of the inverse of the target COMB agree
        let comb_target_vec_of_vec_simple     
            =   VecOfVec::from_iterable_of_iterables( 
                    (0..num_indices_row)
                        .map( |k| comb_target_inv.row(&k) ) 
                ).ok().unwrap();
        for row_index in 0..num_indices_row {
            itertools::assert_equal( 
                    comb_target.column_reverse( & row_index ),
                    (& comb_target_vec_of_vec_simple).column_reverse( & row_index )
                )
        }

        // check that the major and columns of `TargetCombInverse` agree
        let comb_target_inv_vec_of_vec_simple     
            =   VecOfVec::from_iterable_of_iterables( 
                    (0..num_indices_row).map( |k| comb_target_inv.row(&k) ) 
                ).ok().unwrap();
        for row_index in 0..num_indices_row {
            // println!("PRINTING HERE: see below");
            // println!("ROW INDEX = {:?}", row_index );
            // println!("{:?}", comb_target_inv.column_reverse( row_index ).collect_vec());
            // println!("{:?}", (& comb_target_inv_vec_of_vec_simple).column_reverse( row_index ).collect_vec());            
            assert_equal( 
                    comb_target_inv.column_reverse( & row_index ).collect_vec(),
                    (& comb_target_inv_vec_of_vec_simple).column_reverse( & row_index ).collect_vec()
                )
        }

    }



    /// Checks that Umatch decomposition is correct (using a random example matrix, D) in the following sense:
    /// T^{-1} * T = I
    /// S^{-1} * S = I
    /// T^{-1} * D * S = M   
    /// And the rows of T, T^{-1}, S, and S^{-1} appear in strictly ascending order 
    #[test]
    fn test_umatchrowmajor_comprehensive_overall() {   

        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
        use crate::algebra::matrices::types::product::ProductMatrix;
        use crate::algebra::matrices::query::MatrixOracle;               

        let num_indices_row             =   10;
        let num_indices_col             =   20;
        let approximate_density           =   0.2;
        let modulus                     =   17;
        let allow_nonstructural_zero     =   true;

        let ring_operator           =   PrimeOrderField::new( modulus );
        let matrix_to_factor_data       =   VecOfVec::random_mod_p_with_density( num_indices_row, num_indices_col, approximate_density, modulus, allow_nonstructural_zero );
        let matrix_to_factor = MatrixAlgebraPacket::with_default_order( & matrix_to_factor_data, ring_operator );
     

        let umatch 
            =   Umatch::new( 
                    matrix_to_factor, 
                    (0..num_indices_row).rev(), 
                );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch );
        let matching = umatch.generalized_matching_matrix_ref();
  

        let comb_target = umatch.target_comb();
        let comb_target_inv = umatch.target_comb_inverse();        
        let comb_source = umatch.source_comb();        
        let comb_source_inv = umatch.source_comb_inverse();  
        let target_comb_inverse_times_matrix_to_factor_matched_block = umatch.target_comb_inverse_times_matrix_to_factor_matched_block();  
        let target_comb_inverse_times_matrix_to_factor_matched_block_rows_indexed_by_column_index = umatch.target_comb_inverse_times_matrix_to_factor_matched_block_with_rows_indexed_by_matched_column_index();

        let comb_target_ref         =   & comb_target;
        let comb_target_inv_ref         =   & comb_target_inv;
        let comb_source_ref         =   & comb_source;
        let comb_source_inv_ref         =   & comb_source_inv;            
        let target_comb_inverse_times_matrix_to_factor_matched_block_ref     =   & target_comb_inverse_times_matrix_to_factor_matched_block;
        let target_comb_inverse_times_matrix_to_factor_matched_block_rows_indexed_by_column_index_ref = & target_comb_inverse_times_matrix_to_factor_matched_block_rows_indexed_by_column_index;
        

        let product_source = ProductMatrix::new( comb_source_ref, comb_source_inv_ref );
        let product_target = ProductMatrix::new( comb_target_ref, comb_target_inv_ref );        
        let product_target_comb_inv_times_matrix_to_factor = ProductMatrix::new( comb_target_inv_ref, matrix_to_factor );      
        let product_target_comb_inv_times_matrix_to_factor_times_source_comb = ProductMatrix::new( product_target_comb_inv_times_matrix_to_factor, comb_source_ref );                        


        // ----------------------------------------------------------------------------------------------------------------
        // check that both representations of the seed matrix A (= matched part of Tinv * D) are internally valid        

        let matched_row_indices     =   umatch.matched_row_indices_in_ascending_order();
        let matched_column_indices   =   umatch.matched_column_indices_in_ascending_order();

        assert!(
            matrix_oracle_is_internally_consistent(
                target_comb_inverse_times_matrix_to_factor_matched_block_rows_indexed_by_column_index_ref, 
                matched_column_indices.iter().cloned(), // sorted row indices (which happen to be the same as the column indices)
                matched_column_indices.iter().cloned(), // sorted column indices
            )
            &&
            matrix_oracle_is_internally_consistent(
                target_comb_inverse_times_matrix_to_factor_matched_block_ref, 
                matched_row_indices.iter().cloned(), // sorted row indices
                matched_column_indices.iter().cloned(), // sorted column indices
            )            
        );


        // ----------------------------------------------------------------------------------------------------------------
        // check that all four (inverse) COMB's are internally valid
        // see documentation for `matrix_oracle_is_internally_consistent`, for details    

        assert!(
            matrix_oracle_is_internally_consistent(
                comb_source_ref, 
                0..num_indices_col, 
                0..num_indices_col
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_source_inv_ref, 
                0..num_indices_col, 
                0..num_indices_col
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_target_ref, 
                0..num_indices_row, 
                0..num_indices_row
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_target_inv_ref, 
                0..num_indices_row, 
                0..num_indices_row
            )                               
        );
      

        // ----------------------------------------------------------------------------------------------------------------
        // check that all four (inverse) COMB's return entries in the proper order
        // see documentation for `matrix_order_operators_are_internally_consistent`, for details

        assert!(
            matrix_order_operators_are_internally_consistent(
                comb_source_ref, 
                0..num_indices_col, 
                0..num_indices_col
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_source_inv_ref, 
                0..num_indices_col, 
                0..num_indices_col
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_target_ref, 
                0..num_indices_row, 
                0..num_indices_row
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_target_inv_ref, 
                0..num_indices_row, 
                0..num_indices_row
            ).is_ok()                                 
        );        

        // ----------------------------------------------------------------------------------------------------------------



        // println!("matrix_to_factor:");
        // print_indexed_rows( & matrix_to_factor, 0 .. num_indices_row );
        // println!("matching:");
        // print_indexed_rows( & matching, 0 .. num_indices_row );        
        // println!("comb_source:");
        // print_indexed_rows( & comb_source, 0 .. num_indices_col );        
        // println!("comb_source_inv:");
        // print_indexed_rows( & comb_source_inv, 0 .. num_indices_col );     
        // println!("comb_target:");
        // print_indexed_rows( & comb_target, 0 .. num_indices_row );        
        // println!("comb_target_inv:");
        // print_indexed_rows( & comb_target_inv, 0 .. num_indices_row );                        
        // println!("comb_target_inv * matrix_to_factor * comb_source:");
        // print_indexed_rows( & product_target_comb_inv_times_matrix_to_factor_times_source_comb, 0 .. num_indices_row );                                        
        for column_index in 0 .. num_indices_col {
            println!("row: {:?}", column_index );
            println!("row {:?}: {:?}", column_index, product_source.row( & column_index ).collect_vec() );
            itertools::assert_equal( product_source.row( & column_index ), std::iter::once( (column_index, 1) )   );
        }


        // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
        for column_index in 0 .. num_indices_col { 
            assert_eq!(
                product_source.row( & column_index ).collect_vec(),
                vec![ (column_index, 1) ]
            ) 
        }
      

        // check that the product of the target COMB with its inverse is identity T * T^{-1} = I
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target.row( & row_index ).collect_vec(),
                vec![ (row_index, 1) ]
            ) 
        }   
       
        
        // check the factorization T^{-1} * D * S = M
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target_comb_inv_times_matrix_to_factor_times_source_comb.row( &row_index ).collect_vec(),
                matching.row( & row_index ).collect_vec()
            ) 
        }   
       

        // check that the major and columns of `SourceCombInverseMatchedBlockRowsIndexedByColumnIndex` agree (i.e. that, taken all together, they run over the same entries)     
        verify_rows_compatible_with_columns(
                target_comb_inverse_times_matrix_to_factor_matched_block_rows_indexed_by_column_index_ref,
                umatch.generalized_matching_matrix_ref().bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order().iter().cloned(),
                umatch.generalized_matching_matrix_ref().bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order().iter().cloned(),                
            );
           
        
        // check that the major and columns of `SourceCombInverseMatchedBlockRowsIndexedByColumnIndex` are sorted
        for column_index in umatch.generalized_matching_matrix_ref().bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order().iter().cloned() {
            assert!(    target_comb_inverse_times_matrix_to_factor_matched_block_rows_indexed_by_column_index_ref.column_reverse( & column_index ).is_sorted_by( |x, y| x.0 > y.0 )     );
        }


        for column_index in umatch.generalized_matching_matrix_ref().bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order().iter().cloned() {
            assert!(    target_comb_inverse_times_matrix_to_factor_matched_block_rows_indexed_by_column_index_ref.row( & column_index ).is_sorted_by( |x, y| x.0 < y.0 )     );
        }     
        

        // ----------------------------------------------------------------------------------------------------------------
        // check that all four (inverse) COMB's are internally valid
        // see documentation for `matrix_oracle_is_internally_consistent`, for details

        assert!(
            matrix_oracle_is_internally_consistent(
                comb_source_ref, 
                0..num_indices_col, 
                0..num_indices_col
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_source_inv_ref, 
                0..num_indices_col, 
                0..num_indices_col
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_target_ref, 
                0..num_indices_row, 
                0..num_indices_row
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_target_inv_ref, 
                0..num_indices_row, 
                0..num_indices_row
            )                               
        );   

        // ----------------------------------------------------------------------------------------------------------------
        // check that all four (inverse) COMB's return entries in the proper order
        // see documentation for `matrix_order_operators_are_internally_consistent`, for details

        assert!(
            matrix_order_operators_are_internally_consistent(
                comb_source_ref, 
                0..num_indices_col, 
                0..num_indices_col
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_source_inv_ref, 
                0..num_indices_col, 
                0..num_indices_col
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_target_ref, 
                0..num_indices_row, 
                0..num_indices_row
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_target_inv_ref, 
                0..num_indices_row, 
                0..num_indices_row
            ).is_ok()                                 
        );        

        // ----------------------------------------------------------------------------------------------------------------



        // ----------------------------


        // check that `SourceCombInverseMatchedBlockRowsIndexedByColumnIndex` is upper triangular
        for column_index in umatch.generalized_matching_matrix_ref().bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order().iter().cloned() {
            assert!(    target_comb_inverse_times_matrix_to_factor_matched_block_rows_indexed_by_column_index_ref.column_reverse( & column_index ).next().unwrap().0 == column_index     );
        }        


// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:BEGIN        
        // check that columns are sorted in strictly descending order
        //  NOTE: THIS IS UNNECESSARY FOR THE COMBS, SINCE WE TEST THAT THEIR columnS EQUAL THOSE OF VecOfVec objects, WHOSE MINOR DESCENDING VIEWS ARE *ALWAYS* STRICTLY DECREASING IN INDEX
        for row_index in 0 .. num_indices_col { 
            assert!(    matrix_to_factor.column_reverse( & row_index             ).is_sorted_by( |x, y| x.0 > y.0 )     );
            assert!(    comb_target.column_reverse( & row_index       ).is_sorted_by( |x, y| x.0 > y.0 )     );
            // assert!(    comb_target_inv.column_reverse( row_index.clone()   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            // assert!(    comb_source.column_reverse( row_index.clone()         ).is_sorted_by( |x, y| x.0 < y.0 )     );
            // assert!(    comb_source_inv.column_reverse( row_index.clone()     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
        }          
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:END  

        
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG: BEGIN   
        // check that the major and columns of the inverse of the target COMB agree
        let comb_target_vec_of_vec_simple     
            =   VecOfVec::from_iterable_of_iterables(
                    (0..num_indices_row).map( |k| comb_target.row(&k) ) 
                ).ok().unwrap();
        for row_index in 0..num_indices_row {
            // println!("VIEW MAJOR DESCEND IS STARTING FOR THIS ROUND: row_index = {:?}", row_index);
            // println!("VIEW MAJOR DESCEND LAZY CONSTRUCTION: {:?}", comb_target.column_reverse( row_index ).collect_vec());
            // println!("VIEW MAJOR DESCEND FROM ROW: {:?}", (& comb_target_vec_of_vec_simple).column_reverse( row_index ).collect_vec());            
            // println!("VIEW MAJOR DESCEND IS FINISHED FOR THIS ROUND");
            itertools::assert_equal( 
                    comb_target.column_reverse( & row_index ),
                    (& comb_target_vec_of_vec_simple).column_reverse( & row_index )
                )
        }      
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:END      
        
        // check that rows are sorted in strictly ascending order
        for row_index in 0 .. num_indices_row { 
            assert!(    matrix_to_factor.row(   & row_index   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_target.row(        & row_index   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_target_inv.row(    & row_index   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_source.row(        & row_index   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_source_inv.row(    & row_index   ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
        }  
              

    }



    //  ===================================================================================
    //  RECOVERY OF COMBS -- TARGETTED AT SPECIFIC POINTS IN THE PROCESS
    //  ===================================================================================

    
    //  COMB DOMAIN INV (SMALL + LARGE)
    //  ------------------------------------------------------------------------------------
    //  TAG USED TO IDENTIFY THESE TESTS: #SourceCombInvImplementRow    

    #[test]
    fn test_retreival() {
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch, TargetCombInverseTimesMatrixToFactorMatchedBlock};
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;                 
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;

        let matrix_to_factor_data: VecOfVec< usize, usize >   =   
        VecOfVec::new(  
                vec![
                    vec![ (0, 1), (1, 1), (2, 2) ],
                    vec![ (0, 1),         (2, 1) ],
                    vec![         (1, 1), (2, 1) ],
                ]
            ).ok().unwrap();
        let ring_operator       =   PrimeOrderField::new( 13 );
        let matrix_to_factor_packet     =   MatrixAlgebraPacket::with_default_order( &matrix_to_factor_data, ring_operator );

        let umatch  
            =   Umatch::new( 
                        matrix_to_factor_packet, 
                        (0..3).rev(),                    
                    );
        let _umatch_ref = & umatch;
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( umatch_ref);
        // let umatch_with_refs_ref = & umatch_with_refs;
        
        //  the "seed" matrix A (equal to the pivot block of the inverse of the target COMB times the pivot block of the matching array)
        let A   =   TargetCombInverseTimesMatrixToFactorMatchedBlock::new( & umatch );

        //  the source COMB
        let comb_source_inv = umatch.source_comb_inverse();
        let comb_source_inv_ground_truth
                =   VecOfVec::new(
                            vec![
                                vec![ (0, 1),         (2, 1) ],
                                vec![         (1, 1), (2, 1) ],
                                vec![                 (2, 1 )                        ],
                            ]
                    ).ok().unwrap();
        let comb_source_inv_ground_truth_ref = & comb_source_inv_ground_truth;
        for row_index in 0 .. 3 {
            // println!("GROUND TRUTH  : {:?}", comb_source_inv_ground_truth_ref.row( &row_index ).collect_vec() );
            // println!("UNPACKED      : {:?}", comb_source_inv.row( &row_index ).into_iter().collect_vec() );   
            // println!("SCALE FACTORS : {:?}", umatch_with_refs.matching.structural_nonzero_values_in_sequence() );    
            // println!("row_index        : {:?}", row_index );                                
            itertools::assert_equal(
                    comb_source_inv_ground_truth_ref.row( & row_index ),
                    comb_source_inv.row( &row_index ),
                )
        }
    }





    //  COMB DOMAIN (SMALL)
    //  ------------------------------------------------------------------------------------
    //  TAG USED TO IDENTIFY THESE TESTS: #SourceCombImplementRow

    #[test]
    fn test_umatchrowmajor_comb_source_small_example() {

        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;                 
        use crate::algebra::matrices::types::product::ProductMatrix;     

        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;

        // let num_rows = 2; let num_cols = 2; let modulus = 7;
        // let matrix_to_factor = VecOfVec::new( vec![ vec![ (0usize,5usize), (1,5)], vec![ (1,6)]] );
        let num_rows = 1; let num_cols = 4; let modulus = 7;

        let ring_operator       =   PrimeOrderField::new( 13 );        
        let matrix_to_factor_data = VecOfVec::new( vec![ vec![ (2usize, 6usize), (3,1)], ]  ).ok().unwrap();        
        //  NOTE: matrix_to_factor can be regarded as     [  0  0  6  1  ]
        let matrix_to_factor_packet     =   MatrixAlgebraPacket::with_default_order( &matrix_to_factor_data, ring_operator );    

        let umatch_root = 
                Umatch::new( 
                        matrix_to_factor_packet, 
                        (0 .. num_rows).rev(), 
                    );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch_root );
        // let umatch_ref = & umatch_with_refs;

        let comb_source = umatch_root.source_comb();
        let comb_source_inv = umatch_root.source_comb_inverse();     
        let matrix_to_factor_matched_columns_only = umatch_root.matrix_to_factor_matched_columns_only();          
                
        // check that S * S^{-1} = identity
        let c_times_c_inv = 
            ProductMatrix::new( 
                    &comb_source, 
                    &comb_source_inv, 
                );

        // println!("matrix_to_factor:");
        // print_indexed_rows( & matrix_to_factor_ref, 0 .. num_rows );
        // println!("matrix_to_factor_matched_columns_only:");
        // print_indexed_rows( & matrix_to_factor_matched_columns_only, 0 .. num_rows );        
        // println!("matching:");
        // print_indexed_rows( & umatch_root.generalized_matching_matrix_ref(), 0 .. num_rows );    
        // println!("target_comb_inverse_times_matrix_to_factor_matched_block (recall that num_rows = {:?}):", num_rows);        
        // print_indexed_rows( && umatch_root.target_comb_inverse_times_matrix_to_factor_matched_block(), 0 .. num_rows );
        // println!("comb_source (recall that num_cols = {:?}) (THIS FUNCTION CALL SEEMS TO BREAK DOWN INTERNALLY WHERE target_comb_inverse_times_matrix_to_factor_matched_block IS CALLED):", num_cols);
        // print_indexed_rows( & comb_source, 0 .. num_cols );        
        // println!("comb_source_inv:");
        // print_indexed_rows( & comb_source_inv, 0 .. num_cols );                
        for column_index in 0 .. num_cols {
            // println!("{:?}", c_times_c_inv.row( column_index ).collect_vec() );
            itertools::assert_equal( c_times_c_inv.row( & column_index ), std::iter::once( (column_index, 1) )   );
        }

        // check D * S is right-reduced

    }


    //  COMB DOMAIN (LARGER + RANDOM)
    //  ------------------------------------------------------------------------------------
    //  TAG USED TO IDENTIFY THESE TESTS: #SourceCombImplementRow)    


    #[test]
    fn test_umatchrowmajor_comb_source() {

        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::types::product::ProductMatrix;         
        use crate::algebra::matrices::query::MatrixOracle;        
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField; 

        let num_rows = 10; let num_cols = 10; let modulus = 7;
        let matrix_to_factor_data = VecOfVec::random_mod_p(num_rows, num_cols, modulus);
        let ring_operator = PrimeOrderField::new( modulus );
        let matrix_to_factor_packet     =   MatrixAlgebraPacket::with_default_order( &matrix_to_factor_data, ring_operator );     

        let umatch_root = 
                Umatch::new( 
                        matrix_to_factor_packet, 
                        (0 .. num_rows).rev(),                  
                    );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch_root );
        // let umatch_ref = & umatch_with_refs;

        let comb_source = umatch_root.source_comb();
        let comb_source_inv = umatch_root.source_comb_inverse();     
                
        // check that S * S^{-1} = identity
        let c_times_c_inv = 
            ProductMatrix::new( 
                    & comb_source, 
                    & comb_source_inv, 
                );

        // println!("matrix_to_factor:");
        // print_indexed_rows( & matrix_to_factor_ref, 0 .. num_rows );
        // println!("matching:");
        // print_indexed_rows( & umatch_root.generalized_matching_matrix_ref(), 0 .. num_rows );        
        // println!("comb_source:");
        // print_indexed_rows( & comb_source, 0 .. num_cols );        
        // println!("comb_source_inv:");
        // print_indexed_rows( & comb_source_inv, 0 .. num_cols );                
        for column_index in 0 .. num_cols {
            // println!("{:?}", c_times_c_inv.row( column_index ).collect_vec() );
            itertools::assert_equal( c_times_c_inv.row( & column_index ), std::iter::once( (column_index, 1) )   );
        }

        // check D * S is right-reduced

    }


    #[test]
    fn doc_test() {
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, product::ProductMatrix};
        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
        use crate::algebra::matrices::debug::product_is_identity_matrix;
        use crate::algebra::matrices::display::print_indexed_columns;      
        use crate::algebra::matrices::query::MatrixOracle;
        
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;        
        
        use itertools::Itertools;

        // DEFINE INPUTS
        // ===============================

        // define the ring operator and order operator
        let modulus               =   5;
        let ring_operator         =   PrimeOrderField::new( modulus );        

        // define the matrix we wish to factor
        let matrix_to_factor_data          =   VecOfVec::new( 
                                                vec![   
                                                            vec![(0,1), (1,1), (2,1)],
                                                            vec![                   ], 
                                                            vec![              (2,1)], 
                                                ] 
                                            ).ok().unwrap();
        let matrix_to_factor_packet     =   MatrixAlgebraPacket::with_default_order( &matrix_to_factor_data, ring_operator );     
                                        
        // COMPUTE U-MATCH
        // ===============================
                                        
        let umatch
            =   Umatch::new(
                    matrix_to_factor_packet,  // the matrix we wish to factor
                    (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
                );

        println!("matrix_major");
        umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal.matrix_rows_data.print_dense(0);
        println!("{:?}", umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal.matrix_rows_data.inner_vec_of_vec_ref() );
        
        println!("matrix_minor");
        umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal.matrix_columns_data.print_dense(0);        
        println!("{:?}", umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal.matrix_columns_data.inner_vec_of_vec_ref() );        
        
        println!("number of pairs: {:?}", umatch.matching.number_of_structural_nonzeros() );        
        
        println!(
            "result of transpose: {:?}", 
            umatch.matched_block_of_target_comb_inverse_indexed_by_ordinal_of_matched_row_off_diagonal.matrix_rows_data.transpose_deep( umatch.generalized_matching_matrix_ref().number_of_structural_nonzeros() ).unwrap().inner_vec_of_vec_ref()
        );
            
            
        // INSPECT FACTORIZATION
        // ===============================
            
        // extract T, T^{-1}, S, S^{-1}, and M
        let t           =   umatch.target_comb();        // the target COMB
        let tinv        =   umatch.target_comb_inverse();    // inverse of the the target COMB
        let s           =   umatch.source_comb();          // the source COMB
        let sinv        =   umatch.source_comb_inverse();      // inverse of the source COMB
        let m           =   umatch.generalized_matching_matrix_ref();         // the generalized matching matrix
            
            
        println!("\nColumns of the target COMB");   print_indexed_columns( &t, 0..3 ); 
        println!("\nColumns of the   source COMB");   print_indexed_columns( &s, 0..3 ); 
        println!("\nColumns of the generalized matching matrix"); print_indexed_columns( &m, 0..3 ); 
            
        // this will print the following:
        //
        // Columns of the target COMB
        // column 0: [(0, 1)]
        // column 1: [(1, 1)]
        // 
        // Columns of the   source COMB
        // column 0: [(0, 1)]
        // column 1: [(1, 1), (0, 3)]
        // column 2: [(2, 1), (0, 2)]
        // 
        // Columns of the generalized matching matrix
        // column 0: [(0, 1)]
        // column 1: []
        // column 2: [(1, 1)]

        // SOLVE Ax = b FOR x
        // ===============================
        
        let b   =   [ (0,1), (2,1) ]; 
        let x   =   umatch.solve_dx_equals_b( b.clone() ).unwrap();
        let dx  =   umatch.multiply_dx(x).collect_vec();
        assert!( dx.eq( & b ) );
            
            
        // VERIFY THE CALCULATION
        // ===============================
            
        // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
        product_is_identity_matrix( &s, &sinv, 0..3 );
            
        // check that the product of the target COMB with its inverse is identity: T * T^{-1} = I
        product_is_identity_matrix( &t, &tinv, 0..3 );
            
        // check the factorization: T^{-1} * D * S = M
        let rinv_d   = ProductMatrix::new( &tinv,   &matrix_to_factor_packet );      
        let rinv_d_c = ProductMatrix::new( &rinv_d, &s );                
        for row_index in 0 .. 3 { 
            assert_eq!(
                rinv_d_c.row( &row_index ).collect_vec(),
                m.row( &row_index ).collect_vec()
            ) 
        }            
    }


    
 
}    


#[cfg(test)]
mod doc_test_solvers {
    use crate::algebra::matrices::types::packet::MatrixAlgebraPacket;

    

    

    


    #[test]
    fn doc_test_solve_dx_equals_b() {
        
        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;        
        use itertools::Itertools;

        // DEFINE THE MATRIX
        // ===============================
        let matrix          =   VecOfVec::new( 
                                                vec![   
                                                            vec![(0,true), (1,true), (2,true)],
                                                            vec![                            ], 
                                                            vec![                    (2,true)], 
                                                ] 
                                            ).ok().unwrap();
        let matrix_packet       =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( & matrix );
                                        
        // COMPUTE U-MATCH
        // ===============================
                                        
        let umatch
            =   Umatch::new(
                    & matrix_packet,  // the matrix we wish to factor
                    (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
                );        

        // SOLVE Ax = b FOR x
        // ===============================
        
        let b   =   [ (0,true), (2,true) ]; 
        let x   =   umatch.solve_dx_equals_b( b ).unwrap();
        let dx  =   umatch.multiply_dx(x);
        assert!( dx.eq( b ) );        
    }






    #[test]
    fn doc_test_solve_xd_equals_b() {

        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::operations::umatch::row_major::Umatch;     
        
        // define the matrix
        // -----------------
        let d = VecOfVec::new(
                vec![
                                        vec![  (0,true), (1,true),           ],
                                        vec![            (1,true), (2,true), ],
                                    ]
            ).ok().unwrap();

        let d_packet       =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( & d );

        // obtain a u-match factorization
        // ------------------------------
        let umatch  =   Umatch::new( 
            &d_packet, 
            0..2,      
        );
        
        // try solving xd = b
        // ------------------
        
        // Case 1: a solution exists; in this case we are gauaranteed to find one
        let x = umatch.solve_xd_equals_b( vec![ (0,true), (2,true), ] );        
        assert!( x.is_some() );
        assert!( x.unwrap().eq( & vec![ (0,true), (1,true), ] ) );

        // Case 2: no solution exists; in this case we get a certificate that no solution exists
        let x = umatch.solve_xd_equals_b( vec![ (0,true), (1,true), (2,true) ] );        
        assert!( x.is_none() );
    }


    #[test]
    fn doc_test_solve_xd_equals_b__withfloatcoefficients() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
         
        
        let d = VecOfVec::new(
                vec![
                                        vec![  (0,3.), (1,3.),         ],
                                        vec![          (1,3.), (2,3.), ],
                                    ]
            ).ok().unwrap();
        let d_packet       =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( & d );

        let umatch  =   Umatch::new( 
            & d_packet, 
            0..2,             
        );
        
        // instance where there exists a solution
        let x = umatch.solve_xd_equals_b( vec![ (0,3.), (1,6.), (2,3.), ] );        
        assert!( x.is_some() );
        assert!( x.unwrap().eq( &vec![ (0,1.), (1,1.) ] ) );

        // instance where there does not exist a solution
        let x = umatch.solve_xd_equals_b( vec![ (0,1.), (1,-1.) ] );        
        assert!( x.is_none() );
    }    



}    














