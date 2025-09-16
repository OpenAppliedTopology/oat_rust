

// -----------------------------------------------------------------------------------


//! Invert a triangular matrix; see also [`solve`](crate::algebra::matrices::operations::solve).
//! 
//! # How to improve performance
//! 
//! If you only need to solve `Ax = b` for `x`, where `A` is a sparse matrix, then you will probably get
//! better results by using the [`solve`](crate::algebra::matrices::operations::solve) module to calculate `x` via back-substitution than by computing the inverse of `A`,
//! then multiplying `x` with `A^{-1}.`.  See the blog post
//! [don't invert that matrix](https://www.johndcook.com/blog/2010/01/19/dont-invert-that-matrix/), by John Cook,
//! for general thoughts on inverting sparse matrices.
//! 
//! # Examples
//! 
//! Let's invert the matrix
//! 
//! ```text 
//!       |  1    1  |
//!       |  0    1  |
//! ```
//! 
//! with coefficients in the prime field of order 1049.
//! 
//! 
//! ```
//! use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
//! use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket;
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::invert::InverseUpperTriangularMatrix;
//! use oat_rust::utilities::order::OrderOperatorByLessThan;
//! use oat_rust::algebra::matrices::display::print_indexed_rows;
//! 
//! // Define the ring operator (prime field of order 1049)
//! let modulus = 1049;
//! let ring_operator  =   PrimeOrderField::new(modulus);                  
//! 
//! // Define a row-major upper triangular amtrix
//! let matrix  =   VecOfVec::new(
//!                     vec![   
//!                         vec![ (0,1), (1,1), ], 
//!                         vec![        (1,1), ],
//!                     ],
//!                 ).ok().unwrap();
//! let matrix  =   MatrixAlgebraPacket::with_default_order( &matrix, ring_operator.clone() );
//! 
//! // Define the inverse
//! let inverse =   InverseUpperTriangularMatrix::new( &matrix );   
//! 
//! // Print the inverse
//! let row_indices = vec![ 0, 1 ];
//! print_indexed_rows( &inverse, row_indices );
//! ```
//! 
//! This should print the following:
//! 
//! ```text
//! $ row 0: [(0, 1), (1, 1048)]
//! $ row 1: [(1, 1)] 
//! ```
//! 
//! which is the correct solution, since the inverse of `matrix` is
//! ```text 
//!       |  1   -1  |
//!       |  0    1  |
//! ```
//! and -1 = 1048, modulo 1049.


use std::vec::IntoIter;
use std::cmp::Ordering;
use std::fmt::Debug;

use crate::{algebra::{matrices::{query::{MatrixAlgebra, MatrixOracle}, types::transpose::OrderAntiTranspose}, rings::traits::{DivisionRingOperations, RingOperations}, vectors::operations::{Simplify, VectorOperations}}, utilities::{iterators::{general::PeekUnqualified, merge::hit::{hit_bulk_insert, hit_merge_by_predicate}}, order::JudgeOrder}};
use crate::algebra::matrices::operations::combine_rows_and_columns::{LinearCombinationOfRows};
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet};




//  ---------------------------------------------------------------------------
//  INVERT A TRIANGULAR ARRAY


/// An iterator that represents a row of the inverse of an upper triangular matrix.
/// 
/// This struct will only function correctly if the initial value of `next_entry_of_inv` and `entries_to_eliminate_simplified_heap`
/// are set correctly when the iterator is created. 
/// 
/// # Developer notes
/// 
/// This object is returned by `matrix.row( .. )` where `matrix` is a matrix of type [InverseUpperTriangularMatrix].
/// Alternatively, one could implement a different oracle for the same matrix (that is, for the inverse of a
/// triangular array) where `row( .. )` returns an object of type [`TriangularSolve`](crate::algebra::matrices::operations::triangular_solve::TriangularSolve).
/// In that implementation, `row( .. )` would simply return `TriangularSolve::new( x, matrix, ring_operator, order_operator)`,
/// where `x` is an entry iterator representing a standard unit vector.  However, there are some small differences between that approach and the one
/// implemented here, namely:
/// - the alternate approach wraps a (possibly large) number of iterators in wrappers of type [`TwoTypeIterator`](crate::utilities::iterators::general::TwoTypeIterator)
/// - the approach which is actually applied here has a `head-tail` structure which allows peeking
/// 
/// It would be interesting to contrast the performance of these two approaches, to see if there is a meaningful difference.
pub struct RowOfInverseOfUpperTriangularMatrix
                < Matrix > 
    where
        Matrix:     MatrixAlgebra< 
                        RowEntry:   KeyValSet,
                    >    
{
    ring_operator:                          Matrix::RingOperator, // the operator for the coefficient ring_operator    
    matrix:                                 Matrix, // the matrix to be inverted    
    next_entry_of_inv:                      Option< Matrix::RowEntry  >, // the next entry in this row of the inverse matrix     
    entries_to_eliminate_simplified_heap:   LinearCombinationOfRows< Matrix >,
}                                   

impl    < Matrix, Index > 

        Iterator for    

        RowOfInverseOfUpperTriangularMatrix
            < Matrix >

    where
        Index:      PartialEq,
        Matrix:     MatrixAlgebra< 
                        RowEntry:       KeyValSet,
                        ColumnIndex =   Index,
                        RowIndex =      Index, 
                        RingOperator:   DivisionRingOperations,
                    >                   

{

    type Item = Matrix::RowEntry;

    fn next( &mut self ) -> Option<Self::Item> { 

        match self.next_entry_of_inv.take() {
            // RETURN NONE IF `next_entry_of_inv` IS EMPTY
            None => None,  

            // OTHERWISE `next_entry_of_inv` IS NOT EMPTY, AND WE DO THE FOLLOWING
            Some( return_value ) => {

                // IF THE HEAD OF HEADTAIL IS NONEMPTY, IT BECOMES THE NEXT TARGET FOR ELIMINATION
                if let Some( mut entry_to_eliminate )       =   self.entries_to_eliminate_simplified_heap.next() {

                    // TO ELIMINATE THIS ENTRY, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX WE WANT TO INVERT; DENOTE THE UNSCLAED VECTOR BY v
                    let mut seed_of_eliminating_iterator    =   self.matrix.row( &entry_to_eliminate.key() ).into_iter();
                    
                    // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
                    let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

                    // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
                    let scale_factor     =    self.ring_operator.negate(
                                                            self.ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                        );
                    
                    // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
                    // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_eliminate_simplified_heap
                    let eliminating_iterator    =   seed_of_eliminating_iterator.scale_by( scale_factor.clone(), self.ring_operator.clone() );

                    // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_eliminate_simplified_heap`
                    hit_bulk_insert( &mut self.entries_to_eliminate_simplified_heap.unsimplified, vec![eliminating_iterator] ); 

                    // NOW RE-USE `entry_to_eliminate` AS A CONTAINER FOR THE NEXT ENTRY IN THE INVERSE MATRIX
                    entry_to_eliminate.set_val( scale_factor );

                    // THE NEXT NONZERO ENTRY ON THE INVERSE MATRIX (AFTER THE ONE WE ARE ABOUT TO RETURN, NAMELY return_value) HAS INDEX EQUAL TO `entry_to_eliminate.key()` AND COEFFICIENT EQUAL TO `scale_factor`
                    self.next_entry_of_inv = Some( entry_to_eliminate );

                }

                // RETURN THE NEXT ENTRY OF THE INVERSE MATRIX
                // we could have done this much earlier, in principle, but we had to set up the next value of `self.next_entry_of_inv` before we returned `return_value`
                Some( return_value )

            }
        }
    }
}         

impl    < Matrix, Index > 

        PeekUnqualified for    

        RowOfInverseOfUpperTriangularMatrix< 
                Matrix,
            > 

    where 
        Index:      PartialEq,
        Matrix:     MatrixAlgebra< 
                        RowEntry:       KeyValSet,
                        ColumnIndex =   Index,
                        RowIndex =      Index, 
                        RingOperator:   DivisionRingOperations,
                    >      


{
    fn peek_unqualified(&mut self) -> std::option::Option<&<Self as Iterator>::Item> { 
        match &self.next_entry_of_inv {
            None => { None },
            Some(x) => { Some( x ) }
        }
    }
}






/// Given an invertible upper triangular matrix `A`, returns a row of the inverse of `A`.
pub fn row_of_inverse_of_upper_triangular_matrix < 
            Matrix,
        > 
        (   
            row_index: & Matrix::RowIndex, 
            matrix: Matrix, 
        ) 
        ->      
        RowOfInverseOfUpperTriangularMatrix< 
                Matrix,
            > 
    where 
        Matrix:                         MatrixAlgebra< RingOperator: DivisionRingOperations >,
        Matrix::RowEntry:               KeyValSet,    
{
    // DEFINE THE RING AND ORDER OPERATORS    
    let ring_operator   =   matrix.ring_operator();
    let order_operator          =   matrix.order_operator_for_row_entries();

    let mut unscaled_seed_of_entries_to_be_eliminated   =   matrix.row( &row_index );
    
    // CONSTRUCT THE FIRST ENTRY OF THE INVERSE ARRAY; NOTE THAT WE 
    //      (1) TAKE THE FIRST ENTRY OF THE ARRAY TO BE INVERTED
    //      (2) **INVERT THE SCALAR VALUE** OF THAT ENTRY
    let mut diagonal_entry_of_inverse                  =   unscaled_seed_of_entries_to_be_eliminated.next().unwrap();
    diagonal_entry_of_inverse.set_val(
            ring_operator.invert( diagonal_entry_of_inverse.val() )
        );

    // CONSTRUCT THE OBJECT THAT ITERATORS OVER ENTRIES TO ELIMINATE 
    // This is obtained by taking the row/column of the original matrix, removing its leading entry, then scaling by 1/(the value of the leading entry)
    let scaled_seed_of_tail_to_be_eliminated = unscaled_seed_of_entries_to_be_eliminated.scale_by( 
                    diagonal_entry_of_inverse.val(),
                    ring_operator.clone(),   
                );
    // let head_to_be_eliminated   =   scaled_seed_of_tail_to_be_eliminated.next();
    // let tail_to_be_eliminated   =   hit_merge_by_predicate(
    //                                         vec![ scaled_seed_of_tail_to_be_eliminated ],
    //                                         order_operator,
    //                                     );

    let entries_to_eliminate_simplified_heap
            =   Simplify::new(
                        hit_merge_by_predicate(vec![ scaled_seed_of_tail_to_be_eliminated ], order_operator),
                        ring_operator.clone(),                                                                            
                    );
    
    RowOfInverseOfUpperTriangularMatrix {
        ring_operator, // the operator for the coefficient ring_operator    
        matrix, // the matrix to be inverted    
        next_entry_of_inv:              Some( diagonal_entry_of_inverse ), // the next entry in this row of the inverse matrix    
        entries_to_eliminate_simplified_heap,
    }

 }









/// Given an invertible upper triangular matrix `A`, returns a row of the inverse of `A`.
/// 
/// Returns `Err( row_index )` if the matrix does not have a row for `row_index`
pub fn row_of_inverse_of_upper_triangular_matrix_result < 
            Matrix,
        > 
        (   
            row_index: Matrix::RowIndex, 
            matrix: Matrix, 
        ) 
        ->
        Result<
            RowOfInverseOfUpperTriangularMatrix< Matrix >,
            Matrix::RowIndex,
        >
    where 
        Matrix:                         MatrixAlgebra< RingOperator: DivisionRingOperations >,
        Matrix::RowEntry:               KeyValSet,   
{

    // DEFINE THE RING AND ORDER OPERATORS
    let ring_operator           =   matrix.ring_operator();
    let order_operator          =   matrix.order_operator_for_row_entries();    

    let mut unscaled_seed_of_entries_to_be_eliminated   =   matrix.row_result( &row_index )?;
    
    // CONSTRUCT THE FIRST ENTRY OF THE INVERSE ARRAY; NOTE THAT WE 
    //      (1) TAKE THE FIRST ENTRY OF THE ARRAY TO BE INVERTED
    //      (2) **INVERT THE SCALAR VALUE** OF THAT ENTRY
    let mut diagonal_entry_of_inverse                  =   unscaled_seed_of_entries_to_be_eliminated.next().unwrap();
    diagonal_entry_of_inverse.set_val(
            ring_operator.invert( diagonal_entry_of_inverse.val() )
        );

    // CONSTRUCT THE OBJECT THAT ITERATORS OVER ENTRIES TO ELIMINATE 
    // This is obtained by taking the row/column of the original matrix, removing its leading entry, then scaling by 1/(the value of the leading entry)
    let scaled_seed_of_tail_to_be_eliminated = unscaled_seed_of_entries_to_be_eliminated.scale_by( 
                    diagonal_entry_of_inverse.val(),
                    ring_operator.clone(),   
                );
    // let head_to_be_eliminated   =   scaled_seed_of_tail_to_be_eliminated.next();
    // let tail_to_be_eliminated   =   hit_merge_by_predicate(
    //                                         vec![ scaled_seed_of_tail_to_be_eliminated ],
    //                                         order_operator,
    //                                     );

    let entries_to_eliminate_simplified_heap
            =   Simplify::new(
                        hit_merge_by_predicate(vec![ scaled_seed_of_tail_to_be_eliminated ], order_operator),
                        ring_operator.clone(),                                                                            
                    );
    
    let row     =   RowOfInverseOfUpperTriangularMatrix {
        ring_operator, // the operator for the coefficient ring_operator    
        matrix, // the matrix to be inverted    
        next_entry_of_inv:              Some( diagonal_entry_of_inverse ), // the next entry in this row of the inverse matrix    
        entries_to_eliminate_simplified_heap,
    };

    Ok(row)

 }











/// The inverse of an upper triangular matrix `M`, computed in a lazy fashion.
/// 
/// # Requirements
/// 
/// The order operator for row entries (respectively, indices) should be the same as the order operator for column entries (respectively, indices).
/// For a discussion of matrix order operators, see the documenation for the [MatrixAlgebra](crate::algebra::matrices::query::MatrixAlgebra) trait.
/// 
/// # Caution
/// 
/// By design, this object will clone `M` every time the user looks up a row or column of `M` inverse.
/// This is costly if `M` occupies a large amount of memory.
/// To avoid excess memory use, consider constructing the inverse with a reference: `InverseUpperTriangularMatrix::new( & M )`
/// rather than `InverseUpperTriangularMatrix::new( M )`.
/// 
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket;
/// use oat_rust::algebra::matrices::operations::invert::InverseUpperTriangularMatrix;
/// use oat_rust::utilities::order::OrderOperatorByLessThan;
/// use oat_rust::algebra::matrices::display::print_indexed_rows;
/// 
/// // Define the ring operator
/// let modulus = 1049;
/// let ring_operator  =   PrimeOrderField::new(modulus);                  
/// 
/// // Define an upper triangular amtrix
/// let data    =   VecOfVec::new(
///                       vec![   
///                           vec![ (0,1), (1,1), ], 
///                           vec![        (1,1), ],
///                       ],
///                   ).ok().unwrap();
/// 
/// // Wrap the matrix in a packet that contains the ring operator and (default) order operators
/// let matrix  =   MatrixAlgebraPacket::with_default_order( &data, ring_operator );
/// 
/// // Define the inverse
/// let inverse = InverseUpperTriangularMatrix::new( & matrix );  
/// 
/// // Print the inverse
/// let row_indices = vec![ 0, 1 ];
/// print_indexed_rows( &inverse, row_indices );
/// ```
/// 
/// This should print the following:
/// 
/// ```bash
/// $ row 0: [(0, 1), (1, 1048)]
/// $ row 1: [(1, 1)] 
/// ```
/// 
/// which is the correct solution, since the inverse of `matrix` is
/// ```bash 
///       |  1   -1  |
///       |  0    1  |
/// ```
/// and -1 = 1048, modulo 1049.
 #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
 pub struct InverseUpperTriangularMatrix
                < Matrix, > 

{
    matrix_to_invert:             Matrix, // the matrix to be inverted    
}





impl < Matrix > 
 
InverseUpperTriangularMatrix
    < Matrix > 
{
    pub fn new( matrix_to_invert: Matrix ) -> Self {
        InverseUpperTriangularMatrix{ matrix_to_invert }              
    }
}  







impl < Matrix, Index, Entry >

    MatrixOracle for 

    InverseUpperTriangularMatrix
        < Matrix > 

    where 
        Matrix:                 Clone + MatrixAlgebra< 
                                            RowIndex=Index, 
                                            ColumnIndex=Index, 
                                            RowEntry=Entry, 
                                            ColumnEntry=Entry,
                                            RingOperator: DivisionRingOperations,
                                        >,
        Entry:                  Clone + Debug + PartialEq + KeyValSet <  Key = Index,  Val = Matrix::Coefficient   >,
        Index:                  Clone + Debug + Eq,
    
{
    type Coefficient            =   Matrix::Coefficient;
    type Column                 =   IntoIter< Matrix::RowEntry >; // THIS IS THE DATA STRUCTURE OBTAINED BY CALLING `.into_iter()` on a Rust Vec
    type ColumnEntry            =   Matrix::RowEntry;
    type ColumnIndex            =   Matrix::RowIndex;
    type ColumnReverse          =   RowOfInverseOfUpperTriangularMatrix<  // column i of M^{-1} equals row i of (M antitransposed)^{-1}
                                        OrderAntiTranspose< Matrix >, 
                                    >;
    type Row                    =   RowOfInverseOfUpperTriangularMatrix< 
                                        Matrix, 
                                    >;
    type RowEntry               =   Matrix::RowEntry;
    type RowIndex               =   Matrix::RowIndex;
    type RowReverse             =   IntoIter< Matrix::RowEntry >;


    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.matrix_to_invert.has_column_for_index(index)
    }
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        self.matrix_to_invert.has_row_for_index(index)
    }
    /// Unlike [InverseUpperTriangularMatrix::column_reverse], we can't generate the entries of a column
    /// of the inverse in order; therefore to run this computation we generate all entries in reverse order,
    /// store them in a `Vec`, and reverse the `Vec`.
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {        
        let mut vec: Vec<_>     =   self.column_reverse( index ).collect(); // collect entries into a Rust Vec
        vec.reverse(); // reverse the order of entries in the Vec
        vec.into_iter()
    }
    /// Calculates entries in order via back-substitution (this operation is lazy; it does not calculate all entries at once).
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
        let antitransopose      =   OrderAntiTranspose::new( self.matrix_to_invert.clone() );
        row_of_inverse_of_upper_triangular_matrix(  // define an iterator that calculates the row in lazy fasion
            index, 
            antitransopose, 
        )        
    }
    fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
        row_of_inverse_of_upper_triangular_matrix(  // define an iterator that calculates the row in lazy fasion
            index, 
            self.matrix_to_invert.clone(), 
        )           
    }
    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
        let mut row: Vec<_>     =   self.row( index ).collect();
        row.reverse();
        row.into_iter()
    }
  

    /// Find the nonzero entry at row `row` and column `column`
    /// 
    /// This is a relatively expensive calcuation. It works by computing elements in the given row one at a time,
    /// until reaching the desired column index.
    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        let row                 =   self.row( row );
        let order_operator      =   self.matrix_to_invert.order_operator_for_row_indices();

        for entry in row {
            match order_operator.judge_cmp( & entry.key(), column ) {
                Ordering::Less => { continue }
                Ordering::Equal => { return Some( entry.val() ) }
                Ordering::Greater => { return None }
            }
        }
        return None
    }

}






impl < Matrix, Index, Entry >

    MatrixAlgebra for 

    InverseUpperTriangularMatrix
        < Matrix > 

    where 
        Matrix:                 Clone + MatrixAlgebra<
                                            RowIndex=Index, 
                                            ColumnIndex=Index, 
                                            RowEntry=Entry, 
                                            ColumnEntry=Entry,
                                            RingOperator: DivisionRingOperations,
                                        >,
        Entry:                  Clone + Debug + PartialEq + KeyValSet <  Key = Index,  Val = Matrix::Coefficient   >,
        Index:                  Clone + Debug + Eq,                                       
    
{

    type OrderOperatorForColumnEntries      =   Matrix::OrderOperatorForColumnEntries;
    type OrderOperatorForColumnIndices      =   Matrix::OrderOperatorForColumnIndices;
    type OrderOperatorForRowEntries         =   Matrix::OrderOperatorForRowEntries;
    type OrderOperatorForRowIndices         =   Matrix::OrderOperatorForRowIndices;    
    type RingOperator                       =   Matrix::RingOperator;

    /// Returns the order operator for the column entries of the matrix `M` which we are inverting
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.matrix_to_invert.order_operator_for_column_entries()
    }
    /// Returns the order operator for the column indices of the matrix `M` which we are inverting    
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.matrix_to_invert.order_operator_for_column_indices()
    }
    /// Returns the order operator for the row entries of the matrix `M` which we are inverting    
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.matrix_to_invert.order_operator_for_row_entries()        
    }
    /// Returns the order operator for the row indices of the matrix `M` which we are inverting        
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.matrix_to_invert.order_operator_for_row_indices()        
    }    
    /// Returns the ring operator of the matrix `M` which we are inverting
    fn ring_operator( &self ) -> Self::RingOperator {
        self.matrix_to_invert.ring_operator()
    }

}












#[cfg(test)]
mod doc_test_drafts {
    use crate::{algebra::matrices::types::packet::MatrixAlgebraPacket, utilities::order::{OrderOperatorAuto, OrderOperatorByKey}};

    


    #[test]
    fn test_inverse_small() {
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::operations::invert::InverseUpperTriangularMatrix;
        use crate::algebra::matrices::display::print_indexed_rows;

        // Define the ring operator (prime field of order 1049)
        let modulus = 1049;
        let ring_operator  =   PrimeOrderField::new(modulus);  

        // Define a row-major upper triangular amtrix
        let matrix  =   VecOfVec::new(
                                                        vec![   
                                                            vec![ (0,1), (1,1), ], 
                                                            vec![        (1,1), ],
                                                        ],
        ).ok().unwrap();         

        let matrix_packet       =   MatrixAlgebraPacket::with_default_order(&matrix, ring_operator);

        // Define the inverse
        let inverse = InverseUpperTriangularMatrix::new( & matrix_packet, );  

        // Print the inverse
        let row_indices = vec![ 0, 1 ];
        print_indexed_rows( &inverse, row_indices );
        
        // This should print the following:
        //
        // $ row 0: [(0, 1), (1, 1048)]
        // $ row 1: [(1, 1)] 
        // 
        // This is the correct solution, since the inverse of `matrix` is
        //
        //      |  1   -1  |
        //      |  0    1  |
        //
        // and -1 = 1048, modulo 1049
    }

}





















#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use itertools::Itertools;
    use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    use crate::algebra::matrices::types::product::ProductMatrix;
    use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
    use crate::algebra::matrices::query::MatrixOracle;

    // Invert upper triangular matrices
    // =====================================================================================================================

    /// [this is a subroutine / a helper function] 
    /// Computes the inverse of the given matrix, then checks that the product of the two matrices is identity.
    #[cfg(test)] // because this is a helper function, we mark it with the macro shown left to prevent "warning: unused function" from showing up
    fn test_inversion_of_input< Coefficient, RingOperator >( 
                matrix:         & VecOfVec< usize, Coefficient, >, 
                ring_operator:  RingOperator, 
                matrix_size:    usize,
            ) 
        where   Coefficient:        Clone + PartialEq + std::fmt::Debug,
                RingOperator:       DivisionRingOperations< Element = Coefficient > + Clone,
    {
        use crate::{algebra::matrices::{debug::matrix_order_operators_are_internally_consistent, types::packet::MatrixAlgebraPacket}, utilities::order::{OrderOperatorAuto, OrderOperatorByKey}};

        let matrix_packet       =   MatrixAlgebraPacket{
                                        matrix:                             & matrix,
                                        ring_operator,
                                        order_operator_for_column_entries:  OrderOperatorByKey,
                                        order_operator_for_column_indices:  OrderOperatorAuto,
                                        order_operator_for_row_entries:     OrderOperatorByKey,
                                        order_operator_for_row_indices:     OrderOperatorAuto,
                                    };

        // compute the inverse
        let inverse =   InverseUpperTriangularMatrix::new( & matrix_packet ); 

        // check that the inverse is internally consistent
        assert!(
            crate::algebra::matrices::debug::matrix_oracle_is_internally_consistent(
                matrix_packet.clone(),
                0..matrix_size, // iterator that runs over all row indices in order
                0..matrix_size, // iterator that runs over all column indices in order
            )
        );

        // check that matrix entries are properly sorted
        // (in fact this is unnecessary, because the VecOfVec data structure is protected such that it should be impossible to create an improperly ordered matrix)        
        assert!(
            matrix_order_operators_are_internally_consistent(
                matrix_packet.clone(),
                0..matrix_size, 
                0..matrix_size,
            ).is_ok()
        );

        // product of the putative inverse with the original matrix                                
        let product =   ProductMatrix::new( 
                                                                                    matrix_packet.clone(), 
                                                                                    & inverse, 
                                                                                );
        
        // function to extract a row and format as a rust Vec
        let row_vec = | p: usize | 
                                                    product.row( &p )
                                                    .collect_vec();

        // get a copy of the number one                      
        let one = <RingOperator as crate::algebra::rings::traits::SemiringOperations>::one();

        // check that every row of the product looks like the corresponding row of the identity matrix                                            
        for row_index in 0 .. matrix_size {
            assert_eq!( vec![ ( row_index, one.clone() ) ] , row_vec( row_index ) );
        }                                                                
    }


    
    #[test]
    fn test_inversion_of_specific_matrices() {
        use num::rational::Ratio;        
        use crate::algebra::rings::types::native::RingOperatorForNativeRustNumberType;

        // Define the integer type for the rational numbers we'll use
        type IntegerType = isize;
        let modulus = 1049;

        // Define shorthand functions for creating new rational numbers
        let r = Ratio::new;    
        let q = Ratio::from_integer;   
        
        // Define the ring operators
        let ring_operator_q  =   RingOperatorForNativeRustNumberType::< Ratio<IntegerType> >::new();  
        let ring_operator_f  =   RingOperatorForNativeRustNumberType::< f64 >::new();          
        let ring_operator_p  =   PrimeOrderField::new(modulus);                  

        let matrix  =   VecOfVec::new(
                                                        vec![   
                                                            vec![ (0,1.), (1,1.), ], 
                                                            vec![         (1,1.), ],
                                                        ],
        ).ok().unwrap();
        test_inversion_of_input( & matrix, ring_operator_f, 2 );        

        let matrix  =   VecOfVec::new(
                            vec![   
                                    vec![ (0,q(1)),  (1,q(-1)), (2,q(3))       ], 
                                    vec![            (1,q(-1)), (2,q(4))       ],
                                    vec![                       (2,q(5))       ],                                    
                                ],
        ).ok().unwrap();     
        test_inversion_of_input( & matrix, ring_operator_q, 3 );        
        
        let matrix  =   VecOfVec::new(
                            vec![   
                                    vec![ (0,r(4,1)),  (1,r(-6,1))                 ], 
                                    vec![              (1,r(4,1)),   (2,r(-2,1))   ],
                                    vec![                            (2,r(7,1))    ],                                    
                                ],
        ).ok().unwrap();        
        test_inversion_of_input( & matrix, ring_operator_q, 3 );           

        let matrix  =   VecOfVec::new(
                                vec![   
                                    vec![   (0,r(2,1)),     (1,r(6,1)),                     (3,r(8,1)),  ], 
                                    vec![                   (1,r(2,1)),     (2,r(9,1)),                  ],
                                    vec![                                   (2,r(2,1)),     (3,r(-4,1)), ],  
                                    vec![                                                   (3,r(4,1)),  ],                                  
                                ],
        ).ok().unwrap();     
        test_inversion_of_input( & matrix, ring_operator_q, 4 );         
        
        // MUCH LARGER randomized matrix to invert, with nonzero diagonal entries
        // use rand::Rng;        // we use this module to generate random elements
        // let matrix_size =   20;

        // let mut rng = rand::thread_rng(); // this generates random integers
        // let mut vec_of_vec = vec![];
        // for row_index in 0 .. matrix_size {
        //     let coefficient_leading         =   rng.gen_range( 1 .. modulus );
        //     let mut new_vec     =   vec![ (row_index, coefficient_leading) ]; // start with a nonzero diagonal element            
        //     for q in row_index+1 .. matrix_size { // fill out the rest of this row of the matrix
        //         let coefficient   = rng.gen_range( 0 .. modulus );
        //         let flag = rng.gen_range(0usize .. 3usize);
        //         if      flag == 0 { new_vec.push( ( q, 0 )           ) }
        //         else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
        //         else              { continue }
        //     }
        //     vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
        // }

        // let matrix  =   VecOfVec::new( vec_of_vec ); // formally wrap the matrix in a VecOfVec struct

        let matrix_size = 20;
        let matrix = VecOfVec::random_mod_p_upper_unitriangular( matrix_size, modulus );
        test_inversion_of_input( & matrix, ring_operator_p, matrix_size );                 

    }
}

