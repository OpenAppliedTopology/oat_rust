//! Solve `Ax = b` for `x`, where `A` is a upper triangular with nonzero diagonal entries.
//! 
//! Concretely, `A` is a matrix oracle; the major and column indices of `A` must have the same type.
//! 
//! # Alternate solutions
//! 
//! You can also solve the exact same problem, `Ax = b`, using the [echelon](crate::algebra::matrices::operations::solve::echelon) module.  That module doesn't assume that your matrix is triangular, so you have to use a function that "matches" the row and column indices of diagonal elements; for us, that function can just be the [IdentityFunction].



use crate::algebra::matrices::query::{MatrixAlgebra};
use crate::algebra::matrices::types::transpose::OrderAntiTranspose;
use crate::algebra::rings::traits::{ RingOperations, DivisionRingOperations, };
use crate::utilities::iterators::merge::hit::{IteratorsMergedInSortedOrder, hit_bulk_insert, hit_merge_by_predicate};
use crate::algebra::vectors::entries::{KeyValSet, KeyValGet};
use crate::algebra::vectors::operations::{Simplify, Scale, VectorOperations};
use crate::utilities::iterators::general::{TwoTypeIterator};
use crate::utilities::order::is_sorted_strictly;

use std::fmt::Debug;
use std::iter::Iterator;

use derive_getters::{Getters, Dissolve};


// TRIANGULAR SOLUTION -- ASCEND
// ---------------------------------------------------------------------------

/// Solve `xA = b` where `A` is  is upper-triangular and row-major
/// 
/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, **in strictly ascending order**.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::{ vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket };
/// use oat_rust::algebra::matrices::operations::solve::triangle::TriangularSolveForRowVector;       
/// use crate::oat_rust::algebra::matrices::operations::MatrixOracleOperations;
/// use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
/// use oat_rust::utilities::order::{ OrderOperatorAuto, OrderOperatorByKey };
/// 
/// 
/// // Define an upper triangular matrix
/// let matrix_data     =   VecOfVec::new(
///                            vec![                                     
///                                vec![ (0,true),  (1,true), (2,true)       ], 
///                                vec![            (1,true), (2,true)       ],
///                                vec![                      (2,true)       ],                                       
///                            ],
///                         ).ok().unwrap(); 
///                     
/// // Wrap the matrix in a MatrixAlgebraPacket, which determines how to perform ring operations
/// // and compare order of entries. Here we use the two-element field, and our order operators
/// // just impliment the ordinary order on integers.
/// let matrix          =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients(&matrix_data);
///                     
/// // SOLVE xA = b FOR x
/// // --------------------------
///                     
/// // define a sparse vector b
/// let b = vec![ (0,true), (1,true) ];        
///                     
///                     
/// // create a solver to solve xA = b for x
/// let solution        =   TriangularSolveForRowVector::solve(
///                             b.iter().cloned(),   // the solver requires a bona-fide iterator, not an iterable
///                             & matrix,            // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
///                         ).ok().unwrap();         // unwrap the value from the enclosing Result struct
///                 
/// // check the solution, i.e. check that xA = b
/// let product: Vec<_> =   matrix.multiply_with_row_vector( solution ).collect();
/// assert_eq!( product, b );        
/// ```
/// 
#[derive(Debug, Clone, PartialEq, Getters, Dissolve)]
pub struct TriangularSolveForRowVector< 
                    ProblemVector,
                    Matrix,
                > 
    where 
        ProblemVector:                      IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra,
        Matrix::RowEntry:                   KeyValSet,      

{
    matrix:                                 Matrix, // the `A` in `Ax = b`
    entries_to_eliminate_simplified_heap:   Simplify<
                                              IteratorsMergedInSortedOrder< 
                                                      // the IteratorsMergedInSortedOrder struct merges a collection of iterators into a single iterator; the type of those iterators is given below (they are essentially scaled vectors)
                                                      Scale< 
                                                              TwoTypeIterator< // this enum allows us to treat two different iterator types as a single type
                                                                      Matrix::Row,  // the iterators returned by the matrix                                                                
                                                                      ProblemVector::IntoIter,
                                                                      // Matrix::RowEntry,
                                                                  >,
                                                              Matrix::RingOperator, // the ring_operator operator
                                                          >,
                                                      // the thing that declares whether one row index comes before of after another    
                                                      Matrix::OrderOperatorForRowEntries
                                                  >,
                                                  Matrix::RingOperator,
                                            >,   
} 






impl    < 
            ProblemVector,
            Matrix,
            Index,        
        > 

    TriangularSolveForRowVector< 
            ProblemVector,
            Matrix,                        
        > 

    where 
        ProblemVector:                      IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra< 
                                                RowIndex        =   Index, 
                                                ColumnIndex     =   Index,
                                                RingOperator:       RingOperations,
                                            >,
        Matrix::RowEntry:                   KeyValSet,
        Index:                              PartialEq,          
 
{        
    /// Solve `Ax = b` for `x`
    /// 
    /// # Requirements
    /// 
    /// - For all `r`, the leading nonzero entry row `r` of the matrix `A` must form `(r,a)`. That is, the leading entry has index `r`. (The leading entry does not have to be a tuple; it can be any equivalent [key-value pair](algebra::vectors::entries::KeyValGet))
    /// - The matrix `A` must have a row for every index `r` such that either (i) the vector  `b` has an entry of form `(r,a)`, or (ii) `A` contains has a structural nonzero entry whose column index is `x`
    /// - `b` must iterate over entries in strictly ascending order
    /// 
    ///  See the [TriangularSolveForRowVector] documentation for an example.
    pub fn solve(
                    b:                  ProblemVector,                
                    a:                  Matrix,
                )
            ->
            Result<
                TriangularSolveForRowVector< 
                    Vec< Matrix::RowEntry >,
                    Matrix,      
                >,
                (),
            > 

    {

        // verify that the input vector is sorted
        let b: Vec<_> = b.into_iter().collect(); // collect vector elements so that we can check the order
        let order_operator = (&a).order_operator_for_row_entries();
        if ! is_sorted_strictly(&b, & order_operator) {
            println!("\nThe input vector to `TriangularSolveForRowVector::solve( input_vector, matrix ) is not sorted in strictly ascending order by index.");
            println!("\nThe vector is:\n{:?}\n", &b);
            println!("\nRecall that order may be determined by a customized order operator.\n");
            println!("\nThis error message was generated by OAT.\n");            
            return Err(())            
        }



        // solve the problem

        let ring_operator       =   a.ring_operator();
        let order_operator      =   a.order_operator_for_row_entries();

        let b   =   TwoTypeIterator::Version2( b.into_iter() );   // wrap b in an TwoTypeIterator enum so that it can be easily combined with  other vecgors
        let b   =   b.scale_by( ring_operator.minus_one(), ring_operator.clone() ); // reverse the sign of b, because back-substitution acutally solves for Ax = -b.

        let initial_list_of_vectors =std::iter::once( b );
        let entries_to_eliminate_simplified_heap
                =   hit_merge_by_predicate(  // this merge operation merges a list of vectors together without actually adding their entries; it is the first step in addign them; next we will add togehter entries with equal indices, using the simplify method
                            initial_list_of_vectors,
                            order_operator,
                        )
                        .simplify( ring_operator.clone() );

        Ok( TriangularSolveForRowVector{
            matrix:                             a,
            entries_to_eliminate_simplified_heap, 
        } )
    }     
}






impl    < 
            ProblemVector,
            Matrix,
            Index,        
        > 

        Iterator for    

        TriangularSolveForRowVector< 
                ProblemVector,
                Matrix,       
            > 

    where 
        Index:                              PartialEq,
        ProblemVector:                      IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra< 
                                                RowIndex            =   Index, 
                                                ColumnIndex         =   Index,
                                                RingOperator:           DivisionRingOperations,
                                            >,
        Matrix::RowEntry:                   KeyValSet,  
        Index:                              Debug,                  

{
    type Item = Matrix::RowEntry;

    fn next( &mut self ) -> Option<Self::Item> { 

        // THIS ITERATOR PROCEEDS BY ELIMINATING ENTRIES IN `self.entries_to_eliminate_simplified_heap`

        //  IF THIS ITERATOR IS NONEMPTY THEN WE WILL GENERATE THE NEXT ENTRY OF THE SOLUTION VECTOR; 
        //  note that `entry_to_eliminate` will always have a nonzero coefficient, because `self.entries_to_eliminate_simplified_heap` is a struct of type `Simplify< .. >`; structs of this type only return nonzero entries.
        if let Some( mut entry_to_eliminate )       =   self.entries_to_eliminate_simplified_heap.next() {

            let ring_operator   =   self.matrix.ring_operator();            

            // TO ELIMINATE `entry_to_eliminate`, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX; DENOTE THE UNSCLAED VECTOR BY v
            let mut seed_of_eliminating_iterator    =   self.matrix.row( & entry_to_eliminate.key() ).into_iter();
            
            // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
            let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

            // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
            let scale_factor     =    ring_operator.negate(
                                                    ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                );
            
            // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
            // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_eliminate_simplified_heap
            let eliminating_iterator       
                    =   TwoTypeIterator::Version1( seed_of_eliminating_iterator ) 
                            .scale_by( scale_factor.clone(), ring_operator.clone() );

            // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_eliminate_simplified_heap`
            hit_bulk_insert( &mut self.entries_to_eliminate_simplified_heap.unsimplified, vec![eliminating_iterator] ); 

            // NOW RE-USE `entry_to_eliminate` AS A CONTAINER FOR THE NEXT ENTRY IN THE INVERSE MATRIX
            entry_to_eliminate.set_val( scale_factor );

            Some( entry_to_eliminate ) // recall that we have re-purposed `entry_to_eliminate` to hold the output of this function

        } else {
            None
        }
    }
}  








// TRIANGULAR SOLUTION -- DESCEND
// ---------------------------------------------------------------------------


/// Solve `Ax = b` where `A` is upper-triangular and row-major
/// 
/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, **in reverse order**.
/// For example, if `x = [ (0,1), (1,2), (2,3) ]` then the iterator would return `(2,3)`, then `(1,2)`, then `(0,1)`.
/// 
/// # Why `reverse`?
/// 
/// The 'reverse' in [TriangularSolveForColumnVectorReverse] comes from the fact that entries appear in reverse order.
/// This ordering comes from the fact that back-stubstitution for a column vector with an upper triangular matrix
/// naturally generates entries in reverse order (we are not aware of any method that generates them in forward order).
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::{ vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket };
/// use oat_rust::algebra::matrices::operations::solve::triangle::TriangularSolveForColumnVectorReverse;       
/// use crate::oat_rust::algebra::matrices::operations::MatrixOracleOperations;
/// use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
/// use oat_rust::utilities::order::{ OrderOperatorAuto, OrderOperatorByKey };
/// 
/// // Define an upper triangular matrix
/// let matrix_data     =   VecOfVec::new(
///                            vec![                                     
///                                vec![ (0,true),  (1,true), (2,true)       ], 
///                                vec![            (1,true), (2,true)       ],
///                                vec![                      (2,true)       ],                                       
///                            ],
///                         ).ok().unwrap(); 
///                     
/// // Wrap the matrix in a MatrixAlgebraPacket, which determines how to perform ring operations
/// // and compare order of entries. Here we use the two-element field, and our order operators
/// // just impliment the ordinary order on integers.
/// let matrix          =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients(&matrix_data);
///                     
/// // SOLVE Ax = b FOR x
/// // --------------------------
///                     
/// // define a sparse vector b
/// let b = vec![ (2,true), (0,true) ];        
///                     
/// // create a solver to solve xA = b for x
/// let solution        =   TriangularSolveForColumnVectorReverse::solve(
///                             b.iter().cloned(),  // the solver requires a bona-fide iterator, not an iterable
///                             & matrix,           // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
///                         ).ok().unwrap();        // unwrap the value from the wrapper Result struct
///                        
/// // check the solution, i.e. check that xA = b
/// let product: Vec<_> =   matrix.multiply_with_column_vector_reverse(solution).collect(); // the "reverse" means that we get entries in reverse order
/// assert_eq!( product, b );    
/// ```
// #[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Getters, Dissolve)]
pub struct TriangularSolveForColumnVectorReverse< 
                    ProblemVector,
                    Matrix,
                > 
    where 
        ProblemVector:                      IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra,
        Matrix::ColumnEntry:                KeyValSet,                  

{
    antitranspose_solution: TriangularSolveForRowVector<
                                ProblemVector,
                                OrderAntiTranspose< Matrix >,
                            >,
} 

// custom implementation because #[derive(..)] was throwing errors
impl<ProblemVector, Matrix> 
    Clone for     
    TriangularSolveForColumnVectorReverse
        < ProblemVector, Matrix >    
    where
        TriangularSolveForRowVector< ProblemVector, OrderAntiTranspose< Matrix > >:  Clone,
        
        ProblemVector:                      IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra,
        Matrix::ColumnEntry:                KeyValSet,    
{
    fn clone(&self) -> Self {
        TriangularSolveForColumnVectorReverse{ antitranspose_solution: self.antitranspose_solution.clone() }        
    }
}

// custom implementation because #[derive(..)] was throwing errors
impl<ProblemVector, Matrix> 
    PartialEq for     
    TriangularSolveForColumnVectorReverse
        < ProblemVector, Matrix >    
    where
        TriangularSolveForRowVector< ProblemVector, OrderAntiTranspose< Matrix > >:  PartialEq,

        ProblemVector:                      IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra,
        Matrix::ColumnEntry:                KeyValSet,            
{
    fn eq(&self, other: &Self) -> bool {
        self.antitranspose_solution == other.antitranspose_solution
    }
}

// custom implementation because #[derive(..)] was throwing errors
impl<ProblemVector, Matrix> 
    Eq for     
    TriangularSolveForColumnVectorReverse
        < ProblemVector, Matrix >    
    where
        TriangularSolveForRowVector< ProblemVector, OrderAntiTranspose< Matrix > >:  Eq,
        
        ProblemVector:                      IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra,
        Matrix::ColumnEntry:                KeyValSet,            
{}

// custom implementation because #[derive(..)] was throwing errors
impl<ProblemVector, Matrix> 
    PartialOrd for     
    TriangularSolveForColumnVectorReverse
        < ProblemVector, Matrix >    
    where
        TriangularSolveForRowVector< ProblemVector, OrderAntiTranspose< Matrix > >:  PartialOrd,

        ProblemVector:                      IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra,
        Matrix::ColumnEntry:                KeyValSet,            
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.antitranspose_solution.partial_cmp(&other.antitranspose_solution)
    }
}

// custom implementation because #[derive(..)] was throwing errors
impl<ProblemVector, Matrix> 
    Ord for     
    TriangularSolveForColumnVectorReverse
        < ProblemVector, Matrix >    
    where
        TriangularSolveForRowVector< ProblemVector, OrderAntiTranspose< Matrix > >:  Ord,
        
        ProblemVector:                      IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra,
        Matrix::ColumnEntry:                KeyValSet,            
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.antitranspose_solution.cmp(&other.antitranspose_solution)
    }
}




impl    < 
            ProblemVector,
            Matrix,
            Index,        
        > 

    TriangularSolveForColumnVectorReverse< 
            ProblemVector,
            Matrix,       
        > 

    where 
        ProblemVector:                      IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra< 
                                                RowIndex            =   Index, 
                                                ColumnIndex         =   Index,
                                                RingOperator:           RingOperations,
                                            >,
        Matrix::ColumnEntry:                KeyValSet,
        Index:                              PartialEq,
 
{        
    /// Solve `Ax = b` for `x`
    /// 
    /// # Requirements
    /// 
    /// - `A.column_reverse( x )` returns an iterator whose first entry has index `x` and a nonzero coefficient.
    /// - `b` must iterate over entries in strictly **descending** order of row index.
    /// 
    ///  See the [TriangularSolveForColumnVectorReverse] for an example.
    pub fn solve(
                    b:              ProblemVector,                
                    a:              Matrix,
                )
            ->
            Result< 
                TriangularSolveForColumnVectorReverse< 
                    Vec<Matrix::ColumnEntry>, // we convert the ProblemVector into a Vec in order to check that it is sorted
                    Matrix,      
                >,
                (),
            >

    {

        // check that b has entries sorted in reverse order
        let b: Vec<_> = b.into_iter().collect(); // collect vector elements so that we can check the order
        let order_operator = (&a).order_operator_for_column_entries_reverse();

        if ! is_sorted_strictly( &b, &order_operator) {
            println!("Error: the entries of b must be sorted in strictly descending order of row index.");
            return Err( () );
        }


        

        // solve the problem
        let antitranspose = OrderAntiTranspose::new(a);
        match TriangularSolveForRowVector::solve( b, antitranspose, ) {
            Err( () ) => {
                println!("\n\nThe preceding error message originated when a call to `TriangularSolveForColumnVectorReverse::solve( input_vector, matrix )`");
                println!("made a secondary call to `TriangularSolveForRowVector::solve( input_vector, antitranspose_of_matrix )`.");
                println!("So the entries of the vector should appear in the REVERSE of the order that the matrix imposes on its row indices.");
                return Err( () );
            },
            Ok( antitranspose_solution ) => { 
                return Ok( TriangularSolveForColumnVectorReverse{ antitranspose_solution } )
            }
        }
        
    }     
}






impl    < 
            ProblemVector,
            Matrix,  
            Index,                
        > 

        Iterator for    

        TriangularSolveForColumnVectorReverse< 
                ProblemVector,
                Matrix,         
            > 

    where 
        TriangularSolveForRowVector< ProblemVector, OrderAntiTranspose< Matrix >, >: Iterator< Item = Matrix::ColumnEntry >,
    
        ProblemVector:                      IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             MatrixAlgebra< RowIndex=Index, ColumnIndex=Index, RingOperator: DivisionRingOperations >,
        Matrix::ColumnEntry:                KeyValSet,     

{
    type Item = Matrix::ColumnEntry;

    fn next( &mut self ) -> Option<Self::Item> { 
        self.antitranspose_solution.next()
    }
}  









//  DOCSTRING TESTS
//  ===========================================================================






#[cfg(test)]
mod doctring_tests {
    use crate::algebra::matrices::operations::MatrixOracleOperations;



    #[test]
    fn test_docstring_solve_row() {
        use crate::algebra::matrices::types::{ vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket };
        use crate::algebra::matrices::operations::solve::triangle::TriangularSolveForRowVector;       
        use crate::algebra::rings::types::field_prime_order::BooleanField;
        use crate::utilities::order::{ OrderOperatorAuto, OrderOperatorByKey };


        // Define an upper triangular matrix
        let matrix_data     =   VecOfVec::new(
                                   vec![                                     
                                       vec![ (0,true),  (1,true), (2,true)       ], 
                                       vec![            (1,true), (2,true)       ],
                                       vec![                      (2,true)       ],                                       
                                   ],
                                ).ok().unwrap(); 
                            
        // Wrap the matrix in a MatrixAlgebraPacket, which determines how to perform ring operations
        // and compare order of entries. Here we use the two-element field, and our order operators
        // just impliment the ordinary order on integers.
        let matrix          =   MatrixAlgebraPacket{
                                    matrix: & matrix_data,
                                    ring_operator: BooleanField,
                                    order_operator_for_row_entries: OrderOperatorByKey, // this order operator compares (integer, coefficient) pairs according to the default order on integers
                                    order_operator_for_row_indices: OrderOperatorAuto,
                                    order_operator_for_column_entries: OrderOperatorByKey,  // this order operator compares (integer, coefficient) pairs according to the default order on integers
                                    order_operator_for_column_indices: OrderOperatorAuto,
                                };
                            
        // SOLVE xA = b FOR x
        // --------------------------
                            
        // define a sparse vector b
        let b = vec![ (0,true), (1,true) ];        
                            
                            
        // create a solver to solve xA = b for x
        let solution =  TriangularSolveForRowVector::solve(
                                b.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
                                & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                            ).ok().unwrap();
                        
        // check the solution, i.e. check that xA = b
        let product: Vec<_> = matrix.multiply_with_row_vector( solution ).collect();
        assert_eq!( product, b );      
    }


    #[test]
    fn test_docstring_solve_column() {
        use crate::algebra::matrices::types::{ vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket };
        use crate::algebra::matrices::operations::solve::triangle::TriangularSolveForColumnVectorReverse;       
        use crate::algebra::rings::types::field_prime_order::BooleanField;
        use crate::utilities::order::{ OrderOperatorAuto, OrderOperatorByKey };

        // Define an upper triangular matrix
        let matrix_data     =   VecOfVec::new(
                                   vec![                                     
                                       vec![ (0,true),  (1,true), (2,true)       ], 
                                       vec![            (1,true), (2,true)       ],
                                       vec![                      (2,true)       ],                                       
                                   ],
                                ).ok().unwrap(); 
                            
        // Wrap the matrix in a MatrixAlgebraPacket, which determines how to perform ring operations
        // and compare order of entries. Here we use the two-element field, and our order operators
        // just impliment the ordinary order on integers.
        let matrix          =   MatrixAlgebraPacket{
                                    matrix: & matrix_data,
                                    ring_operator: BooleanField,
                                    order_operator_for_row_entries: OrderOperatorByKey, // this order operator compares (integer, coefficient) pairs according to the default order on integers
                                    order_operator_for_row_indices: OrderOperatorAuto,
                                    order_operator_for_column_entries: OrderOperatorByKey,  // this order operator compares (integer, coefficient) pairs according to the default order on integers
                                    order_operator_for_column_indices: OrderOperatorAuto,
                                };
                            
        // SOLVE Ax = b FOR x
        // --------------------------
                            
        // define a sparse vector b
        let b = vec![ (2,true), (0,true) ];        
                            
        // create a solver to solve xA = b for x
        let solution =  TriangularSolveForColumnVectorReverse::solve(
                                b.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
                                & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                            ).ok().unwrap();
                        
        // check the solution, i.e. check that xA = b
        let product: Vec<_> = matrix.multiply_with_column_vector_reverse(solution).collect(); // the "reverse" means that we get entries in reverse order
        assert_eq!( product, b );    
    }

}





















//  TESTS
//  ===========================================================================


#[cfg(test)]
mod tests {

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use itertools::Itertools;
    use crate::{algebra::matrices::types::vec_of_vec::sorted::VecOfVec, algebra::rings::types::field_prime_order::PrimeOrderField};


    // Invert upper triangular matrices
    // =====================================================================================================================

    /// [this is a subroutine / a helper function] 
    /// Computes the inverse of the given matrix, then checks that the product of the two matrices is identity.
    #[cfg(test)] // because this is a helper function, we mark it with the macro shown left to prevent "warning: unused function" from showing up
    fn test_triangular_solve< RingOperator >( 
                matrix:         & VecOfVec< usize, RingOperator::Element >, 
                ring_operator:  RingOperator, 
                matrix_size:    usize 
            ) 
        where   RingOperator::Element:                      PartialEq + Ord + std::fmt::Debug,
                RingOperator:                               DivisionRingOperations + Clone,
    {
        use crate::{algebra::matrices::types::packet::MatrixAlgebraPacket, utilities::order::{OrderOperatorAuto, OrderOperatorByKey}};

        
        let matrix          =   MatrixAlgebraPacket{
                                    matrix,
                                    ring_operator:                      ring_operator.clone(),
                                    order_operator_for_row_entries:     OrderOperatorByKey, // this order operator compares (integer, coefficient) pairs according to the default order on integers
                                    order_operator_for_row_indices:     OrderOperatorAuto,
                                    order_operator_for_column_entries:  OrderOperatorByKey,  // this order operator compares (integer, coefficient) pairs according to the default order on integers
                                    order_operator_for_column_indices:  OrderOperatorAuto,
                                };

        for num_nonzeros in [0, 1, 2, 3, matrix_size].iter().map(|x| std::cmp::min(x, &matrix_size) ) {
            // iterate over all vectors with a given number of nonzeros
            for vector_support in (0 .. matrix_size).combinations( *num_nonzeros ) {                
                let mut vec = vector_support
                                .iter()
                                .cloned()
                                .map(|x| ( x, RingOperator::one() ))
                                .collect_vec();
                vec.sort();

                // solve
                let solution = TriangularSolveForRowVector::solve(
                                        vec.iter().cloned(),                    
                                        & matrix,
                                    )
                                    .ok()
                                    .unwrap()
                                    .peekable()
                                    .simplify( ring_operator.clone() )
                                    .collect_vec();
                assert!(
                    vec.into_iter()
                        .eq( solution.multiply_self_as_a_row_vector_with_matrix( & matrix ) )
                );
            }
        }                                                               
    }


    
    #[test]
    fn test_triangular_solve_on_specific_matrices() {
        use num::rational::Ratio;        
        // use crate::algebra::matrices::query::MajorDimension;
        use crate::algebra::rings::types::native::RingOperatorForNativeRustNumberType;

        // Define the integer type for the rational numbers we'll use
        type IntegerType = isize;
        let modulus = 1049;

        // Define shorthand functions for creating new rational numbers
        let r = Ratio::new;    
        let q = Ratio::from_integer;   
        
        // Define the ring operators
        let ring_operator_q  =   RingOperatorForNativeRustNumberType::< Ratio<IntegerType> >::new();  
        let ring_operator_p  =   PrimeOrderField::new(modulus);                  

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(1)), ], 
                                    vec![            (1,q(1)), ],
                                ],
                        ).ok().unwrap();
        test_triangular_solve( &matrix, ring_operator_q, 2 );

        // these lines result in an ERROR becuase the constructor `VecOfVec::new` panics when it receives a vector of vectors where an internal vector is not sorted in *strictly* asecneding order
        // let matrix  =   VecOfVec::new(
        //                         vec![ 
        //                             vec![ (0,q(1)), (0,q(1)),   (1,q(1)), (1,q(1)), ], 
        //                             vec![                       (1,q(1)), (1,q(1)), ],
        //                         ],
        // );
        // test_triangular_solve( &matrix, ring_operator_q.clone(), 2 );        

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(-1)), (2,q(3))       ], 
                                    vec![            (1,q(-1)), (2,q(4))       ],
                                    vec![                       (2,q(5))       ],                                    
                                ],
                        ).ok().unwrap();     
        test_triangular_solve( &matrix, ring_operator_q, 3 );        
        
        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,r(4,1)),  (1,r(-6,1))                 ], 
                                    vec![              (1,r(4,1)),   (2,r(-2,1))   ],
                                    vec![                            (2,r(7,1))    ],                                    
                                ],
                        ).ok().unwrap();        
        test_triangular_solve( &matrix, ring_operator_q, 3 );           

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![   (0,r(2,1)),     (1,r(6,1)),                     (3,r(8,1)),  ], 
                                    vec![                   (1,r(2,1)),     (2,r(9,1)),                  ],
                                    vec![                                   (2,r(2,1)),     (3,r(-4,1)), ],  
                                    vec![                                                   (3,r(4,1)),  ],                                  
                                ],
                        ).ok().unwrap();     
        test_triangular_solve( &matrix, ring_operator_q, 4 );         
        
        //  MUCH LARGER randomized matrix to invert, with nonzero diagonal entries
        //  ----------------------------------------------------------------------

        let matrix_size =   20;
        let matrix = VecOfVec::random_mod_p_upper_triangular_invertible( matrix_size, modulus );

        test_triangular_solve( &matrix, ring_operator_p, matrix_size );                 

    }
}