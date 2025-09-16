//! Tools used internally to compute a U-match factorization



use serde_json::map::Entry;

use super::*;




use crate::algebra::matrices::operations::transform_vector_wise::PutbackIteratorMatrix;
use crate::algebra::vectors::entries::{KeyValSet, KeyValNew};

use crate::algebra::matrices::types::matching::{GeneralizedMatchingMatrixWithSequentialOrder};


use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
use crate::algebra::matrices::query::{MatrixOracle, MatrixAlgebra};


use crate::utilities::iterators::general::{SkipUntil, IterWrappedVec};
use crate::utilities::iterators::merge::hit::{IteratorsMergedInSortedOrder, hit_bulk_insert};
use crate::algebra::vectors::entries::{KeyValPair};
use crate::algebra::rings::traits::{SemiringOperations, RingOperations, DivisionRingOperations};

use crate::utilities::order::{JudgePartialOrder, ReverseOrder};


use crate::algebra::vectors::operations::{Scale, VectorOperations, Simplify, OnlyIndicesInsideCollection, LinearCombinationSimplified};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;
use std::iter::Peekable;
use std::vec::IntoIter;








//  HELPER FUNCTION 
//  ----------------------------------------------------------------------------------

/// Helper function for `try_writing_next_row_of_target_comb_to_buffer`.  That function
/// involves multiplying a row of a target COMB by a factored matrix; sometimes entries
/// at the beginning of the resulting row vector can be deleted without penalty; the 
/// current function performs this deletion on a (scalar multiple) of an individual row of the 
/// factored matrix.
fn target_comb_entry_to_scaled_truncated_row_of_matrix_to_factor<
        MatrixToFactor,
    > 
    (
        target_comb_inv_entry:              ( usize, MatrixToFactor::Coefficient ),
        scale_factor:                       MatrixToFactor::Coefficient,        
        truncation_limit:                   & MatrixToFactor::RowEntry,
        matrix_to_factor:                    & MatrixToFactor,     
        matching:                           & GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient >,        
    )
    ->
    Peekable<
            Scale < 
                    MatrixToFactor::Row,  // a row of the factored matrix
                    MatrixToFactor::RingOperator, 
                >, 
        >
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
                                    RowIndex:                       Hash, // required for the hashing performed by the generalized matching array
                                    RowEntry:                       KeyValPair, 
                                    ColumnEntry:                    KeyValPair,                                    
                                    RingOperator:                   DivisionRingOperations,           
                                >,                
{
    let ring_operator           =   matrix_to_factor.ring_operator();
    let order_operator_for_row_entries          =   matrix_to_factor.order_operator_for_row_entries();

    let ( target_comb_inv_entry_index, target_comb_inv_entry_coeff ) = target_comb_inv_entry;

    // we multiply the entry coefficient `target_comb_inv_entry.val()` by the `scale_factor`, 
    let scale_factor_aggregate = ring_operator.multiply( scale_factor, target_comb_inv_entry_coeff ); 
    
    // then we multiply `target_comb_inv_entry.val() * scale_factor_aggregate` by the row of the **factored matrix**
    // (not matching array) that is indexed by `target_comb_inv_entry.key()`
    
    matrix_to_factor // <-- note: matrix_to_factor, **not matching**
                    .row( 
                            // this picks out the row index corresponding to the entry; this index might *not* be an integer, though `target_comb_inv_entry.key()` *is* an integer
                            & matching.row_index_for_ordinal( target_comb_inv_entry_index )  
                        )
                    .into_iter()
                    .scale_by( scale_factor_aggregate, ring_operator )
                    // making the iterator peekable will help us prune some entries
                    .peekable() 
                    // we discard the entries of this iterator that have index <= `leading_entry_to_eliminate.key()` (because those entries would eventually be cancelled out, anyway -- this saves possibly a lot of heap insert + pop operations)
                    // in doing so, we must take care not to throw out too many entries by accident; that is why peekable is used                    
                    .skip_until( |x| {
                        let this = &order_operator_for_row_entries; this.judge_partial_cmp(truncation_limit, x) == Some(Ordering::Less) }  ) 
}


//  HELPER FUNCTION 
//  ----------------------------------------------------------------------------------


/// A row of an attempted solution to Ax = -b, calculated via cohomology
/// described in "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// 
/// The output of this function is best explained in terms of an example.  Let M be a matrix,
/// and suppose that we have partially executed the reduction algorithm described in 
/// "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek to 
/// obtain a U-match factorization.  This algorithm proceeds by reducing rows iteratively, 
/// from bottom to top.  To process a new row, the reduction procees by adding a sequence of
/// linear combinations of lower rows.  The addition of each linear combination eliminates a
/// new leading entry in the vector to be reduced.  The output of this reduction algorithm is a
/// sequence of index-coefficient pairs related to the sequence of linear combinations which we use to
/// reduce the vector.  The `RowOfUmatchReducingVector` returns this sequence of 
/// entries in order, as an interator, assuming that 
/// - the `matrix_to_factor` field stores matrix `M` in row-major form
/// - the portions of the matching array and inverse-codomain-COMB are stored in the fields
/// `matching` and `array_target_comb_inv_off_diag`.
/// - the field `order_operator_for_row_entries` correctly determines which column index precedes which 
/// 
/// # How to identify pivot indices
/// 
/// This struct modifies one of its private fields, the object `entries_to_eliminate_simplified_heap`.
/// The modified iterator represents the result of reducing `entries_to_eliminate_simplified_heap` via
/// clearing operations, as described in "U-match factorization, ..." by Hang, Ziegelmeier, Giusti,
/// and Henselman-Petrusek.
/// - if calling `entries_to_eliminate_simplified_heap.next()` returns `Some( column_index, coeff )` after
///  `RowOfUmatchReducingVector` returns its last item,
///   then (assuming the struct is initialized correctly vis-a-vis the reduction algorithm
///   in  "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek) 
///   `column_index` is a pivot index, as is the corresponding `row_index` (that is, the row index of
///   the factored matrix, which indexes the vector that initialized `entries_to_eliminate_simplified_heap`).
/// - otherwise, `entries_to_eliminate_simplified_heap.next()` returns `None`, and we do not have a
///   new pivot index.
/// 


/// Needs new docs.
/// 
/// Note: entries are not inserted into the buffer vector in sorted order, and this function
/// does not sort the vector before returning it.
fn try_writing_next_row_of_target_comb_to_buffer< 
        MatrixToFactor,
    > 
    (
        target_comb_inv_off_diag_row_buffer:      &mut Vec< ( usize, MatrixToFactor::Coefficient ) >,          
        entries_to_eliminate_simplified_heap:    & mut Simplify <
                                                        IteratorsMergedInSortedOrder < 
                                                                Peekable<
                                                                        // the IteratorsMergedInSortedOrder struct merges a collection of iterators into a single iterator; the type of those iterators is given below (they are essentially scaled vectors)
                                                                        Scale < 
                                                                                MatrixToFactor::Row,  // a row of the factored matrix
                                                                                MatrixToFactor::RingOperator, // the ring_operator operator
                                                                            >, 
                                                                    >,
                                                                // the thing that declares whether one row index comes before of after another    
                                                                MatrixToFactor::OrderOperatorForRowEntries                                                                
                                                            >,
                                                        MatrixToFactor::RingOperator,
                                                    >,          
        matrix_to_factor:                      & MatrixToFactor,    
        matching:                     & GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient >,
        array_target_comb_inv_off_diag:   & Vec< Vec< (usize, MatrixToFactor::Coefficient) > >,  // the off-diagonal part of the (block indexed by pivot indices of) inverse of the target COMB       
    ) 
    ->

    Option< MatrixToFactor::RowEntry >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
                                    RowIndex:                       Hash, // required for the hashing performed by the generalized matching array
                                    RowEntry:                       KeyValPair, 
                                    ColumnEntry:                    KeyValPair,                                    
                                    RingOperator:                   DivisionRingOperations,           
                                >,                     
{
    let ring_operator           =   matrix_to_factor.ring_operator();

    // we use this while loop because a for-loop results in a borrowing conflict (the conflict arrises because the for loop 
    // iteratos over `entries_to_eliminate_simplified_heap`, but we modify `entries_to_eliminate_simplified_heap` within the for-loop)
    while let Some( leading_entry_to_eliminate ) = entries_to_eliminate_simplified_heap.next() {

        // if entries_to_eliminate_simplified_heap.unsimplified.len() > 10 { println!("styx[[{:?}]]", entries_to_eliminate_simplified_heap.unsimplified.len()) }
    
        // IF THE MINOR INDEX OF THE ENTRY IS MATCHED, THEN WE CAN ELIMINATE
        if let Some( ordinal )       =   matching.ordinal_for_column_index( & leading_entry_to_eliminate.key() ) {

            let scale_factor                =   ring_operator.negate(
                                                                ring_operator.divide(
                                                                        leading_entry_to_eliminate.val(),                                                                            
                                                                        matching.coefficient_for_ordinal( ordinal ),
                                                                    )
                                                            );
            
            // add the new (scaled and truncated) iterators to `entries_to_eliminate_simplified_heap`
            hit_bulk_insert( 
                    &mut entries_to_eliminate_simplified_heap.unsimplified,  
                    array_target_comb_inv_off_diag[ ordinal ]
                        .iter()
                        .cloned()
                        .chain(  // append a diagonal entry with coefficient 1 to the iterator
                                std::iter::once( 
                                        ( ordinal, MatrixToFactor::RingOperator::one() ) 
                                    ) 
                            )
                        .map(   
                                | target_comb_inv_entry | 
                                {
                                // println!(
                                //         "target_comb_entry_to_scaled_truncated_row_of_matrix_to_factor:      {:?}",
                                //         target_comb_entry_to_scaled_truncated_row_of_matrix_to_factor(
                                //                 target_comb_inv_entry.clone(),
                                //                 scale_factor.clone(),
                                //                 & leading_entry_to_eliminate, // truncation_limit:                   MatrixToFactor::RowEntry,
                                //                 matrix_to_factor, //                     & MatrixToFactor,     
                                //                 matching, //                    & GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient >,        
                                //                 ring_operator.clone(), //                    RingOperator, // the operator for the coefficient ring_operator    
                                //                 order_operator_for_row_entries.clone() //               OrderOperatorForRowEntries,                                      
                                //             )
                                //             .collect_vec()
                                //     );

                                target_comb_entry_to_scaled_truncated_row_of_matrix_to_factor(
                                        target_comb_inv_entry,
                                        scale_factor.clone(),
                                        & leading_entry_to_eliminate, // truncation_limit:                   MatrixToFactor::RowEntry,
                                        matrix_to_factor, //                     & MatrixToFactor,     
                                        matching, //                    & GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient >,        
                                    )
                                }
                            )
                );            
            // if entries_to_eliminate_simplified_heap.unsimplified.len() > 10 { println!("styx2[[{:?}]]", entries_to_eliminate_simplified_heap.unsimplified.len()) }                                         

            // push ( scale_factor * eliminating_row_with_diagonal_entry_re-added ) to the buffer; we do this in two steps
            // step 1: push the diagonal entry == (the scale factor and the index of the target COMB row that we used to perform the elimination) to the buffer
            target_comb_inv_off_diag_row_buffer.push(   ( ordinal, scale_factor.clone() )   );
            // step 1: push off-diagonal entries
            target_comb_inv_off_diag_row_buffer
                    .extend(  
                            array_target_comb_inv_off_diag[ ordinal ]
                                .iter()
                                .cloned() // note: we have to do this b/c `target_comb_inv_off_diag_row_buffer` is a `VecOfVec`, not a `VecOfVec`
                                .scale_by( scale_factor, ring_operator.clone() )
                        );            
            // if entries_to_eliminate_simplified_heap.unsimplified.len() > 10 { println!("styxdone[[{:?}]]", entries_to_eliminate_simplified_heap.unsimplified.len()) }                                         
        } else {

            // REMARK: in this case we cannot clear `leading_entry_to_eliminate` via elementary operations; 
            // therefore `leading_entry_to_eliminate` is a matched (i.e. pivot) entry

            return Some( leading_entry_to_eliminate );
        }   
    }
    None    
}





//  FUNCTION(S) TO EXTRACT Ripiv (the pivot portion of T^{-1}[pivot_indices, pivot_indices])
//  ------------------------------------------------------------------------------------------

/// Returns the block of the target COMB indexed by pivot indices, and a copy of the matching matrix.
/// 
/// For details on this factorization, see [this preprint](https://arxiv.org/pdf/2108.08831.pdf).
/// 
/// - This U-match factorization is computed via the standard "cohomology algorithm." The idea of this algorithm is
///   perform Gauss-Jordan elimination, adding multiples of lower rows to higher rows in order to clear leading entries.
/// 
/// - The `iter_row_index` argument must iterate over rows indices in **REVERSE ORDER**, according to 
///   `matrix_to_factor.order_operator_for_row_indices()`. This is necessary for the low-to-high operations to
///   execute correctly.
/// 
/// This function returns a `GeneralizedMatchingMatrixWithSequentialOrder` data structure, which encodes the generalized matching matrix
/// for the U-match decomposition. This data structure also encodes bijections `{0 .. N} <--> {matched row indices}`
/// and `{0 .. N} <--> {matched column indices }`
///   - These bijections allow us to place the matched row (respectively, column) indices into sequences `r0 .. rN` and `c0 .. cN`.
///   - The pairs `(r0,c0) .. (rN,cN)` encode the locations of the nonzero entries in the generalized matching matrix.
///   - The sequence `r0 < .. < rN` is strictly sorted in ascending order, according to the order operator on row indices
///     provided by `matrix_to_factor.order_operator_on_row_indices()`.
///   - The sequence `c0 .. cN` is **not sorted according to the order operator on column indices** provided by `matrix_to_factor.order_operator_on_row_indices()`.
///     (nor can it be, in general, if we require `r0 < .. < rN` to be sorted and require `(r0,c0) .. (rN,cN)` to be the locations
///     of the nonzero elements of `M`)  
pub fn get_pivot_block_of_target_comb_inverse_with_deleted_diagonal< MatrixToFactor, IterRowIndex >
    ( 
            matrix_to_factor:                       &MatrixToFactor,  // use of a reference poses no problem for future calls to this function, since output has no references        
            iter_row_index:                         IterRowIndex,
    ) 
    -> 
    ( 
        VecOfVec< usize, MatrixToFactor::Coefficient >,          
        GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient >,
    )
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
                                    RowIndex:                       Hash, // required for the hashing performed by the generalized matching array
                                    RowEntry:                       KeyValPair, 
                                    ColumnEntry:                    KeyValPair,                                    
                                    RingOperator:                   DivisionRingOperations,           
                                >,                    
        IterRowIndex:                                   Iterator< Item = MatrixToFactor::RowIndex >,

{
    
    // Initialize some objects
    let ring_operator   =   matrix_to_factor.ring_operator();
    let order_operator_for_row_entries      =   matrix_to_factor.order_operator_for_row_entries();

    let mut entries_to_eliminate_simplified_heap    =   IteratorsMergedInSortedOrder::new( order_operator_for_row_entries.clone() )
                                                                            .simplify( ring_operator.clone() ); // an empty (simplified) merger of iterators
    let mut matching          =   GeneralizedMatchingMatrixWithSequentialOrder::new(); // an all-zero generalized matching array
    let mut array_target_comb_inv_off_diag: Vec< Vec< ( usize, MatrixToFactor::Coefficient ) > >  =   Vec::new();    
    let mut target_comb_inv_off_diag_row_buffer   =   Vec::new();
    let mut target_comb_inv_off_diag_row_buffer_simplified   =   Vec::new();    

    // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE row_index APPEAR IN STRICTLY ASCENDING ORDER
    // set up a variable to track the last row index produced by `iter_row_index`; we will use this to ensure that the row indices returned by `iter_row_index` appear in strictly ascending order
    // let mut prior_row_index_opt = None;

    // let mut counter = 0;
    // build the (pivot block of the) target COMB row by row
    for row_index in iter_row_index {

        // counter +=1;
        // print!("row({:?})", counter);

        // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE row_index APPEAR IN STRICTLY ASCENDING ORDER
        // check that this row index is strictly greater than the last
        // if let Some( ref prior_keyaj ) = prior_row_index_opt {
        //     match order_operator_of_row_index.judge_lt( prior_keyaj, &row_index ) {
        //         true => {  },
        //         false => { panic!("Umatch factorization cannot be completed because the sequence of row indices is not strictly ascending") }
        //     }               
        // }
        // prior_row_index_opt.replace( row_index.clone() );

        // clear the collection of entries to eliminate
        entries_to_eliminate_simplified_heap.unsimplified.clear();        

        // insert the sequence of entries in row `row_index`
        entries_to_eliminate_simplified_heap.unsimplified.insert_one_iter(
                matrix_to_factor.row( & row_index )
                    .scale_by( MatrixToFactor::RingOperator::one(), ring_operator.clone() ) // !!! might be more efficient to use a two-type iter than to perform this multiplication, which only serves to make the iterator compatible with the IteratorsMergedInSortedOrder struct
                    .peekable()
        );

        target_comb_inv_off_diag_row_buffer.clear();
        let leading_entry_uneliminable_opt =
        try_writing_next_row_of_target_comb_to_buffer(
                & mut target_comb_inv_off_diag_row_buffer,
                & mut entries_to_eliminate_simplified_heap,        
                  matrix_to_factor,
                & matching,
                & array_target_comb_inv_off_diag,
            );

        match leading_entry_uneliminable_opt {
            // in this case we do not have a pivot
            Some( leading_entry_uneliminable )    => { 
                let pivot_column_index        =   leading_entry_uneliminable.key();
                let pivot_coeff         =   leading_entry_uneliminable.val();
                // println!("leading_entry_uneliminable {:?}", &leading_entry_uneliminable);              
                // update the matching array
                matching.push( row_index, pivot_column_index, pivot_coeff ); 
                // sort the buffer vector
                target_comb_inv_off_diag_row_buffer.sort_by( |a, b| b.0.cmp( &a.0 ) ); // note this yields DESCENDING ORDER -- which is what we want, because the order of indices is currently inverted (the greatest row_index gets the lowest ordinal, and the least row_index gets the highest ordinal; we will correct this later)
                // println!("THE BUFFER:      {:?}", &target_comb_inv_off_diag_row_buffer);
                // simplify the sequence of entries in the buffer (note that the buffer itself could have entries with duplicate indices)
                // !!! hypothetically there's a opportunity for a small optimization here, where we simultaneously simplify and sort at the same time; this would avoid a few sorting operations
                let iter_to_push     =   target_comb_inv_off_diag_row_buffer
                                                            .iter()
                                                            .cloned()
                                                            .peekable()
                                                            .simplify( ring_operator.clone() );
                
                // write the simplified sequence to a buffer vector (this avoids having to reallocate, later, since we can't know the length of the simplified sequence 'till we iterate over it)
                target_comb_inv_off_diag_row_buffer_simplified.clear();
                target_comb_inv_off_diag_row_buffer_simplified.extend( iter_to_push );
                
                // write the simplified sequence to a new vector of exactly the right length
                let num_entries = target_comb_inv_off_diag_row_buffer_simplified.len();
                let mut vec_to_push     =   Vec::with_capacity( num_entries );
                vec_to_push.extend( target_comb_inv_off_diag_row_buffer_simplified.drain( 0..num_entries ) );
                
                // update the vector that stores the (off diagonal entries of the) "pivot" part of the inverse of the target COMB
                array_target_comb_inv_off_diag.push( vec_to_push );
            }, 

            // in this case we don't have a pivot entry, so there is nothing to do
            None   => {}
        }
    }
    
    // if there are no pairs, return empty matrices (otherwise we will need to perform some transformations)
    if array_target_comb_inv_off_diag.is_empty() { return ( VecOfVec::new(array_target_comb_inv_off_diag).ok().unwrap(), matching ) }

    // remove excess capacity
    array_target_comb_inv_off_diag.shrink_to_fit();
    

    // reverse the order of rows (because they are currently inverted, since we started with the bottom row and worked up)
    ( &mut array_target_comb_inv_off_diag ).reverse();

    // invert the ordinal used to index each entry
    let num_matched_pairs_minus_one = array_target_comb_inv_off_diag.len() - 1;
    for row_vec in array_target_comb_inv_off_diag.iter_mut() {  
        for entry in row_vec.iter_mut() {
            entry.set_key( num_matched_pairs_minus_one - entry.key() )
        }
    }

    // invert the ordinals of entries in the matching array
    matching.reverse_order_of_matches();

    // return the (off-diagonal entries of the pivot block of the) target COMB
    ( VecOfVec::new(array_target_comb_inv_off_diag).ok().unwrap(), matching )

}










/// Returns (i) the block of the target COMB indexed by pivot indices, and
/// (ii)) the matching matrix; output is correct ONLY if pivot column indices and pivot row indices are disjoint.
/// 
/// For details on this factorization, see [this preprint](https://arxiv.org/pdf/2108.08831.pdf).
/// 
/// This U-match factorization is computed via the standard "cohomology algorithm,"
/// using the clear/compress/twist optimization.  This optimization skips over rows of the matrix indexed by keys
/// that already index pivot columns.
/// 
/// 
/// This function returns a `GeneralizedMatchingMatrixWithSequentialOrder` data structure, which encodes the generalized matching matrix
/// for the U-match decomposition. This data structure also encodes bijections `{0 .. N} <--> {matched row indices}`
/// and `{0 .. N} <--> {matched column indices }`
///   - These bijections allow us to place the matched row (respectively, column) indices into sequences `r0 .. rN` and `c0 .. cN`.
///   - The pairs `(r0,c0) .. (rN,cN)` encode the locations of the nonzero entries in the generalized matching matrix.
///   - The sequence `r0 < .. < rN` is strictly sorted in ascending order, according to the order operator on row indices
///     provided by `matrix_to_factor.order_operator_on_row_indices()`.
///   - The sequence `c0 .. cN` is **not sorted according to the order operator on column indices** provided by `matrix_to_factor.order_operator_on_row_indices()`.
///     (nor can it be, in general, if we require `r0 < .. < rN` to be sorted and require `(r0,c0) .. (rN,cN)` to be the locations
///     of the nonzero elements of `M`)  
pub fn target_comb_inv_off_diag_pivot_block_skipmatched< MatrixToFactor, IterRowIndex, IndexForRowsAndColumns, EntryForRowsAndColumns, Coefficient >
    ( 
            matrix_to_factor:                       &MatrixToFactor,  // use of a reference poses no problem for future calls to this function, since output has no references        
            iter_row_index:                         IterRowIndex,
    ) 
    -> 
    ( 
        VecOfVec< usize, MatrixToFactor::Coefficient >, 
        GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient >,
    )
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
        IterRowIndex:               IntoIterator< Item = MatrixToFactor::RowIndex >,                    
        Coefficient:                Clone + Debug + PartialEq,

{
    let ring_operator                                           =   matrix_to_factor.ring_operator();
    let order_operator_for_row_entries                          =   matrix_to_factor.order_operator_for_row_entries();
    let order_operator_for_row_indices                             =   matrix_to_factor.order_operator_for_row_indices();
    let matrix_to_factor_with_putback                                        =   PutbackIteratorMatrix::new( matrix_to_factor );
    
    // Initialize some objects
    let mut entries_to_eliminate_simplified_heap                =   IteratorsMergedInSortedOrder::new( order_operator_for_row_entries.clone() )
                                                                            .simplify( ring_operator.clone() ); // an empty (simplified) merger of iterators
    let mut matching          =   GeneralizedMatchingMatrixWithSequentialOrder::new(); // an all-zero generalized matching array
    let mut array_target_comb_inv_off_diag: Vec< Vec< ( usize, MatrixToFactor::Coefficient ) > >  =   Vec::new();    
    let mut target_comb_inv_off_diag_row_buffer                =   Vec::new();
    let mut target_comb_inv_off_diag_row_buffer_simplified     =   Vec::new();    

    // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE row_index APPEAR IN STRICTLY ASCENDING ORDER
    // set up a variable to track the last row index produced by `iter_row_index`; we will use this to ensure that the row indices returned by `iter_row_index` appear in strictly ascending order
    // let mut prior_row_index_opt = None;

    // let mut counter = 0;
    // let pb = indicatif::ProgressBar::new(100000);
    let mut sc_counter = 0;

    // build the (pivot block of the) target COMB row by row
    for row_index in iter_row_index {

    
        if matching.has_a_match_for_column_index( & row_index ) { 
            continue }      
          

        // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE row_index APPEAR IN STRICTLY ASCENDING ORDER
        // check that this row index is strictly greater than the last
        // !! NB: It's hard to use this with the ordinary order operator on weighted simplices, because
        //        people usually want that operator to put lower dimensional simplices first, but in row
        //        reduction operations we want to visit lower dimensional simplices first, which means putting them
        //        last in the filtration.
        // if let Some( ref prior_row_index ) = prior_row_index_opt {
        //     match order_operator_for_row_indices.judge_le( prior_row_index, &row_index ) {
        //         false => { 
        //             // do nothing; this what we want
        //          },
        //         true => { 
        //             panic!(
        //                 "Umatch factorization cannot be completed because the provided sequence of row indices \
        //                 is not sorted in strictly descending order; one index is {:?} and the next is {:?}.",
        //                 & prior_row_index,
        //                 & row_index,
        //             )                     
        //         }
        //     }               
        // }
        // prior_row_index_opt.replace( row_index.clone() );

        // clear the collection of entries to eliminate
        entries_to_eliminate_simplified_heap.unsimplified.clear();        

        // load the next row of the matrix to factor
        let mut next_iter = matrix_to_factor_with_putback.row( & row_index );        
        
        
        // attempt to short circuit
        if let Some( leading_entry ) = next_iter.next() {
            // if possible, short circuit
            if matrix_to_factor
                .bottom_entry_for_column( &leading_entry.key() )
                .as_ref() 
                == Some(&leading_entry) {

                
                matching.push( row_index, leading_entry.key(), leading_entry.val() ); // push entry to the matching matrix
                array_target_comb_inv_off_diag.push( vec![] ); // record the fact that there are no off-diagonal entries for this row in the inverse target COMB
                sc_counter += 1; // increment the short circuit counter
                continue
                
            // otherwise put the leading entry back
            } else {
                next_iter.put_back( leading_entry );
            }
        }
        


        // insert the sequence of entries in row `row_index`
        entries_to_eliminate_simplified_heap.unsimplified.insert_one_iter(
                next_iter
                    .into_iter()
                    .scale_by( MatrixToFactor::RingOperator::one(), ring_operator.clone() ) // !!! might be more efficient to use a two-type iter than to perform this multiplication, which only serves to make the iterator compatible with the IteratorsMergedInSortedOrder struct
                    .peekable()
        );

        target_comb_inv_off_diag_row_buffer.clear();

      
        let leading_entry_uneliminable_opt =
            try_writing_next_row_of_target_comb_to_buffer(
                    & mut target_comb_inv_off_diag_row_buffer,
                    & mut entries_to_eliminate_simplified_heap,        
                    & matrix_to_factor_with_putback,
                    & matching,
                    & array_target_comb_inv_off_diag,
                );


        match leading_entry_uneliminable_opt {
            // in this case we do not have a pivot
            Some( leading_entry_uneliminable )    => { 
                let pivot_column_index        =   leading_entry_uneliminable.key();
                let pivot_coeff         =   leading_entry_uneliminable.val();
                // println!("leading_entry_uneliminable {:?}", &leading_entry_uneliminable);              
                // update the matching array
                matching.push( row_index, pivot_column_index, pivot_coeff ); 
                // sort the buffer vector
                target_comb_inv_off_diag_row_buffer.sort_by( |a, b| b.0.cmp( &a.0 ) ); // note this yields DESCENDING ORDER -- which is what we want, because the order of indices is currently inverted (the greatest row_index gets the lowest ordinal, and the least row_index gets the highest ordinal; we will correct this later)
                // println!("THE BUFFER:      {:?}", &target_comb_inv_off_diag_row_buffer);
                // simplify the sequence of entries in the buffer (note that the buffer itself could have entries with duplicate indices)
                // !!! hypothetically there's a opportunity for a small optimization here, where we simultaneously simplify and sort at the same time; this would avoid a few sorting operations
                let iter_to_push     =   target_comb_inv_off_diag_row_buffer
                                                            .iter()
                                                            .cloned()
                                                            .peekable()
                                                            .simplify( ring_operator.clone() );
                
                // write the simplified sequence to a buffer vector (this avoids having to reallocate, later, since we can't know the length of the simplified sequence 'till we iterate over it)
                target_comb_inv_off_diag_row_buffer_simplified.clear();
                target_comb_inv_off_diag_row_buffer_simplified.extend( iter_to_push );
                
                // write the simplified sequence to a new vector of exactly the right length
                let num_entries = target_comb_inv_off_diag_row_buffer_simplified.len();
                let mut vec_to_push     =   Vec::with_capacity( num_entries );
                vec_to_push.extend( target_comb_inv_off_diag_row_buffer_simplified.drain( 0..num_entries ) );
                
                // update the vector that stores the (off diagonal entries of the) "pivot" part of the inverse of the target COMB
                array_target_comb_inv_off_diag.push( vec_to_push );
            }, 

            // in this case we don't have a pivot entry, so there is nothing to do
            None   => {}
        }       
    }
    
    // if there are no pairs, return empty matrices (otherwise we will need to perform some transformations)
    if array_target_comb_inv_off_diag.is_empty() { return ( VecOfVec::new(array_target_comb_inv_off_diag).ok().unwrap(), matching ) }

    // remove excess capacity
    array_target_comb_inv_off_diag.shrink_to_fit();
    

    // reverse the order of rows (because they are currently inverted, since we started with the bottom row and worked up)
    ( &mut array_target_comb_inv_off_diag ).reverse();

    // invert the ordinal used to index each entry
    let num_matched_pairs_minus_one = array_target_comb_inv_off_diag.len() - 1;
    for row_vec in array_target_comb_inv_off_diag.iter_mut() {  
        for entry in row_vec.iter_mut() {
            entry.set_key( num_matched_pairs_minus_one - entry.key() )
        }
    }

    // invert the ordinals of entries in the matching array
    matching.reverse_order_of_matches();

    // println!("pareto pairs: {:?}", sc_counter );
    // println!("number of matched pairs: {:?}", matching.number_of_structural_nonzeros() );    

    // return the (off-diagonal entries of the pivot block of the) target COMB
    ( VecOfVec::new(array_target_comb_inv_off_diag).ok().unwrap(), matching )

}











//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- ROWS INDEXED BY COLUMN_INDEX, COLUMNS BY COLUMN_INDEX
//  =========================================================================================================





/// Represents the matrix A =  (R_{ρρ})^{−1}D_{ρκ} defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// Concretely, this matrix equals (the pivot block of the inverse of the target COMB) * (the pivot block of the factored matrix).
/// 
/// This struct equals the matched block of the inverse source COMB up to reindexing and rescaling rows. In particular:
/// 
/// 1) Mathematically, if S_{κκ} is the matched part of the
/// source COMB, M_{ρκ} is the matched part of the matching matrix, M_{ρκ}^{-1} is the inverse of M_{ρκ}
/// obtained by transposing and inverting the nonzero entries, and D_{ρκ} is the matched part of the factored matrix, then `S_{κκ}^{-1} = M_{ρκ}^{-1}A`
/// 
/// 2) Note in particular that the matrix A has rows indexed by matched column indices of the row indices, rather than 
/// by the matched row indices themselves. Correspondingly, the column lookup commands return iterators that iterate over entries in ascending or descending order 
/// of *column* index.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex< 'a, MatrixToFactor > 
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
                                    RowIndex:                       Hash, // required for the hashing performed by the generalized matching array
                                    RowEntry:                       KeyValPair, 
                                    ColumnEntry:                    KeyValPair,                                    
                                    RingOperator:                   DivisionRingOperations,           
                                >,                    
{
    pub umatch:    &'a Umatch< MatrixToFactor > ,  
}

impl < 'a, MatrixToFactor > 
    
    TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex
        < 'a, MatrixToFactor > 
    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
                                    RowIndex:                       Hash, // required for the hashing performed by the generalized matching array
                                    RowEntry:                       KeyValPair, 
                                    ColumnEntry:                    KeyValPair,                                    
                                    RingOperator:                   DivisionRingOperations,           
                                >,                      
{
    // Make a new [`SourceCombInverseMatchedBlockRowsIndexedByColumnIndex`].
    pub fn new( umatch_ref: &'a Umatch< MatrixToFactor >  ) -> Self {
        TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex{ umatch: umatch_ref }
    }
}


// Implement `Eq` if matrix coefficients implement `Eq`.
// (in fact, Eq can be implemented for this struct if and only if matrix coefficients implement Eq)
impl < 'a, MatrixToFactor, > 

    Eq for
    
    TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex
        < 'a, MatrixToFactor > 
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





impl < 'a, MatrixToFactor >
    
    MatrixOracle for 
    
    TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
                                    RowIndex:                       Hash, // required for the hashing performed by the generalized matching array
                                    RowEntry:                       KeyValPair, 
                                    ColumnEntry:                    KeyValPair,                                    
                                    RingOperator:                   DivisionRingOperations,           
                                >,        

{   
    type Coefficient            =   MatrixToFactor::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex               =   MatrixToFactor::ColumnIndex;       // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   MatrixToFactor::ColumnIndex;       // The type of column indices

    type RowEntry               =   MatrixToFactor::RowEntry   ;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   MatrixToFactor::RowEntry   ;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                    =   LinearCombinationSimplified< 
                                        OnlyIndicesInsideCollection< MatrixToFactor::Row, &'a HashMap<MatrixToFactor::ColumnIndex, usize> >, 
                                        MatrixToFactor::RingOperator, 
                                        MatrixToFactor::OrderOperatorForRowEntries 
                                    >;  // What you get when you ask for a row.
    type RowReverse             =   IntoIter< MatrixToFactor::RowEntry >;        // What you get when you ask for a row with the order of entries reversed
    type Column                 =   IntoIter< MatrixToFactor::RowEntry >;            // What you get when you ask for a column   
    type ColumnReverse          =   IterWrappedVec< MatrixToFactor::RowEntry >;     // What you get when you ask for a column with the order of entries reversed                                

    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { 
            // our strategy is just to lookup the entry in another matrix which is essentially the same, up to re-indexing
            let row_ordinal = self.umatch.generalized_matching_matrix_ref().ordinal_for_column_index( row ).unwrap();
            let  matched_blocks_of_target_comb_inverse_times_matrix_to_factor_oc = self.umatch.matched_blocks_of_target_comb_inverse_times_matrix_to_factor_oc();
            return matched_blocks_of_target_comb_inverse_times_matrix_to_factor_oc.structural_nonzero_entry( &row_ordinal, column );        
        }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool 
        { self.umatch.generalized_matching_matrix_ref().has_a_match_for_column_index( index ) }
    fn has_row_for_index(  &   self, index: & Self::RowIndex)   -> bool 
        { self.umatch.generalized_matching_matrix_ref().has_a_match_for_column_index( index ) }
    
    fn row(                     & self,  index: & MatrixToFactor::ColumnIndex    )       -> Self::Row
    {
        // This implementation simply invokes a row of MatchingMatrixTimesInverseSourcCombMatchedBlockIndexedByRowOrdinal
        // This works because the two matrices have the same set of rows -- they simply look up these rows via different indices 
        let ordinal = self.umatch.matching.ordinal_for_column_index( &index ).unwrap();
        let proxy 
            = self.umatch.matched_blocks_of_target_comb_inverse_times_matrix_to_factor_oc();
        proxy.row( & ordinal )
    }
    fn row_reverse(             &self,  index: &Self::RowIndex    )   -> Self::RowReverse   { 
        let mut vec = self.row(index).collect_vec(); 
        (&mut vec).reverse(); 
        vec.into_iter()
    }
    
    fn column(                  &self,  index: &Self::ColumnIndex )   -> Self::Column       { 
        let mut vec = self.column_reverse(index).collect_vec();
        (& mut vec).reverse();
        vec.into_iter() 
    }
    fn column_reverse(          &self,  index: &Self::ColumnIndex )   -> Self::ColumnReverse
    {

        // a column of the matched block of T^{-1}*D; entries of the column are indexed by integers
        let col_indexed_by_ordinal   =   self.umatch.matched_blocks_of_target_comb_inverse_times_matrix_to_factor_oc().column_reverse( & index );

        // reindex the entries with column indices
        let mut col_reindexed      =   col_indexed_by_ordinal.map( 
                                                |x| 
                                                MatrixToFactor::RowEntry::new(
                                                        self.umatch.matching.column_index_for_ordinal( x.0 ),
                                                        x.1
                                                    )
                                            )
                                            .collect_vec();
        col_reindexed.shrink_to_fit();

        // now all we have to do is sort the column according to (column index) index!

        // repackage the order comparator for entries in rows; we do this because the resulting struct has methods to produce std::comp::Ordering, which we will need to sort the column vector
        let order_decider =
                //     InferTotalOrderFromJudgePartialOrder::new( 
                            ReverseOrder::new( // we have to reverse the order because we want our output to be sorted in *descending* order
                                    self.umatch.order_operator_for_row_entries() 
                                );
                        // );
        // wrap the order_decider in a closure
        let order_decider_function = |x: &MatrixToFactor::RowEntry, y: &MatrixToFactor::RowEntry|  order_decider.judge_partial_cmp(x, y).unwrap();
        // use the closure to the sort the vector
        col_reindexed.sort_by( order_decider_function );


        IterWrappedVec::new( col_reindexed )
    }
} 













//  IMPLEMENT MATRIX OPERATORS
//  --------------------------------------------------------------------------------------------------------------






impl < 'a, MatrixToFactor >
    
    MatrixAlgebra for 
    
    TargetCombInverseTimesMatrixToFactorMatchedBlockRowsIndexedByColumnIndex
        < 'a, MatrixToFactor >

    where   
        MatrixToFactor:         MatrixAlgebra<
                                    ColumnIndex:                    Hash, // required for the hashing performed by the generalized matching array
                                    RowIndex:                       Hash, // required for the hashing performed by the generalized matching array
                                    RowEntry:                       KeyValPair, 
                                    ColumnEntry:                    KeyValPair,                                    
                                    RingOperator:                   DivisionRingOperations,           
                                >,          
{
    type OrderOperatorForColumnEntries          =   MatrixToFactor::OrderOperatorForRowEntries;
    
    type RingOperator                           =   MatrixToFactor::RingOperator;
    
    type OrderOperatorForRowEntries             =   MatrixToFactor::OrderOperatorForRowEntries;
    
    type OrderOperatorForRowIndices             =   MatrixToFactor::OrderOperatorForColumnIndices;
    
    type OrderOperatorForColumnIndices          =   MatrixToFactor::OrderOperatorForColumnIndices;
    
    fn ring_operator( &self ) -> Self::RingOperator {
        self.umatch.ring_operator()
    }
    
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.umatch.matrix_to_factor.order_operator_for_row_entries()
    }
    
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.umatch.matrix_to_factor.order_operator_for_column_indices()
    }
    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.umatch.matrix_to_factor.order_operator_for_row_entries()
    }
    
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.umatch.matrix_to_factor.order_operator_for_column_indices()
    }
}
















//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- INDEXED BY MINOR KEY ORDINAL INTEGERS
//  =========================================================================================================


//  !!! UNDER CONSTRUCTION
//  !!! AN IMPORTANT CHALLENGE IS THAT MATCHING MATRICES ARE NOT CURRENTLY DESIGNED TO STORE MINOR KEY ORDINALS



// /// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
// /// Concretely, this matrix equals (the pivot block of the inverse of the target COMB) * (the pivot block of the matching array).
// /// 
// /// This marix is indexed by **integers** (its major and column indices are `usize`).  Specifically, row `i`
// /// corresponds to `\rho_i` and column `j` corresponds to `\kappa_j`, in the notation of "U-match factorization, ...".
// pub struct SourceCombInverseMatchedBlockIndexedByColumnOrdinal< 'a, 'b, MatrixToFactor, MatrixToFactor::Row, MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient, RingOperator, OrderOperatorForRowEntries, OrderOperatorForColumnEntries > 
//     where     
//         MatrixToFactor::ColumnIndex:                     Clone + Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
//         MatrixToFactor::RowIndex:                     Clone + Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
//         MatrixToFactor::Row:            IntoIterator,
//         MatrixToFactor::RowEntry:      KeyValGet < Key = MatrixToFactor::ColumnIndex, Val = MatrixToFactor::Coefficient >,
// {
//     umatch_ref:     Umatch< MatrixToFactor > ,  
// }

// impl < 'a, 'b, MatrixToFactor, MatrixToFactor::Row, MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient, RingOperator, OrderOperatorForRowEntries, OrderOperatorForColumnEntries > 
    
//     SourceCombInverseMatchedBlockIndexedByColumnOrdinal
//         < 'a, 'b, MatrixToFactor, MatrixToFactor::Row, MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient, RingOperator, OrderOperatorForRowEntries, OrderOperatorForColumnEntries > 
//     where   
//         MatrixToFactor::ColumnIndex:                     Clone + Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
//         MatrixToFactor::RowIndex:                     Clone + Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
//         MatrixToFactor::Row:            IntoIterator,
//         MatrixToFactor::RowEntry:      KeyValGet < Key = MatrixToFactor::ColumnIndex, Val = MatrixToFactor::Coefficient >,
//         RingOperator:               Clone,
//         OrderOperatorForRowEntries:            Clone,
// {
//     // Make a new [`SourceCombInverseMatchedBlockIndexedByColumnOrdinal`].
//     fn new( umatch_ref: & Umatch< MatrixToFactor >  ) -> Self {
//         SourceCombInverseMatchedBlockIndexedByColumnOrdinal{ umatch_ref: (*umatch_ref).clone() }
//     }
// }

    


//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- INDEXED BY ROW INDEX ORDINAL INTEGERS
//  =========================================================================================================


//  !!! UNDER CONSTRUCTION
//
//  !!! ISSUE: THE OBJECT matrix_to_factor DOES NOT RETURN ENTRIES IN ASCENDING ORDER OF ROW INDEX ORDINAL; WOULD HAVE TO GENERATE ALL ENTRIES AND RETURN AS A VECTOR
