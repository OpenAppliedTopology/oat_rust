//! Analyze and debug matrices.


use crate::{algebra::rings::traits::SemiringOperations, utilities::iterators::is_sorted::IsSortedBy};
use crate::algebra::vectors::entries::{KeyValGet, KeyValNew, KeyValSet};
use std::{collections::HashMap, fmt::Debug};
use itertools::Itertools;
use crate::algebra::matrices::types::product::ProductMatrix;

use super::operations::solve::echelon::RowEchelonSolver;
use super::{query::{MatrixAlgebra, MatrixOracle, }, types::transpose::OrderAntiTranspose};



//  VERIFY THAT ROWS AND COLUMNS ARE COMPATIBLE
//  -----------------------------------------------------------------------------------------------



/// Ensures that each entry in each row also appears in the corresponding column; panics otherwise.
/// 
/// Only checks the rows returned by `iter_row_index`.
/// 
/// See comments in source code for details.
pub fn verify_rows_compatible_with_columns_helper_function< Matrix, RowIndexIterator >(
                matrix:             Matrix,
                iter_row_index:        RowIndexIterator,              
            )
    where
        Matrix:                             MatrixOracle,
        Matrix::RowEntry:                   KeyValNew,
        Matrix::ColumnEntry:                KeyValNew,
        RowIndexIterator:                   Clone + Iterator< Item = Matrix::RowIndex >,
{
    // Make sure all the entries in the row also appear in the corresponding column
    for row_index in iter_row_index {
        let row = matrix.row( & row_index );
        // for each entry in the row ..
        for row_entry in row {
            // construct the corresponding entry of the corresponding (reverse) column
            let columnreverse_entry 
                    = Matrix::ColumnEntry::new( row_index.clone(), row_entry.val() );
            // check that this entry does in fact lie in the (reverse) column
            let reversed_column = matrix.column( & row_entry.key() );
            let exists_in_column = reversed_column.into_iter().any( |x| x == columnreverse_entry );
            if ! exists_in_column {
                println!("\nerror: entry {:?} appears in row {:?} of the matrix, but entry {:?} does not appear in column {:?}.\n", 
                    (row_entry.key(), row_entry.val()), 
                    row_index.clone(),  
                    (columnreverse_entry.key(), columnreverse_entry.val()),            
                    row_entry.key(),
                );      
                println!("\nrow {:?} of the matrix is:\n{:?}\n", row_index.clone(), matrix.row( & row_index ).collect_vec() );
                println!("\ncolumn {:?} of the matrix is:\n{:?}\n", row_entry.key(),  matrix.column( &  row_entry.key() ).collect_vec(), );                
            }
            assert!( exists_in_column );
        }
    }  
}


/// Ensures that each entry in each row also appears in the corresponding column, and vice versa; panics otherwise.
/// 
/// Specifically checks that 
/// - each entry `(column_index,coefficient)` in `matrix.row(&row_index)` has a corresponding entry `(row_index,coefficient)` in `matrix.column(&column_index)`
/// - each entry `(row_index,coefficient)` in `matrix.column_reverse(&row_index)` has a corresponding entry `(column_index,coefficient)` in `matrix.row_reverse(&row_index)`
/// 
/// # Caveats
/// 
/// - This test only looks at rows and columns specified by `iter_row_index` and `iter_column_index`
/// - Does not check that 
///   - `matrix.row(&index)` returns the same entries (in reverse order) as `matrix.row_reverse(&index)`
///   - `matrix.column(&index)` returns the same entries (in reverse order) as `matrix.column_reverse(&index)`
pub fn verify_rows_compatible_with_columns< Matrix, RowIndexIterator, ColumnIndexIterator >(
                matrix:         Matrix,
                iter_column_index:    ColumnIndexIterator,                
                iter_row_index:    RowIndexIterator,
            )
    where
        Matrix:                             MatrixOracle,
        Matrix::Row:            IntoIterator,
        Matrix::RowEntry:       Debug + std::cmp::PartialEq + KeyValNew < Key = Matrix::ColumnIndex, Val = Matrix::Coefficient >,
        Matrix::ColumnReverse:           IntoIterator,
        Matrix::ColumnEntry:      Debug + std::cmp::PartialEq + KeyValNew < Key = Matrix::RowIndex, Val = Matrix::Coefficient >,
        Matrix::ColumnIndex:                     Clone + Debug,
        Matrix::RowIndex:                     Clone + Debug,
        Matrix::Coefficient:                     Debug,        
        RowIndexIterator:                         Clone + Iterator< Item = Matrix::RowIndex >,
        ColumnIndexIterator:                         Clone + Iterator< Item = Matrix::ColumnIndex >,        
{
    verify_rows_compatible_with_columns_helper_function( 
            & matrix, 
            iter_row_index,                 
        );

    verify_rows_compatible_with_columns_helper_function( 
            OrderAntiTranspose::new( matrix ), 
            iter_column_index,                           
        );        
}





//  CHECK THAT MATRIX ORACLE IS VALID
//  ==============================================================================



/// Checks that the commands to lookup rows, columns, and entries are mutually consistent
/// 
/// Specifically checks that
/// - `.row()` and `.row_reverse()` return the same entries in opposite order
/// - `.column()` and `.column_reverse()` return the same entries in opposite order
/// - `.row(r)` returns the same sequence of entries that one would obtain by iterating over column indices in sorted order, calling `.structural_nonzero_entry( r, c )` for each index `c`
/// - `.column(c)` returns the same sequence of entries that one would obtain by iterating over row indices in sorted order, calling `.structural_nonzero_entry( r, c )` for each index `r`
/// 
/// # Inputs
/// 
/// The function takes as input a vector called `sorted_row_indices`, which should contain the row indices of the matrix in sorted order,
/// and `sorted_column_indices`, which contains the column indices of the matrix in sorted order.
/// It's necessary for the user to provide these inputs because matrix oracles have no way of knowing what the set of all valid row and column indices, a priori;
/// they are simply passive structures that return information about specific rows and columns when requested by the user.
pub fn matrix_oracle_is_internally_consistent< Matrix, I, J >
    (
        matrix:                 Matrix,
        sorted_row_indices:     I,
        sorted_column_indices:  J,
    )
    ->
    bool

    where
        Matrix:     MatrixOracle,
        I:          IntoIterator< Item = Matrix::RowIndex >,
        J:          IntoIterator< Item = Matrix::ColumnIndex >,        
{

    // collect the row and column indices into vectors
    let sorted_row_indices      =   sorted_row_indices.into_iter().collect::<Vec<_>>();
    let sorted_column_indices   =   sorted_column_indices.into_iter().collect::<Vec<_>>();    

    // check that matrix.row() and matrix.structural_nonzero_entry() agree.
    for row_index in sorted_row_indices.iter().cloned() {
        let row_from_entries: Vec<_>     =   sorted_column_indices
                        .iter()
                        .cloned()
                        .filter_map(
                            |column_index| 
                            {
                                let coefficient     =   matrix.structural_nonzero_entry( & row_index, & column_index );
                                // if coefficient is Some(a), then return (column_index, a)
                                coefficient.map( |a|  (column_index, a) )
                            }                            
                        ) // filters out None's
                        .collect();
        let row_from_iter: Vec<_>         =   matrix.row( & row_index ).collect();

        if row_from_iter.len() != row_from_entries.len() {
            println!(
                "\nInvalid matrix: calling\n\n   matrix.row({:?})\n\nreturns the following sequence of structural nonzero entries:\n",
                &row_index
            );
                println!("    < begin list of entries >");              
            for entry in &row_from_iter {                      
                println!("    {:?}", entry);                    
            }
                println!("    < end   list of entries >");                
            println!(
                "\nHowever, the sequence of structural nonzero entries obtained by\n\n    matrix.structural_nonzero_entry({:?}, column_index)\n\nfor all values of `column_index` provided by the user is not equal:",
                &row_index
            );
                println!("    < begin list of entries >");              
            for entry in &row_from_entries {                                      
                println!("    {:?}", entry);                  
            }
                println!("    < end   list of entries >");                  
            return false
        }        

        for ( entry_number, entry ) in row_from_iter.iter().cloned().enumerate() {
            let entry       =   ( entry.key(), entry.val() );
            if entry != row_from_entries[ entry_number ] {
                println!(
                    "\nInvalid matrix: the {:?}th structural nonzero entry of row {:?} is\n\n    {:?}\n\naccording to entry look-up, but it is\n\n    {:?}\n\naccording to the row iterator.", 
                    entry_number, 
                    row_index.clone(), 
                    row_from_entries[ entry_number],
                    entry,                    
                );
                return false
            }
        }   
    }

    // check that matrix.column() and matrix.structural_nonzero_entry() agree.
    for column_index in sorted_column_indices.iter().cloned() {
        let column_from_entries: Vec<_>     =   sorted_row_indices
                        .iter()
                        .cloned()
                        .filter_map(
                            |row_index| 
                            {
                                let coefficient     =   matrix.structural_nonzero_entry( &row_index, & column_index );
                                // if coefficient is Some(a), then return (row_index, a)
                                coefficient.map( |a|  (row_index, a) )
                            }                            
                        )                        
                        .collect();
        let column_from_iter: Vec<_>         =   matrix.column( & column_index ).collect();


        if column_from_iter.len() != column_from_entries.len() {
            println!(
                "\nInvalid matrix: calling\n\n    matrix.column({:?})\n\nreturns the following sequence of structural nonzero entries:\n",
                &column_index
            );
                println!("    < begin list of entries >");              
            for entry in &column_from_iter {              
                println!("    {:?}", entry);                   
            }
                println!("    < end   list of entries >");                 
            println!(
                "\nHowever, the sequence of structural nonzero entries obtained by calling\n\n    matrix.structural_nonzero_entry(row_index, {:?})\n\nfor every value of `row_index` provided by the user is not equal:",
                &column_index
            );
                println!("    < begin list of entries >");            
            for entry in &column_from_entries {
                println!("    {:?}", entry);              
            }
                println!("    < end   list of entries >");              
            return false
        }        

        for ( entry_number, entry ) in column_from_iter.iter().cloned().enumerate() {
            let entry       =   ( entry.key(), entry.val() );
            if entry != column_from_entries[ entry_number ] {
                println!(
                    "Invalid matrix: the {:?}th structural nonzero entry of column\n\n    {:?}\n\nis equal to\n\n    {:?}\n\naccording to entry look-up, but it is\n\n    {:?}\n\naccording to the column iterator\n\n    matrix.column({:?}).", 
                    entry_number, 
                    column_index.clone(), 
                    column_from_entries[ entry_number],
                    entry,                
                    column_index.clone(),                         
                );
                return false
            }
        }        
    }   
    
    // check that row() agrees with row_reverse()
    for row_index in sorted_row_indices.iter().cloned() {    
        let     row_from_forward: Vec<_>    =   matrix.row( & row_index ).collect();
        let mut row_from_reverse: Vec<_>    =   matrix.row_reverse( & row_index ).collect();
        ( &mut row_from_reverse ).reverse();
        if ! row_from_forward.iter().eq( row_from_reverse.iter() ) {
            println!(
                "\nInvalid matrix: the sequence of entries returned by `matrix.row_reverse({:?})` is \n\n{:?}\n\nThis does not equal the reverse of the sequence of entries returned by `matrix.row({:?}); that reversed sequence is\n\n{:?}\n\n`.",
                row_index.clone(),
                row_from_reverse,
                row_index.clone(),
                row_from_forward,
            );
            return false
        }
    }

    // check that column() agrees with column_reverse()
    for column_index in sorted_column_indices.iter().cloned() {    
        let     column_from_forward: Vec<_>             =   matrix.column( & column_index ).collect();
        let mut column_from_reverse: Vec<_>     =   matrix.column_reverse( & column_index ).collect();
        ( &mut column_from_reverse ).reverse();        
        if ! column_from_forward.iter().eq( column_from_reverse.iter() ) {
            println!(
                "\nInvalid matrix: the sequence of entries returned by `matrix.column_reverse({:?})` is \n\n{:?}\n\nThis does not equal the reverse of the sequence of entries returned by `matrix.column({:?}); that reversed sequence is\n\n{:?}\n\n`.",                                
                column_index.clone(),
                column_from_reverse,
                column_index.clone(),
                column_from_forward,                
            );
            return false
        }
    }    

    return true

}






/// Verifies that the entries in each row and column appear in strictly sorted order
/// 
/// Concretely, for each `row_index` in `sorted_row_indices` and each `column_index` in `column_indices`,
/// verifies that
/// - the entries in `matrix.row(row_index)` appear in strictly ascending order, as measured by `matrix.order_operator_for_row_entries()`
/// - the entries in `matrix.row(row_index)` appear in strictly ascending order, as measured by `matrix.order_operator_for_row_indices()`
/// - the entries in `matrix.column(column_index)` appear in strictly ascending order, as measured by `matrix.order_operator_for_column_entries()`
/// - the entries in `matrix.column(column_index)` appear in strictly ascending order, as measured by `matrix.order_operator_for_column_indices()`
pub fn matrix_order_operators_are_internally_consistent< Matrix: MatrixAlgebra, I, J >
    (
        matrix:                 Matrix,
        sorted_row_indices:     I,
        sorted_column_indices:  J,
    )
        ->
    Result< (), String >


    where
        I:  IntoIterator< Item = Matrix::RowIndex >,
        J:  IntoIterator< Item = Matrix::ColumnIndex >,        
{

    // get copies of the order operators    
    let order_operator_for_row_entries                  =   matrix.order_operator_for_row_entries();
    let order_operator_for_column_entries               =   matrix.order_operator_for_column_entries();
    let order_operator_for_row_indices                  =   matrix.order_operator_for_row_indices();
    let order_operator_for_column_indices               =   matrix.order_operator_for_column_indices();

    // collect the row and column indices into vectors
    let sorted_row_indices      =   sorted_row_indices.into_iter().collect::<Vec<_>>();
    let sorted_column_indices   =   sorted_column_indices.into_iter().collect::<Vec<_>>(); 

    // check that the user-provided row and column indices are sorted
    // ------------------------------------------------------------------------    
    if ! sorted_row_indices.iter().cloned().is_sorted_strictly_by_order_operator( order_operator_for_row_indices.clone() ).is_ok() {
        let cap = sorted_row_indices.len().min(100);
        println!( "The user-provided list of row indices is not sorted. Here are the first 100 (or fewer) entries:\n{:?}", &sorted_row_indices[..cap]);
        return Err( "The user-provided list of row indices is not sorted.".to_owned() )
    }
    if ! sorted_column_indices.iter().cloned().is_sorted_strictly_by_order_operator( order_operator_for_column_indices.clone() ).is_ok() {
        let cap = sorted_column_indices.len().min(100);        
        println!( "The user-provided list of column indices is not sorted. Here are the first 100 (or fewer) entries:\n{:?}", &sorted_column_indices[..cap]);        
        return Err( "The user-provided list of column indices is not sorted.".to_owned() )
    }  


    // check that rows are sorted
    // ------------------------------------------------------------------------    
    for index in sorted_row_indices.iter().cloned() {     
        // check sorting by the entry order
        if !    matrix.row( & index )
                    .is_sorted_strictly_by_order_operator( order_operator_for_row_entries.clone() )
                    .is_ok()
        {
            let row_vec = matrix.row( & index ).collect::<Vec<_>>();
            let mut msg = format!("\nThe entries of row {:?} do not appear in strictly sorted (ascending) order. The entries are:\n", index);
            for entry in &row_vec {
                msg.push_str(&format!("    {:?}\n", entry));
            }
            return Err(msg);
        }
        // check sorting by the index order
        if !    matrix.row( & index )
                    .map(
                        |entry|
                        entry.key()
                    )
                    .is_sorted_strictly_by_order_operator( order_operator_for_column_indices.clone() )
                    .is_ok()
        {
            let row_vec = matrix.row( & index ).collect::<Vec<_>>();
            let mut msg = format!("\nThe entries of row {:?} do not appear in strictly sorted (ascending) order. The entries are:\n", index);
            for entry in &row_vec {
                msg.push_str(&format!("    {:?}\n", entry));
            }
            return Err(msg);
        }        
    }  


    // check that columns are sorted
    // ------------------------------------------------------------------------
    for index in sorted_column_indices.iter().cloned() {
        // check sorting by the entry order
        if !    matrix.column( & index )
                    .is_sorted_strictly_by_order_operator( order_operator_for_column_entries.clone() )
                    .is_ok()
        {
            let column_vec = matrix.column( & index ).collect::<Vec<_>>();                
            let mut msg = format!("\nThe entries of column {:?} do not appear in strictly sorted (ascending) order. The entries are:\n", index);
            for entry in &column_vec {
                msg.push_str(&format!("    {:?}\n", entry));
            }
            return Err(msg);
        }
        // check sorting by the index order
        if !    matrix.column( & index )
                    .map(
                        |entry|
                        entry.key()
                    )
                    .is_sorted_strictly_by_order_operator( order_operator_for_row_indices.clone() )
                    .is_ok()
        {
            let column_vec = matrix.column( & index ).collect::<Vec<_>>();            
            let mut msg = format!("\nThe entries of column {:?} do not appear in strictly sorted (ascending) order. The entries are:\n", index);
            for entry in &column_vec {
                msg.push_str(&format!("    {:?}\n", entry));
            }
            return Err(msg);
        }        
    }      


    // check that the information contained in the matrix oracle is internally consistent
    if !    matrix_oracle_is_internally_consistent(
                & matrix, 
                sorted_row_indices.iter().cloned(), 
                sorted_column_indices.iter().cloned()
            )
    {
        return Err( "Underlying matrix oracle failed to pass the `matrix_oracle_is_internally_consistent` test".to_owned() )
    } else {
        return Ok(())
    }
}











//  CHECK THAT TWO MATRICES ARE EQUAL
//  ==============================================================================


/// Determines if both matrices are internally consistent and equal
/// 
/// Returns a `HashMap< String, bool >` which records the following data:
/// - matrix 1 is internally consistent:              true/false 
/// - matrix 2 is internally consistent:              true/false 
/// - matrices are internally consistent and equal:   true/false   
/// 
/// For details about "internal consistency," see the documentation for [matrix_oracle_is_internally_consistent].
fn matrices_are_internally_consistent_and_equal< Matrix1, Matrix2, RowIndexIter, ColumnIndexIter >
    ( 
        matrix_1:                   Matrix1,
        matrix_2:                   Matrix2,
        sorted_row_indices:         RowIndexIter,
        sorted_column_indices:      ColumnIndexIter
    ) 
    -> HashMap< String, bool > 

    where 
        Matrix1:            MatrixOracle,
        Matrix2:            MatrixOracle< RowIndex=Matrix1::RowIndex, ColumnIndex=Matrix1::ColumnIndex, Coefficient=Matrix1::Coefficient >,
        RowIndexIter:       IntoIterator< Item = Matrix1::RowIndex    >,
        ColumnIndexIter:    IntoIterator< Item = Matrix1::ColumnIndex >, 

{
    let sorted_row_indices: Vec<_>              =   sorted_row_indices.into_iter().collect();
    let sorted_column_indices: Vec<_>           =   sorted_column_indices.into_iter().collect();


    let a = matrix_oracle_is_internally_consistent( &matrix_1, sorted_row_indices.clone(), sorted_column_indices.clone() );
    let b = matrix_oracle_is_internally_consistent( &matrix_2, sorted_row_indices.clone(), sorted_column_indices.clone() );





    // Check that "structural nonzero entry" look-up operattions agree.
    let mut c                                   =   true;
    for row_index in sorted_row_indices.iter().cloned() {
        let iter_1  =   matrix_1
                                                        .row( & row_index )
                                                        .map(
                                                            |x| 
                                                            ( x.key(), x.val() )  // convert entry to a tuple for ease of comparison
                                                        );
        let iter_2  =   matrix_2
                                                        .row( & row_index.clone() )
                                                        .map(
                                                            |x| 
                                                            ( x.key(), x.val() )  // convert entry to a tuple for ease of comparison
                                                        );
        if ! iter_1.eq( iter_2 ) { c = false }                                                        

    }

    let mut output = HashMap::new();
    output.insert( String::from("matrix 1 is internally consistent"), a  );
    output.insert( String::from("matrix 2 is internally consistent"), b  );
    output.insert( String::from("matrices are internally consistent and equal"), a && b && c  );   

    output     

}











//  UTILITY FUNCTION FOR TESTING: CHECK THAT A PRODUCT OF TWO MATRICES IS IDENTITY
//  ==============================================================================

/// Checks that a (user-specified) set of rows of the product of two matrices equals the corresponding set of rows of an identityt matrix
/// 
/// Concretely, returns `false` if `product.row( k )` does not equal the `k`th standard unit vector, for any `k` in `iter_row_index`, where
/// `product` is the product of the two matrices.
pub fn product_is_identity_matrix<
            Matrix1, 
            Matrix2,               
            RowIndexIterator,
        > 
        (
            matrix_1:           Matrix1,
            matrix_2:           Matrix2,
            iter_row_index:     RowIndexIterator,
        )
        ->
        bool
    where
        Matrix1:                            MatrixAlgebra<
                                                RowEntry:       KeyValSet,
                                                ColumnEntry:    KeyValSet,
                                            >,
        Matrix2:                            MatrixAlgebra< 
                                                Coefficient                 =   Matrix1::Coefficient, 
                                                RowIndex                    =   Matrix1::ColumnIndex, 
                                                ColumnIndex                 =   Matrix1::RowIndex, 
                                                RingOperator                =   Matrix1::RingOperator,
                                                OrderOperatorForRowIndices  =   Matrix1::OrderOperatorForColumnIndices,
                                                RowEntry:                       KeyValSet,
                                            >, 
        RowIndexIterator:                   Iterator< Item = Matrix1::RowIndex >,                    
{

    let product = ProductMatrix::new( matrix_1, matrix_2, );
    let one = Matrix1::RingOperator::one();

    for row_index in iter_row_index {
        let row = product.row( & row_index );
        itertools::assert_equal(
            row.map( |x| (x.key(), x.val()) ),
            std::iter::once( ( row_index, one.clone() ) )
        );
        // match equals_standard_unit_vector { true => {continue}, false=> {return false } }
    }
    true
}