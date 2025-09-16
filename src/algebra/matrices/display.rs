//! Display data from a matrix oracle.


use itertools::Itertools;

use super::query::MatrixOracle;



/// Print one row of the matrix (represented as a vector of entries) for each item in the iterator.
/// 
/// # Example
/// 
/// ```
/// // import the relevant crates
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; // a particular type of sparse matrix
/// use oat_rust::algebra::matrices::query::MatrixOracle; // the trait for looking up rows and columns
/// use oat_rust::algebra::matrices::display::print_indexed_rows;
///         
/// // define sparse version of the following matrix
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] ).ok().unwrap();
/// 
/// // define an iterator that runs over the row indices (i.e. the row indices)
/// let iter_row_index = 0..2;  
/// 
/// // print the rows specified by the iterator.  this should show the following:
/// // $ row 0: [(1, 5), (2, 5)]
/// // $ row 1: [(2, 7)]
/// print_indexed_rows( &&matrix, iter_row_index ); // we pass `&&matrix` because `&matrix` implements the oracle trait and we have to pass a *reference* to something that implements the oracle trait
/// ```
pub fn print_indexed_rows
            < Matrix, RowIndexIterator >
            ( matrix: & Matrix, iter_row_index: RowIndexIterator ) 
    where   
        Matrix:                         MatrixOracle,     
        RowIndexIterator:               IntoIterator< Item = Matrix::RowIndex >,        
{
    for row_index in iter_row_index {        
        println!("row for index {:?}", row_index.clone());
        for entry in matrix.row( &row_index ).into_iter()  { 
            println!("{:?}", entry) 
        }
    }
}

/// Print one column of the matrix (represented as a vector of entries) for each item in the iterator.
/// 
/// # Example
/// 
/// ```
/// // import the relevant crates
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; // a particular type of sparse matrix
/// use oat_rust::algebra::matrices::query::MatrixOracle; // the trait for looking up rows and columns
/// use oat_rust::algebra::matrices::display::print_indexed_columns;
///         
/// // define sparse version of the following matrix
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] ).ok().unwrap();
/// 
/// // define an iterator that runs over the column indices (i.e. the row indices)
/// let iter_column_index = 0..3;  
/// 
/// // print the columns specified by the iterator.  this should show the following:
/// // $ column 0: []
/// // $ column 1: [(0, 5)]
/// // $ column 2: [(1, 5), (2, 7)]
/// print_indexed_columns( &&matrix, iter_column_index ); // we pass `&&matrix` because `&matrix` implements the oracle trait and we have to pass a *reference* to something that implements the oracle trait
/// ```
pub fn print_indexed_columns
            < Matrix, ColumnIndexIterator >
            ( matrix: & Matrix, iter_column_index: ColumnIndexIterator ) 
    where   
        Matrix:                         MatrixOracle,      
        ColumnIndexIterator:            IntoIterator< Item = Matrix::ColumnIndex >,        
{
    for column_index in iter_column_index {
        println!(
            "column for index {:?}: {:?}", 
            column_index.clone(), 
            matrix.column_reverse( & column_index )
                .into_iter().collect_vec()  );
    }
}









#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    

    #[test]
    fn test_print_indexed_rows() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; // a particular type of sparse matrix
         // the trait for looking up rows and columns
        use crate::algebra::matrices::display::print_indexed_rows;

        // define a row-major sparse matrix representing the array
        //  0   5   5
        //  0   0   7
        let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] ).ok().unwrap();
        let iter_row_index = 0..2;  // iterates over the row indices, i.e. over the row indices

        // print the rows specified by the iterator
        print_indexed_rows( &&matrix, iter_row_index ); // we pass `&&matrix` because `&matrix` implements the oracle trait and we have to pass a *reference* to something that implements the oracle trait
        // this should show the following:
    }


}
