//! A rare-use module to aid iteration over entries in a column.
//! 
//! The main purpose of this module is to provide a convenient (though inefficient!)
//! method to iterate over the entries in a column of a sparse matrix, when
//! it is otherwise hard to do so.


use derive_getters::{Getters, Dissolve};
use derive_new::new;

use crate::algebra::vectors::entries::KeyValGet;

use super::MatrixOracle;








/// Iterates over the entries of a column
/// 
/// This struct contains three pieces of data
/// 
/// - a matrix oracle, M
/// - an iterator that runs over every row index for the matrix, I
/// - a column index, C
/// 
/// When we call `self.next()` on the struct, it uses I to generate a
/// row index, r; it then searches row r for an entry of form (C, v).
/// If it finds one, then the struct returns `Some((C,v))`. Otherwise
/// it takes another row index from I, and repeats.  
/// 
/// # Example
/// 
/// ```
/// use oat_rust::algebra::matrices::query::{MatrixOracle};
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; 
/// 
/// 
/// // data for a matrix:
/// // |  1   0 |
/// // | -1   0 |
/// // |  0   0 |
/// // |  0   0 |           
/// let matrix = vec![ 
///         vec![ (0, 1.)  ], 
///         vec![ (0,-1.)  ], 
///         vec![          ], 
///         vec![          ],                 
///     ];
/// 
/// // wrap the data in a VecOfVec sparse matrix struct
/// let matrix = & VecOfVec::new(matrix).ok().unwrap();
/// 
/// // get some sparse columns
/// let mut column_0        =   matrix.column(&0);
/// let mut column_0_rev    =   matrix.column_reverse(&0);
/// let mut column_1        =   matrix.column(&1);
/// let mut column_1_rev    =   matrix.column_reverse(&1);
/// 
/// // check the columns are correct
/// assert_eq!( column_0.next(),        Some( (0,  1.) )    );
/// assert_eq!( column_0.next(),        Some( (1, -1.) )    );
/// assert_eq!( column_0.next(),        None                );    
/// 
/// assert_eq!( column_0_rev.next(),    Some( (1, -1.) )    );
/// assert_eq!( column_0_rev.next(),    Some( (0,  1.) )    );
/// assert_eq!( column_0_rev.next(),    None                );   
/// 
/// assert_eq!( column_1.next(),        None        );
/// assert_eq!( column_1_rev.next(),    None        );   
/// ```  
#[derive(Debug, Copy, Clone, Eq, PartialEq, Getters, Dissolve, new)]
pub struct SparseColumn< Matrix, RowIndexIterator >
    where
        Matrix:             MatrixOracle,
        RowIndexIterator:   Clone + Iterator< Item = Matrix::RowIndex >,
{
    pub matrix:                 Matrix,                 // oracle for a sparse matrix
    pub row_index_iterator:     RowIndexIterator,       // iterates over the row indices of the matrix, exhaustively
    pub column_index:           Matrix::ColumnIndex,    // we want the entries from this column
} 

 // Iterator
impl < Matrix, RowIndexIterator >

        Iterator for 

        SparseColumn< Matrix, RowIndexIterator >

    where
        Matrix:                 MatrixOracle,
        Matrix::RowIndex:       Clone,
        Matrix::ColumnIndex:    Eq,
        RowIndexIterator:       Clone + Iterator< Item = Matrix::RowIndex >,
{
    type Item = ( Matrix::RowIndex, Matrix::Coefficient );

    fn next( &mut self ) -> Option< Self::Item > {
        
        // this while-loop runs over row indices, until a condition breaks the loop
        while let Some( row_index ) = self.row_index_iterator.next() {

            // given a row index, search the corresponding row for an entry with the correct column index
            for entry in self.matrix.row( & row_index ) {
                if entry.key() == self.column_index {
                    return Some( ( row_index, entry.val() ) )
                }
            }
        }
        return None
    }
}  