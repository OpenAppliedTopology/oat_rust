//! Sparse matrices
//! 
//! - [About](#about)
//! - [Tools](#modules)
//! - [Build your own matrix](#build-your-own-matrix)
//! 
//! 
//! 
//! # About
//!
//! OAT provides a variety of tools to work with sparse matrices.  Key features include
//! - flexible indexing: matrices can be indexed by arbitrary keys; for example, the boundary matrix of a simplicial complex can be indexed by simplices; this feature is powerful in a variety of homology computations, where explit enumeration of row and column indices is expensive
//! - flexible coefficients: you can work over any coeficient ring or field using [ring operators](crate::algebra::rings#terminology).
//! - lazy look up: rows and columns of matrices can be constructed in a lazy fashion, consistent with current state-of-the-art practices in large scale PH computation
//! - extensive unit testing
//! 
//! 
//! # Build your own matrix
//!
//! Users can work with any object *as if it were a matrix*, provided the object implements the right [query traits](crate::algebra::matrices::query).
//! Any object that implements one of these traits is called a **matrix oracle**.
//!
//! Let's illustrate with an example.  We'll define an object that represents an identity matrix, but throw in a twist:
//! - rows are indexed integers
//! - columns are indexed by strings
//! - for example, row 0 has a single nonzero entry, in column "0".  Note the difference: the row index is a number, but the column index is a string.
//! 
//! This object uses no memory at all!
//! 
//! 
//! ```
//! use oat_rust::algebra::matrices::query::MatrixOracle;
//! use std::iter::Once;
//! 
//! 
//! //  -------------------------------------------------------------------
//! //  Define the new struct, and create an instance
//! //  -------------------------------------------------------------------
//! 
//! pub struct MyMatrix; // this struct contains no data, and uses zero memory
//! let matrix = MyMatrix; // here's an instance
//! 
//! 
//! //  -------------------------------------------------------------------
//! //  Implement `MatrixOracle`
//! //  -------------------------------------------------------------------
//! 
//! impl MatrixOracle for MyMatrix {
//!     type Coefficient            =   i32;
//! 
//!     type RowIndex               =   usize;
//! 
//!     type ColumnIndex            =   String;
//! 
//!     type RowEntry               =   ( String, i32 );
//! 
//!     type ColumnEntry            =   ( usize, i32 );
//! 
//!     type Row                    =   Once<( String, i32 )>;
//! 
//!     type RowReverse             =   Once<( String, i32 )>;
//! 
//!     type Column                 =   Once<( usize, i32 )>;
//! 
//!     type ColumnReverse          =   Once<( usize, i32 )>;
//! 
//!     fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
//!         let column_index        =   index.to_string(); // convert the integer to a String
//!         std::iter::once( ( column_index, 1 ) )
//!     }
//! 
//!     fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
//!         let column_index        =   index.to_string(); // convert the integer to a String
//!         std::iter::once( ( column_index, 1 ) )
//!     }
//! 
//!     fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
//!         let row_index = index.parse::<usize>().unwrap(); // convert the string to an integer
//!         std::iter::once( ( row_index, 1 ) )
//!     }
//! 
//!     fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
//!         let row_index = index.parse::<usize>().unwrap(); // convert the string to an integer
//!         std::iter::once( ( row_index, 1 ) )
//!     }
//! 
//!     fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
//!         true // this matrix has a row for every possible integer
//!     }
//! 
//!     fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
//!         index.parse::<u64>().is_ok() // the matrix has a column for index `index` if and only if the `index` can be converted to a usize
//!     }
//! 
//!     fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
//!         if column == &row.to_string() {
//!             // if the row and column indices agree, then return 1.0
//!             Some(1)
//!         } else {
//!             // otherwise return None to indicate that the entry is structurally zero
//!             None
//!         }
//!     }
//! }
//! 
//! 
//! 
//! //  -------------------------------------------------------------------
//! //  Test
//! //  -------------------------------------------------------------------
//! 
//! 
//! assert!( matrix.row( & 0 ).eq( std::iter::once( ( 0.to_string(), 1 ) ) ) );
//! assert!( matrix.row( & 1 ).eq( std::iter::once( ( 1.to_string(), 1 ) ) ) );   
//! assert!( matrix.row_reverse( & 0 ).eq( std::iter::once( ( 0.to_string(), 1 ) ) ) );
//! assert!( matrix.row_reverse( & 1 ).eq( std::iter::once( ( 1.to_string(), 1 ) ) ) ); 
//! 
//! assert!( matrix.column( & 0.to_string() ).eq( std::iter::once( ( 0, 1 ) ) ) );
//! assert!( matrix.column( & 1.to_string() ).eq( std::iter::once( ( 1, 1 ) ) ) );   
//! assert!( matrix.column_reverse( & 0.to_string() ).eq( std::iter::once( ( 0, 1 ) ) ) );
//! assert!( matrix.column_reverse( & 1.to_string() ).eq( std::iter::once( ( 1, 1 ) ) ) );
//! 
//! assert_eq!( matrix.structural_nonzero_entry( & 0, & 0.to_string() ),  Some(1) );
//! assert_eq!( matrix.structural_nonzero_entry( & 0, & 1.to_string() ),  None      );
//! 
//! assert!( matrix.has_row_for_index( & 0 ) );
//! ```
//!
//! 
//! 


pub mod operations;
pub mod query;
pub mod types;
pub mod debug;
pub mod display;








#[cfg(test)]
mod tests {



    // Note this useful idiom: importing names from outer (for mod tests) scope.


    #[test]
    fn doc_test_build_your_own_matrix() {

        use crate::algebra::matrices::query::MatrixOracle;
        use std::iter::Once;


        //  -------------------------------------------------------------------
        //  Define the new struct, and create an instance
        //  -------------------------------------------------------------------

        pub struct MyMatrix; // this struct contains no data, and uses zero memory
        let matrix = MyMatrix; // here's an instance


        //  -------------------------------------------------------------------
        //  Implement `MatrixOracle`
        //  -------------------------------------------------------------------

        impl MatrixOracle for MyMatrix {
            type Coefficient            =   i32;
        
            type RowIndex               =   usize;
        
            type ColumnIndex            =   String;
        
            type RowEntry               =   ( String, i32 );
        
            type ColumnEntry            =   ( usize, i32 );
        
            type Row                    =   Once<( String, i32 )>;
        
            type RowReverse             =   Once<( String, i32 )>;
        
            type Column                 =   Once<( usize, i32 )>;
        
            type ColumnReverse          =   Once<( usize, i32 )>;
        
            fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
                let column_index        =   index.to_string(); // convert the integer to a String
                std::iter::once( ( column_index, 1 ) )
            }
        
            fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
                let column_index        =   index.to_string(); // convert the integer to a String
                std::iter::once( ( column_index, 1 ) )
            }
        
            fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
                let row_index = index.parse::<usize>().unwrap(); // convert the string to an integer
                std::iter::once( ( row_index, 1 ) )
            }
        
            fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
                let row_index = index.parse::<usize>().unwrap(); // convert the string to an integer
                std::iter::once( ( row_index, 1 ) )
            }
        
            fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
                true // this matrix has a row for every possible integer
            }
        
            fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
                index.parse::<u64>().is_ok() // the matrix has a column for index `index` if and only if the `index` can be converted to a usize
            }
        
            fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
                if column == &row.to_string() {
                    // if the row and column indices agree, then return 1.0
                    Some(1)
                } else {
                    // otherwise return None to indicate that the entry is structurally zero
                    None
                }
            }
        }



        //  -------------------------------------------------------------------
        //  Test
        //  -------------------------------------------------------------------


        assert!( matrix.row( & 0 ).eq( std::iter::once( ( 0.to_string(), 1 ) ) ) );
        assert!( matrix.row( & 1 ).eq( std::iter::once( ( 1.to_string(), 1 ) ) ) );   
        assert!( matrix.row_reverse( & 0 ).eq( std::iter::once( ( 0.to_string(), 1 ) ) ) );
        assert!( matrix.row_reverse( & 1 ).eq( std::iter::once( ( 1.to_string(), 1 ) ) ) ); 

        assert!( matrix.column( & 0.to_string() ).eq( std::iter::once( ( 0, 1 ) ) ) );
        assert!( matrix.column( & 1.to_string() ).eq( std::iter::once( ( 1, 1 ) ) ) );   
        assert!( matrix.column_reverse( & 0.to_string() ).eq( std::iter::once( ( 0, 1 ) ) ) );
        assert!( matrix.column_reverse( & 1.to_string() ).eq( std::iter::once( ( 1, 1 ) ) ) );

        assert_eq!( matrix.structural_nonzero_entry( & 0, & 0.to_string() ),  Some(1) );
        assert_eq!( matrix.structural_nonzero_entry( & 0, & 1.to_string() ),  None      );

        assert!( matrix.has_row_for_index( & 0 ) );

    }
}
