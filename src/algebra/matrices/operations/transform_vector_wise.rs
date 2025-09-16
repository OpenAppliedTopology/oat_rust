
//! Create a new matrix by applying a transformation to each row or column of an existing matrix.
//! 
//! See also [`transform_entry_wise`](oat_rust::algebra::matrices::operations::transform_entry_wise)



use derive_new::new;
use itertools::PutBack;
use itertools::put_back;

use crate::{algebra::matrices::query::{MatrixAlgebra, MatrixOracle,}, utilities::sets::MapHasKeyOrSequenceHasElement, };











//  PUT BACK ROWS
//  -------------------------------------------------------------------------------------



/// Wraps every row, column, reverse row, and reverse column of a matrix in a [`PutBack`](itertools::PutBack) iterator.
/// 
/// This is useful when you want to iterate over the rows or columns of a matrix, but you also want to be able to "put back" the last entry that you have already iterated over.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::algebra::matrices::operations::transform_vector_wise::PutbackIteratorMatrix;
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// 
/// let matrix = VecOfVec::new( 
///     vec![ 
///         vec![(1,1), (2,2), (3,3)], 
///         vec![(4,4), (5,5), (6,6)] 
///     ] 
/// ).ok().unwrap();    
/// let put_back_matrix = PutbackIteratorMatrix::new( & matrix );
/// 
/// let mut row = put_back_matrix.row(&0);
/// assert_eq!( row.next(), Some( (1,1) ) );
/// row.put_back( (1,1) );
/// assert_eq!( row.next(), Some( (1,1) ) );
/// 
/// let mut column = put_back_matrix.column(&2);
/// assert_eq!( column.next(), Some( (0,2) ) );
/// column.put_back( (0,2) );
/// assert_eq!( column.next(), Some( (0,2) ) );
/// ```
/// 
#[derive(new,Clone,Copy,Debug,Eq,PartialEq,Ord,PartialOrd)]
pub struct PutbackIteratorMatrix<
        Matrix,
    >
{
    matrix:                         Matrix,
}
  

// Implement MatrixOracle
impl < Matrix >

    MatrixOracle for  

    PutbackIteratorMatrix
        < Matrix >

    where
        Matrix:                             MatrixOracle,    
{
    type Coefficient            =     Matrix::Coefficient;   // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =     Matrix::RowIndex;      // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =     Matrix::ColumnIndex;   // The type of column indices    
    type RowEntry               =     Matrix::RowEntry;      // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =     Matrix::ColumnEntry;   // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =     PutBack< Matrix::Row >;           // What you get when you ask for a row.
    type RowReverse             =     PutBack< Matrix::RowReverse >;    // What you get when you ask for a row with the order of entries reversed
    type Column                 =     PutBack< Matrix::Column >;        // What you get when you ask for a column
    type ColumnReverse          =     PutBack< Matrix::ColumnReverse >; // What you get when you ask for a column with the order of entries reversed                             


    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row
        { put_back( self.matrix.row(index) ) }         
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse
        { put_back( self.matrix.row_reverse(index) ) }   
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column
        { put_back( self.matrix.column(index) ) }
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse
        { put_back( self.matrix.column_reverse(index) ) }

    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool
        { self.matrix.has_row_for_index(index) }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { 
            self.matrix.has_column_for_index(index)
        }    
    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { self.matrix.structural_nonzero_entry( row, column ) }

}



// Implement MatrixAlgebra
impl < Matrix >

    MatrixAlgebra for  

    PutbackIteratorMatrix
        < Matrix >

    where
        Matrix:                             MatrixAlgebra,    
{
    type RingOperator                   =   Matrix::RingOperator;
    type OrderOperatorForRowEntries     =   Matrix::OrderOperatorForRowEntries;
    type OrderOperatorForRowIndices     =   Matrix::OrderOperatorForRowIndices;
    type OrderOperatorForColumnEntries  =   Matrix::OrderOperatorForColumnEntries;
    type OrderOperatorForColumnIndices  =   Matrix::OrderOperatorForColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator { self.matrix.ring_operator() }
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries { self.matrix.order_operator_for_row_entries() }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices { self.matrix.order_operator_for_row_indices() }    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries { self.matrix.order_operator_for_column_entries() }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices { self.matrix.order_operator_for_column_indices() }   
}
















//  MASK: ONLY INDICES OUTSIDE
//  -------------------------------------------------------------------------------------

//  COLUMNS
//  -----------------------------------------

use crate::algebra::vectors::operations::OnlyIndicesOutsideCollection;


/// A matrix oracle whose rows iterate only over entries that have indices **outside**
/// a given collection.
#[derive(Clone,Copy,Debug,Eq,PartialEq,Ord,PartialOrd)]
pub struct OnlyColumnIndicesOutsideCollection<
        Matrix, ColumnIndicesToExclude,
    >
    where
        Matrix:                         MatrixOracle,    
        ColumnIndicesToExclude:         MapHasKeyOrSequenceHasElement< Matrix::ColumnIndex >,                  
{
    matrix:                         Matrix,
    column_indices_to_exclude:       ColumnIndicesToExclude,
}

// implement the object
impl < Matrix, ColumnIndicesToExclude, >

    OnlyColumnIndicesOutsideCollection
        < Matrix, ColumnIndicesToExclude >

    where
        Matrix:                         MatrixOracle,    
        ColumnIndicesToExclude:               MapHasKeyOrSequenceHasElement< Matrix::ColumnIndex >,                          

    {
        /// Create a new [OnlyColumnIndicesOutsideCollection] from a matrix and a collection
        /// of column indices.
        /// 
        /// In order to operate properly, the collection of column indices should implement
        /// the [`MapHasKeyOrSequenceHasElement`](crate::utilities::sets::MapHasKeyOrSequenceHasElement) 
        /// trait.
        pub fn new( matrix: Matrix, column_indices_to_exclude: ColumnIndicesToExclude ) 
        -> 
        OnlyColumnIndicesOutsideCollection
            < Matrix, ColumnIndicesToExclude >
        where
            Matrix:                         MatrixOracle,    
            ColumnIndicesToExclude:               MapHasKeyOrSequenceHasElement< Matrix::ColumnIndex >,              
        { OnlyColumnIndicesOutsideCollection{ matrix, column_indices_to_exclude: column_indices_to_exclude, } }
    }     

// Implement MatrixOracle
impl < Matrix, ColumnIndicesToExclude >

    MatrixOracle for  

    OnlyColumnIndicesOutsideCollection
        < Matrix, ColumnIndicesToExclude >

    where
        Matrix:                             MatrixOracle,    
        ColumnIndicesToExclude:             Copy + MapHasKeyOrSequenceHasElement< Matrix::ColumnIndex >,    
{
    type Coefficient            =     Matrix::Coefficient;   // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =     Matrix::RowIndex;      // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =     Matrix::ColumnIndex;   // The type of column indices    
    type RowEntry               =     Matrix::RowEntry;      // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =     Matrix::ColumnEntry;   // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =     OnlyIndicesOutsideCollection< Matrix::Row, ColumnIndicesToExclude >;           // What you get when you ask for a row.
    type RowReverse             =     OnlyIndicesOutsideCollection< Matrix::RowReverse, ColumnIndicesToExclude >;    // What you get when you ask for a row with the order of entries reversed
    type Column                 =     Matrix::Column;        // What you get when you ask for a column
    type ColumnReverse          =     Matrix::ColumnReverse; // What you get when you ask for a column with the order of entries reversed                             


    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row
        { OnlyIndicesOutsideCollection::new( self.matrix.row(index), self.column_indices_to_exclude  )  }        
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse
        { OnlyIndicesOutsideCollection::new( self.matrix.row_reverse(index), self.column_indices_to_exclude  )  }    
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column
        { self.matrix.column(index) }    
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse
        { self.matrix.column_reverse(index) }

    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool
        { self.matrix.has_row_for_index(index) }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { 
            let matrix_has_column = self.matrix.has_column_for_index(index);
            let column_is_excluded = self.column_indices_to_exclude.map_has_key_or_sequence_has_element(index);
            return matrix_has_column && (! column_is_excluded)
        }    
    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { self.matrix.structural_nonzero_entry( row, column ) }

}



// Implement MatrixAlgebra
impl < Matrix, ColumnIndicesToExclude >

    MatrixAlgebra for  

    OnlyColumnIndicesOutsideCollection
        < Matrix, ColumnIndicesToExclude >

    where
        Matrix:                             MatrixAlgebra,    
        ColumnIndicesToExclude:             Copy + MapHasKeyOrSequenceHasElement< Matrix::ColumnIndex >,   
{
    type RingOperator                   =   Matrix::RingOperator;
    type OrderOperatorForRowEntries     =   Matrix::OrderOperatorForRowEntries;
    type OrderOperatorForRowIndices     =   Matrix::OrderOperatorForRowIndices;
    type OrderOperatorForColumnEntries  =   Matrix::OrderOperatorForColumnEntries;
    type OrderOperatorForColumnIndices  =   Matrix::OrderOperatorForColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator { self.matrix.ring_operator() }
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries { self.matrix.order_operator_for_row_entries() }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices { self.matrix.order_operator_for_row_indices() }    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries { self.matrix.order_operator_for_column_entries() }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices { self.matrix.order_operator_for_column_indices() }   
}




//  ROWS
//  -----------------------------------------


/// A matrix oracle whose columns iterate only over entries that have indices **outside**
/// a given collection.
#[derive(Clone,Copy,Debug,Eq,PartialEq,Ord,PartialOrd)]
pub struct OnlyRowIndicesOutsideCollection<
        Matrix, RowIndicesToExclude,
    >
    where
        Matrix:                             MatrixOracle,    
        RowIndicesToExclude:                MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,                  
{
    matrix:                                 Matrix,
    row_indices_to_exclude:                 RowIndicesToExclude,
}

// implement the object
impl < Matrix, RowIndicesToExclude, >

    OnlyRowIndicesOutsideCollection
        < Matrix, RowIndicesToExclude >

    where
        Matrix:                        MatrixOracle,    
        RowIndicesToExclude:               MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,                          

    {
        /// Create a new [OnlyRowIndicesOutsideCollection] from a matrix and a collection
        /// of column indices.
        /// 
        /// In order to operate properly, the collection of column indices should implement
        /// the [`MapHasKeyOrSequenceHasElement`](crate::utilities::sets::MapHasKeyOrSequenceHasElement) 
        /// trait.
        pub fn new( matrix: Matrix, row_indices_to_exclude: RowIndicesToExclude ) 
        -> 
        OnlyRowIndicesOutsideCollection
            < Matrix, RowIndicesToExclude >
        { OnlyRowIndicesOutsideCollection{ matrix, row_indices_to_exclude: row_indices_to_exclude, } }
    }    


// Implement MatrixOracle
impl < Matrix, RowIndicesToExclude >

    MatrixOracle for  

    OnlyRowIndicesOutsideCollection
        < Matrix, RowIndicesToExclude >

    where
        Matrix:                         MatrixOracle,    
        RowIndicesToExclude:            Copy + MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,    
{
    type Coefficient            =     Matrix::Coefficient;   // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =     Matrix::RowIndex;      // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =     Matrix::ColumnIndex;   // The type of column indices    
    type RowEntry               =     Matrix::RowEntry;      // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =     Matrix::ColumnEntry;   // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =     Matrix::Row;           // What you get when you ask for a row.
    type RowReverse             =     Matrix::RowReverse;    // What you get when you ask for a row with the order of entries reversed
    type Column                 =     OnlyIndicesOutsideCollection< Matrix::Column,        RowIndicesToExclude >;        // What you get when you ask for a column
    type ColumnReverse          =     OnlyIndicesOutsideCollection< Matrix::ColumnReverse, RowIndicesToExclude >; // What you get when you ask for a column with the order of entries reversed                             


    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row
        { self.matrix.row(index) }    
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse
        { self.matrix.row_reverse(index) }
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column
        { OnlyIndicesOutsideCollection::new( self.matrix.column(index), self.row_indices_to_exclude  )  }        
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse
        { OnlyIndicesOutsideCollection::new( self.matrix.column_reverse(index), self.row_indices_to_exclude  )  }    

    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool
        { 
            let matrix_has_row = self.matrix.has_row_for_index(index);
            let row_is_excluded = self.row_indices_to_exclude.map_has_key_or_sequence_has_element(index);
            return matrix_has_row && (! row_is_excluded)
        }    
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { self.matrix.has_column_for_index(index) }    
    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { self.matrix.structural_nonzero_entry( row, column ) }

}



// Implement MatrixAlgebra
impl < Matrix, RowIndicesToExclude >

    MatrixAlgebra for  

    OnlyRowIndicesOutsideCollection
        < Matrix, RowIndicesToExclude >

    where
        Matrix:                         MatrixAlgebra,    
        RowIndicesToExclude:            Copy + MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,    
{
    type RingOperator                   =   Matrix::RingOperator;
    type OrderOperatorForRowEntries     =   Matrix::OrderOperatorForRowEntries;
    type OrderOperatorForRowIndices     =   Matrix::OrderOperatorForRowIndices;
    type OrderOperatorForColumnEntries  =   Matrix::OrderOperatorForColumnEntries;
    type OrderOperatorForColumnIndices  =   Matrix::OrderOperatorForColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator { self.matrix.ring_operator() }
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries { self.matrix.order_operator_for_row_entries() }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices { self.matrix.order_operator_for_row_indices() }    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries { self.matrix.order_operator_for_column_entries() }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices { self.matrix.order_operator_for_column_indices() }   
}



//  MASK: ONLY INDICES INSIDE
//  -------------------------------------------------------------------------------------



//  COLUMNS
//  -----------------------------------------


use crate::algebra::vectors::operations::OnlyIndicesInsideCollection;


/// A matrix oracle whose rows iterate only over entries that have indices **outside**
/// a given collection.
#[derive(Clone,Copy,Debug,Eq,PartialEq,Ord,PartialOrd)]
pub struct OnlyColumnIndicesInsideCollection
                < Matrix, ColumnIndicesToInclude, >
{
    matrix:                             Matrix,
    column_indices_to_include:          ColumnIndicesToInclude,   
}

// implement the object
impl < Matrix, ColumnIndicesToInclude >

    OnlyColumnIndicesInsideCollection
        < Matrix, ColumnIndicesToInclude >

    {
        /// Create a new [OnlyColumnIndicesInsideCollection] from a matrix and a collection
        /// of column indices.
        /// 
        /// In order to operate properly, the collection of column indices should implement
        /// the [`MapHasKeyOrSequenceHasElement`](crate::utilities::sets::MapHasKeyOrSequenceHasElement) 
        /// trait.
        pub fn new( matrix: Matrix, column_indices_to_include: ColumnIndicesToInclude ) 
        -> 
        OnlyColumnIndicesInsideCollection
            < Matrix, ColumnIndicesToInclude, >
        where
            Matrix:                         MatrixOracle,    
            ColumnIndicesToInclude:         MapHasKeyOrSequenceHasElement< Matrix::ColumnIndex >,              
        { OnlyColumnIndicesInsideCollection{ matrix, column_indices_to_include } }
    }    





// Implement MatrixOracle
impl < Matrix, ColumnIndicesToInclude >

    MatrixOracle for  

    OnlyColumnIndicesInsideCollection
        < Matrix, ColumnIndicesToInclude >

    where
        Matrix:                             MatrixOracle,    
        ColumnIndicesToInclude:             Copy + MapHasKeyOrSequenceHasElement< Matrix::ColumnIndex >,    
{
    type Coefficient            =     Matrix::Coefficient;   // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =     Matrix::RowIndex;      // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =     Matrix::ColumnIndex;   // The type of column indices    
    type RowEntry               =     Matrix::RowEntry;      // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =     Matrix::ColumnEntry;   // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =     OnlyIndicesInsideCollection< Matrix::Row,        ColumnIndicesToInclude >;           // What you get when you ask for a row.
    type RowReverse             =     OnlyIndicesInsideCollection< Matrix::RowReverse, ColumnIndicesToInclude >;    // What you get when you ask for a row with the order of entries reversed
    type Column                 =     Matrix::Column;        // What you get when you ask for a column
    type ColumnReverse          =     Matrix::ColumnReverse; // What you get when you ask for a column with the order of entries reversed                             


    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row
        { OnlyIndicesInsideCollection::new( self.matrix.row(index), self.column_indices_to_include  )  }        
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse
        { OnlyIndicesInsideCollection::new( self.matrix.row_reverse(index), self.column_indices_to_include  )  }    
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column
        { self.matrix.column(index) }    
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse
        { self.matrix.column_reverse(index) }

    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool
        { self.matrix.has_row_for_index(index) }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { 
            let matrix_has_column = self.matrix.has_column_for_index(index);
            let column_is_included = self.column_indices_to_include.map_has_key_or_sequence_has_element(index);
            return matrix_has_column && column_is_included
        }    
    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { self.matrix.structural_nonzero_entry( row, column ) }

}


// Implement MatrixAlgebra
impl < Matrix, ColumnIndicesToExclude >

    MatrixAlgebra for  

    OnlyColumnIndicesInsideCollection
        < Matrix, ColumnIndicesToExclude >

    where
        Matrix:                             MatrixAlgebra,    
        ColumnIndicesToExclude:             Copy + MapHasKeyOrSequenceHasElement< Matrix::ColumnIndex >,   
{
    type RingOperator                   =   Matrix::RingOperator;
    type OrderOperatorForRowEntries     =   Matrix::OrderOperatorForRowEntries;
    type OrderOperatorForRowIndices     =   Matrix::OrderOperatorForRowIndices;
    type OrderOperatorForColumnEntries  =   Matrix::OrderOperatorForColumnEntries;
    type OrderOperatorForColumnIndices  =   Matrix::OrderOperatorForColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator { self.matrix.ring_operator() }
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries { self.matrix.order_operator_for_row_entries() }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices { self.matrix.order_operator_for_row_indices() }    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries { self.matrix.order_operator_for_column_entries() }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices { self.matrix.order_operator_for_column_indices() }   
}


//  ROWS
//  -----------------------------------------


/// A matrix oracle whose rows iterate only over entries that have indices **outside**
/// a given collection.
#[derive(Clone,Copy,Debug,Eq,PartialEq,Ord,PartialOrd)]
pub struct OnlyRowIndicesInsideCollection<
        Matrix, RowIndicesToInclude,
    >
    where
        Matrix:                         MatrixOracle,    
        RowIndicesToInclude:            MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,                  
{
    matrix:                             Matrix,
    row_indices_to_include:             RowIndicesToInclude,
}

// implement the object
impl < Matrix, RowIndicesToInclude, >

    OnlyRowIndicesInsideCollection
        < Matrix, RowIndicesToInclude >

    where
        Matrix:                        MatrixOracle,    
        RowIndicesToInclude:           MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,                          

    {
        /// Create a new [OnlyRowIndicesInsideCollection] from a matrix and a collection
        /// of row indices.
        /// 
        /// In order to operate properly, the collection of row indices should implement
        /// the [`MapHasKeyOrSequenceHasElement`](crate::utilities::sets::MapHasKeyOrSequenceHasElement) 
        /// trait.
        pub fn new( matrix: Matrix, row_indices_to_include: RowIndicesToInclude ) 
        -> 
        OnlyRowIndicesInsideCollection
            < Matrix, RowIndicesToInclude >
        where
            RowIndicesToInclude:                MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,              
        { OnlyRowIndicesInsideCollection{ matrix, row_indices_to_include } }
    }    



impl < Matrix, RowIndicesToExclude >

    MatrixOracle for  

    OnlyRowIndicesInsideCollection
        < Matrix, RowIndicesToExclude >

    where
        Matrix:                          MatrixOracle,    
        RowIndicesToExclude:             Copy + MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,    
{
    type Coefficient            =     Matrix::Coefficient;   // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =     Matrix::RowIndex;      // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =     Matrix::ColumnIndex;   // The type of column indices    
    type RowEntry               =     Matrix::RowEntry;      // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =     Matrix::ColumnEntry;   // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =     Matrix::Row;           // What you get when you ask for a row.
    type RowReverse             =     Matrix::RowReverse;    // What you get when you ask for a row with the order of entries reversed
    type Column                 =     OnlyIndicesInsideCollection< Matrix::Column,        RowIndicesToExclude >;        // What you get when you ask for a column
    type ColumnReverse          =     OnlyIndicesInsideCollection< Matrix::ColumnReverse, RowIndicesToExclude >; // What you get when you ask for a column with the order of entries reversed                             


    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row
        { self.matrix.row(index) }    
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse
        { self.matrix.row_reverse(index) }
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column
        { OnlyIndicesInsideCollection::new( self.matrix.column(index), self.row_indices_to_include  )  }        
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse
        { OnlyIndicesInsideCollection::new( self.matrix.column_reverse(index), self.row_indices_to_include  )  }    

    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool
        { 
            let matrix_has_row = self.matrix.has_row_for_index(index);
            let row_is_included = self.row_indices_to_include.map_has_key_or_sequence_has_element(index);
            return matrix_has_row && row_is_included
        }    
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { self.matrix.has_column_for_index(index) }    
    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { self.matrix.structural_nonzero_entry( row, column ) }

}    




// Implement MatrixAlgebra
impl < Matrix, RowIndicesToExclude >

    MatrixAlgebra for  

    OnlyRowIndicesInsideCollection
        < Matrix, RowIndicesToExclude >

    where
        Matrix:                         MatrixAlgebra,    
        RowIndicesToExclude:            Copy + MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,    
{
    type RingOperator                   =   Matrix::RingOperator;
    type OrderOperatorForRowEntries     =   Matrix::OrderOperatorForRowEntries;
    type OrderOperatorForRowIndices     =   Matrix::OrderOperatorForRowIndices;
    type OrderOperatorForColumnEntries  =   Matrix::OrderOperatorForColumnEntries;
    type OrderOperatorForColumnIndices  =   Matrix::OrderOperatorForColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator { self.matrix.ring_operator() }
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries { self.matrix.order_operator_for_row_entries() }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices { self.matrix.order_operator_for_row_indices() }    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries { self.matrix.order_operator_for_column_entries() }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices { self.matrix.order_operator_for_column_indices() }   
}









//  TESTS
//  =========================================================================================================

//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {



} 

