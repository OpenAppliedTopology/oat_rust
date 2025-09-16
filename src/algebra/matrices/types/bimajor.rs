//! Stores two matrices, A and B; returns the rows of B as columns.
//! 
//! For example, calling `self.column(&index)` on this struct will return `B.row(&index)`, while calling
//! `self.row(&index)` will return `A.row(&index)`.
//! 
//! This module offers two data structures:
//! 
//! - `MatrixBimajor` holds a pair of oracles
//! - `MatrixBimajorData` holds a pair of objects that become oracles when references (such as [VecOfVec](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec).)

use derive_getters::{Getters, Dissolve};
use derive_new::new;

use crate::algebra::{matrices::{operations::MatrixOracleOperations, query::{MatrixAlgebra, MatrixOracle}}, vectors::entries::KeyValGet};

use std::fmt::Debug;


//  ==========================================================
//  BIMAJOR (REQUIRES TWO ORACLES)
//  ==========================================================

/// Stores two matrices, A and B; returns the rows of B as columns.
/// 
/// For example, calling `self.column(&index)` on this struct will return `B.row(&index)`, while calling
/// `self.row(&index)` will return `A.row(&index)`.
/// 
/// # See also
/// 
/// Also check out [MatrixBimajorData], which serves a similar function.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::algebra::matrices::operations::MatrixOracleOperations;   
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use oat_rust::algebra::matrices::types::bimajor::MatrixBimajor;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// 
/// use itertools::Itertools;        
///   
/// let matrix_data = vec![ 
///         vec![ (0, 1.),                          ], 
///         vec![ (0,-1.),  (1 ,1.),                ], 
///         vec![           (1,-1.),    (2, 1.),    ], 
///         vec![                       (2,-1.),    ],                 
///     ];
///      
/// // define a matrix of form
/// // |  1   0   0 |
/// // | -1   1   0 |
/// // |  0  -1   1 |
/// // |  0   0  -1 |        
/// let matrix_data = VecOfVec::new(matrix_data).ok().unwrap();
/// let matrix                  =   & matrix_data;
/// let num_rows_of_transpose   =   3;
/// let transpose               =   matrix.transpose_deep(num_rows_of_transpose).unwrap(); // a non-lazy transpose; copies the data
/// let bimajor                 =   MatrixBimajor{ matrix_rows: &matrix, matrix_columns: & transpose };
///    
/// // check that rows agree
/// for index in 0..4 {
///     assert_eq!(     matrix.row(&index).collect_vec(), 
///                     bimajor.row(&index).collect_vec(),     );
///     assert_eq!(     matrix.row_reverse(&index).collect_vec(), 
///                     bimajor.row_reverse(&index).collect_vec(),    );                            
///     assert_eq!(     matrix.row(&index).collect_vec(), 
///                     bimajor.row(&index).collect_vec(),            ); 
/// }
/// 
/// // check that minor  views agree
/// for index in 0..3 {
///     assert_eq!(     matrix.column(&index).collect_vec(), 
///                     bimajor.column(&index).collect_vec(),     );
///     assert_eq!(     matrix.column_reverse(&index).collect_vec(), 
///                     bimajor.column_reverse(&index).collect_vec(),    );                            
///     assert_eq!(     matrix.column(&index).collect_vec(), 
///                     bimajor.column(&index).collect_vec(),            );
/// }   
/// ```        
#[derive(new,Clone,Copy,Debug,Getters,Dissolve,Eq,PartialEq,Ord,PartialOrd)]
pub struct MatrixBimajor< MatrixRows, MatrixColumns >
{
    pub matrix_rows:  MatrixRows,
    pub matrix_columns:  MatrixColumns,
}


//  MATRIX ORACLE
//  ---------------------------------------------------------------------

impl < MatrixRows, MatrixColumns >

    MatrixOracle for 

    MatrixBimajor
        < MatrixRows, MatrixColumns > 

    where 
        MatrixRows:    MatrixOracle,
        MatrixColumns:    MatrixOracle< 
                                RowIndex        =   MatrixRows::ColumnIndex, 
                                ColumnIndex     =   MatrixRows::RowIndex,
                                Coefficient     =   MatrixRows::Coefficient,
                            >,    
{
    type Coefficient    =   MatrixRows::Coefficient;
    
    type RowIndex       =   MatrixRows::RowIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex    =   MatrixRows::ColumnIndex;       // The type of column indices
    
    type RowEntry       =   MatrixRows::RowEntry;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry    =   MatrixColumns::RowEntry;          // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    

    
    type Row =              // What you get when you ask for a row.
                            MatrixRows::Row;
    type RowReverse =       // What you get when you ask for a row with the order of entries reversed
                            MatrixRows::RowReverse;
    type Column =           // What you get when you ask for a column
                            MatrixColumns::Row;
    type ColumnReverse =    // What you get when you ask for a column with the order of entries reversed                             
                            MatrixColumns::RowReverse;

    /// Uses the row-major matrix to look up an entry.
    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { self.matrix_rows.structural_nonzero_entry( row, column ) }

    fn row(                     &self,  index: &Self::RowIndex    )   -> Self::Row
        { self.matrix_rows.row( index ) }
    fn row_reverse(             &self,  index: &Self::RowIndex    )   -> Self::RowReverse
        { self.matrix_rows.row_reverse( index ) }
    
    fn column(                  &self,  index: &Self::ColumnIndex )   -> Self::Column
        { self.matrix_columns.row(index) }
    // fn column_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex >
    //     { self.matrix_columns.row_result(index) }    
    fn column_reverse(          &self,  index: &Self::ColumnIndex )   -> Self::ColumnReverse
        { self.matrix_columns.row_reverse(index) }
    // fn column_reverse_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex >
    //     { self.matrix_columns.row_reverse_result(index) }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool
        { self.matrix_rows.has_row_for_index( index ) }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { self.matrix_columns.has_row_for_index( index ) }    
        
}






//  MATRIX ALGEBRA
//  ---------------------------------------------------------------------

// impl < MatrixRows, MatrixColumns >
impl < MatrixRows, MatrixColumns, RowIndex, ColumnIndex, Coefficient, RowEntry, ColumnEntry >

    MatrixAlgebra for 

    MatrixBimajor
        < MatrixRows, MatrixColumns > 

    // where 
    //      MatrixRows:            MatrixAlgebra,
    //      MatrixColumns:         MatrixAlgebra<
    //                                 RowIndex        =   < MatrixRows as MatrixOracle >::ColumnIndex, 
    //                                 ColumnIndex     =   < MatrixRows as MatrixOracle >::RowIndex,
    //                                 Coefficient     =   < MatrixRows as MatrixOracle >::Coefficient,
    //                             >,    

    where 
        MatrixRows:         MatrixAlgebra<
                                RowIndex        =   RowIndex, 
                                ColumnIndex     =   ColumnIndex, 
                                Coefficient     =   Coefficient, 
                                RowEntry        =   RowEntry,
                            >,
        MatrixColumns:      MatrixAlgebra< 
                                RowIndex        =   ColumnIndex, 
                                ColumnIndex     =   RowIndex,
                                Coefficient     =   Coefficient,
                                RowEntry        =   ColumnEntry,
                            >,    
        Coefficient:            Clone + Debug + PartialEq,
        RowIndex:               Clone + Debug + Eq,
        ColumnIndex:            Clone + Debug + Eq,
        RowEntry:               Clone + Debug + PartialEq + KeyValGet< Key = ColumnIndex, Val = Coefficient >,
        ColumnEntry:            Clone + Debug + PartialEq + KeyValGet< Key = RowIndex, Val = Coefficient >,     
{
    type RingOperator                                   =   < MatrixRows as MatrixAlgebra >::RingOperator;

    type OrderOperatorForRowEntries                     =   < MatrixRows as MatrixAlgebra >::OrderOperatorForRowEntries;

    type OrderOperatorForRowIndices                     =   < MatrixColumns as MatrixAlgebra >::OrderOperatorForColumnIndices;

    type OrderOperatorForColumnEntries                  =   < MatrixColumns as MatrixAlgebra >::OrderOperatorForRowEntries;

    type OrderOperatorForColumnIndices                  =   < MatrixRows as MatrixAlgebra >::OrderOperatorForColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator {
        self.matrix_rows.ring_operator()
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.matrix_rows.order_operator_for_row_entries()
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.matrix_columns.order_operator_for_column_indices()
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.matrix_columns.order_operator_for_row_entries()
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.matrix_rows.order_operator_for_column_indices()
    }
}







//  MATRIX ORACLE OPERATIONS
//  ---------------------------------------------------------------------



impl < MatrixRows, MatrixColumns >

    MatrixOracleOperations for 
    
    MatrixBimajor< MatrixRows, MatrixColumns >

{}    







//  ==========================================================
//  BIMAJOR DATA (REQUIRES **DATA** FOR TWO ORACLES)
//  ==========================================================





//  ==========================================================
//  BIMAJOR (REQUIRES TWO ORACLES)
//  ==========================================================

/// Similar to [MatrixBimajor], but works with objects that only become oracles when referenced.
/// 
/// Stores two matrices, A and B; returns the rows of B as columns.
/// 
/// For example, calling `self.column(&index)` on this struct will return `(&B).row(&index)`, while calling
/// `self.row(&index)` will return `(&A).row(&index)`.
/// 
/// The key distinction with [MatrixBimajor] is as follows:
/// - [MatrixBimajor] 
///   - is itself an oracle
///   - holds two oracles, A and B
/// - [MatrixBimajorData]
///   - is not an oracle, but [& MatrixBimajorData](MatrixBimajorData) is
///   - holds two objects A and B such that &A and &B *are* oracles (even if A and B are not)
/// 
/// # Motivation
/// 
/// We use this matrix type for several operations related to U-match factorization.
/// In this use case, the [MatrixBimajorData] contains two [VecOfVec](crate::algebra::matrices::types::sorted::VecOfVec) 
/// objects; these don't implement the [MatrixOracle] trait themselves (and indeed cannot, due to lifetime considerations),
/// but references to them do!
/// 
/// Were it not for this important application, we would likely delete [MatrixBimajorData] from the library.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::algebra::matrices::query::MatrixOracle;     
/// use oat_rust::algebra::matrices::types::bimajor::MatrixBimajorData;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// 
/// use itertools::Itertools;   
/// 
/// 
/// // data for a matrix:
/// 
/// // |  1   0   0 |
/// // | -1   1   0 |
/// // |  0  -1   1 |
/// // |  0   0  -1 |           
/// let matrix_data = vec![ 
///         vec![ (0,1.),                       ], 
///         vec![ (0,-1.),  (1,1.),             ], 
///         vec![           (1,-1.), (2,1.),    ], 
///         vec![                    (2,-1.),   ],                 
///     ];
/// 
/// // define matrices
/// let matrix_data = VecOfVec::new(matrix_data).ok().unwrap();
/// let transpose_data   = matrix_data.transpose_deep(3).unwrap(); // a non-lazy transpose; copies the data
/// let bimajor_data     =   MatrixBimajorData{ matrix_rows_data: matrix_data.clone(), matrix_columns_data: transpose_data };
/// let matrix  = & matrix_data;
/// 
/// // check that rows agree
/// for index in 0..4 {
///     assert_eq!(     matrix.row(&index).collect_vec(), 
///                     ( & bimajor_data ).row(&index).collect_vec(),     );
///     assert_eq!(     matrix.row_reverse(&index).collect_vec(), 
///                     ( & bimajor_data ).row_reverse(&index).collect_vec(),    );                            
///     assert_eq!(     matrix.row(&index).collect_vec(), 
///                     ( & bimajor_data ).row(&index).collect_vec(),            ); 
/// }
/// 
/// // check that minor  views agree
/// for index in 0..3 {
///     assert_eq!(     matrix.column(&index).collect_vec(), 
///                     ( & bimajor_data ).column(&index).collect_vec(),     );
///     assert_eq!(     matrix.column_reverse(&index).collect_vec(), 
///                     ( & bimajor_data ).column_reverse(&index).collect_vec(),    );                            
///     assert_eq!(     matrix.column(&index).collect_vec(), 
///                     ( & bimajor_data ).column(&index).collect_vec(),            );
/// }                         
/// ```        
#[derive(new,Clone,Copy,Debug,Getters,Dissolve,Eq,PartialEq)]
pub struct MatrixBimajorData< MatrixRows, MatrixColumns >
{
    pub matrix_rows_data:       MatrixRows,
    pub matrix_columns_data:    MatrixColumns,
}


//  MATRIX ORACLE
//  ---------------------------------------------------------------------

impl < 'a, MatrixRows, MatrixColumns, RowIndex, ColumnIndex, Coefficient, RowEntry, ColumnEntry >

    MatrixOracle for 

    &'a MatrixBimajorData
        < MatrixRows, MatrixColumns > 

    where 
        &'a MatrixRows:     MatrixOracle<
                                RowIndex        =   RowIndex, 
                                ColumnIndex     =   ColumnIndex, 
                                Coefficient     =   Coefficient, 
                                RowEntry        =   RowEntry,
                            >,
        &'a MatrixColumns:  MatrixOracle< 
                                RowIndex        =   ColumnIndex, 
                                ColumnIndex     =   RowIndex,
                                Coefficient     =   Coefficient,
                                RowEntry        =   ColumnEntry,
                            >,    
        Coefficient:            Clone + Debug + PartialEq,
        RowIndex:               Clone + Debug + Eq,
        ColumnIndex:            Clone + Debug + Eq,
        RowEntry:               Clone + Debug + PartialEq + KeyValGet< Key = ColumnIndex, Val = Coefficient >,
        ColumnEntry:            Clone + Debug + PartialEq + KeyValGet< Key = RowIndex, Val = Coefficient >,        
{
    type Coefficient        =   Coefficient;    
    type RowIndex           =   RowIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   ColumnIndex;       // The type of column indices    
    type RowEntry           =   RowEntry;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   ColumnEntry;          // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`                            
    type Row                =   <&'a MatrixRows as MatrixOracle >::Row;       // What you get when you ask for a row.                            
    type RowReverse         =   <&'a MatrixRows as MatrixOracle >::RowReverse;// What you get when you ask for a row with the order of entries reversed                            
    type Column             =   <&'a MatrixColumns as MatrixOracle >::Row;// What you get when you ask for a column                            
    type ColumnReverse      =   <&'a MatrixColumns as MatrixOracle >::RowReverse;// What you get when you ask for a column with the order of entries reversed                                                         
                        
    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { ( & self.matrix_rows_data ).structural_nonzero_entry( row, column ) }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool
        { ( & self.matrix_rows_data ).has_row_for_index( index ) }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { ( & self.matrix_columns_data ).has_row_for_index( index ) }      
    fn row(                     &self,  index: &Self::RowIndex    )   -> Self::Row
        { ( & self.matrix_rows_data ).row( index ) }
    fn row_reverse(             &self,  index: &Self::RowIndex    )   -> Self::RowReverse
        { ( & self.matrix_rows_data ).row_reverse( index ) }
    fn column(                  &self,  index: &Self::ColumnIndex )   -> Self::Column
        { ( & self.matrix_columns_data ).row(index) }
    fn column_reverse(          &self,  index: &Self::ColumnIndex )   -> Self::ColumnReverse
        { ( & self.matrix_columns_data ).row_reverse(index) }        
}



//  MATRIX ALGEBRA
//  ---------------------------------------------------------------------

impl < 'a, MatrixRows, MatrixColumns, RowIndex, ColumnIndex, Coefficient, RowEntry, ColumnEntry >

    MatrixAlgebra for 

    &'a MatrixBimajorData
        < MatrixRows, MatrixColumns > 

    where 
        &'a MatrixRows:     MatrixAlgebra<
                                RowIndex        =   RowIndex, 
                                ColumnIndex     =   ColumnIndex, 
                                Coefficient     =   Coefficient, 
                                RowEntry        =   RowEntry,
                            >,
        &'a MatrixColumns:  MatrixAlgebra< 
                                RowIndex        =   ColumnIndex, 
                                ColumnIndex     =   RowIndex,
                                Coefficient     =   Coefficient,
                                RowEntry        =   ColumnEntry,
                            >,    
        Coefficient:            Clone + Debug + PartialEq,
        RowIndex:               Clone + Debug + Eq,
        ColumnIndex:            Clone + Debug + Eq,
        RowEntry:               Clone + Debug + PartialEq + KeyValGet< Key = ColumnIndex, Val = Coefficient >,
        ColumnEntry:            Clone + Debug + PartialEq + KeyValGet< Key = RowIndex, Val = Coefficient >,        
{
    type RingOperator                                   =   < &'a MatrixRows as MatrixAlgebra >::RingOperator;

    type OrderOperatorForRowEntries                     =   < &'a MatrixRows as MatrixAlgebra >::OrderOperatorForRowEntries;

    type OrderOperatorForRowIndices                     =   < &'a MatrixColumns as MatrixAlgebra >::OrderOperatorForColumnIndices;

    type OrderOperatorForColumnEntries                  =   < &'a MatrixColumns as MatrixAlgebra >::OrderOperatorForRowEntries;

    type OrderOperatorForColumnIndices                  =   < &'a MatrixRows as MatrixAlgebra >::OrderOperatorForColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator {
        (& self.matrix_rows_data).ring_operator()
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        (&self.matrix_rows_data).order_operator_for_row_entries()
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        (&self.matrix_columns_data).order_operator_for_column_indices()
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        (&self.matrix_columns_data).order_operator_for_row_entries()
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        (&self.matrix_rows_data).order_operator_for_column_indices()
    }
}







//  MATRIX ORACLE OPERATIONS
//  ---------------------------------------------------------------------



impl < 'a, MatrixRows, MatrixColumns >

    MatrixOracleOperations for 
    
    &'a MatrixBimajor< MatrixRows, MatrixColumns >

{}    








//  ==========================================================
//  TESTS
//  ==========================================================







#[cfg(test)]
mod tests {
    

    



    #[test]
    fn test_matrix_bimajor() {   

        use crate::algebra::matrices::operations::MatrixOracleOperations;   
        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::matrices::types::bimajor::MatrixBimajor;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;

        use itertools::Itertools;   


        // data for a matrix:

        // |  1   0   0 |
        // | -1   1   0 |
        // |  0  -1   1 |
        // |  0   0  -1 |           
        let matrix_data = vec![ 
                vec![ (0,1.),                       ], 
                vec![ (0,-1.),  (1,1.),             ], 
                vec![           (1,-1.), (2,1.),    ], 
                vec![                    (2,-1.),   ],                 
            ];
        
        // define matrices
        let matrix_data = VecOfVec::new(matrix_data);
        let matrix = & matrix_data.unwrap();
        let transpose   = matrix.transpose();
        let bimajor     =   MatrixBimajor{ matrix_rows: &matrix, matrix_columns: & transpose };
        
        // check that rows agree
        for index in 0..4 {
            assert_eq!(     matrix.row(&index).collect_vec(), 
                            bimajor.row(&index).collect_vec(),     );
            assert_eq!(     matrix.row_reverse(&index).collect_vec(), 
                            bimajor.row_reverse(&index).collect_vec(),    );                            
            assert_eq!(     matrix.row(&index).collect_vec(), 
                            bimajor.row(&index).collect_vec(),            ); 
        }

        // check that minor  views agree
        for index in 0..3 {
            let index = &index;
            assert_eq!(     matrix.column(&index).collect_vec(), 
                            bimajor.column(&index).collect_vec(),     );
            assert_eq!(     matrix.column_reverse(&index).collect_vec(), 
                            bimajor.column_reverse(&index).collect_vec(),    );                            
            assert_eq!(     matrix.column(&index).collect_vec(), 
                            bimajor.column(&index).collect_vec(),            );
        }                    
    }  








    #[test]
    fn test_matrix_bimajor_data() {   
 
        use crate::algebra::matrices::query::MatrixOracle;       
        use crate::algebra::matrices::types::bimajor::MatrixBimajorData;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;

        use itertools::Itertools;   


        // data for a matrix:

        // |  1   0   0 |
        // | -1   1   0 |
        // |  0  -1   1 |
        // |  0   0  -1 |           
        let matrix_data = vec![ 
                vec![ (0,1.),                       ], 
                vec![ (0,-1.),  (1,1.),             ], 
                vec![           (1,-1.), (2,1.),    ], 
                vec![                    (2,-1.),   ],                 
            ];
        
        // define matrices
        let matrix_data = VecOfVec::new(matrix_data).unwrap();
        let transpose_data   = matrix_data.transpose_deep(3).unwrap(); // a non-lazy transpose; copies the data
        let bimajor_data     =   MatrixBimajorData{ matrix_rows_data: matrix_data.clone(), matrix_columns_data: transpose_data };
        let matrix  = & matrix_data;

        // check that rows agree
        for index in 0..4 {
            assert_eq!(     matrix.row(&index).collect_vec(), 
                            ( & bimajor_data ).row(&index).collect_vec(),     );
            assert_eq!(     matrix.row_reverse(&index).collect_vec(), 
                            ( & bimajor_data ).row_reverse(&index).collect_vec(),    );                            
            assert_eq!(     matrix.row(&index).collect_vec(), 
                            ( & bimajor_data ).row(&index).collect_vec(),            ); 
        }

        // check that minor  views agree
        for index in 0..3 {
            assert_eq!(     matrix.column(&index).collect_vec(), 
                            ( & bimajor_data ).column(&index).collect_vec(),     );
            assert_eq!(     matrix.column_reverse(&index).collect_vec(), 
                            ( & bimajor_data ).column_reverse(&index).collect_vec(),    );                            
            assert_eq!(     matrix.column(&index).collect_vec(), 
                            ( & bimajor_data ).column(&index).collect_vec(),            );
        }                    
    }  





}

