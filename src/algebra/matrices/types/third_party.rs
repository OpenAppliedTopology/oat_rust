//! Wrappers for third-party matrix types
//! 
//! # Tips on construction
//! 
//! The best reference for third-party matrix structs is their documentation.  However,
//! here are a few common construction methods
//! 
//! Build a dense 2 dimensional ndarray from a vector of vectors
//! 
//! ```
//! let mat = vec![ vec![1;2]; 3 ]; // a vector of vectors
//! let mat = ndarray::Array2::from_shape_vec((3, 2), mat);
//! ```
//! 
//! Build a sparse 2x2 identity matrix two different ways
//! 
//! ```
//! use oat_rust::algebra::matrices::types::third_party::IntoCSR;
//! 
//! let eye =   vec![ vec![(0,1)], vec![(1,1)] ].into_csr(2,2);
//! 
//! let eye =   vec![ (0,0,1), (1,1,1) ].into_csr(2,2);
//! ```

use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::matrices::query::{ MatrixOracle, column_helper::SparseColumn};

use std::iter::Rev;
use std::ops::{Deref, Range};
use std::sync::Arc;
use std::fmt::Debug;

use sprs::{ TriMat, CsMatBase, SpIndex, TriMatBase };
use sprs::vec::{VectorIterator, IntoSparseVecIter};


//  =====================================================================================
//  & sprs::CsMatBase
//  =====================================================================================


//  VIEW TYPE
//  ------------------------

/// A wrapper for the `sprs` struct `VectorIterator<'a, N, I >`, which returns
/// `(usize, N)` instead of `(usize, &N)`
#[derive(Clone)]
pub struct SprsRowIteratorCloned< 'a, N, I > {
    vector_iterator: VectorIterator<'a, N, I >
}

impl < 'a, N, I >

    Iterator for 

    SprsRowIteratorCloned
        < 'a, N, I >

    where
        I: SpIndex,
        N: Clone,          
{
    type Item = ( usize, N );
    fn next(&mut self) -> Option<Self::Item> {
        self.vector_iterator.next().map(|(k,v)| (k,v.clone()))
    }
}  


/// Runs over the entries of an outer view of a `sprs::CsMatBase`
/// 
/// **Note** "Outer views" of a CSR matrix are rows.  Outer views of a
/// CSC matrix are columns.
/// 
/// **In OAT we currently regard every instance of `sprs::CsMatBase` as
/// a CSR matrix, even if the `sprs` module regards it as a CSC matrix
/// (that is, even if `matrix.storage()` returns `sprs::CompressedStorage::CSR`)**.
/// This is equivalent to saying that we always regard "outer views" of a `CsMatBase`
/// as rows.
#[derive(Clone,Debug,Eq,PartialEq)]
pub struct SprsRowIterator< 'a, N, I, > {
    ordinals:           Range< usize >,
    coefficients:       &'a [N],
    column_indices:     &'a [I],
}

impl < 'a, N, I, >

    Iterator for 

    SprsRowIterator< 'a, N, I, >

    where 
        I:      SpIndex,
        N:      Clone,
{
    type Item   =   (usize, N);
    fn next( &mut self ) -> Option< Self::Item > {
        self.ordinals.next().map(
            |i| 
            ( self.column_indices[i].index(), self.coefficients[i].clone() )
        )
    }
}    

/// Runs over the entries of an outer view of a `sprs::CsMatBase` **in reverse order**
/// 
/// **Note** "Outer views" of a CSR matrix are rows.  Outer views of a
/// CSC matrix are columns.
/// 
/// **In OAT we currently regard every instance of `sprs::CsMatBase` as
/// a CSR matrix, even if the `sprs` module regards it as a CSC matrix
/// (that is, even if `matrix.storage()` returns `sprs::CompressedStorage::CSR`)**.
/// This is equivalent to saying that we always regard "outer views" of a `CsMatBase`
/// as rows.
#[derive(Clone,Debug)]
pub struct SprsRowReverseIterator< 'a, N, I > {
    ordinals:           Rev<Range<usize>>,
    coefficients:       &'a [N],
    column_indices:     &'a [I],
}

impl < 'a, N, I >

    Iterator for 

    SprsRowReverseIterator< 'a, N, I >

    where 
        I:      SpIndex,
        N:      Clone,    
{
    type Item   =   (usize, N);
    fn next( &mut self ) -> Option< Self::Item > {
        self.ordinals.next().map(
            |i| {
                let i = i.index(); // convert the index to usize
                ( self.column_indices[i].index(), self.coefficients[i].clone() ) // return the column index and coefficient
            }
        )
    }
}  


//  -------------------------------------------------------------------------------------
//  & sprs::CsMatBase
//  -------------------------------------------------------------------------------------


//  MATRIX ORACLE
//  -------------


/// This implements the [MatrixOracle] trait for [CsMatBase].
/// 
/// # Example
///
/// ```
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use oat_rust::algebra::matrices::types::third_party::IntoCSR;
/// use oat_rust::algebra::matrices::debug::matrix_oracle_is_internally_consistent;
/// 
/// 
/// // data for a matrix:
/// // |  1   0 |
/// // | -1   0 |
/// // |  0   0 |
/// // |  0   0 |           
/// let matrix = vec![ 
///         vec![ (0, 1)   ], 
///         vec![ (0,-1)   ], 
///         vec![          ], 
///         vec![          ],                 
///     ];
/// 
/// // wrap the data in a VecOfVec sparse matrix struct (see documentation for details)
/// let matrix = matrix.into_csr(4,2);
/// 
/// // check that matrix is internally consistent
/// let sorted_row_indices = 0..4;
/// let sorted_column_indices = 0..2;
/// assert!(
///     matrix_oracle_is_internally_consistent( &matrix, sorted_row_indices, sorted_column_indices)
/// );                
/// 
/// // get some sparse columns
/// let mut column_0        =   (& matrix ).column(&0);
/// let mut column_0_rev    =   (& matrix ).column_reverse(&0);
/// let mut column_1        =   (& matrix ).column(&1);
/// let mut column_1_rev    =   (& matrix ).column_reverse(&1);
/// 
/// // check the columns are correct
/// assert_eq!( column_0.next(),        Some( (0,  1) )    );
/// assert_eq!( column_0.next(),        Some( (1, -1) )    );
/// assert_eq!( column_0.next(),        None                );    
/// 
/// assert_eq!( column_0_rev.next(),    Some( (1, -1) )    );
/// assert_eq!( column_0_rev.next(),    Some( (0,  1) )    );
/// assert_eq!( column_0_rev.next(),    None                );   
/// 
/// assert_eq!( column_1.next(),        None        );
/// assert_eq!( column_1_rev.next(),    None        );     
/// ```
impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >

    MatrixOracle for

    &'a CsMatBase
        < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

    where
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,
        N: Clone + Debug + Eq,
{           

    type Coefficient            =   N; // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =   usize; // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   usize; // The type of column indices    
    type RowEntry               =   (usize,N); // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   (usize,N); // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =   SprsRowIterator< 'a, N, I, >;  // What you get when you ask for a row.
    type RowReverse             =   SprsRowReverseIterator< 'a, N, I, >;  // What you get when you ask for a row with the order of entries reversed
    type Column                 =   SparseColumn< // What you get when you ask for a row.
                                        &'a CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr >,
                                        Range< usize >,
                                    >;
    type ColumnReverse          =   SparseColumn<  // What you get when you ask for a column with the order of entries reversed                             
                                        &'a CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr >,
                                        Rev< Range< usize > >,
                                    >;   

    fn structural_nonzero_entry(                   &   self, row: & Self::RowIndex, col: & Self::ColumnIndex ) ->  Option< Self::Coefficient > { 
        if self.is_csc() {
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats.");  
        }
        self.get( * row, * col ).cloned()
    }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { *index < self.cols() }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool 
        { *index < self.rows() }
    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row {   
        if *index >= self.rows() {
            panic!("A program has called `self.row( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows.", index, self.rows() );
        }    
        self.row_result( index ).unwrap()  
    } 
    fn row_result(                 &   self, index: & Self::RowIndex   )   -> Result< Self::Row, Self::RowIndex > {         
        if self.is_csc() { 
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if *index >= self.rows() {
            Err( index.clone() )
        } else {
            let ordinals_raw: Range<Iptr>   =   self.indptr().outer_inds(*index);
            let ordinals: Range< usize >    =   ordinals_raw.start.index() .. ordinals_raw.end.index();
            let column_indices: &[I]        =   self.indices();
            let coefficients: &[N]          =   self.data();
            Ok( SprsRowIterator{ ordinals, coefficients, column_indices } )
        }
    }
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse {  
        if *index >= self.rows() {
            panic!("A program has called `self.row( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows.", index, self.rows() );
        }         
        self.row_reverse_result( index ).unwrap()  
    }     
    fn row_reverse_result(         &   self, index: &Self::RowIndex   )   -> Result< Self::RowReverse, Self::RowIndex >
    { 
        if self.is_csc() { 
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if *index >= self.rows() {
            Err( index.clone() ) // the index is too big
        } else {
            let ordinals_raw: Range<Iptr>       =   self.indptr().outer_inds( *index );
            let ordinals: Rev< Range< usize > > =   ( ordinals_raw.start.index() .. ordinals_raw.end.index() ).rev();
            let column_indices: &[I]            =   self.indices();
            let coefficients: &[N]              =   self.data();
            Ok( SprsRowReverseIterator{ ordinals, coefficients, column_indices } )
        }
    }    
    /// Get a column
    /// 
    /// # Example
    ///
    /// ```
    /// use oat_rust::algebra::matrices::query::MatrixOracle;
    /// use oat_rust::algebra::matrices::types::third_party::IntoCSR;
    /// use oat_rust::algebra::matrices::debug::matrix_oracle_is_internally_consistent;
    /// 
    /// 
    /// // data for a matrix:
    /// // |  1   0 |
    /// // | -1   0 |
    /// // |  0   0 |
    /// // |  0   0 |           
    /// let matrix = vec![ 
    ///         vec![ (0, 1)   ], 
    ///         vec![ (0,-1)   ], 
    ///         vec![          ], 
    ///         vec![          ],                 
    ///     ];
    /// 
    /// // wrap the data in a VecOfVec sparse matrix struct (see documentation for details)
    /// let matrix = matrix.into_csr(4,2);
    /// 
    /// // check that matrix is internally consistent
    /// let sorted_row_indices = 0..4;
    /// let sorted_column_indices = 0..2;
    /// assert!(
    ///     matrix_oracle_is_internally_consistent( &matrix, sorted_row_indices, sorted_column_indices)
    /// );                
    /// 
    /// // get some sparse columns
    /// let mut column_0        =   (& matrix ).column(&0);
    /// let mut column_0_rev    =   (& matrix ).column_reverse(&0);
    /// let mut column_1        =   (& matrix ).column(&1);
    /// let mut column_1_rev    =   (& matrix ).column_reverse(&1);
    /// 
    /// // check the columns are correct
    /// assert_eq!( column_0.next(),        Some( (0,  1) )    );
    /// assert_eq!( column_0.next(),        Some( (1, -1) )    );
    /// assert_eq!( column_0.next(),        None                );    
    /// 
    /// assert_eq!( column_0_rev.next(),    Some( (1, -1) )    );
    /// assert_eq!( column_0_rev.next(),    Some( (0,  1) )    );
    /// assert_eq!( column_0_rev.next(),    None                );   
    /// 
    /// assert_eq!( column_1.next(),        None        );
    /// assert_eq!( column_1_rev.next(),    None        );     
    /// ```
    fn column(                  &   self, index: &Self::ColumnIndex)   -> Self::Column    { 
        if *index >= self.cols() {
            panic!("A program has called `self.col( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} columns", index, self.rows() );
        }
        self.column_result(index).unwrap() 
    }
    fn column_result(              &   self, index: &Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex > {
        if self.is_csc() { 
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if *index >= self.cols() {
            return Err( index.clone() ) // if the index is too great, then return None
        } else {
            let row_index_iterator  =   0 .. self.rows();
            let matrix              =   self;
            let column_index        =   index.clone();
            Ok( SparseColumn{ matrix, row_index_iterator, column_index } )
        }
    }   
    /// Get a column with entries in reverse order
    /// 
    /// # Example
    ///
    /// ```
    /// use oat_rust::algebra::matrices::query::MatrixOracle;
    /// use oat_rust::algebra::matrices::types::third_party::IntoCSR;
    /// use oat_rust::algebra::matrices::debug::matrix_oracle_is_internally_consistent;
    /// 
    /// 
    /// // data for a matrix:
    /// // |  1   0 |
    /// // | -1   0 |
    /// // |  0   0 |
    /// // |  0   0 |           
    /// let matrix = vec![ 
    ///         vec![ (0, 1)   ], 
    ///         vec![ (0,-1)   ], 
    ///         vec![          ], 
    ///         vec![          ],                 
    ///     ];
    /// 
    /// // wrap the data in a VecOfVec sparse matrix struct (see documentation for details)
    /// let matrix = matrix.into_csr(4,2);
    /// 
    /// // check that matrix is internally consistent
    /// let sorted_row_indices = 0..4;
    /// let sorted_column_indices = 0..2;
    /// assert!(
    ///     matrix_oracle_is_internally_consistent( &matrix, sorted_row_indices, sorted_column_indices)
    /// );                
    /// 
    /// // get some sparse columns
    /// let mut column_0        =   (& matrix ).column(&0);
    /// let mut column_0_rev    =   (& matrix ).column_reverse(&0);
    /// let mut column_1        =   (& matrix ).column(&1);
    /// let mut column_1_rev    =   (& matrix ).column_reverse(&1);
    /// 
    /// // check the columns are correct
    /// assert_eq!( column_0.next(),        Some( (0,  1) )    );
    /// assert_eq!( column_0.next(),        Some( (1, -1) )    );
    /// assert_eq!( column_0.next(),        None                );    
    /// 
    /// assert_eq!( column_0_rev.next(),    Some( (1, -1) )    );
    /// assert_eq!( column_0_rev.next(),    Some( (0,  1) )    );
    /// assert_eq!( column_0_rev.next(),    None                );   
    /// 
    /// assert_eq!( column_1.next(),        None        );
    /// assert_eq!( column_1_rev.next(),    None        );     
    /// ```   
    fn column_reverse(          &   self, index: &Self::ColumnIndex)   -> Self::ColumnReverse       {
        if *index > self.cols() {
            panic!("A program has called `self.col( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows", index, self.rows() );
        }
        self.column_reverse_result(index).unwrap() 
    }        
    fn column_reverse_result(      &   self, index: &Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex > {
        if self.is_csc() { 
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if *index >= self.cols() {
            return Err( index.clone() ) // if the index is too great, then return None
        } else {
            let row_index_iterator  =   ( 0 .. self.rows() ).rev();
            let matrix              =   self;
            let column_index        =   index.clone();
            Ok( SparseColumn{ matrix, row_index_iterator, column_index } )        
        }
    }

} 






//  MATRIX ORACLE OPERATIONS
//  -------------------------------------------------------------------------------------


impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >

    MatrixOracleOperations for

    &'a CsMatBase
        < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

    where
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,
        N: Clone + Debug + Eq,        
{}        













//  =====================================================================================
//  Arc< sprs::CsMatBase >
//  =====================================================================================


//  ROW TYPE
//  ------------------------


/// Iterates over some entries of a `sprs::CsMatBase` matrix
/// 
/// # Internal structure
/// 
/// This object stores an arc to a matrix `Arc< CsMatBase >` and a an iterator `linear_indices` 
/// of type `LinearIndices`. Internally, the matrix stores its nonzero entries in vectors
/// 
/// - `matrix.indices()` which stores the row or column indices of the structural nonzero entries.
///   In OAT we typically regard a `sprs::CsMatBase` matrix as a CSR matrix, so these indices
///   are typically column indices.
/// - `matrix.data()` which stores the coefficients of the structural nonzero entries.
/// 
/// When the iterator is called, it first 
/// calls `self.linear_indices.next()` to get the next index, and then uses that index to
/// look up the corresponding entry in `matrix.indices()` and `matrix.data()`.
/// 
/// 
/// # Key features
/// 
/// - the `CsMatBase` matrix is wrapped in an `Arc` (not a reference)
/// - returns items of type `(usize, N)` instead of `(usize, &N)`
/// 
/// # Errors
/// 
/// The constructor for this object doesn't check that each index in the `linear_indices` iterator
/// is less than the number of entries in the matrix.
#[derive(Clone,Debug,Eq,PartialEq)]
pub struct MatrixEntryIteratorArc
            < N, I, IptrStorage, IndStorage, DataStorage, Iptr, LinearIndices > 
    where
        N: Clone,
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,
        LinearIndices: Iterator<Item = usize>,  
{
    linear_indices:  LinearIndices,
    matrix: Arc< CsMatBase< N, I, IptrStorage, IndStorage, DataStorage, Iptr > >
}

impl < N, I, IptrStorage, IndStorage, DataStorage, Iptr, LinearIndices > 

    Iterator for 

    MatrixEntryIteratorArc
        < N, I, IptrStorage, IndStorage, DataStorage, Iptr, LinearIndices > 

    where
        N: Clone,
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,    
        LinearIndices: Iterator<Item = usize>,  
{
    type Item = ( usize, N );
    fn next(&mut self) -> Option<Self::Item> {
        self.linear_indices.next().map( 
                |x| 
                ( self.matrix.indices()[ x ].index(), self.matrix.data()[ x ].clone() ) 
            )
    }
}  








//  MATRIX ORACLE
//  -------------

impl < N, I, IptrStorage, IndStorage, DataStorage, Iptr >

    MatrixOracle for

    Arc<
        CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 
    >

    where
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,
        N: Clone + Debug + Eq,
{           

    type Coefficient            =   N; // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =   usize; // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   usize; // The type of column indices    
    type RowEntry               =   (usize,N); // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   (usize,N); // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =   MatrixEntryIteratorArc < N, I, IptrStorage, IndStorage, DataStorage, Iptr, Range<usize> >;  // What you get when you ask for a row.
    type RowReverse             =   std::vec::IntoIter<(usize,N)>;  // What you get when you ask for a row with the order of entries reversed
    type Column                 =   SparseColumn< 
                                        Arc< CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr > >,
                                        Range< usize >,
                                    >;
    type ColumnReverse          =   SparseColumn<  
                                        Arc< CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr > >,
                                        Rev< Range< usize > >,
                                    >;   

    fn structural_nonzero_entry(                   &   self, row: & Self::RowIndex, col: & Self::ColumnIndex ) ->  Option< Self::Coefficient > { 
        if self.is_csc() {
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats.");  
        }
        self.get( * row, * col ).cloned()
    }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { *index < self.cols() }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool 
        { *index < self.rows() }
    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row {   
        if *index >= self.rows() {
            panic!("A program has called `self.row( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows.", index, self.rows() );
        }    
        self.row_result( index ).unwrap()  
    } 
    fn row_result(                 &   self, index: & Self::RowIndex   )   -> Result< Self::Row, Self::RowIndex > {         
        if self.is_csc() { 
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if *index >= self.rows() {
            Err( index.clone() )
        } else {
            let range = self.indptr().outer_inds_sz( *index );
            Ok( MatrixEntryIteratorArc{ matrix: self.clone(), linear_indices: range } )
            // let ordinals_raw: Range<Iptr>   =   self.indptr().outer_inds( * index );
            // let ordinals: Range< usize >    =   ordinals_raw.start.index() .. ordinals_raw.end.index();
            // let column_indices: &[I]        =   self.indices();
            // let coefficients: &[N]          =   self.data();
            // Ok( SprsRowIterator{ ordinals, coefficients, column_indices } )
        }
    }
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse {  
        if *index >= self.rows() {
            panic!("A program has called `self.row( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows.", index, self.rows() );
        }         
        self.row_reverse_result( index ).unwrap()  
    }     
    fn row_reverse_result(         &   self, index: &Self::RowIndex   )   -> Result< Self::RowReverse, Self::RowIndex >
    { 
        if self.is_csc() { 
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if *index >= self.rows() {
            Err( index.clone() ) // the index is too big
        } else {
            let mut vec = self.row(index).collect::<Vec<_>>();
            ( &mut vec ).reverse();
            Ok( vec.into_iter() )
        }
    }    
    /// Get a column
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
    fn column(                  &   self, index: &Self::ColumnIndex)   -> Self::Column    { 
        if *index >= self.cols() {
            panic!("A program has called `self.col( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} columns", index, self.rows() );
        }
        self.column_result(index).unwrap() 
    }
    fn column_result(              &   self, index: &Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex > {
        if self.is_csc() { 
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if *index >= self.cols() {
            return Err( index.clone() ) // if the index is too great, then return None
        } else {
            let row_index_iterator  =   0 .. self.rows();
            let matrix              =   self;
            let column_index        =   index.clone();
            Ok( SparseColumn{ matrix: matrix.clone(), row_index_iterator, column_index } )
        }
    }   
    /// Get a column with entries in reverse order
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::query::MatrixOracle;
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
    /// let matrix = VecOfVec::new(matrix).ok().unwrap();
    /// 
    /// // get some sparse columns
    /// let mut column_0        =   (& matrix ).column(&0);
    /// let mut column_0_rev    =   (& matrix ).column_reverse(&0);
    /// let mut column_1        =   (& matrix ).column(&1);
    /// let mut column_1_rev    =   (& matrix ).column_reverse(&1);
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
    fn column_reverse(          &   self, index: &Self::ColumnIndex)   -> Self::ColumnReverse       {
        if *index > self.cols() {
            panic!("A program has called `self.col( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows", index, self.rows() );
        }
        self.column_reverse_result(index).unwrap() 
    }        
    fn column_reverse_result(      &   self, index: &Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex > {
        if self.is_csc() { 
            panic!("A program has called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if *index >= self.cols() {
            return Err( index.clone() ) // if the index is too great, then return None
        } else {
            let row_index_iterator  =   ( 0 .. self.rows() ).rev();
            let matrix              =   self.clone();
            let column_index        =   index.clone();
            Ok( SparseColumn{ matrix, row_index_iterator, column_index } )        
        }
    }

}









//  MATRIX ORACLE OPERATIONS
//  ---------------------------------------------------------------


impl < N, I, IptrStorage, IndStorage, DataStorage, Iptr >

    MatrixOracleOperations for

    Arc<
        CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 
    >

    where
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,
        N: Clone + Debug + Eq,
{}     














      


//  -------------------------------------------------------------------------------------
//  Conversion to sparsemat CSR matrices
//  -------------------------------------------------------------------------------------

/// Sparse matrix data structures that implement this trait have a `.into_csr()` method to convert themselves into CSR matrices.
pub trait IntoCSR< I, N > 
    where
        N: Clone + num::Num,
        I: SpIndex,
{
    /// Converts `self` into a CSR matrix of size `nrows` by `ncols`.
    fn into_csr( self, nrows: usize, ncols: usize ) -> CsMatBase< N, I, Vec<I>, Vec<I>, Vec<N> >;
}


//  Implement for Vec< Vec< N > > 
//  ---------------------------------------- 
impl < N >
    
    IntoCSR
        < usize, N > for
    
    Vec< Vec< N > > 

    where
        N: Clone + num::Num,

{
    /// Converts Vec< Vec< N > > into a CSR matrix   of size  `nrows` by `ncols`.
    /// 
    /// If `v` has type `Vec< Vec< N > >` and `v[i][j]=z` then 
    /// `v.into_csr(nrows, ncols)[i, j]=z`.
    /// 
    /// Panics if `nrows` or `ncols` is not large enough.
    fn into_csr( self, nrows: usize, ncols: usize ) 
        -> 
        CsMatBase<N, usize, Vec<usize>, Vec<usize>, Vec<N>>
    {
        let mut mat = TriMat::new( ( nrows, ncols ) );
        for ( rownum, row ) in self.into_iter().enumerate() {
            for ( colnum,val ) in row.into_iter().enumerate() {
                mat.add_triplet( rownum, colnum, val );
            }
        }
        mat.to_csr()
    }
}



//  Implement for Vec< Vec< (usize, N) > > 
//  ---------------------------------------- 
impl < N >
    
    IntoCSR
        < usize, N > for
    
    Vec< Vec< (usize, N) > > 

    where
        N: Clone + num::Num,

{
    /// Converts Vec< Vec< (usize, N) > > into a CSR matrix  of size  `nrows` by `ncols`.
    /// 
    /// If `v` has type `Vec< Vec< (usize, N) > >` and `v[i]` contains
    /// `(j,z)` then `v.into_csr(nrows, ncols)[i, j]=z`.
    /// 
    /// 
    /// Panics if `nrows` or `ncols` is not large enough.
    fn into_csr( self, nrows: usize, ncols: usize ) 
        -> 
        // CsMatBase<N, usize, Vec<_>, Vec<usize>, Vec<N>, _>
        // CsMatI< N, usize, usize >
        CsMatBase<N, usize, Vec<usize>, Vec<usize>, Vec<N>>
    {
        let mut mat = TriMat::new( ( nrows, ncols ) );
        for ( rownum, row ) in self.into_iter().enumerate() {
            for ( colnum,val ) in row {
                mat.add_triplet( rownum, colnum, val );
            }
        }
        mat.to_csr()
    }
}


//  Implement for Vec< (I, I, N) > 
//  ---------------------------------------- 
impl < I, N >
    
    IntoCSR
        < I,N > for
    
    Vec< (I, I, N) >

    where
        N: Clone + num::Num,
        I: SpIndex,

{
    /// Converts Vec< Vec< (I, I, N) > > into a CSR matrix  of size  `nrows` by `ncols`.
    /// 
    /// If `v` has type `Vec< (I, I, N) >` then for each
    /// `(i,j,z)` in `v` we have `v.into_csr(nrows, ncols)[i, j]=z`.  If `v` contains
    /// multiple triplets with the equal values for `i` and `j` then 
    /// one will be overwritten.
    /// 
    /// 
    /// Panics if `nrows` or `ncols` is not large enough.
    /// 
    /// This is just a wrapper for the method `from_triplets` in the 
    /// `sprs` module.
    fn into_csr( self, nrows: usize, ncols: usize ) 
        -> 
        CsMatBase< N, I, Vec<I>, Vec<I>, Vec<N> >
    {
        let capacity = self.len();
        let mut iv = Vec::with_capacity( capacity );
        let mut jv = Vec::with_capacity( capacity );
        let mut zv = Vec::with_capacity( capacity );                
        for (i,j,z) in self {
            iv.push(i); jv.push(j); zv.push(z);
        }        
        TriMatBase::from_triplets( (nrows, ncols), iv, jv, zv ).to_csr()
    }
}








//  Implement for sorted::VecOfVec< usize, OrderedFloat >
//  ---------------------------------------- 
impl < N >
    
    IntoCSR
        < usize, N > for
    
    crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec< usize, N >

    where
        N: Clone + num::Num,

{
    /// Converts self into a CSR matrix of size  `nrows` by `ncols`.
    /// 
    /// Panics if `nrows` or `ncols` is not large enough.
    fn into_csr( self, nrows: usize, ncols: usize ) 
        -> 
        CsMatBase< N, usize, Vec<usize>, Vec<usize>, Vec<N> >
    {
        let capacity = self.number_of_structural_nonzeros();
        let mut iv = Vec::with_capacity( capacity );
        let mut jv = Vec::with_capacity( capacity );
        let mut zv = Vec::with_capacity( capacity );                
        for (i,j,z) in self.triplets() {
            iv.push(i); jv.push(j); zv.push(z);
        }        
        TriMatBase::from_triplets( (nrows, ncols), iv, jv, zv ).to_csr()
    }
}








#[cfg(test)]
mod tests {
    use crate::algebra::matrices::{debug::{matrix_oracle_is_internally_consistent, matrix_order_operators_are_internally_consistent}, types::vec_of_vec::sorted::VecOfVec};


    #[test]
    fn test_csmatbase_small() {   

        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::matrices::types::third_party::IntoCSR;


        // data for a matrix:
        // |  1   0 |
        // | -1   0 |
        // |  0   0 |
        // |  0   0 |           
        let matrix = vec![ 
                vec![ (0, 1)   ], 
                vec![ (0,-1)   ], 
                vec![          ], 
                vec![          ],                 
            ];
        
        // wrap the data in a VecOfVec sparse matrix struct (see documentation for details)
        let matrix = matrix.into_csr(4,2);

        // check that matrix is internally consistent
        let sorted_row_indices = 0..4;
        let sorted_column_indices = 0..2;
        assert!(
            matrix_oracle_is_internally_consistent( &matrix, sorted_row_indices, sorted_column_indices)
        );                

        // get some sparse columns
        let mut column_0        =   (& matrix ).column(&0);
        let mut column_0_rev    =   (& matrix ).column_reverse(&0);
        let mut column_1        =   (& matrix ).column(&1);
        let mut column_1_rev    =   (& matrix ).column_reverse(&1);

        // check the columns are correct
        assert_eq!( column_0.next(),        Some( (0,  1) )    );
        assert_eq!( column_0.next(),        Some( (1, -1) )    );
        assert_eq!( column_0.next(),        None                );    

        assert_eq!( column_0_rev.next(),    Some( (1, -1) )    );
        assert_eq!( column_0_rev.next(),    Some( (0,  1) )    );
        assert_eq!( column_0_rev.next(),    None                );   

        assert_eq!( column_1.next(),        None        );
        assert_eq!( column_1_rev.next(),    None        );                        
    }  



    #[test]
    fn test_csmatbase_random() {   

        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::matrices::types::{third_party::IntoCSR, vec_of_vec::sorted::VecOfVec};

        let num_rows          =   10;
        let num_columns       =   10;
        let modulus             =   11;

        // generate random matrix
        let matrix  =   VecOfVec::random_mod_p(num_rows, num_columns, modulus).inner_vec_of_vec(); // first as vec-of-vec
        let matrix  =   matrix.into_csr(num_rows, num_columns); // convert to CSR

        // check that matrix is internally consistent
        let sorted_row_indices = 0..num_rows;
        let sorted_column_indices = 0..num_columns;
        assert!(
            matrix_oracle_is_internally_consistent( &matrix, sorted_row_indices, sorted_column_indices)
        );                                     
    }      





}

