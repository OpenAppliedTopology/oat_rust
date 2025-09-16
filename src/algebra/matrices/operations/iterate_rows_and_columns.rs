//! Iterators representing sequences of rows or columns
//! 
//! This module provides some objects that are essentially wrappers for more complex objects that represent
//! sequences of rows and columns; their function is to help
//! keep code short and concise.


use crate::algebra::matrices::query::MatrixOracle;
use derive_new::new;




/// Contains a matrix and an iterable; returns one row of the matrix for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new,PartialEq,Eq,Hash)]
pub struct SequenceOfRows
                < T, I > 
{
    matrix:                 T,
    row_index_iterator:     I,
}        

impl < T, I >

    Iterator for 
    
    SequenceOfRows
        < T, I > 

    where
        T:      MatrixOracle,
        I:      Iterator< Item = T::RowIndex >  
{
    type Item = T::Row;
    fn next(&mut self) -> Option<Self::Item> {
        self.row_index_iterator.next().map(|x| self.matrix.row(&x) )
    }
}   



/// Contains a matrix and an iterable; returns one row of the matrix (with order of entries reversed) for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new,PartialEq,Eq,Hash)]
pub struct SequenceOfReverseRows
                < T, I >      
{
    matrix:     T,
    indices:    I,
}        

impl < T, I >

    Iterator for 
    
    SequenceOfReverseRows
        < T, I > 

    where
        T:      MatrixOracle,
        I:      Iterator< Item = T::RowIndex >  
{
    type Item = T::RowReverse;
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|x| self.matrix.row_reverse(&x) )
    }
}   


/// Contains a matrix and an iterable; returns one column of the matrix for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new,PartialEq,Eq,Hash)]
pub struct SequenceOfColumns
                < T, I > 
{
    matrix:     T,
    indices:    I,
}        

impl < T, I >

    Iterator for 
    
    SequenceOfColumns
        < T, I > 

    where
        T:      MatrixOracle,
        I:      Iterator< Item = T::ColumnIndex >  
{
    type Item = T::Column;
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|x| self.matrix.column(&x) )
    }
}   


/// Contains a matrix and an iterable; returns one column of the matrix (with order of entries reversed) for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new,PartialEq,Eq,Hash)]
pub struct SequenceOfReverseColumns
                < T, I >      
{
    matrix:     T,
    indices:    I,
}        

impl < T, I >

    Iterator for 
    
    SequenceOfReverseColumns
        < T, I > 

    where
        T:      MatrixOracle,
        I:      Iterator< Item = T::ColumnIndex >  
{
    type Item = T::ColumnReverse;
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|x| self.matrix.column_reverse(&x) )
    }
}   


