//! Create a new matrix by applying a transformation to each entry of an existing matrix.

use std::{fmt::Debug, marker::PhantomData};

use derive_getters::{Dissolve, Getters};

use crate::{algebra::matrices::query::{ MatrixOracle}, utilities::{functions::evaluate::EvaluateFunction, iterators::general::MapByTransform}};
use crate::algebra::vectors::entries::{KeyValNew, ReindexEntry};









//  REINDEX SQUARE MATRIX
//  ---------------------------------------------------------------------

/// Reindexes a matrix given a pair of functions that map old indices to new indices, and vice versa.
/// 
/// Concretely, suppose `f` is a function that maps new indices to old ones. If the original matrix `M` satisfies 
/// ```text
/// M[f(i),f(j)] = t
/// ```
/// then the new matrix `N` then satisfies 
/// ```text
/// N[ i, j ] = t
/// ```
/// 
/// # Requirements and warnings
/// 
/// The functions `reindex_function_old_to_new` and `reindex_function_new_to_old` must be muutally inverse.
/// That is, `reindex_function_new_to_old(    reindex_function_old_to_new(  old_index   )    )` must equal
/// `old_index`, for all row/column indices of `M`.  Otherwise the matrix `N` may violate some of the conditions
/// required for a valid implementation of the [MatrixOracle] trait.
/// 
/// # Design notes
/// 
/// This construction relies on user-provided functions to map old indices to new and vice versa. 
/// Currently the functions take full values as inputs, rather than references.
/// References would probably give better performance in general, but they would introduce the
/// need for explcit lifetime parameters and bounds, which gives added complexity. 
/// At the time this documenation was written, the [ReindexSquareMatrix] struct is only used 
/// once in the OAT library: to reindex a square, integer-indexed matrix representing the matched
/// part of the target COMB in a U-match decomposition. Since cloning versus referencing has
/// very little performance difference for integers, we elected to stick with the simpler implementation
/// for now (aknowledging that this is only half the battle; if you want to index into the matrix then
/// you will have to clone the index). If users find strong motivation for a more complex, reference-based implementation, this
/// can be revisited.
#[derive(Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct ReindexSquareMatrix< MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 
{
    matrix_unreindexed:             MatrixUnreindexed,
    reindex_function_old_to_new:    ReindexFunctionOldToNew,
    reindex_function_new_to_old:    ReindexFunctionNewToOld,    
    phantom_indexold:               PhantomData< IndexOld >,
    phantom_indexnew:               PhantomData< IndexNew >,
    phantom_entrynew:               PhantomData< EntryNew >,
}

//  Implement the struct
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 
{
    pub fn new( matrix_unreindexed: MatrixUnreindexed, reindex_function_old_to_new: ReindexFunctionOldToNew, reindex_function_new_to_old: ReindexFunctionNewToOld )
        -> 
        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 
    {
        ReindexSquareMatrix{ matrix_unreindexed, reindex_function_old_to_new, reindex_function_new_to_old, phantom_indexnew: PhantomData, phantom_indexold: PhantomData, phantom_entrynew: PhantomData }
    }                    
}




//  MatrixOracle
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        MatrixOracle for 

        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        where
            MatrixUnreindexed:          MatrixOracle< RowIndex = IndexOld, ColumnIndex = IndexOld >,
            ReindexFunctionOldToNew:    Clone + EvaluateFunction< IndexOld, IndexNew >,
            ReindexFunctionNewToOld:    Clone + EvaluateFunction< IndexNew, MatrixUnreindexed::RowIndex >,
            EntryNew:                   Clone + Debug + PartialEq + KeyValNew <  Key = IndexNew,  Val = MatrixUnreindexed::Coefficient   >,
            IndexNew:                   Clone + Debug + Eq,
{   
    type RowEntry               =   EntryNew;    
    type ColumnEntry            =   EntryNew;        
    type ColumnIndex            =   IndexNew;     
    type RowIndex               =   IndexNew;     
    type Coefficient            =   MatrixUnreindexed::Coefficient;

    type Row                    =   MapByTransform< 
                                                         MatrixUnreindexed::Row, 
                                                         EntryNew,
                                                         ReindexEntry< ReindexFunctionOldToNew >
                                                     >;    
    type RowReverse             =   MapByTransform< 
                                                         MatrixUnreindexed::RowReverse, 
                                                         EntryNew,
                                                         ReindexEntry< ReindexFunctionOldToNew >
                                                     >; 
    type Column                 =   MapByTransform< 
                                                         MatrixUnreindexed::Column, 
                                                         EntryNew,
                                                         ReindexEntry< ReindexFunctionOldToNew >
                                                     >;    
    type ColumnReverse          =   MapByTransform< 
                                                         MatrixUnreindexed::ColumnReverse, 
                                                         EntryNew,
                                                         ReindexEntry< ReindexFunctionOldToNew >
                                                     >;                                                                                                              
    
    fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index.clone() );
        let entry_reindexer = ReindexEntry::new( self.reindex_function_old_to_new.clone() );
        
        MapByTransform::new(
            self.matrix_unreindexed.row( & selecting_index_new ),
            entry_reindexer,
        )
    }
    
    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index.clone() );
        let entry_reindexer = ReindexEntry::new( self.reindex_function_old_to_new.clone() );
        
        MapByTransform::new(
            self.matrix_unreindexed.row_reverse( & selecting_index_new ),
            entry_reindexer,
        )
    }
    
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index.clone() );
        let entry_reindexer = ReindexEntry::new( self.reindex_function_old_to_new.clone() );
        
        MapByTransform::new(
            self.matrix_unreindexed.column( & selecting_index_new ),
            entry_reindexer,
        )
    }
    
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index.clone() );
        let entry_reindexer = ReindexEntry::new( self.reindex_function_old_to_new.clone() );
        
        MapByTransform::new(
            self.matrix_unreindexed.column_reverse( & selecting_index_new ),
            entry_reindexer,
        )
    }
    
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index.clone() );
        return self.matrix_unreindexed.has_row_for_index( & selecting_index_new )
    }
    
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index.clone() );
        return self.matrix_unreindexed.has_column_for_index( & selecting_index_new )
    }
    
    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        let row = (self.reindex_function_new_to_old).evaluate_function( row.clone() );
        let column = (self.reindex_function_new_to_old).evaluate_function( column.clone() );

        return self.matrix_unreindexed.structural_nonzero_entry( & row, & column )
    }    
}










//  REINDEX COLUMNS
//  ---------------------------------------------------------------------

/// Reindexes the *columns* of matrix, given a pair of functions that map old indices to new indices, and vice versa.
/// 
/// Concretely, suppose `f` is a function that maps new indices to old ones. If the original matrix `M` satisfies 
/// ```text
/// M[ i,f(j)] = t
/// ```
/// then the new matrix `N` then satisfies 
/// ```text
/// N[ i, j ] = t
/// ```
/// 
/// The functions `reindex_function_old_to_new` and `reindex_function_new_to_old` must be muutally inverse.
/// That is, `reindex_function_new_to_old(    reindex_function_old_to_new(  old_index   )    )` must equal
/// `old_index`, for all column indices of `M`.  Otherwise the matrix `N` may violate some of the conditions
/// required for a valid implementation of the [MatrixOracle] trait.
/// 
/// # Design notes
/// 
/// This construction relies on user-provided functions to map old indices to new and vice versa. 
/// Currently the functions take full values as inputs, rather than references.
/// References would probably give better performance in general, but they would introduce the
/// need for explcit lifetime parameters and bounds, which gives added complexity. 
/// At the time this documenation was written, the [ReindexMatrixColumns] struct is only used 
/// once in the OAT library. If users find strong motivation for a more complex, reference-based implementation, this
/// can be revisited.
#[derive(Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct ReindexMatrixColumns< MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 
{
    matrix_unreindexed:             MatrixUnreindexed,
    reindex_function_old_to_new:    ReindexFunctionOldToNew,
    reindex_function_new_to_old:    ReindexFunctionNewToOld,    
    phantom_indexold:               PhantomData< IndexOld >,
    phantom_indexnew:               PhantomData< IndexNew >,
    phantom_entrynew:               PhantomData< EntryNew >,
}

//  Implement the struct
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

    ReindexMatrixColumns
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 
{
    pub fn new( matrix_unreindexed: MatrixUnreindexed, reindex_function_old_to_new: ReindexFunctionOldToNew, reindex_function_new_to_old: ReindexFunctionNewToOld )
        -> 
        ReindexMatrixColumns
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 
    {
        ReindexMatrixColumns{ matrix_unreindexed, reindex_function_old_to_new, reindex_function_new_to_old, phantom_indexnew: PhantomData, phantom_indexold: PhantomData, phantom_entrynew: PhantomData }
    }                    
}


impl< MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew >

    Clone for

    ReindexMatrixColumns
        < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew >
    where
        MatrixUnreindexed:          Clone,
        ReindexFunctionOldToNew:    Clone,
        ReindexFunctionNewToOld:    Clone,
    {
        fn clone(&self) -> Self {
            Self{ 
                matrix_unreindexed:             self.matrix_unreindexed.clone(), 
                reindex_function_old_to_new:    self.reindex_function_old_to_new.clone(), 
                reindex_function_new_to_old:    self.reindex_function_new_to_old.clone(), 
                phantom_indexold:               PhantomData, 
                phantom_indexnew:               PhantomData, 
                phantom_entrynew:               PhantomData 
            }
        }
    }



//  MatrixOracle
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        MatrixOracle for 

        ReindexMatrixColumns
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        where
            MatrixUnreindexed:          MatrixOracle< ColumnIndex = IndexOld >,
            ReindexFunctionOldToNew:    Clone + EvaluateFunction< IndexOld, IndexNew >,
            ReindexFunctionNewToOld:    Clone + EvaluateFunction< IndexNew, MatrixUnreindexed::ColumnIndex >,
            EntryNew:                   Clone + Debug + PartialEq + KeyValNew <  Key = IndexNew,  Val = MatrixUnreindexed::Coefficient   >,
            IndexNew:                   Clone + Debug + Eq,
{   
    type RowEntry               =   EntryNew;    
    type ColumnEntry            =   MatrixUnreindexed::ColumnEntry;        
    type ColumnIndex            =   IndexNew;     
    type RowIndex               =   MatrixUnreindexed::RowIndex;     
    type Coefficient            =   MatrixUnreindexed::Coefficient;

    type Row                    =   MapByTransform< 
                                                         MatrixUnreindexed::Row, 
                                                         EntryNew,
                                                         ReindexEntry< ReindexFunctionOldToNew >
                                                     >;    
    type RowReverse             =   MapByTransform< 
                                                         MatrixUnreindexed::RowReverse, 
                                                         EntryNew,
                                                         ReindexEntry< ReindexFunctionOldToNew >
                                                     >; 
    type Column                 =   MatrixUnreindexed::Column;    
    type ColumnReverse          =   MatrixUnreindexed::ColumnReverse;                                                                                                              
    
    fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
        let entry_reindexer = ReindexEntry::new( self.reindex_function_old_to_new.clone() );
        
        MapByTransform::new(
            self.matrix_unreindexed.row( & index ),
            entry_reindexer,
        )
    }
    
    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
        let entry_reindexer = ReindexEntry::new( self.reindex_function_old_to_new.clone() );
        
        MapByTransform::new(
            self.matrix_unreindexed.row_reverse( & index ),
            entry_reindexer,
        )
    }
    
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index.clone() );
        self.matrix_unreindexed.column( & selecting_index_new )
    }
    
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index.clone() );
        self.matrix_unreindexed.column_reverse( & selecting_index_new )
    }
    
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        return self.matrix_unreindexed.has_row_for_index( index )
    }
    
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index.clone() );
        return self.matrix_unreindexed.has_column_for_index( & selecting_index_new )
    }
    
    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        let column = (self.reindex_function_new_to_old).evaluate_function( column.clone() );

        return self.matrix_unreindexed.structural_nonzero_entry( row, & column )
    }    
}

















//  TESTS
//  =========================================================================================================

//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {

}        
