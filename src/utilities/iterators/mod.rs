//! General tools for working with iterators.

use crate::algebra::vectors::entries::KeyValGet;
use std::fmt::Debug;

pub mod merge;
pub mod general;
pub mod is_sorted;



pub trait OatIteratorMethods: IntoIterator 
{
    /// Convert an interable of iterables into a Vec-of-Vec matrix
    /// 
    /// The input should have the format of a list of lists: `[ row_1, row_2, .. ]`,
    /// where each row is a sequence of `(column_index, coefficient)` pairs.
    /// The list of column indices in each row should appear in strictly ascending
    /// order.
    fn into_sorted_vec_of_vec< ColumnIndex, Coefficient >( self ) -> 
        Result<
            crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec< ColumnIndex, Coefficient >,
            Vec<Vec<(ColumnIndex, Coefficient)>>,
        >
        where
            Self:                                       Sized,
            Self::Item:                                 IntoIterator,
            < Self::Item as IntoIterator >::Item:       KeyValGet< Key = ColumnIndex, Val = Coefficient  >, 
            ColumnIndex:                                Clone + Debug + Ord,
            Coefficient:                                Clone + Debug,
    {
        let vecvec =    self.into_iter()
                            .map( 
                                |x| 
                                x.into_iter()
                                    .map( |index_coefficient_pair| (index_coefficient_pair.key(), index_coefficient_pair.val()) )
                                    .collect::<Vec<_>>() 
                            )
                            .collect();
        crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::new( vecvec )
    }

}