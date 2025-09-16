//! Merge operations specialized for pairs (versus larger groups) of iterators.

use crate::utilities::order::JudgePartialOrder;
use crate::utilities::iterators::general::PeekUnqualified;
use std::cmp::Ordering;


/// Merges two iterators; if both iterators return entries in ascencing order (according to the order comaparator)
/// then so does the merged iterator.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::merge::two_iter_different_type::MergeTwoIteratorsByOrderOperator;        
/// use oat_rust::utilities::order::OrderOperatorAuto;
/// 
/// let vec1 = vec![ 0, 2, 4 ];
/// let vec2 = vec![ 1, 3, 5 ];
/// let merged 
///     = MergeTwoIteratorsByOrderOperator::new( 
///             vec1.iter().peekable(), 
///             vec2.iter().peekable(), 
///             OrderOperatorAuto
///         );
/// itertools::assert_equal( merged, vec![ 0,1,2,3,4,5 ].iter() );
/// ```
pub struct MergeTwoIteratorsByOrderOperator< Iterator1, Iterator2, OrderOperator >
    where 
        Iterator1:          Iterator + PeekUnqualified,
        Iterator2:          Iterator< Item = Iterator1::Item > + PeekUnqualified,
        OrderOperator:    JudgePartialOrder<  Iterator1::Item >
{
    iter1:              Iterator1,
    iter2:              Iterator2,
    order_comparator:   OrderOperator,
}

impl < Iterator1, Iterator2, OrderOperator >

    MergeTwoIteratorsByOrderOperator
        < Iterator1, Iterator2, OrderOperator >

    where 
        Iterator1:          Iterator + PeekUnqualified,
        Iterator2:          Iterator< Item = Iterator1::Item > + PeekUnqualified,
        OrderOperator:    JudgePartialOrder<  Iterator1::Item >        

{
    pub fn new( iter1: Iterator1, iter2: Iterator2, order_comparator: OrderOperator )
            ->
            MergeTwoIteratorsByOrderOperator
                < Iterator1, Iterator2, OrderOperator >                
    {
        MergeTwoIteratorsByOrderOperator{ iter1, iter2, order_comparator }
    }                
}        




impl < Iterator1, Iterator2, OrderOperator >

    Iterator for

    MergeTwoIteratorsByOrderOperator
        < Iterator1, Iterator2, OrderOperator >

    where 
        Iterator1:          Iterator + PeekUnqualified,
        Iterator2:          Iterator< Item = Iterator1::Item > + PeekUnqualified,
        OrderOperator:    JudgePartialOrder<  Iterator1::Item >        
{
    type Item = Iterator1::Item;
    
    fn next( &mut self ) -> Option< Self::Item > {
        match self.iter1.peek_unqualified() {
            Some( item1 ) => {
                match self.iter2.peek_unqualified() {
                    Some( item2 ) => {
                        match {
                            let ref this = self.order_comparator; 
                            this.judge_partial_cmp(item1, item2) == Some(Ordering::Less) } {
                            true    => { return self.iter1.next() },
                            false   => { return self.iter2.next() },
                        }
                    },
                    None => { return self.iter1.next() }
                }
            },
            None => { return self.iter2.next() }            
        }
    }
}      









//  =========================================================================================================
//  TESTS
//  =========================================================================================================


//  ---------------------------------------------------------------------
//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    


    #[test]
    fn test_merge_by_order_comparator() {
        use crate::utilities::iterators::merge::two_iter_different_type::MergeTwoIteratorsByOrderOperator;  
        use crate::utilities::order::OrderOperatorAuto;      

        let vec1 = vec![ 0, 2, 4 ];
        let vec2 = vec![ 1, 3, 5 ];
        let merged 
            = MergeTwoIteratorsByOrderOperator::new( 
                    vec1.iter().peekable(), 
                    vec2.iter().peekable(), 
                    OrderOperatorAuto
                );
        itertools::assert_equal( merged, vec![ 0,1,2,3,4,5 ].iter() );
    }    
}