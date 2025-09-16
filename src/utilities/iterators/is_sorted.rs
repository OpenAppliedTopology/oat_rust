//! Determines whether an iterator is sorted.
//! 
//! Code for this module is inspired by the crate [`is_sorted`](https://docs.rs/is_sorted/0.1.1/is_sorted/trait.IsSorted.html).
//! See below for the copyright and license.  The OAT developers originally attempted to use the `is_sorted` crate as a dependency but encountered build issues.
//! 
//! 
//! 
//! > Copyright (c) 2018 Gonzalo Brito Gadeschi Copyright (c) 2017 The Rust Project Developers
//! >
//! > Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//! > 
//! > The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//! > 
//! > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

use itertools::Itertools;
use crate::utilities::order::JudgePartialOrder;



pub trait IsSortedBy: Iterator {

    /// Returns `true` if the iterator returns items in ascneding order according to `comparator`.
    /// 
    /// Specifically, returns `true` if for every pair of consecutive elemnts `(a,b)` returned
    /// by the iterator, the `comparator` function satisfies `comparator(a,b) = true`.
    /// 
    /// # Examples 
    /// 
    /// ```
    /// use oat_rust::utilities::iterators::is_sorted::IsSortedBy;
    /// 
    /// let v = vec![ 1, 2, 3];
    /// assert!( v.iter().is_sorted_by(| a, b | a < b) );
    /// 
    /// let u = vec![ 1, 2, 1];
    /// assert!( ! u.iter().is_sorted_by(| a, b | a < b) );
    /// ```
    fn is_sorted_by<F: FnMut(& Self::Item, & Self::Item) -> bool > (&mut self, mut comparator: F) -> bool {
        match self.next() {
            Some( mut item_second_from_last ) => {
                for item_last in self {
                    if ! comparator( & item_second_from_last, & item_last ){ 
                        return false 
                    } else {
                        item_second_from_last = item_last;
                    }
                }
                true
            },
            None => { true }
        }
    }



    /// Checks that the iterator is sorted in strictly ascending order according to the order operator
    /// 
    /// If it is sorted strictly, this function returns `Ok(())`. Otherwise it returns `Err( n, a, b )`
    /// where `(a,b)` is the first pair of consecutive elements such that `order_operator.judge_lt(a,b) = false`,
    /// and where `a` is the `n`th item returned by the iterator.
    fn is_sorted_strictly_by_order_operator< OrderOperator >( &mut self, order_operator: OrderOperator ) 
        -> 
        Result< (), (usize, Self::Item, Self::Item ) >

        where
            OrderOperator:          JudgePartialOrder< Self::Item >
    {
        match self.next() {
            Some( mut item_second_from_last ) => {
                for (counter, item_last) in self.enumerate() {
                    if ! order_operator.judge_lt( & item_second_from_last, & item_last ){ 
                        return Err( (counter, item_second_from_last, item_last )  ) 
                    } else {
                        item_second_from_last = item_last;
                    }
                }
                Ok(())
            },
            None => { return Ok(()) }
        }        
    }


}

impl < I: Iterator > IsSortedBy for I
{}












/// Checks consecutive pairs of elements in an iterator with a provided closure.
///
/// The function iterates over consecutive pairs of elements in the given iterator
/// and applies the provided closure to each pair. If the closure returns `false`
/// for any pair, the function returns `Some((index, item1, item2))` where `index`
/// is the position of the first element of the pair, and `item1` and `item2` are
/// the elements that caused the closure to return `false`. If the closure returns
/// `true` for all pairs, the function returns `None`.
///
/// # Parameters
///
/// - `iter`: An iterator of type `I`.
/// - `f`: A closure of type `F` that takes two references to items from the iterator
///        and returns a `bool`.
///
/// # Returns
///
/// `Option<(usize, I::Item, I::Item)>`:
/// - `Some((index, item1, item2))` if the closure returns `false` for any pair,
///   where `index` is the position of the first element of the pair, and `item1`
///   and `item2` are the elements that did not satisfy the closure.
/// - `None` if the closure returns `true` for all pairs.
///
/// # Examples
///
/// ```
/// use itertools::Itertools;
/// use oat_rust::utilities::iterators::is_sorted::check_pairs;
///
/// let vec = vec![1, 2, 3, 5, 4];
/// let result = check_pairs(vec.iter(), |a, b| a < b);
///
/// match result {
///     Some((index, a, b)) => {
///         println!("Pair at index {}: ({:?}, {:?}) does not satisfy the condition", index, a, b);
///     }
///     None => println!("All pairs satisfy the condition"),
/// }
/// ```
///
/// This example checks if each element in the vector is less than the next element.
/// It will print that the pair `(5, 4)` at index 3 does not satisfy the condition.
pub fn check_pairs<I, F>(iter: I, mut f: F) -> Option<(usize, I::Item, I::Item)>
    where
        I: Iterator,
        F: FnMut(&I::Item, &I::Item) -> bool,
        I::Item: Clone,
{
    for (index, (a, b)) in iter.tuple_windows().enumerate() {
        if !f(&a, &b) {
            return Some((index, a.clone(), b.clone()));
        }
    }
    None
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
    fn test_skip_until() {
        use crate::utilities::iterators::is_sorted::IsSortedBy;
        
        let v = vec![ 1, 2, 3];
        assert!( v.iter().is_sorted_by(| a, b | a < b) );
        
        let u = vec![ 1, 2, 1];
        assert!( ! u.iter().is_sorted_by(| a, b | a < b) );    
    }
}