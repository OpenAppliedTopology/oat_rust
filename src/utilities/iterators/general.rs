
//!  General utilities for work with iterators.

use std::{cmp::Ordering, iter::{Iterator, Peekable}, marker::PhantomData, sync::Arc};


use derive_getters::{Dissolve, Getters};
use derive_new::new;
use itertools::{Merge, Itertools};

use crate::{utilities::{functions::evaluate::{EvaluateFunction, EvaluateFunctionRef, LogicalNot}, order::{JudgeOrder, OrderOperatorAuto, OrderOperatorByKey, JudgePartialOrder}, sets::MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper}, algebra::vectors::entries::KeyValGet};



//  ---------------------------------------------------------------------------
//  TYPE ALIASES


/// Type alias for an iterator that excludes elements of a container
pub type FilterOutMembers< J, Container > = 
        Filter< 
                J, 
                LogicalNot< 
                        MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper<
                                Container,
                            >  
                    > 
            >;  

/// Type alias for an iterator that excludes elements that don't belong to a container
pub type FilterOutNonmembers< J, Container > =
        Filter< 
                J, 
                MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper<
                        Container,
                    >  
            >;  


//  ---------------------------------------------------------------------------
//  ELEMENTWISE EQUALITY OF ITERATORS WRAPPED IN RESULTS


/// Compares two `Result`-wrapped iterators for element-wise equality.
/// 
/// Returns `true` if 
/// - both are `Ok` and their iterators are element-wise equal, or
/// - both are `Err` and the errors are equal
/// Returns `false` otherwise.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::utilities::iterators::general::result_iterators_are_elementwise_equal;
/// 
/// let iter1: Result<std::vec::IntoIter<i32>, &str> = Ok(vec![1, 2, 3].into_iter());
/// let iter2: Result<std::vec::IntoIter<i32>, &str> = Ok(vec![1, 2].into_iter());
/// let iter3: Result<std::vec::IntoIter<i32>, &str> = Err("Error");
/// 
/// // Both are Ok  and equal
/// assert!(    result_iterators_are_elementwise_equal(
///     iter1.clone(), 
///     iter1.clone()
/// )); 
/// 
/// // Both are Err and equal   
/// assert!(    result_iterators_are_elementwise_equal(
///     iter3.clone(), 
///     iter3.clone()
/// ));         
/// 
/// // Both are Ok but not equal
/// assert!(  ! result_iterators_are_elementwise_equal(
///     iter1.clone(), 
///     iter2.clone()
/// )); 
/// 
/// // One is Ok, the other is Err  
/// assert!(  ! result_iterators_are_elementwise_equal(
///     iter1.clone(), 
///     iter3.clone()
/// ));        
/// 
/// // One is Err, the other is Ok               
/// assert!(  ! result_iterators_are_elementwise_equal(
///     iter3.clone(),
///     iter1.clone(), 
/// ));  
/// ``` 
pub fn result_iterators_are_elementwise_equal<I1, I2, T>(result1: Result<I1, T>, result2: Result<I2, T>) -> bool
where
    I1: Iterator,
    I2: Iterator,
    I1::Item: PartialEq<I2::Item>,
    T: PartialEq,
{
    match (result1, result2) {
        (Ok(iter1), Ok(iter2)) => iter1.eq(iter2),
        (Err(err1), Err(err2)) => err1 == err2,
        _ => false,
    }
}            





//  ---------------------------------------------------------------------------
//  TRANSFORMATION TRAIT

pub trait TransformIter{

    /// Wraps `self` in a struct that requires items to appear in strictly sorted order.
    /// 
    /// See [RequireStrictAscent]
    fn require_strict_ascent< OrderOperator > ( self, order_operator: OrderOperator ) 
        ->
        RequireStrictAscent< Self, OrderOperator >

    where
        Self:           Sized + Iterator,
        Self::Item:     Clone,
        OrderOperator:  JudgePartialOrder< Self::Item > 
    {
        RequireStrictAscent { iter: self, last_item: None, order_operator }
    }

    /// Wraps `self` in a struct that requires items to appear in strictly sorted order.
    /// 
    /// See [RequireStrictAscentWithPanic]
    fn require_strict_ascent_with_panic< OrderOperator > ( self, order_operator: OrderOperator ) 
        ->
        RequireStrictAscentWithPanic< Self, OrderOperator >

    where
        Self:           Sized + Iterator,
        Self::Item:     Clone,
        OrderOperator:  JudgePartialOrder< Self::Item > 
    {
        RequireStrictAscentWithPanic { iter: self, last_item: None, order_operator }
    }    
        
}


//  Blanket implementation
impl < T > TransformIter for T {}


//  ---------------------------------------------------------------------------
//  PREPEND ONE ELEMENT


/// Iterator coposed of a head (of type `Option<Item>`), and a tail (an iterator)
/// 
/// This struct is similar to several others
/// - `std::iter::Peekable`
/// - `itertools::PutBack`
/// - [HeadTailHit](crate::utilities::iterators::merge::hit::HeadTailHit)
/// 
/// however it combines functionality from all of them: you can both peek put values back.
/// 
/// Calling `self.next()`
/// 
/// - first: tries to take the value from head, 
/// - then (if head is empty): tries to take a value from tail
/// 
/// Thus, if the head is empty, then it will remain empty **until one of two events occur**:
/// 
/// - the user places a value in it
/// - the user calls `.peek_unqualified()`, in which case a value will be stored in head, and a reference to that value returned
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::{HeadTail, PeekUnqualified};
/// 
/// let mut i = HeadTail{ head: Some(0), tail: 1..2 };
/// let mut j = HeadTail{ head: None,    tail: 1..2 };
/// let mut k = HeadTail{ head: None,    tail: 1..2 };
/// let mut l = HeadTail{ head: None,    tail: std::iter::empty::<i32>() };
/// 
/// assert_eq!( i.peek_unqualified(), Some(&0) );
/// assert_eq!( j.peek_unqualified(), Some(&1) );
/// 
/// assert_eq!( i.head.as_ref(), Some(&0) );
/// assert_eq!( j.head.as_ref(), Some(&1) );
/// assert_eq!( k.head.as_ref(), None     ); // no value placed in head because we haven't peeked
/// 
/// println!("{:#?}", &k );
/// 
/// assert_eq!( i.next(), Some(0) );
/// assert_eq!( j.next(), Some(1) );
/// assert_eq!( k.next(), Some(1) );
/// assert_eq!( l.next(), None    );
/// 
/// assert_eq!( i.head.as_ref(), None );
/// assert_eq!( j.head.as_ref(), None );
/// assert_eq!( k.head.as_ref(), None );
/// assert_eq!( l.head.as_ref(), None );
/// ```
#[derive(new,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct HeadTail< I: Iterator >{
    pub head:   Option< I::Item >,
    pub tail:   I,    
}



impl < I: Iterator >

    Iterator for

    HeadTail< I > 
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        self.head.take().or_else(|| self.tail.next() )
    }
}


impl < I: Iterator >

    PeekUnqualified for

    HeadTail< I > 
{
    /// Peek at the next value
    /// 
    /// If the head is empty, try filling it with a value from tail, then return a reference.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::utilities::iterators::general::{HeadTail, PeekUnqualified};
    /// 
    /// let mut i = HeadTail{ head: Some(0), tail: 1..2 };
    /// let mut j = HeadTail{ head: None,    tail: 1..2 };
    /// let mut k = HeadTail{ head: None,    tail: 1..0 };
    /// 
    /// assert_eq!( i.peek_unqualified(), Some(&0) );
    /// assert_eq!( j.peek_unqualified(), Some(&1) );
    /// assert_eq!( k.peek_unqualified(), None     );
    /// 
    /// assert_eq!( i.head.as_ref(), Some(&0) );
    /// assert_eq!( j.head.as_ref(), Some(&1) );
    /// assert_eq!( k.head.as_ref(), None     );
    /// ```
    fn peek_unqualified( &mut self ) -> Option< & Self::Item > {
        self.head
            .is_none()
            .then(|| self.tail.next().map(|x| self.head.replace(x) ) );
        self.head.as_ref()

        // match self.head.is_some() {
        //     true => { return self.head.as_ref() }
        //     false => {
        //         match self.tail.next() {
        //             Some( val ) => { self.head.replace(val);  }
        //         }
        //         self.head.replace( self.tail.next() );

        //     }
        // }
        
        // if self.head.is_none() {
        //     let val = self.tail.next();

        // }
    }
}

impl < I >

    Clone for

    HeadTail< I > 

    where
        I:          Clone + Iterator,
        I::Item:    Clone,
{
    fn clone(&self) -> Self {
        HeadTail { head: self.head.clone(), tail: self.tail.clone() }
    }
}



//  ---------------------------------------------------------------------------
//  REQUIRE STRICT ASCENT


/// Throws an error if items are not sorted in ascending order
/// 
/// Returns `Ok(y)` if `xy` is not preceded by an element `x` strictly greater than `y`.
/// Otherwise returns `Err(x,y)`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::RequireStrictAscent;
/// use oat_rust::utilities::order::OrderOperatorAuto;
/// 
/// let iter            =   vec![ 0, 1, 1, 0, ];
/// let mut required    =   RequireStrictAscent::new(
///                             iter.into_iter(),
///                             OrderOperatorAuto,
///                         );
/// 
/// assert_eq!( required.next(), Some(Ok(0))     );
/// assert_eq!( required.next(), Some(Ok(1))     );
/// assert_eq!( required.next(), Some(Err((1,1))));
/// assert_eq!( required.next(), Some(Err((1,0))));
/// assert_eq!( required.next(), None            );
/// ```
#[derive(Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct RequireStrictAscent
            < Iter, OrderOperator >
    where
        Iter:               Iterator,
        Iter::Item:         Clone,
        OrderOperator:      JudgePartialOrder< Iter::Item >,
{
    pub iter:               Iter,
    pub last_item:          Option< Iter::Item >,
    pub order_operator:     OrderOperator,
}



impl < Iter, OrderOperator >

    RequireStrictAscent
        < Iter, OrderOperator > where

    Iter:               Iterator,
    Iter::Item:         Clone,    
    OrderOperator:      JudgePartialOrder< Iter::Item >,    
{
    /// Create a new instance from an iterator
    pub fn new( iter: Iter, order_operator: OrderOperator ) -> Self {
        RequireStrictAscent { iter, last_item: None, order_operator }
    }

    /// Returns `Result< Vec<_>, (itema, itemb) >`; if an error is returned, then `(itema, itemb)` are consecutive elements such that `itema >= itemb`.
    pub fn into_vec( self ) 
        -> 
        Result< 
                Vec< Iter::Item >, 
                (Iter::Item, Iter::Item),
            > 
    {
        return self.collect();
    } 
}





impl < Iter, OrderOperator >

    Iterator for 

    RequireStrictAscent
        < Iter, OrderOperator > where

    Iter:               Iterator,
    Iter::Item:         Clone,    
    OrderOperator:      JudgePartialOrder< Iter::Item >,    
{
    type Item = Result< Iter::Item, (Iter::Item, Iter::Item) >;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some( val_y ) = self.iter.next() { // if the iterator has another item
            let x_opt = self.last_item.replace( val_y.clone() ); // store a copy in the "last_item" slot, and pull out the prior value
            if let Some( val_x ) = x_opt { // if there is a prior value
                if self.order_operator.judge_ge( &val_x, &val_y ) { // and it is strictly greater than y
                    return Some( Err( (val_x, val_y) ) ) // return an error
                }
            }
            return Some( Ok(val_y) ) // otherwise return y
        }
        None // if the iterator doesn't have another item, then return None
    }
}


//  ---------------------------------------------------------------------------
//  REQUIRE STRICT ASCENT -- OR PANIC


/// Throws an error if items are not sorted in ascending order
/// 
/// Returns `Ok(y)` if `y` is not preceded by an element `x` strictly greater than `y`.
/// Otherwise returns `Err(x,y)`.
/// 
/// # Examples
/// 
/// Panics:
/// 
/// ```
/// use oat_rust::utilities::iterators::general::RequireStrictAscentWithPanic;
/// use oat_rust::utilities::order::OrderOperatorAuto;
/// use std::panic;
/// 
/// 
/// #[should_panic(expected = "\n\n| ERROR: An iterator placed inside a `RequireStrictAscentWithPanic` struct has returned two consecutive entries, (x,y) where x > y.\n| NB: This message can also appear when using a reversed order operator, indicating a failure to strictly *descend*.\n| This error message is generated by OAT.\n\n")]
/// fn test() {
///     let iter            =   vec![ 0, 1, 1, 0, ];
///     let mut required    =   RequireStrictAscentWithPanic::new(
///                                 iter.into_iter(),
///                                 OrderOperatorAuto,
///                             );
//      for item in required {}
/// }
///
/// test();
/// ```
/// 
/// Doesn't panic:
/// 
/// /// ```
/// use oat_rust::utilities::iterators::general::RequireStrictAscentWithPanic;
/// use oat_rust::utilities::order::OrderOperatorAuto;
/// 
/// let iter            =   vec![ 0, 1, 2, 3, ];
/// let mut required    =   RequireStrictAscentWithPanic::new(
///                             iter.into_iter(),
///                             OrderOperatorAuto,
///                         );
///  for item in required {}
/// ```
#[derive(Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct RequireStrictAscentWithPanic
            < Iter, OrderOperator >
    where
        Iter:               Iterator,
        Iter::Item:         Clone,
        OrderOperator:      JudgePartialOrder< Iter::Item >,
{
    pub iter:               Iter,
    pub last_item:          Option< Iter::Item >,
    pub order_operator:     OrderOperator,
}



impl < Iter, OrderOperator >

    RequireStrictAscentWithPanic
        < Iter, OrderOperator > where

    Iter:               Iterator,
    Iter::Item:         Clone,    
    OrderOperator:      JudgePartialOrder< Iter::Item >,    
{
    pub fn new( iter: Iter, order_operator: OrderOperator ) -> Self {
        RequireStrictAscentWithPanic { iter, last_item: None, order_operator }
    }
}





impl < Iter, OrderOperator >

    Iterator for 

    RequireStrictAscentWithPanic
        < Iter, OrderOperator > where

    Iter:               Iterator,
    Iter::Item:         Clone,    
    OrderOperator:      JudgePartialOrder< Iter::Item >,    
{
    type Item = Iter::Item;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some( val_y ) = self.iter.next() { // if the iterator has another item
            let x_opt = self.last_item.replace( val_y.clone() ); // store a copy in the "last_item" slot, and pull out the prior value
            if let Some( val_x ) = x_opt { // if there is a prior value
                if self.order_operator.judge_ge( &val_x, &val_y ) { // and it is strictly greater than y
                    panic!("\n\n| ERROR: An iterator placed inside a `RequireStrictAscentWithPanic` struct has returned two consecutive entries, (x,y) where x >= y.\n| NB: This message can also appear when using a reversed order operator, indicating a failure to strictly *descend*.\n| This error message is generated by OAT.\n\n");
                }
            }
            return Some( val_y ) // otherwise return y
        }
        None // if the iterator doesn't have another item, then return None
    }
}


impl < Iter, OrderOperator >

    PeekUnqualified for 

    RequireStrictAscentWithPanic
        < Iter, OrderOperator > where

    Iter:               Iterator + PeekUnqualified,
    Iter::Item:         Clone,    
    OrderOperator:      JudgePartialOrder< Iter::Item >,    
{
    fn peek_unqualified( &mut self ) -> Option< & Self::Item > {
        self.iter.peek_unqualified()
    }
}
//  ---------------------------------------------------------------------------
//  GET OR FIND MINIMA / MAXIMA

/// Returns the index of the minimum value of the iterator.
/// 
/// If several elements are equally minimum, the index of the first element is returned. If the iterator is empty, None is returned.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::find_min;
/// 
/// let a               =   vec![ 1, 3, 0, 2, 0 ];
/// let b: Vec<usize>   =   vec![ ];
/// 
/// assert_eq!( find_min(a), Some(2) );
/// assert_eq!( find_min(b), None    );
/// ```
pub fn find_min< I >( iter: I ) -> Option< usize > 
    where 
        I:          IntoIterator,
        I::Item:    Ord,
{
    println!("replace this with the new position_min tool from itertools");
    iter.into_iter()
        .enumerate()
        .min_by(|x,y| x.1.cmp(&y.1) )
        .map(|x| x.0 )
}

/// The minimum of the maxima of a collection of iterables.
/// 
/// Returns `None` if the collection doesn't contain a nonempty iterable.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::minmax;
/// 
/// assert_eq!( Some(4), minmax(vec![ vec![2,4], vec![3,6] ])  );
/// ```
pub fn minmax< T, Inner, Outer >( outer_iter: Outer ) 
        -> 
        Option< T > 
    where
        T:          Ord,
        Inner:      IntoIterator< Item = T     >,
        Outer:      IntoIterator< Item = Inner >
{
    outer_iter.into_iter().filter_map( |x| x.into_iter().max() ).min()
}



//  ---------------------------------------------------------------------------
//  PEEKING


/// Similar to itertools::PeakingNext, but without the `F` parameter
/// (we hypothesize that omitting this closure will help with type
/// inference)
pub trait PeekUnqualified : Iterator {
    fn peek_unqualified( &mut self ) -> Option < & Self::Item >;
}

impl < I : Iterator > PeekUnqualified for Peekable< I >
{
    fn peek_unqualified( &mut self ) -> Option < &<Self as Iterator>::Item > { self.peek() }
}


//  ---------------------------------------------------------------------------
//  SKIP UNTIL


/// Skips each item returned by the iterator, until some item `i` satisfies 
/// `predicate(&i)=true`.  The resulting iterator returns `i` and every item 
/// thereafter.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::SkipUntil;
/// 
/// let vec = vec![ 0, 1, 2, 1, 0];
/// let iter = vec.iter().cloned().peekable().skip_until( |x| x > &1);
/// itertools::assert_equal( iter, vec![ 2, 1, 0 ])  
/// ```
pub trait SkipUntil : Iterator {
    fn skip_until< P >( self, predicate: P ) -> Self
    where
        P: FnMut(&Self::Item) -> bool;
}

impl < I : Iterator + PeekUnqualified > SkipUntil for I 
{
    fn skip_until< P >( mut self, mut predicate: P ) -> Self 
    where
        P: FnMut(&Self::Item) -> bool    
    
    {  
        while let Some( peek ) = self.peek_unqualified() {   if predicate( peek ) { break } else { let _ = self.next(); } }
        self
    }
}






//  ---------------------------------------------------------------------------
//  TRANSFORM ENTRY-WISE


/// Similar to `std::iter::Map` but instead of taking a transformation that implements `FnMut`, 
/// this object takes a transformation that implements `EvaluateFunction`.
/// 
/// See also [`FilterMapByTransform`].
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::MapByTransform;
/// use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
/// 
/// let entrywise_transform_unwrapped = |x : usize| -> i32 { x as i32 };
/// let entrywise_transform_wrapped = EvaluateFunctionFnMutWrapper::new(entrywise_transform_unwrapped);
/// let vec = vec![ 1, 2, 3];
/// let iter = vec.iter().cloned();
/// 
/// let transformed_entrywise = MapByTransform::new( iter, entrywise_transform_wrapped );
/// itertools::assert_equal( transformed_entrywise, vec![ 1i32, 2, 3] );  
/// ```
#[derive(Copy,Debug,Dissolve,Eq,PartialEq,Ord,PartialOrd)]
pub struct MapByTransform
                < Iter, ItemNew, ItemwiseTransform >
    where
        Iter:               Iterator,
        ItemwiseTransform:  EvaluateFunction< Iter::Item, ItemNew >
{
    iter:                   Iter,
    itemwise_transform:     ItemwiseTransform,
    phantom_itemnew:        PhantomData< ItemNew >,
}                

impl < Iter, ItemNew, ItemwiseTransform >

    MapByTransform
        < Iter, ItemNew, ItemwiseTransform >

    where
        Iter:               Iterator,
        ItemwiseTransform:  EvaluateFunction< Iter::Item, ItemNew > 
{
    pub fn new( iter: Iter, itemwise_transform: ItemwiseTransform ) -> MapByTransform < Iter, ItemNew, ItemwiseTransform > {  
        MapByTransform{ iter, itemwise_transform, phantom_itemnew: PhantomData }
    }
}          


impl < Iter, ItemNew, ItemwiseTransform >

    Iterator for

    MapByTransform
        < Iter, ItemNew, ItemwiseTransform >

    where
        Iter:               Iterator,
        ItemwiseTransform:  EvaluateFunction< Iter::Item, ItemNew >        
{
    type Item = ItemNew;
    fn next( &mut self ) -> Option< Self::Item > {
        match self.iter.next() {
            Some( item_old ) => { Some( self.itemwise_transform.evaluate_function( item_old ) ) },
            None => { None }
        }
    }
}


impl < Iter, ItemNew, ItemwiseTransform >

    Clone for

    MapByTransform
        < Iter, ItemNew, ItemwiseTransform >

    where
        Iter:               Clone + Iterator,
        ItemwiseTransform:  Clone + EvaluateFunction< Iter::Item, ItemNew >        
{
    fn clone( & self ) -> Self { MapByTransform::new( self.iter.clone(), self.itemwise_transform.clone() )}
}


pub trait MapByTransformTrait< ItemwiseTransform > {
    fn map_by_transform< ItemNew >( self, itemwise_transform: ItemwiseTransform ) 
        -> 
        MapByTransform
            < Self, ItemNew, ItemwiseTransform >
    where
        Self:               Iterator + Sized,
        ItemwiseTransform:  EvaluateFunction< Self::Item, ItemNew >
    {
        MapByTransform { iter: self, itemwise_transform, phantom_itemnew: PhantomData }
    }
}

impl < ItemwiseTransform, Iter >

    MapByTransformTrait
        < ItemwiseTransform > for
    
    Iter

    where 
        Iter:   Iterator
{}        


//  ---------------------------------------------------------------------------
//  FILTER MAP -- FILTERTING WITH OBJECTS THAT IMPLMENT EvaluateFunction INSTEAD OF FnMut




/// Similar to `std::iter::FilterMap`, but takes an object that satisfies `FilterMapTransform< Input, Output >` rather than an object
/// that implements `FnMut( Input ) -> Output`.
/// 
/// The `std::iter::filter_map` method for iterators is highly versatile, however it requires an
/// argument that implements `FnMut`.  It is often hard to provide such an argument without
/// using closure operators which have `opaque type`.  As an alternative, the OAT library
/// offers this struct, which takes any object which imlements `FilterMapTransform`.
/// 
/// See also [`MapByTransform`].
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::{FilterMapByTransform};
/// use oat_rust::utilities::functions::evaluate::EvaluateFunction;
/// 
/// // define a struct to execute the filter-map function
/// struct ShiftDown;
/// 
/// // implement the `EvaluateFunction` trait
/// impl EvaluateFunction< isize, Option< isize > > for ShiftDown{
///     fn evaluate_function( &self, input: isize ) -> Option< isize > {
///         match input > 0 {
///             true    => { return None },
///             false   => { return Some( input * input * input ) }
///         }
///     }
/// }
/// 
/// // define an iterator
/// let iter = [-2,-1,0,1,2].iter().cloned();
/// 
/// // create a `filter-mapped` object
/// let iter_transformed = FilterMapByTransform::new( iter, ShiftDown );
/// 
/// // verify this object is correct        
/// itertools::assert_equal( iter_transformed, vec![ -8, -1, 0 ]) 
/// ```
#[derive(Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct FilterMapByTransform< I, FilteringObject, ItemNew > 
    where
        I:                  Iterator,
        FilteringObject:    EvaluateFunction< I::Item, Option< ItemNew > >,
{
    iter:               I,
    filtering_object:   FilteringObject,
    phantom_itemnew:    PhantomData< ItemNew >
}

impl < I, FilteringObject, ItemNew >  

    FilterMapByTransform
        < I, FilteringObject, ItemNew > 

    where
        I:                  Iterator,
        FilteringObject:    EvaluateFunction< I::Item, Option< ItemNew > >,
{
    pub fn new< Iterable: IntoIterator< IntoIter = I > >( iterable: Iterable, filtering_object: FilteringObject ) -> FilterMapByTransform< I, FilteringObject, ItemNew > 
        { FilterMapByTransform{ iter: iterable.into_iter(), filtering_object, phantom_itemnew: PhantomData } }
}        


impl < I, FilteringObject, ItemNew > 

    Iterator for 
    
    FilterMapByTransform< I, FilteringObject, ItemNew > 

    where
        I:                  Iterator,
        FilteringObject:    EvaluateFunction< I::Item, Option< ItemNew > >

{
    type Item = ItemNew;

    fn next( &mut self ) -> Option< Self::Item >  {
        for next_item in self.iter.by_ref() {
            match  self.filtering_object.evaluate_function( next_item ) {
                None    => continue,
                Some( new_item ) => { return Some( new_item ) }
            }
        }
        None
    }
}

impl < I, FilteringObject, ItemNew > 

    Clone for 
    
    FilterMapByTransform< I, FilteringObject, ItemNew > 

    where
        I:                  Clone + Iterator,
        FilteringObject:    Clone + EvaluateFunction< I::Item, Option< ItemNew > >    
{
    fn clone( & self ) -> Self { FilterMapByTransform::new( self.iter.clone(), self.filtering_object.clone() )}
}


// Converts 
pub trait FilterMapByTransformTrait< ItemNew >
    where
        Self:               Iterator + Sized,
        // FilteringObject:    FilterMapTransform< Self::Item, ItemNew >
{
    /// Wraps an iterator in a [`FilterMapByTransform`] struct.  
    /// 
    /// **Note** this method is convenient for chaining transformations of an iterator, however it requires the
    /// iterator to implement `Sized`.  To create a `FilterMapTransform` when the iterator does not implement
    /// `Sized`, one can use [`FilterMapByTransform::new`]
    fn filter_map_by_transform< 
                FilteringObject: EvaluateFunction< Self::Item, Option< ItemNew > > 
            >
        ( self, filtering_object: FilteringObject ) 
        -> 
        FilterMapByTransform< Self, FilteringObject, ItemNew > ;
}

impl < Iter, ItemNew > 

    FilterMapByTransformTrait
        < ItemNew > for
    
    Iter

    where
        Iter:       Iterator
{
    fn filter_map_by_transform< FilteringObject: EvaluateFunction< Self::Item, Option< ItemNew > > >
        ( self, filtering_object: FilteringObject ) -> FilterMapByTransform< Self, FilteringObject, ItemNew > {
        FilterMapByTransform::new(self, filtering_object)
    }
}




//  ---------------------------------------------------------------------------
//  FILTER -- FILTERTING WITH OBJECTS THAT IMPLMENT EvaluateFunction INSTEAD OF FnMut




/// Returns only elements where the internal predicate returns `true`.
/// 
/// See also [`FilterMapByTransform`].
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::{Filter, FilterTrait};
/// use oat_rust::utilities::functions::evaluate::EvaluateFunctionRef;
/// 
/// // define a struct to execute the filter-map function
/// struct Positive;
/// 
/// // implement the `EvaluateFunction` trait
/// impl EvaluateFunctionRef< isize, bool > for Positive {
///     fn evaluate_function_ref( &self, input: &isize ) -> bool { input > &0 }
/// }
/// 
/// // define an iterator
/// let iter = [-2,-1,0,1,2];
/// 
/// // create a `filtered` object
/// let iter_transformed = iter.filter( Positive );
/// 
/// // verify this object is correct        
/// itertools::assert_equal( iter_transformed, vec![ 1, 2 ]) 
/// ```
#[derive(Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct Filter< I, FilteringObject > 
    where
        I:                  Iterator,
        FilteringObject:    EvaluateFunctionRef< I::Item, bool >,
{
    iter:               I,
    filtering_object:   FilteringObject,
}

impl < 'a, I, FilteringObject >  

    Filter
        < I, FilteringObject > 

    where
        I:                  Iterator,        
        FilteringObject:    EvaluateFunctionRef< I::Item, bool >,     
{
    pub fn new< J: IntoIterator< IntoIter = I > >( iter: J, filtering_object: FilteringObject ) -> Filter< I, FilteringObject > 
        { Filter{ iter: iter.into_iter(), filtering_object, } }
}        


impl < I, FilteringObject > 

    Iterator for 
    
    Filter< I, FilteringObject, > 

    where
        I:                  Iterator,     
        FilteringObject:    EvaluateFunctionRef< I::Item, bool >

{
    type Item = I::Item;

    fn next( &mut self ) -> Option< Self::Item >  {
        for next_item in self.iter.by_ref() {
            match  self.filtering_object.evaluate_function_ref( &next_item ) {
                false    => continue,
                true => { return Some( next_item ) }
            }
        }
        None
    }
}

impl < I, FilteringObject > 

    Clone for 
    
    Filter< I, FilteringObject > 

    where
        I:                  Clone + Iterator,
        FilteringObject:    Clone + EvaluateFunctionRef< I::Item, bool >    
{
    fn clone( & self ) -> Self { Filter::new( self.iter.clone(), self.filtering_object.clone() )}
}


// Converts 
pub trait FilterTrait
    where
        Self:               IntoIterator + Sized,
        // FilteringObject:    FilterMapTransform< Self::Item, ItemNew >
{
    /// Wraps an iterator in a [`Filter`] struct.  
    /// 
    /// **Note** this method is convenient for chaining transformations of an iterator, however it requires the
    /// iterator to implement `Sized`.  To create a `Filter` when the iterator does not implement
    /// `Sized`, one can use [`Filter::new`]
    fn filter< 
                FilteringObject: EvaluateFunctionRef< Self::Item, bool > 
            >
        ( self, filtering_object: FilteringObject ) 
        -> 
        Filter< Self::IntoIter, FilteringObject >
    {
        Filter::new(self.into_iter(), filtering_object)
    }            
}

impl < Iter > 

    FilterTrait
        for
    
    Iter

    where
        Iter:           IntoIterator,
{}





// 2-TYPE INTO-ITERATOR
// ---------------------------------------------------------------------------

/// An enum that combines two different into-iterator types into a single into-iterator type.
///
/// There are many operations available in Rust for combining iterators (chaining, merging, etc.).
/// Quite often, the function that performs the combination requires every iterator to have the same
/// type.  However, sometimes we need to combine iterators with different types.
/// This enum allows us to do so.  
/// Given an iterator `iter1` of type `I1` and an iterator `iter2` of type `I2`, we can create two
/// new iterators of the same type: `TwoTypeIterator::Version1( iter 1 )` and `TwoTypeIterator::Version2( iter 2 )`.
/// The new wrappers will iterator over exactly the same items as the iterators they contain.
pub enum IntoTwoTypeIterator < I1, I2 > 
where   I1:  IntoIterator,
        I2:  IntoIterator < Item = I1::Item >,
{
    Version1( I1 ),
    Version2( I2 ),
}

impl < I1, I2 > 
        
    IntoIterator       
    for IntoTwoTypeIterator < I1, I2 > 
    where   I1:  IntoIterator,
            I2:  IntoIterator < Item = I1::Item >,    
    
{
    type Item = I1::Item;
    type IntoIter = TwoTypeIterator< I1::IntoIter, I2::IntoIter >;
    
    fn into_iter(self) -> Self::IntoIter {
        match self {
            Self::Version1( iter1 ) => TwoTypeIterator::Version1( iter1.into_iter() ),
            Self::Version2( iter2 ) => TwoTypeIterator::Version2( iter2.into_iter() ),            
        }
    }
}

impl < I1, I2 > 
        
    Clone       
    for IntoTwoTypeIterator < I1, I2 > 
    where   I1:  Clone + Iterator,
            I2:  Clone + Iterator < Item = I1::Item >,    
    
{
    fn clone(&self) -> Self {
        match self {
            Self::Version1( iter1 ) => { Self::Version1( iter1.clone() ) },
            Self::Version2( iter2 ) => { Self::Version2( iter2.clone() ) },            
        }
    }
}





// 2-TYPE ITERATOR
// ---------------------------------------------------------------------------

/// An enum that combines two different iterator types into a single iterator type.
///
/// There are many operations available in Rust for combining iterators (chaining, merging, etc.).
/// Quite often, the function that performs the combination requires every iterator to have the same
/// type.  However, sometimes we need to combine iterators with different types.
/// This enum allows us to do so.  
/// Given an iterator `iter1` of type `I1` and an iterator `iter2` of type `I2`, we can create two
/// new iterators of the same type: `TwoTypeIterator::Version1( iter 1 )` and `TwoTypeIterator::Version2( iter 2 )`.
/// The new wrappers will iterator over exactly the same items as the iterators they contain.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum TwoTypeIterator < I1, I2 > 
where   I1:  Iterator,
        I2:  Iterator < Item = I1::Item >,
{
    Version1( I1 ),
    Version2( I2 ),
}

impl < I1, I2 > 
        
    Iterator       
    for TwoTypeIterator < I1, I2 > 
    where   I1:  Iterator,
            I2:  Iterator < Item = I1::Item >,    
    
{
    type Item = I1::Item;
    
    fn next( & mut self ) -> Option< Self::Item > {
        use TwoTypeIterator::*; // not sure why we have to include it but this seems to be necessary, c.f. https://stackoverflow.com/questions/33925232/how-to-match-over-self-in-an-enum

        match self {
            Version1( iter1 ) =>   {
                iter1.next()
            },
            Version2( iter2 ) =>   {
                iter2.next()
            }
        }
    }
}

impl < I1, I2 > 
        
    Clone for 
    
    TwoTypeIterator < I1, I2 > 
    
    where   I1:  Clone + Iterator,
            I2:  Clone + Iterator < Item = I1::Item >,    
    
{
    fn clone(&self) -> Self {
        match self {
            Self::Version1( iter1 ) => { Self::Version1( iter1.clone() ) },
            Self::Version2( iter2 ) => { Self::Version2( iter2.clone() ) },            
        }
    }
}


impl < I1, I2 > 
        
    PeekUnqualified for 
    
    TwoTypeIterator < I1, I2 > 
    
    where   I1:  PeekUnqualified + Iterator,
            I2:  PeekUnqualified + Iterator < Item = I1::Item >,    
{   
    fn peek_unqualified( &mut self ) -> Option< & Self::Item > {
        match self {
            Self::Version1( iter ) => { iter.peek_unqualified() }
            Self::Version2( iter ) => { iter.peek_unqualified() }            
        }
    }
}            


// OnlyDuplicates
// ---------------------------------------------------------------------------


/// Returns one copy of each item that repeats one or more times, and skips non-repeated items.
/// 
/// This struct wraps around another iterator.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::OnlyDuplicates;
/// 
/// let a   =   vec![ 1, 2, 2, 3, 3, 3, 4, 5, 5 ];
/// let b   =   vec![ 2, 3, 5 ];
/// 
/// let only_dups       =   OnlyDuplicates::new( a.iter() );
/// 
/// // Check that `only_dups` returns the same sequence of items as `b.iter()`
/// assert!(    only_dups.eq( b.iter() )    );
/// 
/// ```
#[derive(Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct OnlyDuplicates< I: Iterator >{
    iter_unwrapped: I,
    next_candidate: Option< I::Item >,
}

impl < I > 

    OnlyDuplicates< I > 

    where 
        I:          Iterator,
        I::Item:    std::cmp::PartialEq,
{
    pub fn new< J: IntoIterator< IntoIter = I > >( iterable_unwrapped: J ) -> OnlyDuplicates< I > { 
        let mut iter_unwrapped  =   iterable_unwrapped.into_iter();
        let next_candidate  =   iter_unwrapped.next();
        OnlyDuplicates { iter_unwrapped, next_candidate } 
    }
}

impl    < I > 
        
    Iterator for

    OnlyDuplicates< I >

    where 
        I:          Iterator,
        I::Item:    std::cmp::PartialEq,        
{
    type Item = I::Item;

    fn next( &mut self ) -> Option< Self::Item > {

        // initialize two unasigned variables
        let mut is_duplicate; 
        let mut following_candidate;

        loop {
            self.next_candidate.as_ref()?;

            following_candidate     =   self.iter_unwrapped.next();
            is_duplicate            =   self.next_candidate == following_candidate;

            //  iterate over all duplicates
            while following_candidate == self.next_candidate {
                following_candidate = self.iter_unwrapped.next();
            }

            if is_duplicate{
                // swap bindings of our two variables
                std::mem::swap( &mut self.next_candidate, &mut following_candidate );
                // return the old candidate (whose value is now bound to `following_candidate`)
                return following_candidate
            } else {
                self.next_candidate     =   following_candidate;
            }
        }
    }
}    


// SkipDuplicates
// ---------------------------------------------------------------------------


/// Skips items that repeat any number of times.
/// 
/// This struct wraps around another iterator.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::SkipDuplicates;
/// use itertools;
/// 
/// // EXAMPLE 1
/// let a   =   vec![ 1, 2, 2, 3, 3, 3, 4, 5 ];
/// let b   =   vec![ 1, 4, 5 ];
/// 
/// let skip_dups       =   SkipDuplicates::new( a.iter() );
/// 
/// // Check that `skip_dups` returns the same sequence of items as `b.iter()`
/// itertools::assert_equal(    skip_dups,  b.iter()     );
/// 
/// // EXAMPLE 2
/// let a   =   vec![ 1, 2, 2, 3, 3, 3, 4, 5, 5 ];
/// let b   =   vec![ 1, 4 ];
/// 
/// let skip_dups       =   SkipDuplicates::new( a.iter() );
/// 
/// // Check that `skip_dups` returns the same sequence of items as `b.iter()`
/// itertools::assert_equal(    skip_dups,  b.iter()     );
/// 
/// // EXAMPLE 3
/// let a: Vec<usize>   =   vec![  ];
/// let b: Vec<usize>   =   vec![  ];
/// 
/// let skip_dups       =   SkipDuplicates::new( a.iter() );
/// 
/// // Check that `skip_dups` returns the same sequence of items as `b.iter()`
/// itertools::assert_equal(    skip_dups,  b.iter()     );
/// 
/// // EXAMPLE 4
/// let a: Vec<usize>   =   vec![ 1, 1, 2, 2 ];
/// let b: Vec<usize>   =   vec![  ];
/// 
/// let skip_dups       =   SkipDuplicates::new( a.iter() );
/// 
/// // Check that `skip_dups` returns the same sequence of items as `b.iter()`
/// itertools::assert_equal(    skip_dups,  b.iter()     );
/// 
/// ```
#[derive(Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct SkipDuplicates< I: Iterator >{
    iter_unwrapped: I,
    next_candidate: Option< I::Item >,
}

impl < I > 

    SkipDuplicates< I > 

    where 
        I:          Iterator,
        I::Item:    std::cmp::PartialEq,
{
    pub fn new< J: IntoIterator< IntoIter = I > >( iter_unwrapped: J ) -> SkipDuplicates< I > { 
        let mut iter_unwrapped  =   iter_unwrapped.into_iter();
        let next_candidate  =   iter_unwrapped.next();
        SkipDuplicates { iter_unwrapped, next_candidate } 
    }
}

impl    < I > 
        
    Iterator for

    SkipDuplicates< I >

    where 
        I:          Iterator,
        I::Item:    std::cmp::PartialEq,        
{
    type Item = I::Item;

    fn next( &mut self ) -> Option< Self::Item > {

        // initialize two unasigned variables
        let mut is_duplicate; 
        let mut following_candidate;

        loop {
            self.next_candidate.as_ref()?;

            following_candidate     =   self.iter_unwrapped.next();
            is_duplicate            =   self.next_candidate == following_candidate;

            //  iterate over all duplicates
            while following_candidate == self.next_candidate {
                following_candidate = self.iter_unwrapped.next();
            }

            if is_duplicate{
                self.next_candidate     =   following_candidate;
            } else {
                // swap bindings of our two variables
                std::mem::swap( &mut self.next_candidate, &mut following_candidate );
                // return the old candidate (whose value is now bound to `following_candidate`)
                return following_candidate
            }
        }
    }
}    


//  SET OPERTIONS ON ORDERED ITERATORS
//  ==================================================================

/// Determine if one ordered iterator contains another.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::ordered_iter_contains;
/// 
/// let a   =   Vec::<usize>::new();
/// let b   =   vec![ 0, 1,     ];
/// let c   =   vec![    1, 2,  ];
/// let d   =   vec![ 0, 1, 2   ];
/// 
/// assert!(   ordered_iter_contains( d.iter(), a.iter() ) );
/// assert!(   ordered_iter_contains( d.iter(), b.iter() ) );
/// assert!(   ordered_iter_contains( d.iter(), c.iter() ) );
/// assert!( ! ordered_iter_contains( c.iter(), b.iter() ) );
/// assert!( ! ordered_iter_contains( b.iter(), c.iter() ) );
/// assert!( ! ordered_iter_contains( a.iter(), d.iter() ) );
/// ```
pub fn ordered_iter_contains<I>( mut iter_large: I, iter_small: I ) -> bool 
    where 
        I:          Iterator,
        I::Item:    Ord + PartialEq + std::fmt::Debug,
{

    // use size_hints;
    for item in iter_small{
        loop{
            match iter_large.next(){
                None => { return false }
                Some( x ) =>{
                    println!("{:?}", (&x, &item));
                    match x.cmp( &item ) {
                        std::cmp::Ordering::Less    =>  { continue }
                        std::cmp::Ordering::Equal   =>  { break },
                        std::cmp::Ordering::Greater =>  { return false }
                    }
                }
            }   
        }
    }
    true
}



/// Intersect a collection or ordered operators
/// 
/// See the documentation for [IntersectOrderedIterators] for further details.
pub fn intersect_n_ordered_iterators_of_keyval_pairs< T, I >( 
            iterators: T,
        ) 
        -> 
        IntersectOrderedIterators
        < I, OrderOperatorByKey > 
    where
        T:                                              IntoIterator<Item=I>,
        I:                                              Iterator,
        I::Item:                                        KeyValGet,
        <<I as Iterator>::Item as KeyValGet>::Key:      Ord,  // this is the same as I::Item::Key: Ord, Rust just wants us to use fully qualified syntax
{
    let mut iterators: Vec<_>   =   iterators.into_iter().collect();
    iterators.shrink_to_fit();
    IntersectOrderedIterators{ iterators, order_operator: OrderOperatorByKey::new(), confirmed_empty: false }
}

/// Iterates over the intersection of `n` ordered iterators
/// 
/// Each iterator must return elements in strictly ascending order.
/// 
/// Only of the first iterator are returned (this is relevan, for example,
/// with items are (key,val) pairs, and two items are considered equal if
/// the have equal keys).
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::IntersectOrderedIterators;
/// 
/// //  TEST 1: intersection of 3 iterators
/// 
/// let ordered_sets    =   vec![
///                             vec![ 0, 1, 2, 3, 4,    ],
///                             vec![       2, 3, 4,    ],
///                             vec![ 0,    2,    4, 5, ],
///                         ];
/// let intersection    =   IntersectOrderedIterators::new(
///                             ordered_sets.iter().map(|x| x.iter() )
///                         );
/// let v               =   vec![ 2, 4, ];
/// 
/// itertools::assert_equal( v.iter(), intersection );
/// 
/// //  TEST 2: intersection of 0 iterators
/// 
/// let ordered_sets: Vec<Vec<usize>>    =   vec![];
/// let intersection    =   IntersectOrderedIterators::new(
///                             ordered_sets.iter().map(|x| x.iter() )
///                         );
/// let v: Vec<usize>   =   vec![];
/// 
/// itertools::assert_equal( v.iter(), intersection );
/// ```
#[derive(Clone,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct IntersectOrderedIterators< I, OrderOperator > 
    where
        I:                  Iterator,
{
    iterators:              Vec< I >,
    order_operator:         OrderOperator,
    confirmed_empty:        bool, // this is used to block future calls to `next` if the intersection is empty
}

//  Implement struct for custom orders
//  ----------------------------------
impl < I, OrderOperator >

    IntersectOrderedIterators< I, OrderOperator >

    where
        I:                  Iterator,
        OrderOperator:      JudgeOrder< I::Item >,
{
    pub fn with_custom_order< T: IntoIterator<Item=J>, J: IntoIterator< IntoIter = I > >( 
                iterators: T,
                order_operator: OrderOperator
            ) 
            -> 
            IntersectOrderedIterators
                < I, OrderOperator > 
    {
        let mut iterators: Vec<_>   =   iterators.into_iter().map(|iterable| iterable.into_iter()).collect();
        iterators.shrink_to_fit();
        IntersectOrderedIterators{ iterators, order_operator, confirmed_empty: false }
    }
}  

//  Implement struct for default orders
//  -----------------------------------

impl < I >

    IntersectOrderedIterators< I, OrderOperatorAuto >

    where
        I:                  Iterator,
        I::Item:            Ord,
{
    pub fn new< T: IntoIterator<Item=I> >( 
                iterators: T,
            ) 
            -> 
            IntersectOrderedIterators
                < I, OrderOperatorAuto > 
    {
        let mut iterators: Vec<_>   =   iterators.into_iter().collect();
        iterators.shrink_to_fit();        
        IntersectOrderedIterators{ iterators, order_operator: OrderOperatorAuto, confirmed_empty: false }
    }
}  


impl < I, OrderOperator > 

    Iterator for 

    IntersectOrderedIterators< I, OrderOperator >

    where
        I:                  Iterator,
        OrderOperator:      JudgeOrder< I::Item >,
{
    type Item       =   I::Item;

    fn next( &mut self ) -> Option< Self::Item > {

        if self.confirmed_empty { return None } // if the intersection is empty then we do not need to search for more items        

        if self.iterators.is_empty() { return None } // the intersection of zero iterators is empty

        let mut source_old = 0;
        let mut source_new = 0;
        let mut item_old;
        
        match self.iterators[source_old].next() {
            None => { 
                self.confirmed_empty = true; // if the first iterator is empty then the intersection is empty
                return None
            } Some( item_new ) => { 
                item_old = item_new;

                loop {
                    source_new = ( source_new + 1 ) % self.iterators.len();
                    if source_new == source_old { return Some( item_old ) }

                    loop {
                        match self.iterators[source_new].next() {
                            None => { 
                                self.confirmed_empty = true; // if the next iterator is empty then the intersection is empty
                                return None 
                            } Some( item_new ) => {
                                match self.order_operator.judge_cmp( &item_new, &item_old ) {
                                    Ordering::Less => { 
                                        // increment upward until we meet or exceed item_old
                                        continue; 
                                    } Ordering::Equal => { 
                                        // stop incrementing if item_new equals item_old
                                        if source_new == 0 { item_old = item_new } // this change ensures we only return items in the first iterator
                                        break; // break the loop
                                    } Ordering::Greater => {
                                        // if item_new exceeds item_old then item_old is not in the intersection;
                                        // therefore re-initialize a search                                
                                        item_old = item_new;
                                        source_old = source_new;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // // The intersection of 0 iterators is emtpy
            // if self.iterators.is_empty() { return None }

            // let mut min_index;
            // let mut min_val;        
            // let mut all_min;
            // let mut indices_to_check;

            // loop {
            //     // find the index of the minimum value
            //     min_index   =   find_min( self.next_elts.iter() ).unwrap();

            //     // pull out the next element from the iterator at index min_index, and swap its bound value with self.next_elts[min_index]
            //     min_val     =   self.iterators[ min_index ].next(); // this will become the true min val on the next line
            //     std::mem::swap(&mut min_val, &mut self.next_elts[min_index]); // now min_val equals the minimum value of next_elts, before we modified next_elts

            //     // If one of the iterators is empty then (i) min_val = None, and (ii) there are no more elements in the intersection
            //     if min_val.is_none() { return None }
        
            //     // Check whether every element equals min, and replace the min vals
            //     all_min     =   true;
            //     indices_to_check = (0..min_index).chain( min_index+1 .. self.next_elts.len()); // iterate over every index except min_index
            //     for p in indices_to_check {
            //         match self.next_elts[p] == min_val {
            //             true    =>  { self.next_elts[p] = self.iterators[p].next()  },
            //             false   =>  { all_min = false }
            //         }
            //     }

            //     if all_min { return min_val }
            // }
        }
    }
}


/// Returns an iterator that runs over elments in the inersection of `iter1` and `iter2`.
/// 
/// Each iterator must be sorted in **strictly** ascending order.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::utilities::iterators::general::intersect_ordered_iterators;
/// 
/// let a = [0,1,2].iter().cloned();
/// let b = [1,2,3].iter().cloned();
/// 
/// assert!( intersect_ordered_iterators(a,b).eq( [1,2] ) );
/// ```
/// 
/// # Opportunities for improvement
/// 
/// Terminate search after one iterator terminates
pub fn intersect_ordered_iterators< I >( iter1: I, iter2: I )
    ->
    OnlyDuplicates< Merge< I, I > >

    where
        I:          Iterator,
        I::Item:    PartialOrd,
{
    OnlyDuplicates::new( iter1.merge( iter2 ))
}

/// Returns an iterator that runs over elments in the symetric difference of `iter1` and `iter2`
/// 
/// Each iterator must be sorted in **strictly** ascending order.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::utilities::iterators::general::symmetric_difference_of_ordered_iterators;
/// 
/// let a = [0,1,2].iter().cloned();
/// let b = [1,2,3].iter().cloned();
/// 
/// assert!( symmetric_difference_of_ordered_iterators(a,b).eq( [0,3] ) );
/// ```
pub fn symmetric_difference_of_ordered_iterators< I >( iter1: I, iter2: I )
    ->
    SkipDuplicates< Merge< I, I > >

    where
        I:          Iterator,
        I::Item:    PartialOrd,
{
    SkipDuplicates::new( iter1.merge( iter2 ))
}





// OncePeekable
// ---------------------------------------------------------------------------

#[derive(Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct OncePeekable< T > { item_opt: Option< T > }

impl < T > OncePeekable< T > {
    pub fn new( item: T ) -> OncePeekable< T > { OncePeekable{ item_opt: Some(item) } }
}

impl < T > Iterator for OncePeekable< T > {
    type Item = T;
    fn next( &mut self ) -> Option< Self::Item > {
        self.item_opt.take() // replace the internally stored item with None, and return whatever was stored there previously
    }
}

impl < T > 

    PeekUnqualified for 

    OncePeekable< T >
{
    fn peek_unqualified( &mut self ) -> Option < & Self::Item > { self.item_opt.as_ref() }
}







// VECTOR WRAPPED IN AN ITERATOR
// ---------------------------------------------------------------------------

/// Wrapper containing an `Arc<Vec>` and a pointer `usize` to an elment of the `Arc<Vec>`; iterates over the elements of
/// the vector, incrementing the pointer as it goes.
/// 
/// This struct is essentially a vector that implements `Iterator`.  It differs from 
/// `Arc<Vec<T>>` because `Arc<Vec<T>>` does not implement `Iterator`.  Moreover, it differs 
/// from the iterator we usually associate with a vector, `Iter< '_, T >` in that 
/// `IterWrappedArcVec< T >` has no lifetime parameter.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::IterWrappedArcVec;
/// use std::sync::Arc;
/// 
/// // error[E0106]: missing lifetime specifier
/// // fn iter_from_nothing() -> std::slice::Iter<'_, usize> { ( vec![ 0, 0, 0 ] ).iter() } 
/// 
/// // no error
/// fn arc_vec_wrapped_in_iter_from_nothing() -> IterWrappedArcVec< usize > { 
///     IterWrappedArcVec::new( Arc::new(vec![2,3,4]) ) 
/// } 
/// 
/// let iter = arc_vec_wrapped_in_iter_from_nothing();
/// itertools::assert_equal(iter, vec![2,3,4])
/// ```
#[derive(Clone,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct IterWrappedArcVec< T > { unwrapped_vec: Arc<Vec< T >>, pointer_to_next_val: usize }

impl < T > IterWrappedArcVec < T > {
    /// Create a new `IterWrappedVec< T >` from a `Vec< T >`.  
    pub fn new( unwrapped_vec: Arc<Vec< T >> ) -> IterWrappedArcVec < T > { IterWrappedArcVec{ unwrapped_vec, pointer_to_next_val: 0  } }

    /// Returns an iterature that runs backward across all elements of the vector
    /// 
    /// The result is the same as calling `IterWrappedVecReverse::new( vec )`, where `vec` is the underlying vector.
    /// 
    /// **This behavoir differs from the `.rev()` method for two-sided iterators** because the new iterator
    /// will run over all the elements in the vector. By contrast, if we have a Rust Range, iterate forwrad through a few
    /// elements, and then iterate backward, the backward iterator will not return the elements that were already returned
    /// by the forward iterations.
    pub fn reset_and_run_backward( self ) -> IterWrappedArcVecReverse< T > { IterWrappedArcVecReverse::new( self.unwrapped_vec )  }    
}

impl < T: Clone > Iterator for IterWrappedArcVec< T > { 
    type Item = T;

    fn next( &mut self ) -> Option< Self::Item > {
        match self.pointer_to_next_val < self.unwrapped_vec.len() {
            true    =>  { 
                let a = self.unwrapped_vec[ self.pointer_to_next_val ].clone(); // retreive the next item
                self.pointer_to_next_val += 1; // increment the pointer
                Some( a ) // return the item
            }
            false   =>  { None }
        }
    }
}






/// Wrapper containing an `Arc<Vec>` and a pointer `usize`; iterates over the elements of
/// the vector *in reverse order**, decrementing the pointer as it goes.
/// 
/// This struct is essentially a vector that implements `Iterator`.  It differs from 
/// `Vec<T>` because `Vec<T>` does not implement `Iterator`.  Moreover, it differs 
/// from the iterator we usually associate with a vector, `Iter< '_, T >` in that 
/// `IterWrappedArcVecReverse< T >` owns its own data.  Thus, for example, we can write
/// a function which takes no arguments, but which returns a nonempty 
/// `IterWrappedArcVecReverse< T >`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::IterWrappedArcVecReverse;
/// use std::sync::Arc;
/// 
/// // error[E0106]: missing lifetime specifier
/// // fn iter_from_nothing() -> Rev<std::slice::Iter<'_, usize>> { ( vec![ 0, 0, 0 ] ).iter() } 
/// 
/// // no error
/// fn vec_wrapped_in_iter_from_nothing() -> IterWrappedArcVecReverse< usize > { 
///     IterWrappedArcVecReverse::new( Arc::new(vec![2,3,4]) ) 
/// } 
/// 
/// // notice that the iterator returns items in reverse order
/// let iter = vec_wrapped_in_iter_from_nothing();
/// itertools::assert_equal(iter, vec![4,3,2])
/// ```
#[derive(Clone,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct IterWrappedArcVecReverse< T > { unwrapped_vec: Arc<Vec< T >>, steps_back_to_next_val: usize }

impl < T > IterWrappedArcVecReverse < T > {
    /// Create a new `IterWrappedVec< T >` from a `Vec< T >`.  
    pub fn new( unwrapped_vec: Arc<Vec< T >> ) -> IterWrappedArcVecReverse < T > { IterWrappedArcVecReverse{ unwrapped_vec, steps_back_to_next_val: 1  } }

    /// Returns an iterature that runs forward across all elements of the vector
    /// 
    /// The result is the same as calling `IterWrappedVec::new( vec )`, where `vec` is the underlying vector.
    /// 
    /// **This behavoir differs from the `.rev()` method for two-sided iterators** because the new iterator
    /// will run over all the elements in the vector. By contrast, if we have a Rust Range, iterate backward through a few
    /// elements, and then iterate forward, the forward iterator will not return the elements that were already returned
    /// by the backward iterations.
    pub fn reset_and_run_forward( self ) -> IterWrappedArcVec< T > { IterWrappedArcVec::new( self.unwrapped_vec )  }
}

impl < T: Clone > Iterator for IterWrappedArcVecReverse< T > { 
    type Item = T;

    fn next( &mut self ) -> Option< Self::Item > {
        let len     =   self.unwrapped_vec.len();
        match self.steps_back_to_next_val <= len {
            true    =>  { 
                let a = self.unwrapped_vec[ len - self.steps_back_to_next_val ].clone(); // retreive the next item
                self.steps_back_to_next_val += 1; // increment the pointer
                Some( a ) // return the item
            }
            false   =>  { None }
        }
    }
 }











// VECTOR WRAPPED IN AN ITERATOR
// ---------------------------------------------------------------------------

/// Wrapper containing a `Vec` and a pointer `usize`; iterates over the elements of
/// the vector, incrementing the pointer as it goes.
/// 
/// This struct is essentially a vector that implements `Iterator`.  It differs from 
/// `Vec<T>` because `Vec<T>` does not implement `Iterator`.  Moreover, it differs 
/// from the iterator we usually associate with a vector, `Iter< '_, T >` in that 
/// `IterWrappedVec< T >` owns its own data.  Thus, for example, we can write
/// a function which takes no arguments, but which returns a nonempty 
/// `IterWrappedVec< T >`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::IterWrappedVec;
/// 
/// // error[E0106]: missing lifetime specifier
/// // fn iter_from_nothing() -> std::slice::Iter<'_, usize> { ( vec![ 0, 0, 0 ] ).iter() } 
/// 
/// // no error
/// fn vec_wrapped_in_iter_from_nothing() -> IterWrappedVec< usize > { 
///     IterWrappedVec::new( vec![2,3,4] ) 
/// } 
/// 
/// let iter = vec_wrapped_in_iter_from_nothing();
/// itertools::assert_equal(iter, vec![2,3,4])
/// ```
#[derive(Clone,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct IterWrappedVec< T > { unwrapped_vec: Vec< T >, pointer_to_next_val: usize }

impl < T > IterWrappedVec < T > {
    /// Create a new `IterWrappedVec< T >` from a `Vec< T >`.  
    pub fn new( unwrapped_vec: Vec< T > ) -> IterWrappedVec < T > { IterWrappedVec{ unwrapped_vec, pointer_to_next_val: 0  } }

    /// Returns an iterature that runs backward across all elements of the vector
    /// 
    /// The result is the same as calling `IterWrappedVecReverse::new( vec )`, where `vec` is the underlying vector.
    /// 
    /// **This behavoir differs from the `.rev()` method for two-sided iterators** because the new iterator
    /// will run over all the elements in the vector. By contrast, if we have a Rust Range, iterate forwrad through a few
    /// elements, and then iterate backward, the backward iterator will not return the elements that were already returned
    /// by the forward iterations.
    pub fn reset_and_run_backward( self ) -> IterWrappedVecReverse< T > { IterWrappedVecReverse::new( self.unwrapped_vec )  }    
}

impl < T: Clone > Iterator for IterWrappedVec< T > { 
    type Item = T;

    fn next( &mut self ) -> Option< Self::Item > {
        match self.pointer_to_next_val < self.unwrapped_vec.len() {
            true    =>  { 
                let a = self.unwrapped_vec[ self.pointer_to_next_val ].clone(); // retreive the next item
                self.pointer_to_next_val += 1; // increment the pointer
                Some( a ) // return the item
            }
            false   =>  { None }
        }
    }
}





/// Wrapper containing a `Vec` and a pointer `usize`; iterates over the elements of
/// the vector *in reverse order**, decrementing the pointer as it goes.
/// 
/// This struct is essentially a vector that implements `Iterator`.  It differs from 
/// `Vec<T>` because `Vec<T>` does not implement `Iterator`.  Moreover, it differs 
/// from the iterator we usually associate with a vector, `Iter< '_, T >` in that 
/// `IterWrappedVecReverse< T >` owns its own data.  Thus, for example, we can write
/// a function which takes no arguments, but which returns a nonempty 
/// `IterWrappedVecReverse< T >`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::IterWrappedVecReverse;
/// 
/// // error[E0106]: missing lifetime specifier
/// // fn iter_from_nothing() -> Rev<std::slice::Iter<'_, usize>> { ( vec![ 0, 0, 0 ] ).iter() } 
/// 
/// // no error
/// fn vec_wrapped_in_iter_from_nothing() -> IterWrappedVecReverse< usize > { 
///     IterWrappedVecReverse::new( vec![2,3,4] ) 
/// } 
/// 
/// // notice that the iterator returns items in reverse order
/// let iter = vec_wrapped_in_iter_from_nothing();
/// itertools::assert_equal(iter, vec![4,3,2])
/// ```
#[derive(Clone,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct IterWrappedVecReverse< T > { unwrapped_vec: Vec< T >, steps_back_to_next_val: usize }

impl < T > IterWrappedVecReverse < T > {
    /// Create a new `IterWrappedVec< T >` from a `Vec< T >`.  
    pub fn new( unwrapped_vec: Vec< T > ) -> IterWrappedVecReverse < T > { IterWrappedVecReverse{ unwrapped_vec, steps_back_to_next_val: 1  } }

    /// Returns an iterature that runs forward across all elements of the vector
    /// 
    /// The result is the same as calling `IterWrappedVec::new( vec )`, where `vec` is the underlying vector.
    /// 
    /// **This behavoir differs from the `.rev()` method for two-sided iterators** because the new iterator
    /// will run over all the elements in the vector. By contrast, if we have a Rust Range, iterate backward through a few
    /// elements, and then iterate forward, the forward iterator will not return the elements that were already returned
    /// by the backward iterations.
    pub fn reset_and_run_forward( self ) -> IterWrappedVec< T > { IterWrappedVec::new( self.unwrapped_vec )  }
}

impl < T: Clone > Iterator for IterWrappedVecReverse< T > { 
    type Item = T;

    fn next( &mut self ) -> Option< Self::Item > {
        let len     =   self.unwrapped_vec.len();
        match self.steps_back_to_next_val <= len {
            true    =>  { 
                let a = self.unwrapped_vec[ len - self.steps_back_to_next_val ].clone(); // retreive the next item
                self.steps_back_to_next_val += 1; // increment the pointer
                Some( a ) // return the item
            }
            false   =>  { None }
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
    
    use crate::utilities::functions::evaluate::EvaluateFunction;


    #[test]
    fn test_skip_until() {
        use crate::utilities::iterators::general::SkipUntil;

        let vec = vec![ 0, 1, 2, 1, 0];
        let iter = vec.iter().cloned().peekable().skip_until( |x| x > &1);
        itertools::assert_equal( iter, vec![ 2, 1, 0 ])        
    }


    #[test]
    fn test_filter_map_by_transform() {
        use crate::utilities::iterators::general::FilterMapByTransform;

        // define a struct to execute the filter-map function
        struct ShiftDown;

        // implement the `FilterMapTransform` trait
        impl EvaluateFunction< usize, Option< i32 > > for ShiftDown{
            fn evaluate_function( &self, input: usize ) -> Option< i32 > {
                match input > 3 {
                    true    => { None },
                    false   => { Some( input as i32 - 3 ) }
                }
            }
        }

        // define an iterator
        let iter = 0..9;

        // create a `filter-mapped` object
        let iter_transformed = FilterMapByTransform::new( iter, ShiftDown );

        // verify this object is correct        
        itertools::assert_equal( iter_transformed, vec![ -3, -2, -1, 0 ])        
    }    



    #[test]
    fn test_filter_map_by_transform_2() {
        use crate::utilities::iterators::general::{FilterMapByTransformTrait};
        use std::collections::HashMap;
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        hash.insert( 2usize, -2);
        
        let iter_data = vec![0, 1, 2, 3];
        let iter = iter_data.iter().cloned();
        
        let iter_filtermapped = iter.filter_map_by_transform( & hash );
        itertools::assert_equal( iter_filtermapped.cloned(), vec![-1, -2] );  
    }   

    #[test]
    fn test_implementation_of_filtermapobjectmethod_on_hashmap() {
        use crate::utilities::iterators::general::{FilterMapByTransformTrait};
        use std::collections::HashMap;
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        hash.insert( 2usize, -2);
        
        let iter_data = vec![0, 1, 2, 3];
        let iter = iter_data.iter().cloned();
        
        let iter_filtermapped = iter.filter_map_by_transform( & hash );
        itertools::assert_equal( iter_filtermapped.cloned(), vec![-1, -2] );       
    }

    #[test]
    fn test_map_by_transform() {
        use crate::utilities::iterators::general::MapByTransform;
        use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;

        let entrywise_transform_unwrapped = |x : usize| -> i32 { x as i32 };
        let entrywise_transform_wrapped = EvaluateFunctionFnMutWrapper::new(entrywise_transform_unwrapped);
        let vec = vec![ 1, 2, 3];
        let iter = vec.iter().cloned();

        let transformed_entrywise = MapByTransform::new( iter, entrywise_transform_wrapped );
        itertools::assert_equal( transformed_entrywise, vec![ 1i32, 2, 3] );        
    }
  
}    


//  ---------------------------------------------------------------------
//  Unit tests
//  ---------------------------------------------------------------------

#[cfg(test)]
mod unit_tests {
}