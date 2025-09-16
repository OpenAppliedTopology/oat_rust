//! Partial and total orders
//! 
//! See also [`indexing_and_bijection`](crate::utilities::indexing_and_bijection).
//! 
//! Rust offers two traits to compare the order of elements: [PartialOrd] and [Ord]. For example, if you want to know if
//! `a` is less than `b`, you can run
//! 
//! ```
//! let a = 0; 
//! let b = 1;
//! 
//! // This will assign a value of `true` to `is_less`.
//! let is_less =   a.lt(&b);
//! ```
//! 
//! However, [PartialOrd] and [Ord] offer limited flexibility,
//! especially in situations where multiple different orders might apply to the same type of object (for example, we might want to
//! sort simplices according to dimension, diameter, or lexicographic order, or some combination of the three).
//! 
//! To provide greater flexibility, OAT provides two similar but different traits, [`JudgePartialOrder`] and [`JudgeOrder`].
//! The functions defined by these traits are similar to [PartialOrd] and [Ord], with one important difference: rather than
//! compare to objects directly, you define a third object which performs the comparison for you.
//! 
//! ```
//! use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorAutoReverse};
//! use oat_rust::utilities::order::JudgePartialOrder;
//! 
//! // Define a and b
//! let a = 0;
//! let b = 1;
//! 
//! // Define an order operator:
//! // this one impelements the default order on objects
//! let order_operator = OrderOperatorAuto;
//! 
//! // Use the order operator to compare a and b
//! assert!( order_operator.judge_lt( &a, &b ) );
//! 
//! // Define a DIFFERENT order operator
//! let order_operator_reverse = OrderOperatorAutoReverse;
//! 
//! // Use the new operator to compare a and b with a DIFFERENT order:
//! // this particular operator happens to reverse the default order on objects,
//! // so it will assert that b is less than a
//! assert!( order_operator_reverse.judge_lt( &b, &a ) );
//! ```
//! 
//! This gives greater
//! flexibility, since you can define many different order operators to compare elements of the same type.
//! 
//! # Terminology
//! 
//! We call an object that implements [JudgePartialOrder] and/or [JudgeOrder] an **order operator**.
//! It is analogous to the idea of a [RingOperator](crate::algebra::rings), in the sense that you use
//! order operators to get order information about other objects, just like you use ring operators to get
//! algebraic information about other structures.
//! 
//! # Design notes
//! 
//! In many situations [OrderOperatorAuto] is this simplest to use.  [OrderOperatorByKey] is also common for comparing vector entries; it has a smaller type signature than [OrderOperatorByKeyCustom], but for many operations, e.g. [U-match factorization](crate::algebra::matrices::operations::umatch), the larger signature allows users to pass customized order operators.
// //! object that "externalizes" that implements traits for
// //! comparing two elements -- traits such as [`JudgePartialOrder`], [`ComparePartialOrder`], and [`CompareOrder`].
// //! 
// //! Order plays an important role in the OAT for several reasons.  Some examples include
// //!
// //!   - it is easier to to write a lazy iterator for the sum of two sparse vectors (represented as lists of nonzero entires) if the lists are already sorted.
// //!   - min-heaps rely on order comparisons to function
// //!
// //! Often we find that a single application requires us to order a set of elements in one way at one time, and another way at another time.
// //! 
// //! To achieve this without changing the type of elements being ordered, we use a "helper" object that implements traits for
// //! comparing two elements -- traits such as [`JudgePartialOrder`], [`ComparePartialOrder`], and [`CompareOrder`].


use std::cmp::Ordering;


//  ---------------------------------------------------------------------------
//  SORTING
//  ---------------------------------------------------------------------------

//  SEE ALSO: the sequences_and_ordinals module
//  -------------------------------------------


/// Returns true iff `order_operator.judge_lt( a, b ) == true` for every consequtive pair of
/// elements `a`, `b`, in the slice.  Returns `true` by default for slices of length `<= 1`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::order::{ is_sorted_strictly, OrderOperatorAuto };
/// 
/// assert!( is_sorted_strictly(   &vec![ 0, 1, 2], & OrderOperatorAuto ) );
/// assert!( ! is_sorted_strictly( &vec![ 0, 1, 1], & OrderOperatorAuto ) );
/// ```
pub fn is_sorted_strictly < T, OrderOperator, > ( data: &[T], order_operator: & OrderOperator ) -> bool
    where OrderOperator: JudgePartialOrder<  T >,
{
    data.windows(2).all(|w| order_operator.judge_lt( &w[0], &w[1]) )
}



//  ---------------------------------------------------------------------------
//  ORDER OPERATOR TRAITS
//  ---------------------------------------------------------------------------


//  JudgePartialOrder + JudgeOrder
//  ---------------------------------------------------------------------------


/// Modeled off the Rust trait `std::cmp::PartialOrd`.
/// 
/// It is commonly the case that multiple *different* partial orders are
/// relevant to struct.  For example, one might want to order simplices in
/// ascending/descending lexicographic order, or first by dimension and second
/// by lexicographic order, etc.  The Rust traits for `Ord` and `PartialOrd`
/// only accomodate a single order per struct.  Therefore we create third-party
/// objects or "operators" to operationalize different orders.
/// 
///  REMARK: This trait is similar to (and extends) one copied from the itertools 
///  library for `kmerge`.  The OAT developers
///  suspect that the trait was introduced in order to supply type information 
///  that ordinary closures might miss.  At least, we tried removing the trait
///  from the IteratorsMergedInSortedOrder code we adapted from Itertools, and ran into some type 
///  erros at compilation.
/// 
/// # Examples
/// 
/// For general information about order operators see the documentation page for [order](crate::utilities::order).
/// 
/// The following example shows how to implement [JudgePartialOrder] on a new struct in order to make a new order operator.
/// 
/// ```
/// use oat_rust::utilities::order::{JudgeOrder, JudgePartialOrder};
/// use std::cmp::Ordering;
/// 
/// //  Define a new type of struct
/// //  ---------------------------        
/// pub struct MyOrderOperator;
/// 
/// //  Implement the JudgePartialOrder trait on our new struct
/// //  -------------------------------------------------------
/// impl JudgePartialOrder< usize > for MyOrderOperator {
/// 
///     // this function should return true iff left >= right
///     fn judge_ge(&self,left: &usize, right: &usize) -> bool {
///         left >= right 
///     }
///     // this function should return true iff left > right
///     fn judge_gt(&self,left: &usize, right: &usize) -> bool {
///         left > right
///     }
///     // this function should return true iff left <= right
///     fn judge_le(&self,left: &usize, right: &usize) -> bool {
///         left <= right
///     }
///     // this function should return true iff left < right
///     fn judge_lt(&self, left: &usize, right: &usize ) -> bool {
///         left < right
///     }
///     fn judge_partial_cmp( & self, left: & usize, right: & usize )  -> Option<Ordering> {
///         Some( left.cmp( right ) )
///     }
/// }
/// 
/// //  Implement the JudgeOrder trait on our new struct
/// //  ------------------------------------------------
/// impl JudgeOrder< usize > for MyOrderOperator {
///     // the judge_clamp function in OAT is modeled off of the clamp function here: https://doc.rust-lang.org/std/cmp/trait.Ord.html
///     fn judge_clamp(&self, clampee: usize, min: usize, max: usize) -> usize {
///         clampee.clamp(min,max) // here we just use Rust's built-in clamp function
///     }
///     // the judge_cmp function in OAT is modeled off of the clamp function here: https://doc.rust-lang.org/std/cmp/trait.Ord.html            
///     fn judge_cmp(&self, left: &usize, right: &usize ) -> Ordering {
///         left.cmp(right)
///     }
///     // should return the maximum of left and right
///     fn judge_max(&self, left: usize, right: usize) -> usize {
///         left.max(right)
///     }
///     // should return the minimum of left and right
///     fn judge_min(&self, left: usize, right: usize) -> usize {
///         left.min(right)
///     }
/// }
/// 
/// //  Use the trait
/// fn compare_order() {
///     let a = 1usize;
///     let b = 2usize;
///     let order_operator = MyOrderOperator;
/// 
///     if order_operator.judge_lt( &a, &b ) { println!("a is strictly less than b") }
///     else { println!("a is greater than or equal to b")}
/// }
/// ```
pub trait JudgePartialOrder< T > {
    /// Similar to the `.cmp` method in `std::cmp::PartialOrd`
    fn judge_partial_cmp
        ( & self, left: & T, right: & T ) 
        -> 
        Option<Ordering>;    
    // {
    //     let is_le = self.le(left, right);
    //     let is_ge = self.ge(left, right);
    //     if is_le {
    //         if is_ge { Some( Ordering::Equal ) }
    //         else { Some( Ordering::Less ) }
    //     }
    //     else {
    //         if is_ge { Some( Ordering::Greater ) }
    //         else { None }
    //     }
    // }
    /// Returns true iff left < right
    fn judge_lt(&self, left: &T, right: &T ) -> bool { self.judge_partial_cmp(left, right) == Some(Ordering::Less) }// { self.le( left, right ) && ! self.ge( left, right ) }
    /// Returns true iff left ≤ right    
    fn judge_le(&self,left: &T, right: &T) -> bool { match self.judge_partial_cmp(left,right){ Some(Ordering::Less) => true, Some(Ordering::Equal) => true, _ => false } }
    /// Returns true iff left > right    
    fn judge_gt(&self,left: &T, right: &T) -> bool { self.judge_partial_cmp(left, right) == Some(Ordering::Greater) } // { self.judge_lt( right, left ) }
    /// Returns true iff left ≥ right    
    fn judge_ge(&self,left: &T, right: &T) -> bool { match self.judge_partial_cmp(left,right){ Some(Ordering::Greater) => true, Some(Ordering::Equal) => true, _ => false } }

    /// Returns the reverse order operator
    fn reverse_order_judgements( self ) -> ReverseOrder< Self > 
        where Self: Sized
    {
        ReverseOrder::new( self )
    }
}

/// Modeled off the Rust trait `std::cmp::Ord`.
/// 
/// It is commonly the case that multiple *different* partial orders are
/// relevant to struct.  For example, one might want to order simplices in
/// ascending/descending lexicographic order, or first by dimension and second
/// by lexicographic order, etc.  The Rust traits for `Ord` and `PartialOrd`
/// only accomodate a single order per struct.  Therefore we create third-party
/// objects or "operators" to operationalize different orders.
/// 
/// # Examples
/// 
/// For general information about order operators see the documentation page for [order](crate::utilities::order).
/// 
/// The following example shows how to implement [JudgePartialOrder] on a new struct in order to make a new order operator.
/// 
/// ```
/// use oat_rust::utilities::order::{JudgePartialOrder, JudgeOrder};
/// 
/// use std::cmp::Ordering;
/// 
/// //  Define a new type of struct
/// //  ---------------------------        
/// pub struct MyOrderOperator;
/// 
/// //  Implement the JudgePartialOrder trait on our new struct
/// //  -------------------------------------------------------
/// impl JudgePartialOrder< usize > for MyOrderOperator {
/// 
///     // this function should return true iff left >= right
///     fn judge_ge(&self,left: &usize, right: &usize) -> bool {
///         left >= right 
///     }
///     // this function should return true iff left > right
///     fn judge_gt(&self,left: &usize, right: &usize) -> bool {
///         left > right
///     }
///     // this function should return true iff left <= right
///     fn judge_le(&self,left: &usize, right: &usize) -> bool {
///         left <= right
///     }
///     // this function should return true iff left < right
///     fn judge_lt(&self, left: &usize, right: &usize ) -> bool {
///         left < right
///     }
///     fn judge_partial_cmp( & self, left: & usize, right: & usize )  -> Option<Ordering> {
///         Some( left.cmp( right ) )
///     }
/// }
/// 
/// //  Implement the JudgeOrder trait on our new struct
/// //  ------------------------------------------------
/// impl JudgeOrder< usize > for MyOrderOperator {
///     // the judge_clamp function in OAT is modeled off of the clamp function here: https://doc.rust-lang.org/std/cmp/trait.Ord.html
///     fn judge_clamp(&self, clampee: usize, min: usize, max: usize) -> usize {
///         clampee.clamp(min,max) // here we just use Rust's built-in clamp function
///     }
///     // the judge_cmp function in OAT is modeled off of the clamp function here: https://doc.rust-lang.org/std/cmp/trait.Ord.html            
///     fn judge_cmp(&self, left: &usize, right: &usize ) -> Ordering {
///         left.cmp(right)
///     }
///     // should return the maximum of left and right
///     fn judge_max(&self, left: usize, right: usize) -> usize {
///         left.max(right)
///     }
///     // should return the minimum of left and right
///     fn judge_min(&self, left: usize, right: usize) -> usize {
///         left.min(right)
///     }
/// }
/// 
/// //  Use the trait
/// fn compare_order() {
///     let a = 1usize;
///     let b = 2usize;
///     let order_operator = MyOrderOperator;
/// 
///     if order_operator.judge_lt( &a, &b ) { println!("a is strictly less than b") }
///     else { println!("a is greater than or equal to b")}
/// }
/// ```
pub trait JudgeOrder< T > : JudgePartialOrder< T > 
{
    /// Similar to the `.cmp` method in `std::cmp::Ord`
    fn judge_cmp(&self, left: &T, right: &T ) -> Ordering {
        self.judge_partial_cmp(left, right).unwrap()
    }
    /// Return the greater of left, right; if they are equal, reutrh left
    fn judge_max(&self, left: T, right: T) -> T { match self.judge_lt(&left, &right) { true => { right }, false => { left } } }
    /// Return the lesser of left, rhx; if they are equal, return left
    fn judge_min(&self, left: T, right: T) -> T { match self.judge_le(&left, &right) { true => { left }, false => { right } } }
    /// Similar to Rust's `num::clamp`    
    fn judge_clamp(&self, clampee: T, min: T, max: T) -> T { 
        if self.judge_lt( &max, &min ) { panic!("Cannot call `self.clamp( clampee, min, max )` when max < min") };

        if self.judge_lt( &clampee, &min ) { min }
        else if self.judge_lt( &max, &clampee ) { return max }
        else { return clampee }
    }
}


//  ---------------------------------------------------------------------------
//  NONE GREATER
//  ---------------------------------------------------------------------------


/// Wrapper for `Option<T>`; implements `Ord` and `PartialOrd`, where `None` is the greatest (as opposed to the least) element
/// 
/// Compare this with `Option<T>`, where `None` is the least element
/// 
/// # Example
/// 
/// ```
/// use oat_rust::utilities::order::MakeNoneMaximum;
/// use std::cmp::Ordering;
/// 
/// let a   =   MakeNoneMaximum::from_val(0usize);
/// let b   =   MakeNoneMaximum::from_val(1usize);
/// let c   =   MakeNoneMaximum::from_opt(None);
/// 
/// assert!( a <  b );
/// assert!( a == a );
/// assert!( a <  c );
/// assert!( c == c );
/// 
/// 
/// assert_eq!( a.partial_cmp(&a), Some(Ordering::Equal) );
/// assert_eq!( a.partial_cmp(&b), Some(Ordering::Less ) );
/// assert_eq!( a.partial_cmp(&c), Some(Ordering::Less ) );
/// assert_eq!( c.partial_cmp(&c), Some(Ordering::Equal) );
/// 
/// ```
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MakeNoneMaximum< T > { pub opt: Option< T > }

impl < T > MakeNoneMaximum< T > {
    pub fn new( opt: Option<T> ) -> Self { MakeNoneMaximum{opt} }    
    
    /// Wraps `opt` in a `MakeNoneMaximum{ opt }` struct
    pub fn from_opt( opt: Option<T> ) -> Self { MakeNoneMaximum{opt} }

    /// Wraps `val` in a `MakeNoneMaximum{ opt: Some(val) }` struct
    pub fn from_val( val: T ) -> Self{ MakeNoneMaximum{ opt: Some(val) } }

    /// Returns a clone of the value wrapped in `MakeNoneMaximum`, if it exists
    /// 
    /// # Panics
    /// 
    /// Panics if `self.opt` is `None`
    pub fn val( &self ) -> T where T: Clone { self.opt.clone().unwrap() }

    /// Returns the internally stored `Option<T>`
    pub fn into_inner( self ) -> Option<T> { self.opt }
}

impl < T: Ord + PartialOrd + Eq > Ord for MakeNoneMaximum< T > {
    

    fn cmp(&self, other: &Self) -> Ordering {
        match (&self.opt, &other.opt) {
            (Some(a),Some(b)) => { a.cmp(b) }            
            (None,   None   ) => { Ordering::Equal   }
            (None,   Some(_)) => { Ordering::Greater }
            (Some(_),None   ) => { Ordering::Less    }             
        }
        // match self.opt {
        //     Some( self_val ) => {
        //         match other {
        //             Some( other_val ) => {

        //             }
        //         }
        //     }
        // }
    }
}

impl < T: Ord + PartialOrd > PartialOrd for MakeNoneMaximum< T > {
    
    /// ```
    /// use oat_rust::utilities::order::MakeNoneMaximum;
    /// 
    /// let a   =   MakeNoneMaximum::from_val(0usize);
    /// let b   =   MakeNoneMaximum::from_val(1usize);
    /// let c   =   MakeNoneMaximum::from_opt(None);
    /// 
    /// assert!( a <  b );
    /// assert!( a == a );
    /// assert!( a <  c );
    /// assert!( c == c );
    /// 
    /// ```
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (&self.opt, &other.opt) {
            (Some(a),Some(b)) => { a.partial_cmp(b) }            
            (None,   None   ) => { Some(Ordering::Equal  ) }
            (None,   Some(_)) => { Some(Ordering::Greater) }
            (Some(_),None   ) => { Some(Ordering::Less   ) }            
        }        
    }
}



//  ---------------------------------------------------------------------------
//  DEFAULT IMPLEMENTATION ON REFERENCES
//  ---------------------------------------------------------------------------



//  Auto-implement on references
//  ----------------------------

impl < 'a, T, OrderOperator: JudgePartialOrder< T > >

    JudgePartialOrder< T > for
    
    &'a OrderOperator
{
    fn judge_partial_cmp
        ( & self, left: & T, right: & T ) -> Option<Ordering> { (*self).judge_partial_cmp(left, right) }
}


//  Auto-implement on references
//  ----------------------------

impl < 'a, T, OrderOperator: JudgeOrder< T > >

    JudgeOrder< T > for
    
    &'a OrderOperator
{
    fn judge_cmp(&self, left: &T, right: &T ) -> Ordering {
        (*self).judge_cmp(left, right)
    }
}






//  ---------------------------------------------------------------------------
//  IMPLEMENTORS OF `JudgeOrder`
//  ---------------------------------------------------------------------------

//  WORKS FINE, BUT WE HAVE RESTURCTURED THE ORDER OPERATOR TRAITS, AND THIS STRUCT NO LONGER PROVIDES MUCH FUNCTIONALITY; OK TO DELETE

// /// Wrapper around any struct `ComparatorJudgePartialOrder` that implements the trait `JudgePartialOrder`; the wrapper then implements `JudgePartialOrder` and `JudgeOrder`.
// /// 
// /// # Examples
// /// 
// /// ```
// /// use oat_rust::utilities::order::{InferTotalOrderFromJudgePartialOrder, OrderOperatorAuto, JudgePartialOrder, JudgeOrder};
// /// 
// /// let order_operator      =   InferTotalOrderFromJudgePartialOrder::new( OrderOperatorAuto );
// /// assert_eq!( order_operator.le( &2, &2), true );
// /// assert_eq!( order_operator.judge_lt( &2, &2), false);
// /// assert_eq!( order_operator.le( &1, &2), true );
// /// assert_eq!( order_operator.judge_lt( &1, &2), true );
// /// assert_eq!( order_operator.le( &2, &1), false);
// /// assert_eq!( order_operator.judge_lt( &2, &1), false);
// /// ```
// #[derive(Clone, Copy, Debug, Eq, PartialEq)]
// pub struct InferTotalOrderFromJudgePartialOrder
//             < ComparatorJudgePartialOrder >
//             { comparator_lt: ComparatorJudgePartialOrder }

// impl < ComparatorJudgePartialOrder > 

//     InferTotalOrderFromJudgePartialOrder
//         < ComparatorJudgePartialOrder > 
// {
//     pub fn new( comparator_lt: ComparatorJudgePartialOrder ) -> InferTotalOrderFromJudgePartialOrder< ComparatorJudgePartialOrder > 
//         { InferTotalOrderFromJudgePartialOrder{ comparator_lt  } }
// }

// impl < ComparatorJudgePartialOrder, T > 

//     JudgePartialOrder
//         < T > for
//     InferTotalOrderFromJudgePartialOrder
//         < ComparatorJudgePartialOrder >

//     where 
//         ComparatorJudgePartialOrder:     JudgePartialOrder< T >,
// {
//     fn le(&self,left: &T, right: &T) -> bool { ! self.comparator_lt.judge_lt(right, left) }
//     fn lt(&self,left: &T, right: &T) -> bool { self.comparator_lt.judge_lt(left, right) }  
// }

// impl < ComparatorJudgePartialOrder, T > 

//     JudgeOrder
//         < T > for
//     InferTotalOrderFromJudgePartialOrder
//         < ComparatorJudgePartialOrder >

//     where 
//         ComparatorJudgePartialOrder:     JudgePartialOrder< T >,
// {}


//  ---------------------------------------------------------------------------
//  IMPLEMENTORS OF `JudgePartialOrder`
//  ---------------------------------------------------------------------------










//  Merged type
//  -------------


/// A type that can be either of two different order operators, `OrderOperator1` or `OrderOperator2`.
/// 
/// This enum will judge the order of elements according to whichever order operator it contains internally.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::order::{ ReverseOrder, OrderOperatorAuto, TwoTypeOrderOperator, JudgePartialOrder };
/// 
/// // a function to get either a forward or reverse order operator
/// let get_forward_operator = |forward| {
///     if forward { 
///         TwoTypeOrderOperator::Version1( OrderOperatorAuto ) 
///     } else { 
///         TwoTypeOrderOperator::Version2( ReverseOrder::new( OrderOperatorAuto ) )
///     }
/// };
/// 
/// // get some operators
/// let forward = get_forward_operator( true );
/// let reverse = get_forward_operator( false );
/// 
/// // check that they order elements correctly
/// assert!( forward.judge_lt( &1, &2) );
/// assert!( reverse.judge_lt( &2, &1) );
/// ```
#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum TwoTypeOrderOperator< OrderOperator1, OrderOperator2 >
{
    Version1( OrderOperator1 ),
    Version2( OrderOperator2 ),
}


impl < OrderOperator1, OrderOperator2, T > 

    JudgePartialOrder
        < T > for 
        
    TwoTypeOrderOperator
        < OrderOperator1, OrderOperator2 >

where 
    OrderOperator1: JudgePartialOrder< T >,
    OrderOperator2: JudgePartialOrder< T >,
{
    fn judge_partial_cmp( & self, left: & T, right: & T ) -> Option<Ordering> {
        match self {
            TwoTypeOrderOperator::Version1( op ) => op.judge_partial_cmp( left, right ),
            TwoTypeOrderOperator::Version2( op ) => op.judge_partial_cmp( left, right ),
        }
    }
}




impl < OrderOperator1, OrderOperator2, T > 

    JudgeOrder
        < T > for 
        
    TwoTypeOrderOperator
        < OrderOperator1, OrderOperator2 >

where 
    OrderOperator1: JudgeOrder< T >,
    OrderOperator2: JudgeOrder< T >,
{
    fn judge_cmp( & self, left: & T, right: & T ) -> Ordering {
        match self {
            TwoTypeOrderOperator::Version1( op ) => op.judge_cmp( left, right ),
            TwoTypeOrderOperator::Version2( op ) => op.judge_cmp( left, right ),
        }
    }
}















//  Reverse order
//  -------------

/// Represents the reverse of a total order
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::order::{ ReverseOrder, OrderOperatorAuto, JudgePartialOrder };
/// 
/// let mut order_operator_reverse = ReverseOrder::new( OrderOperatorAuto );
/// assert!( order_operator_reverse.judge_lt( &2, &1) )
/// ```
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReverseOrder
            < UnreversedOrderOperator > 
{ unreversed_order_operator: UnreversedOrderOperator }

impl < UnreversedOrderOperator > 

    ReverseOrder
        < UnreversedOrderOperator > 
{ 
    pub fn new( unreversed_order_operator: UnreversedOrderOperator ) ->  ReverseOrder< UnreversedOrderOperator > { ReverseOrder{ unreversed_order_operator } }
}

impl < T, UnreversedOrderOperator: JudgePartialOrder< T > > 

    JudgePartialOrder
        < T > for 

    ReverseOrder
        < UnreversedOrderOperator >

{
    fn judge_partial_cmp( & self, left: & T, right: & T ) -> Option<Ordering> { self.unreversed_order_operator.judge_partial_cmp( right, left ) }
}

impl < T, UnreversedOrderOperator: JudgePartialOrder< T > + JudgeOrder< T > > 

    JudgeOrder
        < T > for 

    ReverseOrder
        < UnreversedOrderOperator >

{
    fn judge_cmp( & self, left: & T, right: & T ) -> Ordering { self.unreversed_order_operator.judge_cmp( right, left ) }
}

/// A convenience trait which allows the user to reverse an order operator via `operator.reverse_order()`
pub trait IntoReverseOrder 
    where Self: Sized
{
    /// Returns an inverted order operator, consuming `self`
    fn into_reverse_order( self ) -> ReverseOrder< Self > { ReverseOrder::new(self) } // for reasons I don't understand, this auto-implementation doesn't seem to work, but the implementation directly below seems to do the trick
}

impl < T: Sized> IntoReverseOrder for T 
{
    /// Returns an inverted order operator, consuming `self`
    fn into_reverse_order( self ) -> ReverseOrder< Self > { ReverseOrder::new(self) }  
}


//  Less-than (auto implementor, any type)
//  ----------------------------

/// Applies the default order to any object that has a default order
/// 
/// That is, you can use this struct to compare the order of `a` and `b`, whenever
/// `a` and `b` implement [PartialOrd]
/// or [Ord]; 
/// 
/// This struct uses zero memory.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize, Ord, PartialOrd)]
pub struct OrderOperatorAuto;

impl OrderOperatorAuto {
    pub fn new() -> OrderOperatorAuto { OrderOperatorAuto }
}

impl<T: PartialOrd> JudgePartialOrder< T > for OrderOperatorAuto {
    
    fn judge_partial_cmp
            ( & self, left: & T, right: & T ) 
            -> 
            Option<Ordering> {
        left.partial_cmp( right )
    }
}

impl<T: Ord> 

    JudgeOrder
        < T > for 
        
    OrderOperatorAuto {} // all methods in this trait are auto-implemented



//  Greater-than (auto implementor)
//  -------------------------------

/// Applies the REVERSE of the default order to any object that has a default order
/// 
/// If `a < b`, then this struct will say that `b < a`.
/// 
/// You can use this struct to compare the order of `a` and `b`, whenever
/// `a` and `b` implement [PartialOrd]
/// or [Ord].
/// 
/// This struct uses zero memory.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OrderOperatorAutoReverse;

impl OrderOperatorAutoReverse {
    pub fn new() -> OrderOperatorAutoReverse { OrderOperatorAutoReverse }
}

impl<T: PartialOrd> JudgePartialOrder< T > for OrderOperatorAutoReverse {
    fn judge_partial_cmp
            ( & self, left: & T, right: & T ) 
            -> 
            Option<Ordering> {
        right.partial_cmp(left)
    }
}

impl<T: Ord> JudgeOrder< T > for OrderOperatorAutoReverse {} // all methods in this trait are auto-implemented



//  Less-than-by-key (auto implementor)
//  -----------------------------------

/// Determines the order of key-value pairs according to their keys
/// 
/// You can use this struct to compare elements of type `T`, whenever 
/// - `T` implements [`KeyValGet<Key, Val>`](`KeyValGet`), and
/// - `Key` implements [PartialOrd] or [Ord]
/// 
/// This struct will indicate that `x < y` whenever `x.key() < y.key()`
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy, new, Serialize, Deserialize)]
pub struct OrderOperatorByKey;

    

// JudgePartialOrder
impl < KeyValuePairStruct >

    JudgePartialOrder 
        < KeyValuePairStruct > for 
    
    OrderOperatorByKey 

    where   KeyValuePairStruct::Key:      PartialOrd,
            KeyValuePairStruct:           KeyValGet,
    
    {

    fn judge_partial_cmp
            ( & self, left: & KeyValuePairStruct, right: & KeyValuePairStruct ) 
            -> 
            Option<Ordering> {
        ( left.key() ).partial_cmp( &right.key() )
    }
}

// JudgeOrder
impl< KeyValuePairStruct >

    JudgeOrder 
        < KeyValuePairStruct > for 
    
    OrderOperatorByKey 

    where   KeyValuePairStruct::Key:      Ord,
            KeyValuePairStruct:           KeyValGet,
    
    {

    fn judge_cmp
            ( & self, left: & KeyValuePairStruct, right: & KeyValuePairStruct ) 
            -> 
            Ordering {
        ( left.key() ).cmp( &right.key() )
    }
}


//  Less-than-by-key (auto implementor for 2-tuples)
//  -----------------------------------


/// Determines the order of key-value pairs of type `(Key,Val)` according to `Key`
/// 
/// You can use this struct to compare elements of type `(Key,Val)`, whenever `Key` implements
/// [PartialOrd] or [Ord].
/// 
/// This struct will indicate that `(a,x) < (b,y)` whenever `a < b`. It is a little less general
/// than `OrderOperatorByKey`, but it may be more efficient. This is because `OrderOperatorByKey`
/// will clone the first entry of a tuple each time it makes a comparison, whereas `OrderOperatorByKeyTuple`
/// will only make a reference.
/// 
/// This struct uses zero memory.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy, new)]
pub struct OrderOperatorByKeyTuple;



// JudgePartialOrder
impl< Key, Val >

    JudgePartialOrder 
        < (Key, Val) > for 
    
    OrderOperatorByKeyTuple //< Key, Val >

    where   
        Key:            PartialOrd,
    {
    fn judge_partial_cmp
            ( & self, left: & (Key, Val), right: & (Key, Val) ) 
            -> 
            Option<Ordering> {
        left.0.partial_cmp( &right.0)
    }
}

// JudgeOrder
impl< Key, Val >

    JudgeOrder 
        < (Key, Val) > for 
    
    OrderOperatorByKeyTuple 

    where Key:          Ord,
    
    {

    fn judge_cmp
            ( & self, left: & (Key, Val), right: & (Key, Val) ) 
            -> 
            Ordering {
        ( left.0 ).cmp( &right.0 )
    }
}


//  Greater-than-by-key (auto implementor)
//  --------------------------------------


/// Determines the order of key-value pairs according to their keys (but reversed!)
/// 
/// You can use this struct to compare elements of type `T`, whenever 
/// - `T` implements [`KeyValGet<Key, Val>`](`KeyValGet`), and
/// - `Key` implements [PartialOrd] or [Ord]
/// 
/// This struct will indicate that `x < y` whenever `x.key() > y.key()`
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy, new)]
pub struct OrderOperatorByKeyReverse;


// JudgePartialOrder
impl< KeyValuePairStruct >

    JudgePartialOrder 
        < KeyValuePairStruct >
    
    for OrderOperatorByKeyReverse 

    where   KeyValuePairStruct::Key:      PartialOrd,
            KeyValuePairStruct:           KeyValGet,
    
    {

    fn judge_partial_cmp
            ( & self, left: & KeyValuePairStruct, right: & KeyValuePairStruct ) 
            -> 
            Option<Ordering> {
        right.key().partial_cmp( &left.key() )
    }
}

// JudgeOrder
impl< KeyValuePairStruct >

    JudgeOrder 
        < KeyValuePairStruct >
    
    for OrderOperatorByKeyReverse 

    where   KeyValuePairStruct::Key:      Ord,
            KeyValuePairStruct:           KeyValGet,
    
    {

    fn judge_cmp
            ( & self, left: & KeyValuePairStruct, right: & KeyValuePairStruct ) 
            -> 
            Ordering {
        ( right.key() ).cmp( &left.key() )
    }
}



//  Less-than-by-key-WITH-custom-comparator (auto implementor)
//  --------------------------------------


/// Determines the order of key-value pairs by applying a custom order operator to keys.
/// 
/// This struct is a wrapper around a customized order operator, which the user can provide.
/// The struct will indicate that `(a,x) < (b,y)` whenever the order operator says that `a < b`.
/// 
/// Any object can be treated as a key-value pair as long as it implements [`KeyValGet<Key, Val>`](`KeyValGet`).
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy, new)]
pub struct OrderOperatorByKeyCustom < OrderOperatorForKeys > 
{
    order_operator_for_keys: OrderOperatorForKeys,    
}


// JudgePartialOrder
impl< KeyValuePairStruct, OrderOperatorForKeys >

    JudgePartialOrder 
        < KeyValuePairStruct > for 
    
    OrderOperatorByKeyCustom 
        < OrderOperatorForKeys > 
        
    where
        KeyValuePairStruct:         KeyValGet,
        OrderOperatorForKeys:       JudgePartialOrder <  KeyValuePairStruct::Key >,
    
    {
    fn judge_partial_cmp
            ( & self, left: & KeyValuePairStruct, right: & KeyValuePairStruct ) 
            -> 
            Option<Ordering> {
        self.order_operator_for_keys.judge_partial_cmp( & left.key(), & right.key())
    }
}

// JudgeOrder
impl< KeyValuePairStruct, OrderOperatorForKeys >

    JudgeOrder 
        < KeyValuePairStruct > for 
    
    OrderOperatorByKeyCustom 
        < OrderOperatorForKeys > 
        
    where
        KeyValuePairStruct:         KeyValGet,
        OrderOperatorForKeys:       JudgeOrder <  KeyValuePairStruct::Key >,
    {}


//  Mutable closure 
//  ---------------------------------------------

use std::marker::PhantomData;

use derive_new::new;
use serde::{Deserialize, Serialize};

use crate::algebra::vectors::entries::KeyValGet;


/// Determines order based on a user-supplied funcion that returns `Option<Ordering>`
/// 
/// Concretely, this is a wrapper for objects that implement `FnMut(&T, &T)->Option<Ordering>`.
/// If `x` is such an object, then we say `a < b` if `x( &a, &b ) = Some(Ordering::Less)`.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq)]
pub struct OrderOperatorByFn
                < F: Fn(&T, &T)->Option<Ordering>, T > 
{ 
    order_operator_unwrapped:         F, 
    phantom_comparatee:               PhantomData< T > 
}


// Constructor
impl < F: Fn(&T, &T)->Option<Ordering>, T > OrderOperatorByFn< F, T > { 
    pub fn new( order_operator_unwrapped: F ) 
            -> 
            OrderOperatorByFn< F, T > 
    { OrderOperatorByFn{ order_operator_unwrapped, phantom_comparatee: PhantomData } } 
}

// Copy
impl <T, F: Fn(&T, &T)->Option<Ordering> > 

    Clone for 
    
    OrderOperatorByFn
        < F, T > 

    where 
        F:  Clone,
{
    fn clone(&self) -> Self { OrderOperatorByFn::new( self.order_operator_unwrapped.clone() ) }
} 

// Copy
impl <T, F: Fn(&T, &T)->Option<Ordering> > 

    Copy for 
    
    OrderOperatorByFn
        < F, T > 

    where 
        F:  Clone + Copy,      
{ }  


// partial order
impl <T, F: Fn(&T, &T)->Option<Ordering> > 

    JudgePartialOrder
        <T> for 
    OrderOperatorByFn
        < F, T > 
{
    // fn lt( & self, a: &T, b: &T) -> bool {
    //     (self.order_operator_unwrapped)(a, b)
    // }
    fn judge_partial_cmp
            ( & self, left: & T, right: & T ) 
            -> 
            Option<Ordering> {
        (self.order_operator_unwrapped)( left, right )
    }
}



//  Mutable closure -- generatoed from a "less than" comparison
//  ---------------------------------------------

/// Determines order based on a user-supplied funcion that returns `bool`
/// 
/// Concretely, this is a wrapper for objects that implement `FnMut(&T, &T)->bool`.
/// If `x` is such an object, then we say `a < b` if `x( &a, &b ) = true`.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq)]
pub struct OrderOperatorByLessThan
                < F: Fn(&T, &T)->bool, T > 
{ 
    order_operator_unwrapped:         F, 
    phantom_comparatee:                 PhantomData< T > 
}

impl < F: Fn(&T, &T)->bool, T > OrderOperatorByLessThan< F, T > { 
    pub fn new( order_operator_unwrapped: F ) 
            -> 
            OrderOperatorByLessThan< F, T > 
    { OrderOperatorByLessThan{ order_operator_unwrapped, phantom_comparatee: PhantomData } } 
}


// Clone
impl <T, F: Fn(&T, &T)->bool > 

    Clone for 
    
    OrderOperatorByLessThan
        < F, T > 

    where 
        F:  Clone,
{
    fn clone(&self) -> Self { OrderOperatorByLessThan::new( self.order_operator_unwrapped.clone() ) }
} 

// Copy
impl <T, F: Fn(&T, &T)->bool > 

    Copy for 
    
    OrderOperatorByLessThan
        < F, T > 

    where 
        F:  Clone + Copy,      
{ }  


// Judge partial order
impl <T, F: Fn(&T, &T)->bool > 

    JudgePartialOrder
        <T> for 

    OrderOperatorByLessThan
        < F, T > 
{
    fn judge_partial_cmp
            ( & self, left: & T, right: & T ) 
            -> 
            Option<Ordering> {
        if (self.order_operator_unwrapped)( left, right ) { return Some(Ordering::Less) }
        if (self.order_operator_unwrapped)( right, left ) { return Some(Ordering::Greater) }
        Some(Ordering::Equal)
    }

    fn judge_lt( &self, left: & T, right: & T ) -> bool { (self.order_operator_unwrapped)( left, right ) }
}

impl <T, F: Fn(&T, &T)->bool > 

    JudgeOrder
        <T> for 

    OrderOperatorByLessThan
        < F, T > 
{
    fn judge_cmp
            ( & self, left: & T, right: & T ) 
            -> 
            Ordering {
        if (self.order_operator_unwrapped)( left, right ) { return Ordering::Less }
        if (self.order_operator_unwrapped)( right, left ) { return Ordering::Greater }
        Ordering::Equal
    }
}







//  Order first by LENGTH, then LEXICOGRAPHIC
//  -------------------------------------------------------------------------------------------------


/// Orders vectors first by length, then by lexicographic order
/// 
/// # Example
/// 
/// ```
/// use oat_rust::utilities::order::{ JudgePartialOrder, JudgeOrder, LexicographicOrderDominatedBylength};
/// use std::cmp::Ordering;
/// 
/// // create an instance of the order operator
/// let order_operator = LexicographicOrderDominatedBylength::new();
/// 
/// // partial order comparisons
/// // ---------------------------------
/// 
/// assert!( 
///     Some( Ordering::Less ) ==
///     order_operator.judge_partial_cmp( 
///         &vec![0],
///         &vec![1]
///     ) 
/// );
/// 
/// assert!( 
///     Some( Ordering::Greater ) ==
///     order_operator.judge_partial_cmp( 
///         &vec![0,1],
///         &vec![1]
///     )
/// );
/// 
/// assert!( 
///     Some( Ordering::Less ) ==
///     order_operator.judge_partial_cmp( 
///         &vec![0],
///         &vec![1,2]
///     )
/// );
/// 
/// assert!( 
///     Some( Ordering::Equal ) ==
///     order_operator.judge_partial_cmp( 
///         &vec![0],
///         &vec![0]
///     )
/// );
/// 
/// // order comparisons
/// // ---------------------------------
/// 
/// assert!( 
///     Ordering::Less ==
///     order_operator.judge_cmp( 
///         &vec![0],
///         &vec![1]
///     )
/// );
/// 
/// assert!( 
///     Ordering::Greater ==
///     order_operator.judge_cmp( 
///         &vec![0,1],
///         &vec![1]
///     )
/// );
/// 
/// assert!( 
///     Ordering::Less ==
///     order_operator.judge_cmp( 
///         &vec![0],
///         &vec![1,2]
///     )
/// );
/// 
/// assert!( 
///     Ordering::Equal ==
///     order_operator.judge_cmp( 
///         &vec![0],
///         &vec![0]
///     )
/// );    
/// ```
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy, new)]
pub struct LexicographicOrderDominatedBylength{}

impl < T >

    JudgePartialOrder< Vec<T> > for 
    
    LexicographicOrderDominatedBylength

    where
        T:  PartialOrd
{
    fn judge_partial_cmp
        ( & self, left: & Vec<T>, right: & Vec<T> ) 
        -> 
        Option<Ordering>
    {
        let left_len = left.len();
        let right_len = right.len();
        match left_len.cmp( & right_len ) {
            Ordering::Equal => { left.partial_cmp(right)  },
            Ordering::Less => { Some( Ordering::Less ) },
            Ordering::Greater => { Some( Ordering::Greater ) }
        }
    }             
}    


impl < T >

    JudgeOrder< Vec<T> > for 
    
    LexicographicOrderDominatedBylength

    where
        T:  Ord
{
    fn judge_cmp
        ( & self, left: & Vec<T>, right: & Vec<T> ) 
        -> 
        Ordering
    {
        let left_len = left.len();
        let right_len = right.len();
        match left_len.cmp( & right_len ) {
            Ordering::Equal => { left.cmp(right)  },
            Ordering::Less => { Ordering::Less },
            Ordering::Greater => { Ordering::Greater }
        }
    }      
}   
















//  Order first by REVERSE LENGTH, then LEXICOGRAPHIC
//  -------------------------------------------------------------------------------------------------


/// Orders vectors first by length, then by lexicographic order
/// 
/// # Example
/// 
/// ```
/// use oat_rust::utilities::order::{ JudgePartialOrder, JudgeOrder, LexicographicOrderDominatedByReverselength};
/// use std::cmp::Ordering;
/// 
/// // create an instance of the order operator
/// let order_operator = LexicographicOrderDominatedByReverselength::new();
/// 
/// // partial order comparisons
/// // ---------------------------------
/// 
/// assert!( 
///     Some( Ordering::Less ) ==
///     order_operator.judge_partial_cmp( 
///         &vec![0],
///         &vec![1]
///     ) 
/// );
/// 
/// assert!( 
///     Some( Ordering::Less ) ==
///     order_operator.judge_partial_cmp( 
///         &vec![0,1],
///         &vec![1]
///     )
/// );
/// 
/// assert!( 
///     Some( Ordering::Greater ) ==
///     order_operator.judge_partial_cmp( 
///         &vec![0],
///         &vec![1,2]
///     )
/// );
/// 
/// assert!( 
///     Some( Ordering::Equal ) ==
///     order_operator.judge_partial_cmp( 
///         &vec![0],
///         &vec![0]
///     )
/// );
/// 
/// // order comparisons
/// // ---------------------------------
/// 
/// assert!( 
///     Ordering::Less ==
///     order_operator.judge_cmp( 
///         &vec![0],
///         &vec![1]
///     )
/// );
/// 
/// assert!( 
///     Ordering::Less ==
///     order_operator.judge_cmp( 
///         &vec![0,1],
///         &vec![1]
///     )
/// );
/// 
/// assert!( 
///     Ordering::Greater ==
///     order_operator.judge_cmp( 
///         &vec![0],
///         &vec![1,2]
///     )
/// );
/// 
/// assert!( 
///     Ordering::Equal ==
///     order_operator.judge_cmp( 
///         &vec![0],
///         &vec![0]
///     )
/// );    
/// ```
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy, new)]
pub struct LexicographicOrderDominatedByReverselength{}

impl < T >

    JudgePartialOrder< Vec<T> > for 
    
    LexicographicOrderDominatedByReverselength

    where
        T:  PartialOrd
{
    fn judge_partial_cmp
        ( & self, left: & Vec<T>, right: & Vec<T> ) 
        -> 
        Option<Ordering>
    {
        let left_len = left.len();
        let right_len = right.len();
        match left_len.cmp( & right_len ) {
            Ordering::Equal => { left.partial_cmp(right)  },
            Ordering::Less => { Some( Ordering::Greater ) },
            Ordering::Greater => { Some( Ordering::Less ) }
        }
    }             
}    


impl < T >

    JudgeOrder< Vec<T> > for 
    
    LexicographicOrderDominatedByReverselength

    where
        T:  Ord
{
    fn judge_cmp
        ( & self, left: & Vec<T>, right: & Vec<T> ) 
        -> 
        Ordering
    {
        let left_len = left.len();
        let right_len = right.len();
        match left_len.cmp( & right_len ) {
            Ordering::Equal => { left.cmp(right)  },
            Ordering::Less => { Ordering::Greater },
            Ordering::Greater => { Ordering::Less }
        }
    }      
}   



























//  TESTS
//  ========================================================



#[cfg(test)]
mod tests {

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn define_order_operator() {     
        use crate::utilities::order::JudgePartialOrder;

        //  Define a new type of struct
        //  ---------------------------        
        pub struct MyOrderOperator;

        //  Implement the JudgePartialOrder trait on our new struct
        //  -------------------------------------------------------
        impl JudgePartialOrder< usize > for MyOrderOperator {

            // this function should return true iff left >= right
            fn judge_ge(&self,left: &usize, right: &usize) -> bool {
                left >= right 
            }
            // this function should return true iff left > right
            fn judge_gt(&self,left: &usize, right: &usize) -> bool {
                left > right
            }
            // this function should return true iff left <= right
            fn judge_le(&self,left: &usize, right: &usize) -> bool {
                left <= right
            }
            // this function should return true iff left < right
            fn judge_lt(&self, left: &usize, right: &usize ) -> bool {
                left < right
            }
            fn judge_partial_cmp( & self, left: & usize, right: & usize )  -> Option<Ordering> {
                Some( left.cmp( right ) )
            }
        }

        //  Implement the JudgeOrder trait on our new struct
        //  ------------------------------------------------
        impl JudgeOrder< usize > for MyOrderOperator {
            // the judge_clamp function in OAT is modeled off of the clamp function here: https://doc.rust-lang.org/std/cmp/trait.Ord.html
            fn judge_clamp(&self, clampee: usize, min: usize, max: usize) -> usize {
                clampee.clamp(min,max) // here we just use Rust's built-in clamp function
            }
            // the judge_cmp function in OAT is modeled off of the clamp function here: https://doc.rust-lang.org/std/cmp/trait.Ord.html            
            fn judge_cmp(&self, left: &usize, right: &usize ) -> Ordering {
                left.cmp(right)
            }
            // should return the maximum of left and right
            fn judge_max(&self, left: usize, right: usize) -> usize {
                left.max(right)
            }
            // should return the minimum of left and right
            fn judge_min(&self, left: usize, right: usize) -> usize {
                left.min(right)
            }
        }

        //  Use the trait
        fn compare_order() {
            let a = 1usize;
            let b = 2usize;
            let order_operator = MyOrderOperator;

            if order_operator.judge_lt( &a, &b ) { println!("a is strictly less than b") }
            else { println!("a is greater than or equal to b")}
        }

    }

}