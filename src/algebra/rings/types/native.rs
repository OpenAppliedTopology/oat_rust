//! Ring operators for rings that already exist in Rust; packaged in convenient, zero-memory wrappers.
//! 
//! To use a semiring (or ring, or division ring) in OAT, the standard approach is to use a
//! [ring operator](crate::algebra::rings) object that implements the [SemiringOperations], [RingOperations], or [DivisionRingOperations] traits.
//! You can then use the ring operator to
//! perform elementary operations like addition, multiplication, etc.
//! 
//! Rust already has a number of rings "built in."  This module provides a 
//! convenient way to generate a ring operator object for one of these built-in rings, using the [RingOperatorForNativeRustNumberType]
//! struct. 
//! This struct uses zero memory!


use num::rational::Ratio;

use crate::algebra::rings::traits::{SemiringOperations, RingOperations, DivisionRingOperations};
use std::{fmt::Debug, marker::PhantomData};



//----------------------------------------------------------
//  DIVISION RINGS NATIVE TO RUST
//----------------------------------------------------------

/// Zero-memory struct encoding the ring operations for native Rust number types.
///
/// # Examples
/// 
/// 
/// ```
/// use oat_rust::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
/// use oat_rust::algebra::rings::traits::SemiringOperations;
///
/// let ring  =  < RingOperatorForNativeRustNumberType::<usize> >::new();
///
/// assert_eq!( 3, ring.add( 1, 2 ) ); 
/// assert_eq!( 2, ring.multiply( 1, 2 ) );
/// assert_eq!( 0, RingOperatorForNativeRustNumberType::<usize>::zero() );
/// assert_eq!( 1, RingOperatorForNativeRustNumberType::<usize>::one()  );
/// assert!( ! ring.is_0( 1 ) );
/// assert!(   ring.is_0( 0 ) );
/// assert!(   ring.is_1( 1 ) );
/// assert!( ! ring.is_1( 0 ) );
/// ```
/// 
/// 
/// ```
/// use oat_rust::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
/// use oat_rust::algebra::rings::traits::RingOperations;
///
/// let ring = RingOperatorForNativeRustNumberType::<i64>::new();
/// let a : i64 = 1;
/// let b : i64 = 2;
///
/// assert_eq!( -1, ring.subtract( a, b ) );
/// assert_eq!( -1, ring.negate( a ) );
/// ```
/// 
/// ```
/// use oat_rust::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
/// use oat_rust::algebra::rings::traits::DivisionRingOperations;
/// use num::rational::Ratio;
/// 
/// 
/// let ring : RingOperatorForNativeRustNumberType< f64 > =  RingOperatorForNativeRustNumberType::<f64>::new();
/// let a = 2.0 ;
/// let b = 4.0 ;
///
/// assert_eq!( 0.5, ring.divide( a, b ) );
/// assert_eq!( 0.5, ring.invert( a ) );
///
/// ```
///
/// ```
/// use oat_rust::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
/// use oat_rust::algebra::rings::traits::DivisionRingOperations;
/// use num::rational::Ratio;
///
///
/// // The `< .. >` brackets around `RingOperatorForNativeRustNumberType::<Ratio<i64>>`
/// // are used to disambiguate which `new` function should be used
/// // (Rust throws an error if these aren't used)
/// let ring  =     < 
///                     RingOperatorForNativeRustNumberType::<Ratio<i64>> 
///                 >
///                 ::new();
/// let a = Ratio::new( 2, 3 );
/// let b = Ratio::new( 3, 1 );
/// let c = Ratio::new( 2, 9 );
/// let d = Ratio::new( 3, 2 );
///
/// assert_eq!( c, ring.divide( a, b ) );
/// assert_eq!( d, ring.invert( a ) );
/// ```
///
#[derive(Debug, Clone, Copy)]
pub struct RingOperatorForNativeRustNumberType< Element >
{ 
    // This phantom field uses zero memory; it is here only 
    // because rust otherwise complains that `Element` is
    // unused.  See the documentation on `PhantomData` for
    // more details.  
    phantom: PhantomData< Element> 
}

impl    < Element >
        
        RingOperatorForNativeRustNumberType 
            < Element > 
{
    // Generate a [RingOperatorForNativeRustNumberType]
    pub fn new( ) -> Self  
    {
        RingOperatorForNativeRustNumberType { phantom: PhantomData }
    }
}


impl    < Element > 
        
        SemiringOperations for 
        
        RingOperatorForNativeRustNumberType 
            < Element >  
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Mul < Output = Element >  +
                    Clone + 
                    Debug + 
                    PartialEq
{
    /// Element type 
    type Element = Element;

    /// Identity elements
    fn is_0( &self, x: Element ) -> bool { x.is_zero() }
    fn is_1( &self, x: Element ) -> bool { x.is_one() }
    fn zero() -> Element { Element::zero() }
    fn one()  -> Element { Element::one() }

    /// Add
    fn add( &self, x: Element, y: Element ) -> Element { x + y }

    /// Multiply
    fn multiply( &self, x: Element, y: Element ) -> Element { x * y }
}

impl    < Element > 
        
        RingOperations for 
        
        RingOperatorForNativeRustNumberType 
            < Element >  
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    std::ops::Neg  < Output = Element > +
                    Clone + 
                    Debug + 
                    PartialEq
{
    /// Subtract y from x.
    fn subtract( &self, x: Element, y: Element ) -> Element { x - y }

    /// Additive inverse `-x`. 
    fn negate( &self, x: Element ) -> Element { - x }
}

impl    < Element > 

        DivisionRingOperations for 
        
        RingOperatorForNativeRustNumberType 
            < Element >  
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    Clone + 
                    Debug + 
                    PartialEq
{
    /// `x/y` if `y` is nonzero.  
    fn divide( &self, x: Element, y: Element ) -> Element { x / y }

    /// `1/x` if `x` is nonzero.  
    fn invert( &self, x: Element ) -> Element { Element::one() / x }
}


//----------------------------------------------------------
//  TYPE ALIASES
//----------------------------------------------------------

/// A ring operator for the field of real numbers
/// 
/// Operats on elements of type `f64`.
/// 
/// Note: ring operations are performed by floating point arithmetic and do not have infinite precision.
pub type FieldFloat64 = RingOperatorForNativeRustNumberType < f64 >;

/// A ring operator for the field of rational numbers
/// 
/// Operates on elements of type `Ratio<i64>`.
pub type FieldRational64 = RingOperatorForNativeRustNumberType < Ratio<i64> >;

/// A ring operator for the field of rational numbers
/// 
/// Operates on elements of type `Ratio<isize>`.
pub type FieldRationalSize = RingOperatorForNativeRustNumberType < Ratio<isize> >;


/// A ring operator for the ring of integers
/// 
/// Operates on elements of type `i64`.
pub type RingIsize = RingOperatorForNativeRustNumberType < isize >;




#[cfg(test)]
mod tests {


    #[test]
    fn test_misc() {
    }

}



