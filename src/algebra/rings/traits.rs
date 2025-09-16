//! Traits that define what an operator on a semiring, ring, or division ring can do.


//  ---------------------------------------------------------------------------
//  DESIGN NOTES
//  ---------------------------------------------------------------------------

//  * Advantage of this nested structure: makes it straightforward to define matrix multipication
//  over semirings.
//
//  * Reason for deprecating the function "field name" that tells you the underlying mathematical
//  field: 
//  in general, you always know what struct you're working with; so it suffices to describe the
//  mathematical object underlying the struct in the struct's documentation

use std::fmt::Debug;

use auto_impl::auto_impl;





//  ---------------------------------------------------------------------------
//  THE SEMIRING TRAIT
//  ---------------------------------------------------------------------------

/// Arithmetic operations for semirings.
/// 
/// The associated type `Element` is the type of the elements of the semiring. 
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use crate::oat_rust::algebra::rings::traits::SemiringOperations;
/// 
/// // Define a struct that can perform the operations of a semiring.
/// let ring_operator = PrimeOrderField::new(3); 
/// 
/// // This struct implements the `SemiringOperations` trait, so we can use the methods of that trait.
/// // Here is an example of using the method for multiplication:
/// let a = ring_operator.multiply( 1, 2 ); // multiply 1 and 2
/// assert_eq!( a, 2 ); // check that the answer is correct
/// ```
#[auto_impl(&)] // auto-implement this trait on references to objects that implement the trait
pub trait SemiringOperations {

    type Element:     Clone + PartialEq + Debug;


    // IDENTITY ELEMENTS

    /// Return the additive identity.
    fn is_0( &self, x : Self::Element ) -> bool;

    /// Return the multiplicative identity.
    fn is_1( &self, x : Self::Element ) -> bool;

    /// Return the additive identity.
    fn zero() -> Self::Element;

    /// Return the multiplicative identity.
    fn one() -> Self::Element;


    // OPERATIONS 
    
    // DESIGN NOTE: if we changed these functions to take
    // non-references as input, we would have to think about
    // wether to require ring elements to implement the
    // copy trait.

    /// Add
    fn add( &self, x : Self::Element, y : Self::Element ) -> Self::Element;

    /// Multiply
    fn multiply( &self, x : Self::Element, y: Self::Element ) -> Self::Element;

}


//  ---------------------------------------------------------------------------
//  THE RING TRAIT
//  ---------------------------------------------------------------------------


/// Arithmetic operations for **unital** rings.
pub trait RingOperations: SemiringOperations {

    /// Subtract y from x.
    fn subtract( &self, x : Self::Element, y: Self::Element ) -> Self::Element;

    /// Reverse the sign of x.
    fn negate( &self, x : Self::Element ) -> Self::Element;

    /// Return `-1`
    fn minus_one( &self ) -> Self::Element {
        self.negate( Self::one() )
    }    

    /// Return `(-1)^k`
    fn minus_one_to_power( &self, k: usize ) -> Self::Element {
        if k.is_even() { Self::one() }
        else { self.minus_one() }
    }    

}


//----------------------------------------------------------
//  THE DIVISION RING TRAIT 
//----------------------------------------------------------

/// Arithmetic operations for division rings.
pub trait DivisionRingOperations: RingOperations {
    
    /// Divide 
    fn divide( &self, x : Self::Element, y: Self::Element ) -> Self::Element;

    /// Invert 
    fn invert( &self, x : Self::Element ) -> Self::Element;

}







//----------------------------------------------------------
//  TRAITS FO CONVENIENCE
//----------------------------------------------------------

use num::integer::Integer;

pub trait MinusOne: RingOperations {
    fn minus_one( &self ) -> Self::Element;
}

impl    < RingOperator >

        MinusOne for 
        
        RingOperator
        
        where   
            RingOperator: RingOperations
{
    fn minus_one( &self ) -> Self::Element {
        self.negate( RingOperator::one() )
    }
}  

pub trait MinusOneToPower: RingOperations {
    fn minus_one_to_power( &self, k: usize ) -> Self::Element;
}

impl    < RingOperator > 

        MinusOneToPower for 
        
        RingOperator
        
        where   RingOperator: RingOperations
{
    fn minus_one_to_power( &self, k: usize ) -> Self::Element {
        if k.is_even() { RingOperator::one() }
        else { self.minus_one() }
    }
}    