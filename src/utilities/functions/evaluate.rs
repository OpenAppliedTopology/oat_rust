//! Evaluate functions.


use std::collections::HashMap;
use std::hash::Hash;
use std::marker::PhantomData;

use derive_new::new;
 




/// Evaluate a function
/// 
/// This trait is similar in spirit to `Fn`, `FnOnce` and `FnMut`.  A key problem with `Fn`, `FnOnce` and `FnMut` 
/// is that manual implementation is experimental [E0183](https://doc.rust-lang.org/error-index.html#E0183).  The
/// `EvaluateFunction` trait can act as a stand-in for one of these traits when manual implementation is required.
pub trait EvaluateFunction< Input, Output >{
    fn evaluate_function( &self, input: Input ) -> Output;
}

/// Evaluate a function on a reference
/// 
/// This trait is similar in spirit to `Fn`, `FnOnce` and `FnMut`.  A key problem with `Fn`, `FnOnce` and `FnMut` 
/// is that manual implementation is experimental [E0183](https://doc.rust-lang.org/error-index.html#E0183).  The
/// `EvaluateFunction` trait can act as a stand-in for one of these traits when manual implementation is required.
/// 
/// This is also similar to [EvaluateFunction], but the difference in function signature makes subtle (and pragmatically significant) differences in lifetime handling
pub trait EvaluateFunctionRef< Input, Output >{
    fn evaluate_function_ref( &self, input: &Input ) -> Output;
}


//  IDENTITY FUNCTION
//  -------------------------------------------------------------------------------------------------

/// A struct that implements the identify function via the `EvaluateFunction` trait
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::functions::evaluate::IdentityFunction;
/// use oat_rust::utilities::functions::evaluate::EvaluateFunction;
/// 
/// let identity_function = IdentityFunction;
/// assert_eq!( 1, identity_function.evaluate_function(1) );
/// ```
#[derive(Clone,Copy,Debug,new)]
pub struct IdentityFunction;

impl < Input > 

    EvaluateFunction
        < Input, Input > for 
    
    IdentityFunction
{
    fn evaluate_function( &self, input: Input ) -> Input { input }
}


//  NEGATION
//  -------------------------------------------------------------------------------------------------

/// A struct that implements the "not" function
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::functions::evaluate::{IdentityFunction, LogicalNot};
/// use oat_rust::utilities::functions::evaluate::EvaluateFunction;
/// 
/// let negation = LogicalNot::new(IdentityFunction);
/// 
/// assert!( ! negation.evaluate_function(true)  );
/// assert!(   negation.evaluate_function(false) );
/// ```
#[derive(Clone,Copy,Debug,new)]
pub struct LogicalNot< T >{ predicate: T }

impl < T, Input > 

    EvaluateFunction
        < Input, bool > for 
    
    LogicalNot
        < T > where 
    
    T:  EvaluateFunction< Input, bool >

{
    fn evaluate_function( &self, input: Input ) -> bool { ! self.predicate.evaluate_function(input) }
}

impl < T, Input > 

    EvaluateFunctionRef
        < Input, bool > for 
    
    LogicalNot
        < T > where 
    
    T:  EvaluateFunctionRef< Input, bool >

{
    fn evaluate_function_ref( &self, input: &Input ) -> bool { ! self.predicate.evaluate_function_ref(input) }
}


//  REVERSE VECTOR
//  -------------------------------------------------------------------------------------------------


/// Implements the `EvaluateFunction` trait, for the function that consumes a vector and returns the reversed vector.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::functions::evaluate::ReverseVector;
/// use oat_rust::utilities::functions::evaluate::EvaluateFunction;
/// 
/// let reverser = ReverseVector::new();
/// assert_eq!( vec![2,1], reverser.evaluate_function(vec![1,2]) );
/// ```
#[derive(Clone,Copy,Debug,new)]
pub struct ReverseVector;

impl < T > 

    EvaluateFunction
        < Vec<T>, Vec<T> > for 
    
    ReverseVector
{
    fn evaluate_function( &self, input: Vec<T> ) -> Vec<T> { let mut a = input; a.reverse(); a }
}


//  IMPLEMENT ON WRAPPERS -- OF IMPLEMENTORS OF EvaluateFunction
//  -------------------------------------------------------------------------------------------------


/// `EvaluateFunction` does not auto-implement on references by default; however, you can wrap your reference in this struct.
/// 
/// Concretely, if you have a type `T` that implements `EvaluateFunction<Input,Output>`, OAT will not automatically
/// implement `EvaluateFunction<Input,Output>` for `&'a T`.  However, it *will* automatically implement
/// `EvaluateFunction<Input,Output>` for `ReferencedEvaluator< 'a T >`
#[derive(Clone,Copy,Debug,new)]
pub struct ReferencedEvaluator< 'a, T >{
    pub evaluator:  &'a T,
}

impl< 'a, T, Input, Output >

    EvaluateFunction
        < Input, Output > for 
    
    ReferencedEvaluator
        < 'a, T >

    where
        T:  EvaluateFunction< Input, Output >
{
    fn evaluate_function( &self, input: Input ) -> Output { self.evaluator.evaluate_function(input) }
}




//  IMPLEMENT ON WRAPPERS -- OF FNMUT
//  -------------------------------------------------------------------------------------------------


/// Wrapper for objects that implement `Fn( Input ) -> Output`; the wrapper implements `EvaluateFunction< Input, Output >`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionFnMutWrapper};
/// 
/// let fun = | x: usize | -> usize { 2 * x };
/// let wrapped_fun = EvaluateFunctionFnMutWrapper::new( fun );
/// assert_eq!( wrapped_fun.evaluate_function(1), 2 );
/// ```
/// 
/// # Design notes
/// 
/// We found that implementing `EvaluateFunction` on this wrapper helped the compiler to produce helpful error messages when code breaks;
/// however this took place at early stages of development.  It may be worth reviewing whether these wrappers are truly necessary.
#[derive(Debug)]
pub struct EvaluateFunctionFnMutWrapper< Function > { unwrapped_function: Function }

impl < Function > EvaluateFunctionFnMutWrapper< Function > { 
    /// Create new function wrapper.
    pub fn new( unwrapped_function: Function ) -> EvaluateFunctionFnMutWrapper< Function > { EvaluateFunctionFnMutWrapper {unwrapped_function} } 
}

impl < Function, Input, Output >

    EvaluateFunction
        < Input, Output > for

    EvaluateFunctionFnMutWrapper
        < Function >

    where 
        Function:   Fn( Input ) -> Output
{
    fn evaluate_function( & self, input: Input ) -> Output {
        (self.unwrapped_function)( input )
    }
}       

impl < Function > 

    Clone for

    EvaluateFunctionFnMutWrapper
        < Function >   
        
    where
        Function:   Clone
{
    fn clone( &self ) -> EvaluateFunctionFnMutWrapper<Function>
    { EvaluateFunctionFnMutWrapper::new( self.unwrapped_function.clone() ) }
}


//  IMPLEMENT ON WRAPPERS -- OF IMPLEMENTORS
//  -------------------------------------------------------------------------------------------------


/// Wrapper for types that implement the [`EvaluateFunction`] trait in multiple different ways.
/// 
/// Some times, such as `HashMap` implement the [`EvaluateFunction`] trait in multiple different ways.  
/// This can lead to errors where the compiler is uncertain which function the object represents.
/// By wrapping the object in an `EvaluateFunctionDisambiguator` struct, we produce an object that
/// only implements one instance of the [`EvaluateFunction`] trait, thus removing any ambiguity.
#[derive(Clone,Copy,Debug)]
pub struct EvaluateFunctionDisambiguator< Function: EvaluateFunction< Input, Output >, Input, Output > 
{ 
    unwrapped_function: Function,
    phantom_input:      PhantomData< Input >,
    phantom_output:     PhantomData< Output >,
}

impl < Input, Output, Function: EvaluateFunction< Input, Output > > 

    EvaluateFunctionDisambiguator
        < Function, Input, Output, >
{
    pub fn new( unwrapped_function: Function ) -> Self { EvaluateFunctionDisambiguator{ unwrapped_function, phantom_input: PhantomData, phantom_output: PhantomData } }
}        

impl < Input, Output, Function: EvaluateFunction< Input, Output > > 

    EvaluateFunction
        < Input, Output > for

    EvaluateFunctionDisambiguator
        < Function, Input, Output, >
{
    fn evaluate_function( &self, input: Input ) -> Output { self.unwrapped_function.evaluate_function(input) }
}    


//  WRAPPER ANALOGOUS TO RUST'S "MAP" METHOD FOR OPTIONS
//  -------------------------------------------------------------------------------------------------


/// Maps `Option< Input >` to `Option< Output >`, analogous to Rusts's `map` function for options
/// 
/// # Examples
/// 
/// You can combine this creatively with objects that implement multiple instances of `EvaluateFunction`, to get different results.
/// 
/// ```
/// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, Map};
/// use assert_panic::assert_panic;
/// 
/// // Map usize -> Option< usize >
/// let f: Vec<usize> = vec![7,6];
/// let a: Option<usize> = f.evaluate_function(0usize); assert_eq!( a, Some(7usize) );
/// let a: Option<usize> = f.evaluate_function(4usize); assert_eq!( a, None         );
/// 
/// // Map Option< usize > -> Option< usize >
/// let g = Map{ mapping_rule: f };
/// let a: Option<usize> = g.evaluate_function( Some(0usize) ); assert_eq!( a, Some(7usize) );
/// let a: Option<usize> = g.evaluate_function( None::<usize>); assert_eq!( a, None         ); 
///
/// // The following will panic, `Map` calls its inner function on a value to return a value, and a vector of length 2 has no way to map the integer 4 to another integer
/// // let b: Option<usize> = g.evaluate_function(4usize); assert_eq!( b, None         ); 
/// 
/// // Map Option< usize > ->  Option< Option< usize > > with wrapper
/// let a: Option<Option<usize>> = g.evaluate_function( Some(0usize) ); assert_eq!( a, Some( Some(7usize) ) );
/// let a: Option<Option<usize>> = g.evaluate_function( Some(4usize) ); assert_eq!( a, Some( None         ) );
/// let a: Option<Option<usize>> = g.evaluate_function( None::<usize>); assert_eq!( a, None                 );
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Map< MappingRule >{ 
    pub mapping_rule: MappingRule 
}

impl < Input, Output, MappingRule > 

    EvaluateFunction
        < Option< Input >, Option< Output > > for 
    
    Map
        < MappingRule > where 

    MappingRule:    EvaluateFunction< Input, Output >
{
    fn evaluate_function( &self, input: Option< Input > ) -> Option< Output > { 
        input.map(|x| self.mapping_rule.evaluate_function(x)) 
    }
}


//  WRAPPER ANALOGOUS TO RUST'S "MAP" METHOD FOR OPTIONS
//  -------------------------------------------------------------------------------------------------


/// Maps `Option< Input >` to `Option< Output >`, analogous to Rusts's `and_then` method function for options
/// 
/// Note: the distinction between `map` versus `and_then`, is that the closure used in `map` takes a value as input, whereas the closure in `and_then` takes an option as input.
/// 
/// # Examples
/// 
/// # Examples
/// 
/// You can combine this creatively with objects that implement multiple instances of `EvaluateFunction`, to get different results.
/// 
/// ```
/// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, AndThen};
/// use assert_panic::assert_panic;
/// 
/// // Map usize -> Option< usize >
/// let f: Vec<usize> = vec![7,6];
/// let a: Option<usize> = f.evaluate_function(0usize); assert_eq!( a, Some(7usize) );
/// let a: Option<usize> = f.evaluate_function(4usize); assert_eq!( a, None         );
/// 
/// // Map Option< usize > -> Option< usize >
/// let g = AndThen{ mapping_rule: f };
/// let a: Option<usize> = g.evaluate_function( Some(0usize) ); assert_eq!( a, Some(7usize) );
/// let a: Option<usize> = g.evaluate_function( Some(4usize) ); assert_eq!( a, None         );
/// let a: Option<usize> = g.evaluate_function( None::<usize>); assert_eq!( a, None         ); 
/// ```
#[derive(Debug, Clone, Copy)]
pub struct AndThen< MappingRule >{ 
    pub mapping_rule: MappingRule 
}

impl < Input, Output, MappingRule > 

    EvaluateFunction
        < Option< Input >, Option< Output > > for 
    
    AndThen
        < MappingRule > where 

    MappingRule:    EvaluateFunction< Input, Option< Output > >
{
    fn evaluate_function( &self, input: Option< Input > ) -> Option< Output > { 
        input.and_then(|x| self.mapping_rule.evaluate_function(x)) 
    }
}



//  IMPLEMENT ON VECTORS
//  -------------------------------------------------------------------------------------------------


//  UNWRAPPED VALUES
//  --------------------


impl < Output: Clone > 

    EvaluateFunction
        < usize, Output > for

    Vec< Output >
{
    fn evaluate_function( & self, input: usize ) -> Output {
        self[ input ].clone()
    }    
}



impl < 'a, Output: Clone > 

    EvaluateFunction
        < usize, Output > for

    &'a Vec< Output >
{
    fn evaluate_function( & self, input: usize ) -> Output {
        self[ input ].clone()
    }    
}



impl < 'a, Output > 

    EvaluateFunction
        < usize, &'a Output > for

    &'a Vec< Output >
{
    fn evaluate_function( & self, input: usize ) -> &'a Output {
        & self[ input ]
    }    
}


//  WRAPPED IN OPTIONS
//  --------------------


impl < Output: Clone > 

    EvaluateFunction
        < usize, Option< Output > > for

    Vec< Output >
{
    fn evaluate_function( & self, input: usize ) -> Option< Output > {
        match input < self.len() {
            true    =>  { Some( self[ input ].clone() ) },
            false   =>  { None }
        }        
    }    
}



impl < 'a, Output: Clone > 

    EvaluateFunction
        < usize, Option< Output > > for

    &'a Vec< Output >
{
    /// Return `self[input].clone()`.
    fn evaluate_function( & self, input: usize ) -> Option< Output > {
        match input < self.len() {
            true    =>  { Some( self[ input ].clone() ) },
            false   =>  { None }
        }        
    }    
}



impl < 'a, Output > 

    EvaluateFunction
        < usize, Option< &'a Output > > for

    &'a Vec< Output >
{
    /// Return `Some( &self[ input ] )` if `input < self.len()`, and `None` otherwise.
    fn evaluate_function( & self, input: usize ) -> Option< &'a Output > {
        match input < self.len() {
            true    =>  { Some( & self[ input ] ) },
            false   =>  { None }
        }        
    }     
}




//  REFERENCED INPUTS
//----------------------

//  SUBGRUOP: UNWRAPPED VALUES


impl < 'b, Output: Clone > 

    EvaluateFunction
        < &'b usize, Output > for

    Vec< Output >
{
    fn evaluate_function( & self, input: &'b usize ) -> Output {
        self[ *input ].clone()
    }    
}



impl < 'a, 'b, Output: Clone > 

    EvaluateFunction
        < &'b usize, Output > for

    &'a Vec< Output >
{
    fn evaluate_function( & self, input: &'b usize ) -> Output {
        self[ *input ].clone()
    }    
}



impl < 'a, 'b, Output > 

    EvaluateFunction
        < &'b usize, &'a Output > for

    &'a Vec< Output >
{
    fn evaluate_function( & self, input: &'b usize ) -> &'a Output {
        & self[ *input ]
    }    
}


//  SUBGRUOP: WRAPPED VALUES


impl < 'b, Output: Clone > 

    EvaluateFunction
        < &'b usize, Option< Output > > for

    Vec< Output >
{
    fn evaluate_function( & self, input: &'b usize ) -> Option< Output > {
        match *input < self.len() {
            true    =>  { Some( self[ *input ].clone() ) },
            false   =>  { None }
        }        
    }    
}



impl < 'a, 'b, Output: Clone > 

    EvaluateFunction
        < &'b usize, Option< Output > > for

    &'a Vec< Output >
{
    /// Return `self[input].clone()`.
    fn evaluate_function( & self, input: &'b usize ) -> Option< Output > {
        match *input < self.len() {
            true    =>  { Some( self[ *input ].clone() ) },
            false   =>  { None }
        }        
    }    
}



impl < 'a, 'b, Output > 

    EvaluateFunction
        < &'b usize, Option< &'a Output > > for

    &'a Vec< Output >
{
    /// Return `Some( &self[ input ] )` if `input < self.len()`, and `None` otherwise.
    fn evaluate_function( & self, input: &'b usize ) -> Option< &'a Output > {
        match *input < self.len() {
            true    =>  { Some( & self[ *input ] ) },
            false   =>  { None }
        }        
    }     
}






//  IMPLEMENT ON HASHMAPS
//  -------------------------------------------------------------------------------------------------


// 0 - Implement the trait for < Input, Output > on HashMap
impl < Input, Output > 

    EvaluateFunction
        < Input, Output > for 
        
    HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator::< HashMap< usize, i32 >, usize, i32 >
    ///     = EvaluateFunctionDisambiguator::new( hash );     
    /// 
    /// assert_eq!( disambiguated.evaluate_function( 1usize ), -1  );  
    /// ```
    /// 
    /// # Design note
    /// 
    /// This implementation clones the output of the hashmap by default; to do otherwise 
    /// may ultimately require the involvement of lifetime parameters, and a more complicated implementation.
    fn evaluate_function( & self, input: Input ) -> Output {
        return self.get( &input ).unwrap().clone()
    }
}

// 1 - Implement the trait for < &Input, Output > on HashMap
impl < Input, Output > 

    EvaluateFunction
        < &Input, Option< Output > > for 
        
    HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::EvaluateFunction;
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// assert_eq!( hash.evaluate_function( &1 ), Some( -1 )  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// This implementation clones the output of the hashmap by default; to do otherwise 
    /// may ultimately require the involvement of lifetime parameters, and a more complicated implementation.
    fn evaluate_function( & self, input: &Input ) -> Option< Output > {
        return self.get( input ).cloned()
    }
}



// 2 - Implement the trait for < Input, Output > on HashMap
impl < Input, Output > 

    EvaluateFunction
        < Input, Option< Output > > for 
        
    HashMap
        < Input, Output > 
    
    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;    
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator::< HashMap< usize, i32 >, usize, Option< i32 > >
    ///     = EvaluateFunctionDisambiguator::new( hash );
    ///
    /// assert_eq!( disambiguated.evaluate_function( 1usize ), Some( -1i32 )  );
    /// ```
    fn evaluate_function( & self, input: Input ) -> Option< Output > {
        return self.get( &input ).cloned()
    }
}



// 3 - Implement the trait for < &Input, Option< &'a Output > > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < &Input, Option< &'a Output > > for 
        
    &'a HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, &usize, Option<&i32>, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );
    /// 
    /// assert_eq!( disambiguated.evaluate_function( &1 ), Some( & -1 )  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// It would be difficult to implement a version of `EvaluateFunction<_, _>` where `evaulate_function` returns a reference
    /// on `HashMap` directly, because a lifetime parameter is likely needed on the output.  Hence we implement `EvaluateFunction<_, _>`
    /// on `& HashMap`.
    fn evaluate_function( & self, input: &Input ) -> Option< &'a Output > {
        return self.get( input )
    }
}



// 3.1 - Implement the trait for < &Input, Option< &'a Output > > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < &Input, &'a Output > for 
        
    &'a HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, &usize, &i32, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );
    /// 
    /// assert_eq!( disambiguated.evaluate_function( &1 ), & -1  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// It would be difficult to implement a version of `EvaluateFunction<_, _>` where `evaulate_function` returns a reference
    /// on `HashMap` directly, because a lifetime parameter is likely needed on the output.  Hence we implement `EvaluateFunction<_, _>`
    /// on `& HashMap`.
    fn evaluate_function( & self, input: &Input ) -> &'a Output {
        return self.get( input ).unwrap()
    }
}



// 3.2 - Implement the trait for < &Input, Option< &'a Output > > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < Input, &'a Output > for 
        
    &'a HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, &i32, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );
    /// 
    /// assert_eq!( disambiguated.evaluate_function( 1 ), & -1  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// It would be difficult to implement a version of `EvaluateFunction<_, _>` where `evaulate_function` returns a reference
    /// on `HashMap` directly, because a lifetime parameter is likely needed on the output.  Hence we implement `EvaluateFunction<_, _>`
    /// on `& HashMap`.
    fn evaluate_function( & self, input: Input ) -> &'a Output {
        return self.get( &input ).unwrap()
    }
}






// 4 - Implement the trait for < Input, &'a Output > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < Input, Option< &'a Output > > for 
        
    &'a HashMap
        < Input, Output > 
    
    where
        Input:  std::cmp::Eq + Hash,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, Option<&i32>, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );        
    /// 
    /// assert_eq!( disambiguated.evaluate_function( 1 ), Some( &-1 )  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// It would be difficult to implement a version of `EvaluateFunction<_, _>` where `evaulate_function` returns a reference
    /// on `HashMap` directly, because a lifetime parameter is likely needed on the output.  Hence we implement `EvaluateFunction<_, _>`
    /// on `& HashMap`.
    fn evaluate_function( & self, input: Input ) -> Option< &'a Output > {
        return self.get( &input )
    }
}



// 5 - Implement the trait for < Input, &'a Output > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < Input, Option< Output > > for 
        
    &'a HashMap
        < Input, Output > 
    
    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;                
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, Option<i32>, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );  
    /// 
    /// assert_eq!( disambiguated.evaluate_function(1), Some( -1 )  );
    /// ```
    fn evaluate_function( & self, input: Input ) -> Option< Output > {
        return self.get( &input ).cloned()
    }
}



// 6 - Implement the trait for < Input, Output > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < Input, Output > for 
        
    &'a HashMap
        < Input, Output > 
    
    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;        
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, i32, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );      
    /// 
    /// assert_eq!( disambiguated.evaluate_function( 1 ), -1  );
    /// ```
    fn evaluate_function( & self, input: Input ) -> Output {
        return self.get( &input ).cloned().unwrap()
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

    // 0
    // EvaluateFunction
    //     < Input, Output > for 
        
    // HashMap
    //     < Input, Output > 
    #[test]
    fn test_0() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        
        // we must sometimes disambiguate the types of the input and output        
        let disambiguated: EvaluateFunctionDisambiguator::< HashMap< usize, i32 >, usize, i32 >
            = EvaluateFunctionDisambiguator::new( hash );     
        
        assert_eq!( disambiguated.evaluate_function( 1usize ), -1  );   
    } 
    

    // 1
    // EvaluateFunction
    //     < &Input, Option< Output > > for 
        
    // HashMap
    //     < Input, Output >         
    #[test]
    fn test_1() {
        use crate::utilities::functions::evaluate::EvaluateFunction;
        use std::collections::HashMap;
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        
        assert_eq!( hash.evaluate_function( &1 ), Some( -1 )  );
    }
        
    // 2
    // EvaluateFunction
    //     < Input, Option< Output > > for 
        
    // HashMap
    //     < Input, Output > 
    #[test]
    fn test_2() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;    

        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 

        // we must sometimes disambiguate the types of the input and output        
        let disambiguated: EvaluateFunctionDisambiguator::< HashMap< usize, i32 >, usize, Option< i32 > >
            = EvaluateFunctionDisambiguator::new( hash );
        
        assert_eq!( disambiguated.evaluate_function( 1usize ), Some( -1i32 )  );
    }
        

    // 3
    // EvaluateFunction
    //     < &Input, Option< &'a Output > > for 
        
    // &'a HashMap
    //     < Input, Output > 
    #[test]
    fn test_3() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;

        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 

        // we must sometimes disambiguate the types of the input and output        
        let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, &usize, Option<&i32>, > 
            = EvaluateFunctionDisambiguator::new( &hash );
        
        assert_eq!( disambiguated.evaluate_function( &1 ), Some( & -1 )  );
    }


    // 4
    // EvaluateFunction
    //     < Input, Option< &'a Output > > for 
        
    // &'a HashMap
    //     < Input, Output > 
    #[test]
    fn test_4() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;

        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 

        // we must sometimes disambiguate the types of the input and output
        let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, Option<&i32>, > 
            = EvaluateFunctionDisambiguator::new( &hash );        
        
        assert_eq!( disambiguated.evaluate_function( 1 ), Some( &-1 )  );
    }
        

    // 5
    // EvaluateFunction
    //     < Input, Option< Output > > for 
        
    // &'a HashMap
    //     < Input, Output > 
    #[test]
    fn test_5() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;                
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 

        // we must sometimes disambiguate the types of the input and output
        let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, Option<i32>, > 
            = EvaluateFunctionDisambiguator::new( &hash );  
        
        assert_eq!( disambiguated.evaluate_function(1), Some( -1 )  );
    }
            
    // 6
    // EvaluateFunction
    //     < Input, Output > for 
        
    // &'a HashMap
    //     < Input, Output > 
    #[test]
    fn test_6() {        

        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;        

        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        
        // we must sometimes disambiguate the types of the input and output
        let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, i32, > 
            = EvaluateFunctionDisambiguator::new( &hash );      

        assert_eq!( disambiguated.evaluate_function( 1 ), -1  );
    
    }    

    #[test]
    fn test_7() {     
        use std::sync::Arc;
        use crate::utilities::functions::evaluate::EvaluateFunction;

        let a = Arc::new( vec![0,1,2] );
        let b: usize = a.evaluate_function( 0 );
    }
}