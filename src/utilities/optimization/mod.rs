//! Mathematical optimization
//! 
//! 
//! This module provides tools for mathematical optimization, which can be used to [tighten cycle representatives](https://www.frontiersin.org/articles/10.3389/frai.2021.681117/full).
//! 
//! 


use crate::algebra::vectors::entries::KeyValGet;
use crate::utilities::optimization::minimize_l1::SolutionL1;
use std::fmt::Debug;
use std::hash::Hash;
use num::ToPrimitive;


pub mod minimize_l1;

/// This conditionally includes a module which uses the Gurobi solver
///  
/// The Gurobi solver is a powerful optimization engine that can be used to solve linear programming problems efficiently.  It requires a license and some setup, but it can significantly speed up optimization tasks.
/// 
#[cfg(feature = "gurobi")]
pub mod gurobi;








#[cfg(feature = "gurobi")]
pub fn minimize_l1_try_gurobi< 
            Key,
            ConstraintMatrix,
            ConstraintVector, 
            BoundVector,
            ColumnIndexIterable,            
            CostVector,
            Coefficient,
        >
    ( 
        mut a:                      ConstraintMatrix,
        b:                          BoundVector,        
        mut c:                      CostVector,
        column_indices:             ColumnIndexIterable,         
        verbose:                    bool, // if true, then print the progress of the optimization
    ) 
    -> 
    Result< SolutionL1<Key, Coefficient>, String >
    where
        ConstraintMatrix:           FnMut( & Key )->ConstraintVector,
        ConstraintVector:           IntoIterator,        
        ConstraintVector::Item:     KeyValGet < Key = Key, Val = Coefficient >,
        BoundVector:                IntoIterator,        
        BoundVector::Item:          KeyValGet < Key = Key, Val = Coefficient >,        
        Coefficient:                Clone + ToPrimitive, // enables casting to f64
        Key:                        Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array 
        ColumnIndexIterable:        Clone + IntoIterator< Item = Key >, 
        CostVector:                 FnMut( & Key) -> f64,  
{
    if verbose{
        println!("Using Gurobi to minimize l1 norm.\n\n");
    }
    
    gurobi::l1_norm::minimize_l1( a, b, c, column_indices, verbose)
        .map_err( |e| format!("Gurobi error: {}", e) )
}



/// Minimize `c' y` subject to `y = |Ax + b|`
/// 
/// This function solves a problem of form
/// 
/// ```text
/// minimize        c' |y|
/// subject to      y = Ax + b
///                 |y| is the entry-wise absolute value of y
///                 y and x are unconstrained continuous variables
/// ```
/// 
/// We solve this problem by solving the following equivalent linear programming problem:
/// 
/// ```text
/// minimize        c' y
/// subject to      y >=  + Ax + b
///                 y >=  - Ax - b 
///                 y and x are unconstrained continuous variables
/// ```
/// 
/// Arguments:
/// 
/// - `a`: function that consumes `k` and returns the `k`th column  of `A`; order of entries does not matter
/// - `b`: the vector `b`, represented as an iterator over the nonzero entries; entries can by any type that implements `KeyValGet< Key = Key, Val = Coefficient >`, where `Coefficient` implements `ToPrimitive`; order does not matter, but $b$ should contain no repeat entries
/// - `c`: the vector `c`, represented as a function that consumes `i` and returns `c[i]`.  We require the coefficients of `c` to be nonnegative.
/// - `column_indices`: iterator that runs over the column indices of `A`; order does not matter
/// - `verbose`: if true, then print the progress of the optimization
///
/// Returns:
/// 
/// An [SolutionL1]
/// 
/// # Examples
/// 
/// ```
/// // Minimize the l1 norm of
/// // | 1 |       | 2 |
/// // | 1 | * x + | 2 |
/// 
/// use oat_rust::utilities::optimization::minimize_l1::minimize_l1;               
/// 
/// let constraint_data = vec![ vec![ (0,1.), (1,1.) ] ];
/// let a = |i: & usize | -> Vec<(usize,f64)> { constraint_data[*i].clone() }; // the matrix 2x1 matrix [1; 1]
/// let b = vec![ (0usize,2f64), (1usize,2f64) ]; // the vector [2, 2]
/// let c = |x: & usize| -> f64 { 1.0 }; // the vector [1, 1]
/// let column_indices = vec![0];  
/// let verbose = false;      
/// 
/// let solution = minimize_l1( a.clone(),b.clone(),c.clone(), column_indices, verbose).unwrap();
/// 
/// println!("{:?}", solution);
/// assert_eq!( solution.x(), &vec![(0,-2.0)] );
/// assert_eq!( solution.b(), &b );
/// assert_eq!( solution.y(), &Vec::<(usize,f64)>::new() );
/// assert_eq!( solution.cost_b(), &4.0 );
/// assert_eq!( solution.cost_y(), &0.0 ); 
/// ```          
///
/// ```
/// /// Minimize the l1 norm of
/// /// ```
/// /// | 1 |       | 1 |
/// /// | 1 | * x + | 1 |
/// /// | 1 |       | 0 |
/// 
/// use oat_rust::utilities::optimization::minimize_l1::minimize_l1;        
/// 
/// let constraint_data = vec![ vec![ (0,1.), (1,1.), (2,1.), ] ];
/// let a = |i: & usize | -> Vec<(usize,f64)> { constraint_data[*i].clone() };
/// let b = vec![ (0,1.), (1,1.), ];
/// let c = |x: & usize| -> f64 { 1.0 };
/// let column_indices = vec![0]; 
/// let verbose = false;       
/// 
/// let solution = minimize_l1( a.clone(),b.clone(),c.clone(), column_indices, verbose).unwrap();
/// 
/// println!("{:?}", solution);
/// assert_eq!( solution.x(), &vec![(0,-1.0)] );
/// assert_eq!( solution.b(), &b );
/// assert_eq!( solution.y(), &vec![(2,-1.0)] );
/// assert_eq!( solution.cost_b(), &2.0 );
/// assert_eq!( solution.cost_y(), &1.0 );    
/// ```    
#[cfg(not(feature = "gurobi"))]
pub fn minimize_l1_try_gurobi< 
            Key,
            ConstraintMatrix,
            ConstraintVector, 
            BoundVector,
            ColumnIndexIterable,            
            CostVector,
            Coefficient,
        >
    ( 
        mut a:                      ConstraintMatrix,
        b:                          BoundVector,        
        mut c:                      CostVector,
        column_indices:             ColumnIndexIterable,         
        verbose:                    bool, // if true, then print the progress of the optimization
    ) 
    -> 
    Result< SolutionL1<Key, Coefficient>, String >
    where
        ConstraintMatrix:           FnMut( & Key )->ConstraintVector,
        ConstraintVector:           IntoIterator,        
        ConstraintVector::Item:     KeyValGet < Key = Key, Val = Coefficient >,
        BoundVector:                IntoIterator,        
        BoundVector::Item:          KeyValGet < Key = Key, Val = Coefficient >,        
        Coefficient:                Clone + ToPrimitive, // enables casting to f64
        Key:                        Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array 
        ColumnIndexIterable:        Clone + IntoIterator< Item = Key >, 
        CostVector:                 FnMut( & Key) -> f64,  
{
    if verbose {
        println!("Gurobi feature is not enabled. Using the default l1 minimization method for Rust goodlp.\n\n");
    }    
    minimize_l1::minimize_l1( a, b, c, column_indices, verbose)
}