//! Minimize `c' y` subject to `y = |Ax + b|`
//! 
//! 
//! 
//! This module provides a [function](minimize_l1) to solve problems of form 
//! 
//! ```
//! minimize        c' |y|
//! subject to      y = Ax + b
//! where           y and x are unconstrained continuous variables
//!                 |y| is the entry-wise absolute value of y
//! ```
//! 
//! We solve this problem by solving the following equivalent linear programming problem:
//! 
//! ```text
//! minimize        c' y
//! subject to      y >=  + Ax + b
//!                 y >=  - Ax - b 
//!                 y and x are unconstrained continuous variables
//! ```

use gag::Gag;
use grb::attribute::ModelIntAttr::LicenseExpiration;
use grb::constr::IneqExpr;
use grb::prelude::*;
use grb::expr::LinExpr;
use itertools::Itertools;
use num::ToPrimitive;
use num::rational::Ratio;
use ordered_float::OrderedFloat;

use crate::utilities::iterators::merge::hit::HitMergeByPredicateTrait;
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet, KeyValNew};
use crate::algebra::matrices::query::MatrixOracle;
use crate::algebra::matrices::operations::umatch::row_major::Umatch;
use crate::algebra::rings::traits::{SemiringOperations, RingOperations, DivisionRingOperations};
use crate::algebra::rings::types::native::{RingOperatorForNativeRustNumberType, FieldFloat64};
use crate::utilities::order::{JudgePartialOrder, OrderOperatorByKeyCustom, IntoReverseOrder};
use crate::utilities::optimization::minimize_l1::SolutionL1;

use crate::algebra::matrices::operations::umatch::differential::DifferentialComb;
use crate::algebra::vectors::operations::VectorOperations;

use std::fmt::Debug;
use std::hash::Hash;
use std::collections::HashMap;

use derive_getters::{Getters, Dissolve};
use derive_new::new;

type RationalRing = RingOperatorForNativeRustNumberType< Ratio< isize > >;



/// Minimize `c' y` subject to `y = |Ax + b|`
/// 
/// This function solves a problem of form
/// 
/// ```
/// minimize        c' |y|
/// subject to      y = Ax + b
/// where           y and x are unconstrained continuous variables
///                 |y| is the entry-wise absolute value of y
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
/// - `a`: function that consumes `k` and returns the `k`th column  of `A`; a single column should not contain to different entries with the same row indes, but order of entries does not
/// - `b`: the vector `b`, represented as an iterator over the nonzero entries; entries can by any type that implements `KeyValGet< Key = Key, Val = Coefficient >`, where `Coefficient` implements `ToPrimitive`
/// - `c`: the vector `c`, represented as a function that consumes `i` and returns `c[i]`.  We require the coefficients of `c` to be nonnegative.
/// - `column_indices`: iterator that runs over the column indices of `A`; order does not matter
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
/// use oat_rust::utilities::optimization::l1::minimize_l1;               
/// 
/// let constraint_data = vec![ vec![ (0,1.), (1,1.) ] ];
/// let a = |i: & usize | -> Vec<(usize,f64)> { constraint_data[*i].clone() }; // the matrix 2x1 matrix [1; 1]
/// let b = vec![ (0usize,2f64), (1usize,2f64) ]; // the vector [2, 2]
/// let c = |x: & usize| -> f64 { 1.0 }; // the vector [1, 1]
/// let column_indices = vec![0];        
/// 
/// let solution = minimize_l1( a.clone(),b.clone(),c.clone(), column_indices).unwrap();
/// 
/// println!("{:?}", solution);
/// assert_eq!( solution.x(), &vec![(0,-2.0)] );
/// assert_eq!( solution.b(), &b );
/// assert_eq!( solution.y(), &Vec::<(usize,f64)>::new() );
/// assert_eq!( solution.cost_b(), &4.0 );
/// assert_eq!( solution.cost_axb(), &0.0 ); 
/// ```          
///
/// ```
/// /// Minimize the l1 norm of
/// /// ```
/// /// | 1 |       | 1 |
/// /// | 1 | * x + | 1 |
/// /// | 1 |       | 0 |
/// 
/// use oat_rust::utilities::optimization::l1::minimize_l1;        
/// 
/// let constraint_data = vec![ vec![ (0,1.), (1,1.), (2,1.), ] ];
/// let a = |i: & usize | -> Vec<(usize,f64)> { constraint_data[*i].clone() };
/// let b = vec![ (0,1.), (1,1.), ];
/// let c = |x: & usize| -> f64 { 1.0 };
/// let column_indices = vec![0];        
/// 
/// let solution = minimize_l1( a.clone(),b.clone(),c.clone(), column_indices).unwrap();
/// 
/// println!("{:?}", solution);
/// assert_eq!( solution.x(), &vec![(0,-1.0)] );
/// assert_eq!( solution.b(), &b );
/// assert_eq!( solution.y(), &vec![(2,-1.0)] );
/// assert_eq!( solution.cost_b(), &2.0 );
/// assert_eq!( solution.cost_axb(), &1.0 );    
/// ```    
pub fn minimize_l1< 
            Key,
            ConstraintMatrix,
            ConstraintVector, 
            ConstantVector,
            ConstantVectorEntry,
            ColumnIndexIterable,            
            CostVector,
            Coefficient,
        >
    ( 
        mut a:                      ConstraintMatrix,
        b:                          ConstantVector,        
        mut c:                      CostVector,
        column_indices:             ColumnIndexIterable,   
        verbose:                    bool,      
    ) 
    -> 
    grb::Result< SolutionL1<Key, Coefficient> >
    where
        ConstraintMatrix:           FnMut( & Key )->ConstraintVector,
        ConstraintVector:           IntoIterator,
        ConstraintVector::Item:     KeyValGet< Key = Key, Val = Coefficient >,
        ConstantVector:             IntoIterator< Item = ConstantVectorEntry >,
        ConstantVectorEntry:        KeyValGet< Key = Key, Val = Coefficient >,
        Coefficient:                Clone + ToPrimitive, // enables casting to f64
        Key:                        Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   
        ColumnIndexIterable:        Clone + IntoIterator< Item = Key >, 
        CostVector:                 FnMut( & Key) -> f64,                      
{


    let start = std::time::Instant::now();



    //  initialize the model
    let env:            grb::Env;
    let mut model:      grb::Model;


    if verbose {
        env                         =   Env::new("mip.log")?; 
        model                       =   Model::new("model1")?;        
        model.set_param(grb::param::OutputFlag, 1)?;        
    } else {
        let print_gag = Gag::stdout().unwrap();  // this will suppress the output of Gurobi from printing
        env                         =   Env::new("mip.log")?;     
        model                       =   Model::new("model1")?;        
        model.set_param(grb::param::OutputFlag, 0)?;        
    }


    //  initialize lhs, which is a hashmap composed of LinExpr's that collectively represent the vector b
    let mut rowkey_to_rownum                =   HashMap::<Key,usize>::new();
    
    let mut rownum_to_str: Vec<&str> = Vec::new();
    let mut rownum_to_axz: Vec<LinExpr>            = Vec::new();
    let mut rownum_to_y                     =   Vec::new();
    let mut rownum_to_rowkey                =   Vec::new();

    let mut cost_b                          =   0.0;
    let mut b                               =   b.into_iter()
                                                    .map( |entry| (entry.key(), entry.val() ) )
                                                    .collect_vec();

    // add the coefficients of `b` as constants in the constraints
    for entry in b.iter().cloned() {

        // make a new linear expression
        let mut linexpr         =   LinExpr::new();

        // calculate the objective value of the original cycle (incrementally)        
        let coeff               =   entry.val().to_f64().unwrap();
        linexpr.add_constant( coeff ); // note we don't put a variable here  
        rownum_to_axz.push( linexpr ); // push the new constraint to our vector of constraints        

        // make a record of the cost of b (incrementaly)
        cost_b += ( coeff * c( &entry.key() ) ).abs();

        // introduce a y variable for this constriant
        let variable_name=   format!( "y{}", rownum_to_y.len() );
        let variable        =   add_ctsvar!(model, name: &variable_name, bounds: ..)?;      
        rownum_to_y.push( variable );
        rownum_to_rowkey.push(   entry.key() );  
        rowkey_to_rownum.insert( entry.key(), rowkey_to_rownum.len() );    
    }


    // for entry in jordan.column_reverse( birth_column.clone() ) {
    //     let mut linexpr         =   LinExpr::new();
    //     let coeff                       =   to_float(entry.val());
    //     linexpr.add_constant( coeff ); // note we don't put a variable here
        
    //     // book keeping
    //     rowkey_to_rownum.insert( entry.key(), rowkey_to_rownum.len() );
    //     rownum_to_axz.push( linexpr  );
    //     rownum_to_str.push( coeff.to_string() );

    //     // calculate the objective value of the original cycle (incrementally)        
    //     cost_b += ( coeff * c( &entry.key() ) ).abs();

    //     // make a record of the cycle (incrementaly)
    //     b.push( (entry.key(), coeff) );

    //     // introduce a y variable for this constriant
    //     let variable_name=   format!( "y{}", rownum_to_y.len() );
    //     let variable        =   add_ctsvar!(model, name: &variable_name, bounds: ..)?;      
    //     rownum_to_y.push( variable );
    //     rownum_to_rowkey.push( entry.key() );
    // }

    //  extend lhs to a family of expressions for Ax + b
    let mut x_kvpairs              =   Vec::new();
    for column_index in column_indices.into_iter() {

        //  if we've made it here, then we will add a column to the matrix A
        //  ----------------------------------------------------------------
        
        //  define a new x variable, and add it to our indexed set of x's
        let var_x_name=   format!( "x{}", x_kvpairs.len() );
        let var_x        =   add_ctsvar!(model, name: &var_x_name, bounds: ..)?;

        //  generate a column of A
        let column          =   a( & column_index );

        //  upudate the expression Ax + b
        for entry in column.into_iter() {
            let key             =   entry.key();
            let coeff           =   entry.val().to_f64().unwrap();
            if let Some(rownum) = rowkey_to_rownum.get( &key ) {
                // if we already have row for this variable, then add an entry to the row
                rownum_to_axz[rownum.clone()].add_term( coeff, var_x );
                // rownum_to_str[rownum.clone()].push_str( &format!(" + {:?}x{:?}", coeff, x_kvpairs.len()) );
            } else {
                // otherwise start a new row
                let mut linexpr         =   LinExpr::new();
                linexpr.add_term( coeff, var_x );

                rowkey_to_rownum.insert( key.clone(), rowkey_to_rownum.len() );   // update the hashmap rowkey --> rownum   
                rownum_to_rowkey.push( key ); //  update the vector rownum --> rowkey
                rownum_to_axz.push( linexpr ); // push the new constraint to our vector of constraints

                let var_y_name=   format!( "y{}", rownum_to_y.len() ); // name the row
                let var_y_self        =   add_ctsvar!(model, name: &var_y_name, bounds: ..)?;   // add the row variable to the model
                rownum_to_y.push( var_y_self );  // hang on to a copy of the variable                                                   (!!!! not sure if this is necessary)
            }
        }    
        x_kvpairs.push( ( column_index, var_x.clone()  ) );      // keep a copy of the index and the associated variable            
    }

    // add constraints to the model
    // let mut constraints         =   Vec::new(); 
    for (rownum, lhs) in rownum_to_axz.iter().cloned().enumerate() {// (rownum, (key, axz_row)) in rownum_to_axz.iter().enumerate() {
        
        let y                   =   rownum_to_y[ rownum ].clone();        
        let mut rhs = LinExpr::new();
        rhs.add_term( 1.0, y );

            // "positive" constraints: y ≥ Ax + b
            let constraint_name_pos =   format!( "p{}", rownum );                    
            let constr = IneqExpr{ 
                                            lhs: grb::Expr::Linear(lhs.clone()), 
                                            sense: ConstrSense::Less, 
                                            rhs: grb::Expr::Linear(rhs.clone()), 
                                        };                    
            model.add_constr( &constraint_name_pos, constr );

            // "positive" constraints: y ≥ - Ax - b
            let constraint_name_neg =   format!( "n{}", rownum );  
            let mut lhs_neg     =   lhs.clone();
            lhs_neg.mul_scalar( -1.0 );
            let constr = IneqExpr{ 
                                lhs: grb::Expr::Linear(lhs_neg), 
                                sense: ConstrSense::Less, 
                                rhs: grb::Expr::Linear(rhs.clone()), 
                            };                    
            model.add_constr( &constraint_name_neg, constr );          
    }

    // add objective function to the model
    let mut linexpr         =   LinExpr::new();
    for (rowkey, var) in rownum_to_rowkey.iter().zip( rownum_to_y.iter().cloned() ) {
        let cost_coefficient = c( rowkey );
        match OrderedFloat(cost_coefficient) >= OrderedFloat(0f64) {
            true => { linexpr.add_term( cost_coefficient, var ); }
            false => { return Err( grb::Error::NotYetSupported("Error in OAT's optimization utility: the user provided a cost vector `c` with a negative entry; see the documentation for `minimize_l1` for further details".to_string()) ) }
        }        
    }
    model.set_objective( linexpr, Minimize );


    let construction_duration = start.elapsed();
    if verbose {
        println!("Finished constructing L1 optimization program.");
        println!("Construction took {:?} seconds.", construction_duration.as_secs_f64());
        println!("Passing program to solver.\n");
    }    


    
    // solve
    let start = std::time::Instant::now();

    model.write("model.lp");
    model.optimize();

    assert_eq!(model.status()?, Status::Optimal);

    model.write("mip.lp")?;
    model.write("mip.sol")?;


    let solve_duration = start.elapsed();
    if verbose {
        println!("\nFinished solving.\n");
        println!("Solver took {:?} seconds in walltime.", solve_duration.as_secs_f64());

    }    

    // reformat the solution
    let cost_axb     =   model.get_attr(attr::ObjVal)?;

    let x       
        =   x_kvpairs.iter().map(
                    | ( key, var ) |  
                    ( 
                        key.clone(), 
                        model
                            .get_obj_attr(attr::X, var)                                        
                            .expect( &format!("Error retreiving the coefficient of `x` for key {:?}", key.clone()) ) 
                    )
                )
            .filter(|x| x.1 != 0.0 )
            .collect_vec();

    let y       =   rownum_to_axz.iter().enumerate().map(
                                    |(rownum, linexpr)| -> (Key,f64)
                                    {
                                        let key     =   rownum_to_rowkey[rownum].clone();
                                        ( 
                                            key.clone(), 
                                            linexpr
                                                .get_value(&model)
                                                .expect( &format!("Error retreiving the coefficient of `x` for key {:?}", key.clone()) ) 
                                        )    
                                    }                                
                                )
                                .filter( |x: &(Key,f64)| x.1 != 0.0 )
                                .collect_vec();
    
    let answer = SolutionL1{ 
        cost_b, 
        cost_y: cost_axb, 
        b, 
        y, 
        x,
        construction_time: construction_duration.as_secs_f64(),
        solve_time: solve_duration.as_secs_f64(),
    };
    Ok(answer)
}






// // add decision variables with no bounds
// let x1 = add_ctsvar!(model, name: "x1", bounds: ..)?;
// let x2 = add_intvar!(model, name: "x2", bounds: ..)?;

// // add linear constraints
// let c0 = model.add_constr("c0", c!(x1 + 2*x2 >= -14))?;
// let c1 = model.add_constr("c1", c!(-4 * x1 - x2 <= -33))?;
// let c2 = model.add_constr("c2", c!(2* x1 <= 20 - x2))?;

// // model is lazily updated by default
// assert_eq!(model.get_obj_attr(attr::VarName, &x1).unwrap_err(), grb::Error::ModelObjectPending);
// assert_eq!(model.get_attr(attr::IsMIP)?, 0);

// // set the objective function, which updates the model objects (variables and constraints).
// // One could also call `model.update()`
// model.set_objective(8*x1 + x2, Minimize)?;
// assert_eq!(model.get_obj_attr(attr::VarName, &x1)?, "x1");
// assert_eq!(model.get_attr(attr::IsMIP)?, 1);

// // write model to the file.
// model.write("model.lp")?;

// // optimize the model
// model.optimize()?;
// assert_eq!(model.status()?, Status::Optimal);

// // Querying a model attribute
// assert_eq!(model.get_attr(attr::ObjVal)? , 59.0);

// // Querying a model object attributes
// assert_eq!(model.get_obj_attr(attr::Slack, &c0)?, -34.5);
// let x1_name = model.get_obj_attr(attr::VarName, &x1)?;

// // Querying an attribute for multiple model objects
// let val = model.get_obj_attr_batch(attr::X, vec![x1, x2])?;
// assert_eq!(val, [6.5, 7.0]);

// // Querying variables by name
// assert_eq!(model.get_var_by_name(&x1_name)?, Some(x1));




#[cfg(test)]
mod tests {

    /// Minimize the l1 norm of
    /// ```
    /// | 1 |       | 2 |
    /// | 1 | * x + | 2 |
    /// ```
    #[test]
    fn test_l1_a() {       
        use crate::utilities::optimization::l1::minimize_l1;               

        let constraint_data = vec![ vec![ (0,1.), (1,1.) ] ];
        let a = |i: & usize | -> Vec<(usize,f64)> { constraint_data[*i].clone() }; // the matrix 2x1 matrix [1; 1]
        let b = vec![ (0usize,2f64), (1usize,2f64) ]; // the vector [2, 2]
        let c = |x: & usize| -> f64 { 1.0 }; // the vector [1, 1]
        let column_indices = vec![0];        

        let solution = minimize_l1( a.clone(),b.clone(),c.clone(), column_indices).unwrap();

        println!("{:?}", solution);
        assert_eq!( solution.x(), &vec![(0,-2.0)] );
        assert_eq!( solution.b(), &b );
        assert_eq!( solution.y(), &Vec::<(usize,f64)>::new() );
        assert_eq!( solution.cost_b(), &4.0 );
        assert_eq!( solution.cost_axb(), &0.0 );           
        
    }   

    #[test]
    /// Minimize the l1 norm of
    /// ```
    /// | 1 |       | 1 |
    /// | 1 | * x + | 1 |
    /// | 1 |       | 0 |
    /// ```
    fn test_l1_b() {  
        use crate::utilities::optimization::l1::minimize_l1;        

        let constraint_data = vec![ vec![ (0,1.), (1,1.), (2,1.), ] ];
        let a = |i: & usize | -> Vec<(usize,f64)> { constraint_data[*i].clone() };
        let b = vec![ (0,1.), (1,1.), ];
        let c = |x: & usize| -> f64 { 1.0 };
        let column_indices = vec![0];        

        let solution = minimize_l1( a.clone(),b.clone(),c.clone(), column_indices).unwrap();

        println!("{:?}", solution);
        assert_eq!( solution.x(), &vec![(0,-1.0)] );
        assert_eq!( solution.b(), &b );
        assert_eq!( solution.y(), &vec![(2,-1.0)] );
        assert_eq!( solution.cost_b(), &2.0 );
        assert_eq!( solution.cost_axb(), &1.0 );        
    }
}

