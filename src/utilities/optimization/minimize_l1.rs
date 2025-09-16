//! Minimize `c' y` subject to `y = |Ax + b|`
//! 
//! 
//! 
//! This module provides a [function](minimize_l1) to solve problems of form 
//! 
//! ```text
//! minimize        c' |y|
//! subject to      y = Ax + b
//!                 |y| is the entry-wise absolute value of y
//!                 y and x are unconstrained continuous variables
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

use good_lp::{SolverModel, Solution, Expression, Variable, ProblemVariables, default_solver};

use itertools::Itertools;
use num::ToPrimitive;
use sprs::linalg::ordering::start;

use crate::algebra::vectors::entries::{KeyValGet};

use std::fmt::Debug;
use std::hash::Hash;
use std::collections::HashMap;

use derive_getters::{Getters, Dissolve};
use derive_new::new;


/// Solution to an L1 optimization problem
/// 
/// This struct represents an optimal solution to a problem of form
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
/// The struct has methods to return `b`, `x`, `y`, as well as the objective values `cost(b) = c'|b|` and `cost(y) = c`|y|`.
#[derive(new,Clone,Debug,Getters,Dissolve)]
pub struct SolutionL1< Key, Coefficient >{
    pub x:                  Vec< (Key, f64) >,    
    pub b:                  Vec< (Key, Coefficient) >,
    pub y:                  Vec< (Key, f64) >,
    pub cost_b:             f64,
    pub cost_y:             f64,
    pub construction_time:  f64, // time it took to construct the problem
    pub solve_time:         f64, // time it took to solve the problem
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
pub fn minimize_l1< 
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
    let start = std::time::Instant::now();

    //  holds information about the variables
    let mut variables                       =   ProblemVariables::new();


    //  holds a pair (Ax+b, y), which represents a constraint Ax + b ≤ y (we will add the other row, - Ax - b ≤ y, later)
    pub struct Row{ lhs: Expression, rhs: Variable  }

    // hashmaps
    let mut key_to_row                      =   HashMap::<Key,Row>::new(); // hashmap sendin a key to a pair (Ax+b, y), which represents a constraint Ax + b ≤ y (we will add the other row, - Ax - b ≤ y, later)
    let mut key_to_colvar                   =   HashMap::<Key,Variable>::new();  

    let b                   =   b.into_iter()
                                    .map( |entry| (entry.key(), entry.val() ) )
                                    .collect_vec();
    let cost_b                                  
        =   b.iter().map(
                |(k,v)| 
                c(k) * v.to_f64().unwrap().abs() // cost of the entry times the absolute value of the coefficient
            ).sum();

    // add the coefficients of `b` as constants in the constraints
    for entry in b.iter() {

        let key                             =   entry.key();
        if key_to_row.contains_key( &key ) 
            { return Err("Error in OAT's optimization utility: the user provided a cost vector `b` that has multiple entries with the same index.".to_string()) }

        // make a new linear expression
        let mut expression                  =   Expression::with_capacity(1);
        expression                          +=  entry.val().to_f64().unwrap();

        // make a new row variable
        let y                               =   variables.add_variable();
        key_to_row.insert( key, Row{ lhs: expression, rhs: y } );   
    }

    //  extend lhs to a family of expressions for Ax + b
    let mut num_entries         =   0;
    for column_index in column_indices.into_iter() {
        
        //  define a new x variable, and add it to our indexed set of x's
        let x                   =   variables.add_variable();
        key_to_colvar.insert( column_index.clone(), x );

        //  generate a column of A
        let column_entries      =   a( & column_index );

        //  upudate the expression Ax + b
        for entry in column_entries.into_iter() {
            num_entries         +=  1;
            let key             =   entry.key();
            let coeff           =   entry.val().to_f64().unwrap();
            if let Some(row)    =   key_to_row.get_mut( &key ) {
                // if we already have row for this variable, then add an entry to the row
                row.lhs.add_mul( coeff, x );
            } else {
                // otherwise start a new row

                // make a new linear expression
                let mut expression                  =   Expression::with_capacity(1);
                expression.add_mul( coeff, x );

                // make a new row variable
                let y                               =   variables.add_variable();
                key_to_row.insert( key, Row{ lhs: expression, rhs: y } );
            }
        }    
    }

    // form the objective function as an expression
    let mut objective = Expression::with_capacity( key_to_row.len() );
    for (key,row) in key_to_row.iter() {
        let coeff = c(key).to_f64().unwrap();
        objective.add_mul( coeff, row.rhs  );
    }

    // form the model
    let mut problem     =   variables
                            .minimise( objective.clone() )
                            .using( default_solver );

    // add constraints
    for row in key_to_row.values() {
        problem         =   problem.with(   row.lhs.clone() << row.rhs );
        problem         =   problem.with( - row.lhs.clone() << row.rhs );
    }



    let construction_duration = start.elapsed();
    if verbose {        
        println!("Finished construcing L1 optimization program.");
        println!("Construction took {:?} seconds.", construction_duration.as_secs_f64());
        println!("Constraint matrix has {:?} nonzero entries.", num_entries);
        println!("Passing program to solver.\n");
    }



    // solve
    let start = std::time::Instant::now();

    let solution        =   problem.solve().unwrap();

    let solve_duration = start.elapsed();
    if verbose {
        println!("\nFinished solving.\n");
        println!("Solver took {:?} seconds in walltime.", solve_duration.as_secs_f64());

    }


    // reformat the solution
    let cost_y          =   solution.eval(objective);

    let x               =   key_to_colvar.iter().map(
                                        | ( key, var ) |  
                                        ( 
                                            key.clone(), 
                                            solution.value(*var)
                                        )
                                    )
                                .filter(|x| x.1 != 0.0 )
                                .collect_vec();

    let y               =   key_to_row.iter().map(
                                    |(key, row)| -> (Key,f64)
                                        {
                                            ( 
                                                key.clone(), 
                                                row.lhs.eval_with(&solution)
                                            )    
                                        }                                
                                )
                                .filter( |x: &(Key,f64)| x.1 != 0.0 )
                                .collect_vec();


    if verbose {
        println!("GOODLP solution: {:?}", solution.into_inner() );
    }


    let answer = SolutionL1{ 
        cost_b, 
        cost_y, 
        b, 
        y, 
        x,
        construction_time: construction_duration.as_secs_f64(),
        solve_time: solve_duration.as_secs_f64(),
    };
    Ok(answer)
}




/// Solution to an L1 optimization problem
/// 
/// This struct represents an optimal solution to a problem of form
/// 
/// ```text
/// minimize        c ' | b + x |
/// subject to      Ax = 0
///                 | b + x | is the entry-wise absolute value of y
///                 x is an unconstrained continuous variable
/// ```
/// 
/// We solve this problem by solving the following equivalent linear programming problem:
/// 
/// ```text
/// minimize        c' y
/// subject to      y >=  + b + x
///                 y >=  - b - x
///                 Ax =  0
///                 y and x are unconstrained continuous variables
/// ```
/// 
/// The struct has methods to return `b`, `x`, `y`, as well as the objective values `cost(b) = c'|b|` and `cost(y) = c`|y|`.
pub fn minimize_l1_kernel< 
            Key,
            ConstraintMatrix,
            ConstraintMatrixColumn,
            ConstraintVector, 
            ColumnIndexIterable,            
            CostVector,
            Coefficient,
        >
    ( 
        mut a:                      ConstraintMatrix,
        b:                          ConstraintVector,        
        mut c:                      CostVector,
        column_indices:             ColumnIndexIterable,         
    ) 
    -> 
    Result< SolutionL1<Key, Coefficient>, String >
    where
        ConstraintMatrix:               FnMut( & Key )->ConstraintMatrixColumn,
        ConstraintVector:               IntoIterator,
        ConstraintVector::Item:         KeyValGet< Key = Key, Val = Coefficient >,
        ConstraintMatrixColumn:         IntoIterator,
        ConstraintMatrixColumn::Item:   KeyValGet< Key = Key, Val = Coefficient >,        
        Coefficient:                    Clone + ToPrimitive, // enables casting to f64
        Key:                            Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
        ColumnIndexIterable:            Clone + IntoIterator< Item = Key >, 
        CostVector:                     FnMut( & Key) -> f64,                      
{

    let start = std::time::Instant::now();


    //  holds information about the variables
    let mut variables                       =   ProblemVariables::new();

    //  holds a pair (x+b, y), which represents a constraint x + b ≤ y (we will add the other constraint, - x - b ≤ y, later)
    pub struct Col{ x: Variable, y: Variable, x_plus_b: Expression, }

    // hashmaps indexed by row/col indices
    let mut key_to_col                      =   HashMap::<Key,Col>::new(); // hashmap sending a key to an entry of x + b
    let mut key_to_row_expr                 =   HashMap::<Key,Expression>::new(); // hashmap sending a key to a row of Ax

    // drain b into a hashmap
    let b: HashMap<Key, Coefficient>        =   b.into_iter()
                                                    .map( |entry| (entry.key(), entry.val() ) )
                                                    .collect();
    // precompute the cost of b
    let cost_b                                  
        =   b.iter().map(
                |(k,v)| 
                c(k) * v.to_f64().unwrap().abs() // cost of the entry times the absolute value of the coefficient
            ).sum();

    //  update the expression b + x, and form the expressions in Ax
    for column_index in column_indices.into_iter() {

        // generate variables x, y, and an expression x + b for this column
        let x                   =   variables.add_variable();
        let y                   =   variables.add_variable();
        let mut x_plus_b        =   Expression::with_capacity(1);
        
        x_plus_b.add_mul( 1.0, x );
        if let Some(coeff)      =   b.get( & column_index ) {
            x_plus_b            +=  coeff.to_f64().unwrap();
        }

        // record x, y, and x + b in a hashmap
        key_to_col.insert(
            column_index.clone(),
            Col{ x, y, x_plus_b }
        );

        //  generate a column of A
        let column_entries      =   a( & column_index );

        //  upudate the expression Ax
        for entry in column_entries.into_iter() {
            let key             =   entry.key();
            let coeff           =   entry.val().to_f64().unwrap();
            if let Some(row)    =   key_to_row_expr.get_mut( &key ) {
                // if we already have row for this variable, then add an entry to the row
                row.add_mul( coeff, x );
            } else {
                // otherwise start a new row

                // make a new linear expression
                let mut expression                  =   Expression::with_capacity(1);
                expression.add_mul( coeff, x );
                key_to_row_expr.insert( key, expression );
            }
        }    
    }

    // form the objective function as an expression
    let mut objective = Expression::with_capacity( key_to_col.len() );
    for (key, col) in key_to_col.iter() {
        let coeff = c(key).to_f64().unwrap();
        objective.add_mul( coeff, col.y  );
    }

    // form the model
    let mut problem     =   variables
                            .minimise( objective.clone() )
                            .using( default_solver );

    // add constraints
    for col in key_to_col.values() {
        problem         =   problem.with(   col.x_plus_b.clone() << col.y );
        problem         =   problem.with( - col.x_plus_b.clone() << col.y );
    }


    let construction_duration = start.elapsed();



    // solve
    let start = std::time::Instant::now();    
    let solution        =   problem.solve().unwrap();
    let solve_duration = start.elapsed();


    // reformat the solution
    let cost_y          =   solution.eval(objective);

    let x               =   key_to_col.iter().map(
                                        | ( key, col ) |  
                                        ( 
                                            key.clone(), 
                                            solution.value(col.x)
                                        )
                                    )
                                .filter(|x| x.1 != 0.0 )
                                .collect_vec();

    let y               =   key_to_col.iter().map(
                                    |(key, col)| -> (Key,f64)
                                        {
                                            ( 
                                                key.clone(), 
                                                solution.value(col.y)
                                            )    
                                        }                                
                                )
                                .filter( |x: &(Key,f64)| x.1 != 0.0 )
                                .collect_vec();

    let b: Vec<(Key, Coefficient)>  =   b.into_iter().collect(); // convert b from hashmap to vec
    
    let answer = SolutionL1{ 
        cost_b, 
        cost_y, 
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
    fn test_l1_good_a() {       
        use crate::utilities::optimization::minimize_l1::minimize_l1;               

        let constraint_data = vec![ vec![ (0,1.), (1,1.) ] ];
        let a = |i: & usize | -> Vec<(usize,f64)> { constraint_data[*i].clone() }; // the matrix 2x1 matrix [1; 1]
        let b = vec![ (0usize,2f64), (1usize,2f64) ]; // the vector [2, 2]
        let c = |_x: & usize| -> f64 { 1.0 }; // the vector [1, 1]
        let column_indices = vec![0]; 
        let verbose = true;       

        let solution = minimize_l1( a,b.clone(),c, column_indices, verbose).unwrap();

        println!("{:?}", solution);
        assert_eq!( solution.x(), &vec![(0,-2.0)] );
        assert_eq!( solution.b(), &b );
        assert_eq!( solution.y(), &Vec::<(usize,f64)>::new() );
        assert_eq!( solution.cost_b(), &4.0 );
        assert_eq!( solution.cost_y(), &0.0 );           
        
    }   

    #[test]
    /// Minimize the l1 norm of
    /// ```
    /// | 1 |       | 1 |
    /// | 1 | * x + | 1 |
    /// | 1 |       | 0 |
    /// ```
    fn test_l1_good_b() {  
        use crate::utilities::optimization::minimize_l1::minimize_l1;        

        let constraint_data = vec![ vec![ (0,1.), (1,1.), (2,1.), ] ];
        let a = |i: & usize | -> Vec<(usize,f64)> { constraint_data[*i].clone() };
        let b = vec![ (0,1.), (1,1.), ];
        let c = |_x: & usize| -> f64 { 1.0 };
        let column_indices = vec![0];        
        let verbose = true;

        let solution = minimize_l1( a,b.clone(),c, column_indices, verbose).unwrap();

        println!("{:?}", solution);
        assert_eq!( solution.x(), &vec![(0,-1.0)] );
        assert_eq!( solution.b(), &b );
        assert_eq!( solution.y(), &vec![(2,-1.0)] );
        assert_eq!( solution.cost_b(), &2.0 );
        assert_eq!( solution.cost_y(), &1.0 );        
    }


    // THIS PROBLEM MAKES MORE SENSE FOR L2 NORM THAN L1 NORM
    // #[test]
    // /// Minimize the l1 norm of
    // /// ```
    // /// |  1   0   0 |       | 1 |
    // /// | -1   1   0 | * x + | 0 |
    // /// |  0  -1   1 |       | 0 |
    // /// |  0   0  -1 |       | 0 |
    // /// ```
    // fn test_l1_good_c() {  
    //     use crate::utilities::optimization::minimize_l1::minimize_l1;        

    //     let constraint_data = vec![ 
    //             vec![ (0,1.),           ], 
    //             vec![ (0,-1.), (1,1.),  ], 
    //             vec![ (1,-1.), (2,1.),  ], 
    //             vec![ (2,-1.),          ],                 
    //         ];
    //     let a = |i: & usize | -> Vec<(usize,f64)> { constraint_data[*i].clone() };
    //     let b = vec![ (0,1.) ];
    //     let c = |_x: & usize| -> f64 { 1.0 };
    //     let column_indices = vec![0];        

    //     let solution = minimize_l1( a,b.clone(),c, column_indices).unwrap();

    //     println!("{:?}", solution);
    //     assert_eq!( solution.x(), &vec![(0,-1.0)] );
    //     assert_eq!( solution.b(), &b );
    //     assert_eq!( solution.y(), &vec![(2,-1.0)] );
    //     assert_eq!( solution.cost_b(), &2.0 );
    //     assert_eq!( solution.cost_y(), &1.0 );        
    // }    
}

