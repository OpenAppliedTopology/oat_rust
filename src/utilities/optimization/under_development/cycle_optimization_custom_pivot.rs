

use std::hash::HashMap;

use oat_rust::algebra::rings::types::native::FieldRational64;
use crate::utilities::iterators::general::find_min;





// struct OptimizationBase< 'a, MatrixToFactor > 
//     where
//         MatrixToFactor:                           MatrixOracle,
//         MatrixToFactor::ColumnIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
//         MatrixToFactor::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
//         // OrderOperatorByKey<usize, MatrixToFactor::Coefficient, (usize, MatrixToFactor::Coefficient)>: JudgePartialOrder< (usize, MatrixToFactor::Coefficient)>, // this seems extraneous but the compiler seems to want it
//         MatrixToFactor::Row:          IntoIterator,
//         MatrixToFactor::RowEntry:     KeyValGet < Key = MatrixToFactor::ColumnIndex, Val = MatrixToFactor::Coefficient >,
//         RingOperator:                           Clone + Semiring< MatrixToFactor::Coefficient > + Ring< MatrixToFactor::Coefficient > + DivisionRingOperations< MatrixToFactor::Coefficient >,

// {
//     umatch:                         &'a Umatch<MatrixToFactor, RingOperator, OrderOperatorForRowEntries, OrderOperatorForColumnEntries >,
//     order_operator_all_entries:   OrderOperatorForRowEntries,   
//     ring_operator:                  RingOperator,
// }

type Coeff = Ratio< i64 >;

struct OptimizationBase //< Coeff, RingOperator, > 
    // where
    //     RingOperator:                           Clone + Semiring< MatrixToFactor::Coefficient > + Ring< MatrixToFactor::Coefficient > + DivisionRingOperations< MatrixToFactor::Coefficient >,

{
    reduced_head:   Vec< Vec< (usize,Coeff) > >,   
    reduced_tail:   Vec< Vec< (usize,Coeff) > >,   
    ring_operator:  RingOperator,
}







 // =========================

// struct Optimizer< 'a, > // Coeff, RingOperator, > 
//     // where
//     //     RingOperator:                           Clone + Semiring< MatrixToFactor::Coefficient > + Ring< MatrixToFactor::Coefficient > + DivisionRingOperations< MatrixToFactor::Coefficient >,

// {
//     optimization_base:  &'a OptimizationBase, //< Coeff, RingOperator >,
//     row_to_col:         Vec< usize >, // paired row to paired column
//     col_to_row:         Vec< usize >, // paired column to paired row
//     cycle_initial:      HashMap< usize, Coeff >,  
//     cycle_slice:        HashMap< usize, Coeff >,   // = cycle_initial[col_to_row]
//     cob_matrix:         Vec<Vec< (usize, Coeff) >>, 
//     solution_cost:      Ratio<i64>,
//     solution:           Vec< (usize, Coeff) >,    
//     // solution_bdry:      Vec< (usize, Coeff) >,    
//     loose_constraints:  HashMap< usize, Coeff >,
//     loose_rows:         HashMap<    usize, 
//                                     Vec< (usize,Coeff) > 
//                             >,
// }

// impl < 'a > //, Coeff, RingOperator, > 

//     Iterator for 
//     Optimizer
//         // < 'a, Coeff, RingOperator, > 
//     // where
//     //     RingOperator:  Clone + Semiring< MatrixToFactor::Coefficient > + Ring< MatrixToFactor::Coefficient > + DivisionRingOperations< MatrixToFactor::Coefficient >,        
// {
//     type Item = ();

//     fn next( &mut self ) -> Option< Self::Item > {

//         let num_cols                =   self.row_to_col.len();
//         let ring_operator               =   FieldRational64::new();

//         let mut new_cost        =   self.solution_cost;
//         let mut new_pivot_row   =   0;
//         let mut new_pivot_col   =   0;
//         let mut new_constraint_val  =   0;

//         // check to see if there's a direction with negative cost
//         let mut delta_cost      =   Ratio::new(0,1);

//         for ( row_indx, row_coeff) in self.loose_constraints.iter() {
//             let row             =   self.loose_rows.get( & row_indx ).unwrap();
//             for (col_indx, col_coeff) in row.iter() {
//                 prior_zvalue    =   self.cycle_slice.insert( col_indx, self.cycle_initial.get( row_indx ) ).unwrap(); // update z_I temporarily
//                 new_seed        =   veca_vecp_index__vecpivoted( &row, &self.cycle_slice );
//                 _               =   self.cycle_slice.insert( col_indx, prior_zvalue ); // change z_I back                                                            
//                 difference_new  =   multiply_row_vector_with_matrix( 
//                                             new_seed.into_iterator(), 
//                                             & self.cob_matrix,  
//                                             ring_operator.clone(),
//                                             OrderOperatorByKey::new(),
//                                         );
//                 for (key,val1) in self.cycle_initial.iter() {
//                     val2        =   difference_new.get(key).unwrap_or(zero);
//                     difference_new.insert( key, val1-val2 );                    
//                 }
//                 candidate_cost  =   difference_new.into_values().map(|x| x.abs() ).sum();

//                 if candidate_cost < new_cost {
//                     new_cost            =   candidate_cost;
//                     new_pivot_row       =   row_indx;
//                     new_pivot_col       =   col_indx;
//                 }
//             }
//         }

//         if new_cost < self.solution_cost {
//             let row             =   self.loose_rows.get( & new_pivot_row ).unwrap().clone();
//             _       =   self.cycle_slice.insert( col_indx, self.cycle_initial.get( row_indx ) ).unwrap();
//             row_of_inverse = self.cycle_slice.clone();

//         }

//         return Some(())
//     }
// }



// /// Returns the result of pivoting on index `pivot_index` of the active vector, effected on the passive vector.
// pub fn veca_vecp_index__vecpivoted(  
//                 active:         & HashMap< usize, Coeff >, 
//                 passive:        & HashMap< usize, Coeff >,                 
//                 pivot_index:    usize ) 
//             -> HashMap< usize, Coeff > 
// {
//     let zero = Ration::new(0,1);
    
//     let mut pivot_coeff_active = active.get( index ).unwrap();
//     if pivot_coeff_active == zero { panic!("Error: cannot pivot on entry {:?} because this entry is zero", pivot_index) }

//     let mut passive_clone = passive.clone();
//     if let Some( pivot_coeff_passive )  = passive.get( pivot_index ) {
//         if pivot_coeff_passive == zero { return passive_clone }
//         else {
//             for (key,val) in active.iter() {
//                 if key == pivot_index { continue }

//                 new_val = passive_clone.get(key).unwrap_or(zero) - pivot_coeff_passive * val / pivot_coeff_active
//                 passive_clone.insert( key, new_val );
//             }
//         }
//         return passive_clone

//     } else {
//         return passive_clone()
//     }


// }







// =========================

#[allow(non_snake_case)]
struct Optimizer // Coeff, RingOperator, > 
    // where
    //     RingOperator:                           Clone + Semiring< MatrixToFactor::Coefficient > + Ring< MatrixToFactor::Coefficient > + DivisionRingOperations< MatrixToFactor::Coefficient >,

{
    A:                  Vec< SparseVecUsize< Coeff > >,
    C:                  Vec< SparseVecUsize< Coeff > >,
    B:                  Vec< usize >, 
    N:                  Vec< usize >,       
    Ns:                 Vec< bool >,
    b:                  SparseVecUsize< Coeff >,
    c:                  SparseVecUsize< Coeff >,
    x:                  SparseVecUsize< Coeff >,
}

impl Optimizer {
    pub fn m( &self )-> usize { self.A.len() }
    pub fn n( &self )-> usize { self.I.len() }    
}

impl 

    Iterator for 
    Optimizer
        // < 'a, Coeff, RingOperator, > 
    // where
    //     RingOperator:  Clone + Semiring< MatrixToFactor::Coefficient > + Ring< MatrixToFactor::Coefficient > + DivisionRingOperations< MatrixToFactor::Coefficient >,        
{
    type Item = ();

    fn next( &mut self ) -> Option< Self::Item > {

        // check to see if we are at optimum
        let bbar = product( self.C, self.c.iter().map(|(k,v)| (self.B[k],v) ) );
        let iloc = bbar.iter_vals().find(|x| x < zero );
        if iloc.is_none() { return Some(()) } // we are at optimum

        // try to perform a pivot        
        let i    = i.unwrap();
        let Ci   = C.column_reverse(i);
        let score = |(j,sgn)|  {
            constraint_val_plusminus = inner_product( & self.A[j], & self.x );
            cbar =  sgn * ( constraint_val_plusminus - self.c[j] );
            Aij = inner_product( & self.A[j], & Ci );
        };

        let minval;
        let mut piv;
        let mut sgn;
        let m               =   self.m();
        let n               =   self.n();
        for ( j, sgn) in N.iter().cloned().zip( Ns.iter().cloned() ) {

            Abar  = product_with_col( self.A[i], self.C, i );
            cbar  = product( self.A[i], self.x ) - self.b[i];
            if Abar.len() == 0 { continue }
            if cbar == zero { piv_candidate = (i,Abar.next().unwrap().0, zero) }
            j, minval = find_min( Abar );
            piv_candidate = (i,j,minval);
            if piv_candidate.2 < piv { piv = piv_candidate }
        }


        let num_cols                =   self.row_to_col.len();
        let ring_operator               =   FieldRational64::new();

        let mut new_cost        =   self.solution_cost;
        let mut new_pivot_row   =   0;
        let mut new_pivot_col   =   0;
        let mut new_constraint_val  =   0;

        // check to see if there's a direction with negative cost
        let mut delta_cost      =   Ratio::new(0,1);

        for ( row_indx, row_coeff) in self.loose_constraints.iter() {
            let row             =   self.loose_rows.get( & row_indx ).unwrap();
            for (col_indx, col_coeff) in row.iter() {
                prior_zvalue    =   self.cycle_slice.insert( col_indx, self.cycle_initial.get( row_indx ) ).unwrap(); // update z_I temporarily
                new_seed        =   veca_vecp_index__vecpivoted( &row, &self.cycle_slice );
                _               =   self.cycle_slice.insert( col_indx, prior_zvalue ); // change z_I back                                                            
                difference_new  =   multiply_row_vector_with_matrix( 
                                            new_seed.into_iterator(), 
                                            & self.cob_matrix,  
                                            ring_operator.clone(),
                                            OrderOperatorByKey::new(),
                                        );
                for (key,val1) in self.cycle_initial.iter() {
                    val2        =   difference_new.get(key).unwrap_or(zero);
                    difference_new.insert( key, val1-val2 );                    
                }
                candidate_cost  =   difference_new.into_values().map(|x| x.abs() ).sum();

                if candidate_cost < new_cost {
                    new_cost            =   candidate_cost;
                    new_pivot_row       =   row_indx;
                    new_pivot_col       =   col_indx;
                }
            }
        }

        if new_cost < self.solution_cost {
            let row             =   self.loose_rows.get( & new_pivot_row ).unwrap().clone();
            _       =   self.cycle_slice.insert( col_indx, self.cycle_initial.get( row_indx ) ).unwrap();
            row_of_inverse = self.cycle_slice.clone();

        }

        return Some(())
    }
}



/// Returns the result of pivoting on index `pivot_index` of the active vector, effected on the passive vector.
pub fn veca_vecp_index__vecpivoted(  
                active:         & HashMap< usize, Coeff >, 
                passive:        & HashMap< usize, Coeff >,                 
                pivot_index:    usize ) 
            -> HashMap< usize, Coeff > 
{
    let zero = Ration::new(0,1);
    
    let mut pivot_coeff_active = active.get( index ).unwrap();
    if pivot_coeff_active == zero { panic!("Error: cannot pivot on entry {:?} because this entry is zero", pivot_index) }

    let mut passive_clone = passive.clone();
    if let Some( pivot_coeff_passive )  = passive.get( pivot_index ) {
        if pivot_coeff_passive == zero { return passive_clone }
        else {
            for (key,val) in active.iter() {
                if key == pivot_index { continue }

                new_val = passive_clone.get(key).unwrap_or(zero) - pivot_coeff_passive * val / pivot_coeff_active
                passive_clone.insert( key, new_val );
            }
        }
        return passive_clone

    } else {
        return passive_clone()
    }


}