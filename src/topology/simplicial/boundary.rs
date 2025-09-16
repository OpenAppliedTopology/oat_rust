//! Boundaries of filtered and unfiltered simplices


use crate::algebra::rings::traits::RingOperations;




//  ---------------------------------------------------------------------------
//  BOUNDARY VECTORS
//  ---------------------------------------------------------------------------



//  BOUNDARY VECTOR DESCENDING
//  ---------------------------------------------------------------------------


/// Iterates over the terms of the boundary of a simplex, in descending lexicographic order.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::boundary::SimplexBoundaryDescend;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// 
/// let ring_operator = PrimeOrderField::new(3);
/// let simplex = vec![0,1,2];
/// let boundary = SimplexBoundaryDescend::new( simplex, ring_operator );
/// 
/// itertools::assert_equal( boundary, vec![ (vec![1,2], 1), (vec![0,2], 2), (vec![0,1], 1) ]);
/// ```
pub struct SimplexBoundaryDescend
                < Vertex, RingOperator, >   
    where
        Vertex:         Clone + Ord,
        RingOperator:   RingOperations,
{
    facet_opt:              Option< Vec< Vertex > >, // a facet that will be returned when `.next()` is called
    vertex_removed:         Vertex, // the vertex that has been removed
    removal_locus:          usize,  // the slot from which that vertex has been removed
    coefficient:            RingOperator::Element,
    ring_operator:          RingOperator,
}


impl < Vertex, RingOperator, >
    
    SimplexBoundaryDescend
        < Vertex, RingOperator, >   

    where
        Vertex:         Clone + Ord,
        RingOperator:   RingOperations,
{
    /// Throws an error if `simplex` is an empty vector.
    pub fn new( mut simplex: Vec< Vertex >, ring_operator: RingOperator ) -> Self {

        // each vertex has empty boundary; handle this edge case separately
        if simplex.len() == 1 {
            return SimplexBoundaryDescend
            {
                facet_opt:              None,                               
                vertex_removed:         simplex[0].clone(),      // arbitrary values           
                removal_locus:          1,                       // arbitrary values           
                coefficient:            RingOperator::one(),     // arbitrary values      
                ring_operator,           // arbitrary values
            }                   
        }        

        let removal_locus = 0;
        let coefficient = RingOperator::one();                
        let vertex_removed = simplex.remove( removal_locus );
        simplex.shrink_to_fit();

        SimplexBoundaryDescend
        {
            facet_opt:              Some( simplex ), // a facet that will be returned when `.next()` is called
            vertex_removed, // the vertex that has been removed
            removal_locus,  // the slot from which that vertex has been removed
            coefficient,
            ring_operator,
        }        
    }
}        


//  IMPLEMENT ITERATOR
//  ------------------------------------------

impl < Vertex, RingOperator, >

    Iterator for
    
    SimplexBoundaryDescend
        < Vertex, RingOperator, >   

    where
        Vertex:         Clone + Ord,
        RingOperator:   RingOperations,  
{
    type Item = ( Vec< Vertex >, RingOperator::Element );

    fn next( &mut self ) -> Option< Self::Item > {

        match &mut self.facet_opt{
            Some( facet ) => {
                // make a copy of the term to return, before we update the internal state
                let return_value = (facet.clone(), self.coefficient.clone());  
                
                // if there is another vertex to remove, then remove it
                if self.removal_locus < facet.len() {
                    let new_removed = facet[ self.removal_locus ].clone();                    
                    facet[ self.removal_locus ] = self.vertex_removed.clone();
                    self.removal_locus += 1;
                    self.vertex_removed = new_removed;
                    self.coefficient = self.ring_operator.negate( self.coefficient.clone() );
                } 
                // otherwise mark that we have extracted the last facet
                else {
                    self.facet_opt = None;
                }

                Some( return_value )
            }
            None => { None }
        }
    }
}


// pub struct SimplexBoundaryDescend
//                 < Vertex, RingOperator, >   
//     where
//         Vertex:         Clone + Ord,
//         RingOperator:   Ring,
//         RingElement:    Clone,
// {
//     facet:                  Vec< Vertex >,
//     vertex_removed:         Vertex,
//     prior_removal_locus:    usize,
//     coefficient:            RingOperator::ElementOfRing,
//     ring_operator:          RingOperator,
// }


// impl < Vertex, RingOperator, >

//     Iterator for
    
//     SimplexBoundaryDescend
//         < Vertex, RingOperator, >   

//     where
//         Vertex:         Clone + Ord,
//         RingOperator:   Ring,
//         RingElement:    Clone,        
// {
//     type Item = ( Vec< Vertex >, RingElement );

//     fn next( &mut self ) -> Option< Self::Item > {

//         // if this case the internal state does not hold a valid facet
//         if self.prior_removal_locus > self.facet.len() { return None } 

//         // in this case the internal state holds the last facet; we return this value and do not update the internal state, except to mark that we have already returned the last facet (this is achieved by incrementing `prior_removal_locus`)
//         if self.prior_removal_locus == self.facet.len() { 
            
//             self.prior_removal_locus += 1;

//             return  Some( (self.facet.clone(), self.coefficient.clone()) )
//         } 

//         // in this case the internal state contains the next term of the boundary, and this is not the last term; therefore we will return a term, and update the internal state
        
//         // make a copy of the term we with to return, before updating the internal state
//         let return_value = (self.facet.clone(), self.coefficient.clone());

//         // update the internal state
//         let new_removed = self.facet[ self.prior_removal_locus ].clone(); // save a copy of the vertex we're about to remove
//         self.facet[ self.prior_removal_locus ] = self.vertex_removed.clone(); // place the last vertex that we had removed back into the simplex
//         self.coefficient = self.ring_operator.negate( self.coefficient.clone() ); // switch sign on the coefficient
//         self.prior_removal_locus += 1; // increment the locus of the next removal
//         self.vertex_removed = new_removed; // stash the value of the removed vertex in an internal state

//         return Some( return_value )
//     }
// }



// impl < Vertex, RingOperator, >
    
//     SimplexBoundaryDescend
//         < Vertex, RingOperator, >   

//     where
//         Vertex:         Clone + Ord,
//         RingOperator:   Ring,
//         RingElement:    Clone,  
// {
//     /// Throws an error if `simplex` is an empty vector.
//     pub fn new( mut simplex: Vec< Vertex >, ring_operator: RingOperator ) -> Self {

//         // each vertex has empty boundary; handle this edge case separately
//         if simplex.len() == 1 {
//             return SimplexBoundaryDescend
//             {
//                 facet:                  Vec::with_capacity(0), // a facet that will be returned when `.next()` is called
//                 vertex_removed:         simplex[0].clone(), // the vertex that has been removed
//                 prior_removal_locus:    1,  // the slot from which that vertex has been removed
//                 coefficient:            RingOperator::one(),
//                 ring_operator:          ring_operator,
//             }                   
//         }

//         // now handle simplices of dimension ≥ 1        
//         let removal_locus = 0;
//         let coefficient = RingOperator::one();                
//         let vertex_removed = simplex.remove(0);
//         simplex.shrink_to_fit();

//         SimplexBoundaryDescend
//         {
//             facet:                  simplex, // a facet that will be returned when `.next()` is called
//             vertex_removed:         vertex_removed, // the vertex that has been removed
//             prior_removal_locus:    removal_locus,  // the slot from which that vertex has been removed
//             coefficient:            coefficient,
//             ring_operator:          ring_operator,
//         }               
//     }
// }  





//  BOUNDARY VECTOR ASCENDING
//  ---------------------------------------------------------------------------

/// Iterates of the terms of the boundary of a simplex, in ascending lexicographic order.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::boundary::SimplexBoundaryAscend;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// 
/// let ring_operator = PrimeOrderField::new(3);
/// let simplex = vec![0,1,2];
/// let boundary = SimplexBoundaryAscend::new( simplex, ring_operator );
/// 
/// itertools::assert_equal( boundary, vec![ (vec![0,1], 1), (vec![0,2], 2), (vec![1,2], 1) ]);
/// ```
pub struct SimplexBoundaryAscend
                < Vertex, RingOperator, >   
    where
        Vertex:         Clone + Ord,
        RingOperator:   RingOperations,
{
    facet_opt:              Option< Vec< Vertex > >, // a facet that will be returned when `.next()` is called
    vertex_removed:         Vertex, // the vertex that has been removed
    removal_locus:          usize,  // the slot from which that vertex has been removed
    coefficient:            RingOperator::Element,
    ring_operator:          RingOperator,
}


impl < Vertex, RingOperator, >
    
    SimplexBoundaryAscend
        < Vertex, RingOperator, >   

    where
        Vertex:         Clone + Ord,
        RingOperator:   RingOperations,
{
    /// Throws an error if `simplex` is an empty vector.
    pub fn new( mut simplex: Vec< Vertex >, ring_operator: RingOperator ) -> Self {

        // each vertex has empty boundary; handle this edge case separately
        if simplex.len() == 1 {
            return SimplexBoundaryAscend
            {
                facet_opt:              None, // a facet that will be returned when `.next()` is called
                vertex_removed:         simplex[0].clone(), // the vertex that has been removed
                removal_locus:          1,  // the slot from which that vertex has been removed
                coefficient:            RingOperator::one(),
                ring_operator,
            }                   
        }        

        let removal_locus = simplex.len()-1;
        let coefficient = ring_operator.minus_one_to_power( removal_locus );                
        let vertex_removed = simplex.pop().unwrap();
        simplex.shrink_to_fit();

        SimplexBoundaryAscend
        {
            facet_opt:              Some( simplex ), // a facet that will be returned when `.next()` is called
            vertex_removed, // the vertex that has been removed
            removal_locus,  // the slot from which that vertex has been removed
            coefficient,
            ring_operator,
        }        
    }
}        


//  IMPLEMENT ITERATOR
//  ------------------------------------------

impl < Vertex, RingOperator, >

    Iterator for
    
    SimplexBoundaryAscend
        < Vertex, RingOperator, >   

    where
        Vertex:         Clone + Ord,
        RingOperator:   RingOperations,
{
    type Item = ( Vec< Vertex >, RingOperator::Element );

    fn next( &mut self ) -> Option< Self::Item > {

        match &mut self.facet_opt{
            Some( facet ) => {
                // make a copy of the term to return, before we update the internal state
                let return_value = (facet.clone(), self.coefficient.clone());  
                
                // if there is another vertex to remove, then remove it
                if self.removal_locus > 0 {
                    self.removal_locus -= 1;
                    let new_removed = facet[ self.removal_locus ].clone();
                    facet[ self.removal_locus ] = self.vertex_removed.clone();
                    self.vertex_removed = new_removed;
                    self.coefficient = self.ring_operator.negate( self.coefficient.clone() );
                }
                // otherwise mark that we have returned the last facet
                else {
                    self.facet_opt = None;
                }

                Some( return_value )
            }
            None => { None }
        }

        // match self.next_exchange_locus {
        //     Some( exchange_locus ) => {

        //         // make a copy of the term to return, before we update the internal state
        //         let return_value = (self.facet.clone(), self.coefficient.clone()); 



        //         // update the internal state
        //         let new_removed = self.facet[ self.prior_removal_locus ].clone(); // save a copy of the vertex we're about to remove
        //         self.facet[ self.prior_removal_locus ] = self.vertex_removed; // place the last vertex that we had removed back into the simplex
        //         self.coefficient = self.ring_operator.negate( self.coefficient ); // switch sign on the coefficient
        //         self.prior_removal_locus += 1; // increment the locus of the next removal
        //         self.vertex_removed = new_removed; // stash the value of the removed vertex in an internal state

        //         match exchange_locus == 0 {
                    
        //             true => {
        //                 //return value no update
        //             }
        //             false => {
        //                 // return value and update
        //             }
        //         }
        //     }
        //     None => {
        //         // don't return value
        //     }
        // }


        // // if this case the internal state does not hold a valid facet, because we have already returned the last facet
        // if self.prior_removal_locus > self.facet.len() { return None } 

        // // in this case the internal state holds the last facet; we return this value and do not update the internal state, except to mark that we have already returned the last facet
        // if self.prior_removal_locus == self.facet.len() { 
            
        //     self.prior_removal_locus += 1;

        //     return  Some( (self.facet.clone(), self.coefficient.clone()) )
        // } 

        // // in this case the internal state contains the next term of the boundary, and this is not the last term; therefore we will return a term, and update the internal state
        
        // // make a copy of the term we with to return, before updating the internal state
        // let return_value = (self.facet.clone(), self.coefficient.clone());

        // // update the internal state
        // let new_removed = self.facet[ self.prior_removal_locus ].clone(); // save a copy of the vertex we're about to remove
        // self.facet[ self.prior_removal_locus ] = self.vertex_removed; // place the last vertex that we had removed back into the simplex
        // self.coefficient = self.ring_operator.negate( self.coefficient ); // switch sign on the coefficient
        // self.prior_removal_locus += 1; // increment the locus of the next removal
        // self.vertex_removed = new_removed; // stash the value of the removed vertex in an internal state

        // return Some( return_value )
    }
}



//  ===========================================================================
//  ===========================================================================
//  BOUNDARY MATRIX INDEXED VIA A BIMAP
//  ===========================================================================
//  ===========================================================================


