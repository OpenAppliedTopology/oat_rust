//! The nerve of a poset
//! 
//! This module is under construction; see the source code for current progress.

// use std::collections::HashSet;
// use std::marker::PhantomData;

// use itertools::Merge;
// use crate::utilities::iterators::general::{TwoTypeIterator, OnlyDuplicates};






// //  NERVE OF THE POSET OF AN INDEXED FAMILY OF SETS
// //  ==================================================
// //  (ONLY PARTIALLY COMPLETED)

// // pub struct NerveCoboundaryIter< 'a, Vertex, RingOperator, RingElement > {
// //     chain:                  Vec< usize >,
// //     cover:                  &'a Vec< HashSet< Vertex > >,  
// //     next_openid_to_try:     usize,  
// //     ring_operator:          RingOperator,
// //     phantom_ringelement:    RingElement,
// // }

// // impl < 'a, Vertex, RingOperator, RingElement >

// //     NerveCoboundaryIter
// //         < 'a, Vertex, RingOperator, RingElement >
// // {
// //     /// Create an iterator that represents the coboundary of a simplex in the nerve complex.
// //     /// 
// //     /// The argument `cover` represents a family of sets.  The corresponding nerve is the 
// //     /// vertex-ordered combinatorial simplicial complex whose vertices are elements of `cover`,
// //     /// ordered by inclusion.  Since `cover` is a vector of hash-sets, we can represent each
// //     /// chain `u1 < .. < um` as a sequence of integers `n1, .., nm` such that `ui = cover[ni]`
// //     /// for all `i`.
// //     fn new( chain: Vec< usize >, cover: & Vec< HashSet< Vertex > >, ring_operator: RingOperator ) -> Self {
// //         chain.reserve_exact(1); // ensure there is space for one more element
// //         chain.shrink_to( chain.len() + 1); // ensure there is no excess capacity
// //         NerveCoboundaryIter{ chain, cover, ring_operator, next_openid_to_try: 0, phantom_ringelement: PhantomData }
// //     }
// // }            

// // impl < 'a, Vertex, RingOperator, RingElement >

// //     Iterator for

// //     NerveCoboundaryIter
// //         < 'a, Vertex, RingOperator, RingElement >
// // {
// //     type Item = ( Vec< usize >, RingElement );

// //     fn next( &mut self ) -> Option< Self::Item > {

// //         println!("NOTE: THIS NERVE ITERATOR DOES NOT RETURN ELEMENTS IN ASCENDING ORDER; CAN UPDATE THE ITERATOR ALGORITHM TO ACHIEVE THIS, HOWEVER, BY PRECOMPUTING THE POSET AND STORING ANCESTORS/CHILDREN IN ORDER");

// //         let cover_cardinality   =   self.cover.len();
// //         let facet_cardinality   =   self.chain.len();        
// //         let mut insertion_locus =   0;

// //         while self.next_openid_to_try < cover_cardinality {

// //             // skip elements that already belong to the chain
// //             if self.chain.iter().any(|x| x == self.next_openid_to_try ) {
// //                 self.next_openid_to_try += 1;
// //                 continue;
// //             }
            
// //             insertion_locus     =   0;

// //             let next_open_to_try =   self.cover[ self.next_openid_to_try ];


// //             while   ( insertion_locus < facet_cardinality )
// //                     &&
// //                     next_open_to_try.contains( self.cover[ self.facet[ insertion_locus ] ] ) {  
// //                         insertion_locus += 1
// //                     }
            
// //             if  insertion_locus  ==  facet_cardinality
// //                 ||
// //                 ( insertion_locus < facet_cardinality ) && self.cover[ self.facet[ insertion_locus ] ].contains( next_open_to_try ) 
// //                 {
// //                     let mut cofacet = facet.clone();
// //                     cofacet.insert( insertion_locus, self.next_openid_to_try );
// //                     let coefficient     =   self.ring_operator.minus_one_to_power( insertion_locus );
// //                     let return_value = Some( cofacet, coefficient );

// //                     self.next_openid_to_try += 1;
// //                     return return_value
// //                 }
// //             self.next_openid_to_try += 1;

// //         }

// //         return None
// //     }
// // }






// //  POSET NERVES
// //  ==================================================

// pub struct NerveCoboundaryIter< 
//                     Vertex, 
//                     UpwardClosureFn,
//                     DnwardClosureFn,
//                     ClosureIter,                    
//                     RingOperator, 
//                     RingElement 
//                 > {
//     chain:                  Vec< Vertex >,
//     close_up:               UpwardClosureFn,
//     close_dn:               DnwardClosureFn,
//     insertion_locus:        usize,
//     verts_to_insert:        TwoTypeIterator<
//                                     OnlyDuplicates< Merge< ClosureIter > >,
//                                     close_iter,
//                                 >,
//     ring_operator:          RingOperator,
//     phantom_ringelement:    RingElement,
// }

// impl < 'a, Vertex, RingOperator, RingElement >

//     NerveCoboundaryIter
//         < 'a, Vertex, RingOperator, RingElement >
// {
//     /// Create an iterator that represents the coboundary of a simplex in the nerve complex.
//     /// 
//     /// The argument `cover` represents a family of sets.  The corresponding nerve is the 
//     /// vertex-ordered combinatorial simplicial complex whose vertices are elements of `cover`,
//     /// ordered by inclusion.  Since `cover` is a vector of hash-sets, we can represent each
//     /// chain `u1 < .. < um` as a sequence of integers `n1, .., nm` such that `ui = cover[ni]`
//     /// for all `i`.
//     fn new( chain: Vec< usize >, cover: & Vec< HashSet< Vertex > >, ring_operator: RingOperator ) -> Self {
//         chain.reserve_exact(1); // ensure there is space for one more element
//         chain.shrink_to( chain.len() + 1); // ensure there is no excess capacity
//         NerveCoboundaryIter{ chain, cover, ring_operator, next_openid_to_try: 0, phantom_ringelement: PhantomData }
//     }
// }            

// impl < 'a, Vertex, RingOperator, RingElement >

//     Iterator for

//     PosetNerveCoboundaryIter
//         < Vertex, RingOperator, RingElement >
// {
//     type Item = ( Vec< Vertex >, RingElement );

//     fn next( &mut self ) -> Option< Self::Item > {

//         //  A closure operator that we call each time we need to increase the
//         //  insertion locus.  The operator both increases the insertion locus
//         //  and updates the iterator of vertices to insert.
//         let increment_insertion_locus = || -> () {
//             // increment the insertion locus
//             self.insertion_locus += 1;
//             // update the iterator of vertices to insert
//             match self.insertion_locus.cmp( &self.chain.len() ) {
//                 // in this case 1 <= insertion_locus < chain.len()
//                 // in this case any vertex we insert has to be greater than the vertex that precedes it and less than the vertex that follows it
//                 Less => { 
//                     self.verts_to_insert = 
//                         TwoTypeIterator::Version1(
//                             OnlyDuplicates::new(
//                                 self.close_up( self.insertion_locus - 1 )
//                                     .merge( self.close_dn( self.insertion_locus ) )
//                             )
//                         );
//                 },
//                 // in this case insertion_locus == chain.len()
//                 // in this case we are appending an element to the very end of the chain, so it only has to be greater than the last element of the chain              
//                 Equal => { 
//                     self.verts_to_insert = 
//                         TwoTypeIterator::Version2(
//                             self.close_up( self.chain[ self.chain.len() -1 ] )
//                         ); 
//                 }
//                 // in this case insertion_locus > chain.len(); no vertex can be inserted at this index                              
//                 Greater => {}
//             }
//         };

//         //  If there exists a new vertex to insert in the insertion locus, insert it and
//         //  return the resulting entry.
//         //  Otherwise, increase the insertion locus and try again.
//         while self.insertion_locus <= self.chain.len() { // we can insert a vertex in any slot numbered in 0, .., self.chain.len()

//             match self.verts_to_insert{
//                 None => { increment_insertion_locus(); continue }
//                 Some( x ) => {
//                     let mut entry_index     =   self.chain.clone();
//                     entry_index.insert( x, self.insertion_locus );
//                     let entry_coefficient   =   self.ring_operator.minus_one_to_power( self.insertion_locus );
//                     return Some( (entry_index, entry_coefficient) )
//                 }
//             }            

//         }

//         // We reach this point only if the insertion locus is numbered too high; in that case we are done.
//         return None;

//     }
// }







// //  TESTS
// //  -----------------------------------------------------------


// #[cfg(test)]
// mod doc_test_drafts {
    

//     #[test] 
//     fn test_nerve_coboundary() {


//         let open_sets   =   vec![
//                                 vec![                   ],  // 0
//                                 vec![   0,              ],  // 1
//                                 vec![   0,  1,          ],  // 2
//                                 vec![       1,  2,      ],  // 3
//                                 vec![           2,  3,  ],  // 4
//                                 vec![   0,  1,  2,      ],  // 5
//                                 vec![           2,      ],  // 6                             
//                             ];
        
//         // set relations
//         //              <   0   <   1,2,3,4,5,6
//         //         0    <   1   <   2,5
//         //       0,1    <   2   <   5
//         //       0,6    <   3   <   5
//         //         6    <   4   < 
//         //   0,1,2,3    <   5   < 
//         //         0    <   6   <   3,4,5
        
//         let chains      =   vec![
//                                 vec![   0,      ],  // A
//                                 vec![   1,      ],  // B
//                                 vec![   2,      ],  // C
//                                 vec![   3,      ],  // D  
//                                 vec![   4,      ],  // E                                                                                          
//                                 vec![   0,  1,  ],  // F
//                                 vec![   1,  2,  ],  // G
//                                 vec![           ],  // H                                
//                             ];

//          let coboundaries    =  vec![
//                                     // coboundary of A
//                                     vec![ (vec![0,1],2), (vec![0,2],2), (vec![0,3],2), (vec![0,4],2), (vec![0,5],2), (vec![0,6],2), ],
//                                     // coboundary of B
//                                     vec![ (vec![0,1],1), (vec![1,2],2), (vec![1,5],2), ],
//                                     // coboundary of C
//                                     vec![ (vec![0,2],1), (vec![1,2],1), (vec![2,5],2), ],   
//                                     // coboundary of D
//                                     vec![ (vec![0,3],1), (vec![6,3],1), (vec![3,5],2), ],
//                                     // coboundary of E
//                                     vec![ (vec![6,4],1), ],                                    
//                                     // coboundary of F
//                                     vec![ (vec![0,1,2],1), (vec![0,1,5],1), ],   
//                                     // coboundary of G
//                                     vec![ (vec![0,1,2],1), (vec![1,2,5],1), ],                                                                      
//                                     // coboundary of H
//                                     vec![ (vec![0],1), (vec![1],1), (vec![2],1), (vec![3],1), (vec![4],1), (vec![5],1), (vec![6],1), ],                                    
//                                 ];
//     }  
// }