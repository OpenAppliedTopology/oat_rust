//  TEMPORARILY DEPRECATED; HOWEVER, IT SEEMS TO BE WORTHWHILE TO BRING THIS UP TO DATE AND USE IT, FOR UNFILTERED GRAPHS

//! The Clique complex (aka Flag complex) of a simple combinatorial graph

use itertools::Itertools;

use crate::algebra::matrices::query::MatrixOracle;
use crate::algebra::rings::traits::{Semiring, Ring};
use crate::utilities::iterators::general::IntersectOrderedIteratorsUnsafe;

use std::collections::HashSet;
use std::hash::Hash;
use std::iter::FromIterator;
use std::marker::PhantomData;

use crate::topology::simplicial_complex::boundary::{SimplexBoundaryAscend, SimplexBoundaryDescend};

// use super::boundary_matrices::{SimplexBoundaryAscend, SimplexBoundaryDescend};




// pub struct FlagComplexAgileSimplexIteratorLexicographicOrderAscend{ 
//     proximity_matrix:   & Vec<HashSet<usize>>,
//     cardinality:        usize,
//     current_face:       Vec<usize>,
// }

// impl    Iterator for 
//         FlagComplexAgileSimplexIteratorLexicographicOrderAscend {

//     type Item = Vec< usize >;

//     fn next(&mut self) -> Option<Self::Item> {
//         let last_index = self.cardinality - 1;
//         let last_vertex     =   self.current_face[ last_index ];
//         for next_option in last_vertex + 1 .. self.proximity_matrix.len() {
//             if self.current_face.iter().map( |x| self.proximity_matrix.contains(x) ).take( last_index ).all() {
//                 self.current_face[ last_index ] = next_option;
//                 return Some(self.current_face.clone());
//             }
//         }
        
//     }
// }





//  COUBOUNDARY VECTOR DESCENDING
//  -----------------------------------------------------------------------------------



/// Iterates over the terms of the coboundary of a simplex in descending lexicographic order.
/// 
/// The coefficient is calculated as `(-1)^k`, where `xk` is the vertex which appears twice in the template.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::from::relation::DowkerBoundaryMatrixRowReverse;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use std::collections::HashSet;
/// use std::iter::FromIterator;
/// 
/// let ring_operator = PrimeOrderField::new(3);
/// let simplex = vec![1,3];
/// let dowker_simplices = vec![ HashSet::from_iter( vec![0,1,2,3,4] ) ];
/// let coboundary = DowkerBoundaryMatrixRowReverse::from_proximity_matrix( simplex, &dowker_simplices, ring_operator );
/// 
/// itertools::assert_equal( coboundary, vec![ (vec![1,3,4], 1), (vec![1,2,3], 2), (vec![0,1,3], 1) ]);
/// ```
#[derive(Clone, Debug)]
pub struct DowkerBoundaryMatrixRowReverse< Vertex, RingOperator, RingElement >
    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,
{
    next_cofacet_opt:           Option< Vec< Vertex > >,
    next_coefficient:           RingElement,
    vertices_to_insert:         Vec< Vertex >,              // should be sorted in ascending order
    retrieval_locus:            usize,                      // the place where the vertex inserted into `next_cofacet` was drawn from 
    insertion_locus:            usize,                      // this points to the place in `next_cofacet` where we inserted a vertex    
    ring_operator:              RingOperator,
}   


impl < RingOperator, RingElement >

DowkerBoundaryMatrixRowReverse
        < usize, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        // Vertex:             Clone + Hash + Ord,
{
    /// Generates a [DowkerBoundaryMatrixRow] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_proximity_matrix( facet: Vec< usize >, proximity_matrix: &Vec< Vec< usize > >, ring_operator: RingOperator ) -> Self {

        println!("NOTE: we assume that the proximity graph is just a sparsity pattern and DOES NOT CONTAIN the diagonal");
        println!("NOTE: THIS IS INEFFICIENT; CAN REPLACE vertices-to-insert with an iterator");
        println!("NOTE: consider merging the constructors for the ascending and descending ... or maybe not that; consider merging the constructor with the actual Dowker version ");

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if facet.is_empty() {
            // this iterator will be identified as empty, because `vertices_to_insert` is empty
            return DowkerBoundaryMatrixRowReverse{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),   // these are arbitrary values            
                vertices_to_insert:         vec![],                // these are arbitrary values  
                retrieval_locus:            0,                     // these are arbitrary values  
                insertion_locus:            0,                     // these are arbitrary values  
                ring_operator:              ring_operator,         // these are arbitrary values                    
            }                      
        }

        //  obtain a list of vertices to insert
        //  ---------------------------------------------  
        
        let iters = facet.iter().map(|x| proximity_matrix[*x].iter().cloned() );
        let vertices_to_insert = IntersectOrderedIteratorsUnsafe::new( iters ).collect_vec();
        

        //  if there are no vertices to insert, then the coboundary is zero
        //  ---------------------------------------------
        if vertices_to_insert.is_empty() {
            // this iterator is empty, because `next_cofacet_opt` is None
            return DowkerBoundaryMatrixRowReverse{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),     // these are abitrary values            
                vertices_to_insert:         vec![],                  // these are abitrary values  
                retrieval_locus:            0,                       // these are abitrary values  
                insertion_locus:            0,                       // these are abitrary values  
                ring_operator:              ring_operator,           // these are abitrary values                    
            }                      
        }


        //  generate the rest of the initial data
        //  ---------------------------------------------                

        let mut coefficient = ring_operator.minus_one_to_power( facet.len() );
        let mut next_cofacet = facet.clone();
        let mut insertion_locus = facet.len();
        let retrieval_locus = vertices_to_insert.len() - 1;        
        let inserted_vertex = vertices_to_insert[retrieval_locus].clone();
        while   ( insertion_locus > 0 ) 
                && 
                ( inserted_vertex < next_cofacet[ insertion_locus - 1 ] ) {
            insertion_locus -= 1;
            coefficient = ring_operator.negate( coefficient );
        }
        next_cofacet.insert( insertion_locus, inserted_vertex );

        DowkerBoundaryMatrixRowReverse{
            next_cofacet_opt:           Some( next_cofacet ),
            next_coefficient:           coefficient,
            vertices_to_insert:         vertices_to_insert,         // should be sorted in ascending order
            retrieval_locus:            retrieval_locus,            // the place where the inserted vertex was drawn from 
            insertion_locus:            insertion_locus,            // this points to the place in `next_cofacet` where we inserted a vertex
            ring_operator:              ring_operator,             
        }

    }

    // /// Generates a [DowkerBoundaryMatrixRow] for a simplex, given CSR and CSC representations of the binary dowker matrix.
    // fn from_csr_and_csc_matrices( facet: Vec< usize >, dowker_matrix_csr: Vec< HashSet< usize > >, dowker_matrix_csc: Vec< HashSet< usize > > ) {        
    // }
}


impl < Vertex, RingOperator, RingElement >

    Iterator for

    DowkerBoundaryMatrixRowReverse
        < Vertex, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,   
{
    type Item = (Vec<Vertex>, RingElement);

    fn next( &mut self ) -> Option< Self::Item >{

        // println!("{:?} -- DELETE THIS AND ALL DEBUG REQUIREMENTS ON THIS IMPLEMENTATION OF ITER", &self);

        match self.next_cofacet_opt {
            None => { return None }
            Some( ref mut next_cofacet ) => {

                // make a copy of the return value, before we change the internal state of the iterator
                let return_value    =   ( next_cofacet.clone(), self.next_coefficient.clone() );

                // if there are more vertices to insert, update the internal state with data for the next entry
                if self.retrieval_locus > 0 {
                    
                    // grab the next vertex
                    self.retrieval_locus -= 1; 
                    let inserted_vertex     =   self.vertices_to_insert[ self.retrieval_locus ].clone();
                    
                    // update pointers and coefficients
                    while   ( self.insertion_locus > 0 ) 
                            && 
                            ( inserted_vertex < next_cofacet[ self.insertion_locus - 1 ] ) {
                        next_cofacet[ self.insertion_locus ] = next_cofacet[ self.insertion_locus - 1 ].clone();
                        self.insertion_locus -= 1;
                        self.next_coefficient = self.ring_operator.negate( self.next_coefficient.clone() );
                    }
                    // update the cofacet
                    next_cofacet[ self.insertion_locus ] = inserted_vertex;
                } 
                // otherwise flag the fact that there will be no more entries, after we return the one that we have already extracted
                else {
                    self.next_cofacet_opt = None;
                }

                return Some( return_value )
            }
        }

    }
}       



//  COUBOUNDARY VECTOR ASCENDING
//  -----------------------------------------------------------------------------------

/// Iterates over the terms of the coboundary of a simplex in ascending lexicographic order.
/// 
/// The coefficient is calculated as `(-1)^k`, where `xk` is the vertex which appears twice in the template.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::simplicial::from::relation::DowkerBoundaryMatrixRow;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use std::collections::HashSet;
/// use std::iter::FromIterator;
/// 
/// let ring_operator = PrimeOrderField::new(3);
/// let simplex = vec![1,3];
/// let dowker_simplices = vec![ HashSet::from_iter( vec![0,1,2,3,4] ) ];
/// let coboundary = DowkerBoundaryMatrixRow::from_proximity_matrix( simplex, &dowker_simplices, ring_operator );
/// 
/// itertools::assert_equal( coboundary, vec![ (vec![0,1,3], 1), (vec![1,2,3], 2), (vec![1,3,4], 1) ]);
/// ```
#[derive(Clone, Debug)]
pub struct DowkerBoundaryMatrixRow< Vertex, RingOperator, RingElement >
    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,
{
    next_cofacet_opt:           Option< Vec< Vertex > >,
    next_coefficient:           RingElement,
    vertices_to_insert:         Vec< Vertex >,              // should be sorted in ascending order
    retrieval_locus:            usize,                      // the place where the vertex inserted into `next_cofacet` was drawn from 
    insertion_locus:            usize,                      // this points to the place in `next_cofacet` where we inserted a vertex    
    ring_operator:              RingOperator,
}   


impl < RingOperator, RingElement >

    DowkerBoundaryMatrixRow
        < usize, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    /// Generates a [DowkerBoundaryMatrixRow] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_proximity_matrix( facet: Vec< usize >, proximity_matrix: &Vec< Vec< usize > >, ring_operator: RingOperator ) -> Self {

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if facet.is_empty() {
            // this iterator will be identified as empty, because `vertices_to_insert` is empty
            return DowkerBoundaryMatrixRow{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),   // these are arbitrary values            
                vertices_to_insert:         vec![],                // these are arbitrary values  
                retrieval_locus:            0,                     // these are arbitrary values  
                insertion_locus:            0,                     // these are arbitrary values  
                ring_operator:              ring_operator,         // these are arbitrary values                    
            }                      
        }

        //  obtain a list of vertices to insert
        //  ---------------------------------------------        
        let iters = facet.iter().map(|x| proximity_matrix[*x].iter().cloned() );
        let vertices_to_insert = IntersectOrderedIteratorsUnsafe::new( iters ).collect_vec();
        

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if vertices_to_insert.is_empty() {
            // this iterator is empty, because `vertices_to_insert` is empty
            return DowkerBoundaryMatrixRow{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),
                vertices_to_insert:         vec![],     
                retrieval_locus:            0,          
                insertion_locus:            0,          
                ring_operator:              ring_operator,                
            }                      
        }


        //  generate the rest of the initial data
        //  ---------------------------------------------                

        let mut coefficient = RingOperator::one();
        let mut next_cofacet = facet.clone();
        let mut insertion_locus = 0;
        let retrieval_locus = 0;        
        let inserted_vertex = vertices_to_insert[retrieval_locus].clone();
        while   ( insertion_locus < next_cofacet.len() ) 
                && 
                ( inserted_vertex > next_cofacet[ insertion_locus ] ) {
            insertion_locus += 1;
            coefficient = ring_operator.negate( coefficient );
        }
        next_cofacet.insert( insertion_locus, inserted_vertex );

        DowkerBoundaryMatrixRow{
            next_cofacet_opt:           Some( next_cofacet ),
            next_coefficient:           coefficient,
            vertices_to_insert:         vertices_to_insert,         // should be sorted in ascending order
            retrieval_locus:            retrieval_locus,            // the place where the inserted vertex was drawn from 
            insertion_locus:            insertion_locus,            // this points to the place in `next_cofacet` where we inserted a vertex
            ring_operator:              ring_operator,             
        }

    }

    // /// Generates a [DowkerBoundaryMatrixRow] for a simplex, given CSR and CSC representations of the binary dowker matrix.
    // fn from_csr_and_csc_matrices( facet: Vec< usize >, dowker_matrix_csr: Vec< HashSet< usize > >, dowker_matrix_csc: Vec< HashSet< usize > > ) {        
    // }
}


impl < Vertex, RingOperator, RingElement >

    Iterator for

    DowkerBoundaryMatrixRow
        < Vertex, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,   
{
    type Item = (Vec<Vertex>, RingElement);

    fn next( &mut self ) -> Option< Self::Item >{

        // println!("{:?} -- DELETE THIS AND ALL DEBUG REQUIREMENTS ON THIS IMPLEMENTATION OF ITER", &self);

        match self.next_cofacet_opt {
            None => { return None }
            Some( ref mut next_cofacet ) => {

                // make a copy of the return value, before we change the internal state of the iterator
                let return_value    =   ( next_cofacet.clone(), self.next_coefficient.clone() );

                // if there are more vertices to insert, update the internal state with data for the next entry
                if self.retrieval_locus + 1 < self.vertices_to_insert.len() {
                    
                    // grab the next vertex
                    self.retrieval_locus += 1; 
                    let inserted_vertex     =   self.vertices_to_insert[ self.retrieval_locus ].clone();
                    
                    // update pointers and coefficients
                    while   ( self.insertion_locus + 1 < next_cofacet.len() ) 
                            && 
                            ( inserted_vertex > next_cofacet[ self.insertion_locus + 1 ] ) {
                        next_cofacet[ self.insertion_locus ] = next_cofacet[ self.insertion_locus +1 ].clone();
                        self.insertion_locus += 1;
                        self.next_coefficient = self.ring_operator.negate( self.next_coefficient.clone() );
                    }
                    // update the cofacet
                    next_cofacet[ self.insertion_locus ] = inserted_vertex;
                } 
                // otherwise flag the fact that there will be no more entries, after we return the one that we have already extracted
                else {
                    self.next_cofacet_opt = None;
                }

                return Some( return_value )
            }
        }

    }
}    












//  ---------------------------------------------------------------------------
//  DOWKER BOUNDARY MATRIX
//  ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct FlagComplexBoundaryMatrixRowMajor
            < Vertex, RingOperator, RingElement >
    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    proximity_matrix:       Vec< Vec< Vertex > >,
    ring_operator:          RingOperator,
    phantom_ringelement:    PhantomData< RingElement >,
}


impl < Vertex, RingOperator, RingElement >
    
    FlagComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >
    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    pub fn new( proximity_matrix: Vec< Vec< Vertex > >, ring_operator: RingOperator ) -> Self {
        FlagComplexBoundaryMatrixRowMajor{ proximity_matrix, ring_operator, phantom_ringelement: PhantomData }
    }
}


//  INDICES AND COEFFICIENTS
//  ------------------------------------------


impl < Vertex, RingOperator, RingElement >
    
    IndicesAndCoefficients for

    FlagComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type RowIndex = Vec< Vertex >; type ColumnIndex = Vec< Vertex >; type Coefficient = RingElement;
}


//  ORACLE MAJOR ASCEND
//  ------------------------------------------

impl < RingOperator, RingElement >
    
    MatrixOracle for

    FlagComplexBoundaryMatrixRowMajor
        < usize, RingOperator, RingElement >

    where
        // Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type Row            =   DowkerBoundaryMatrixRow< usize, RingOperator, RingElement >;
    type RowIntoIter    =   Self::Row;
    type RowEntry =   ( Self::RowIndex, Self::Coefficient );

    fn row( &self, row_index: Self::RowIndex ) -> Self::Row {
        DowkerBoundaryMatrixRow::from_proximity_matrix( row_index, & self.proximity_matrix, self.ring_operator.clone() )
    }
}

//  ORACLE MAJOR DESCEND
//  ------------------------------------------

impl < RingOperator, RingElement >
    
ViewRowDescend for

    FlagComplexBoundaryMatrixRowMajor
        < usize, RingOperator, RingElement >

    where
        // Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type RowReverse           =   DowkerBoundaryMatrixRowReverse< usize, RingOperator, RingElement >;
    type RowReverseIntoIter   =   Self::RowReverse;
    type RowEntry      =   ( Self::RowIndex, Self::Coefficient );

    fn row_reverse( &self, row_index: Self::RowIndex ) -> Self::RowReverse {

        println!("GOT THROUGH THIS FILE ANE REMOVE THE DEBUG REQUIREMENTS");
        DowkerBoundaryMatrixRowReverse::from_proximity_matrix( row_index, & self.proximity_matrix, self.ring_operator.clone() )
    }
}



//  ORACLE MINOR ASCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
    ViewColAscend for

    FlagComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type Column            =   SimplexBoundaryAscend< Vertex, RingOperator, RingElement >;
    type ColumnIntoIter    =   Self::Column;
    type ColumnEntry       =   ( Self::RowIndex, Self::Coefficient );

    fn row( &self, row_index: Self::RowIndex ) -> Self::Column {
        SimplexBoundaryAscend::new( row_index, self.ring_operator.clone() )
    }
}

//  ORACLE MINOR DESCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
    MatrixOracle for

    FlagComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ColumnReverse           =   SimplexBoundaryDescend< Vertex, RingOperator, RingElement >;
    type ColumnReverseIntoIter   =   Self::ColumnReverse;
    type ColumnEntry      =   ( Self::RowIndex, Self::Coefficient );

    fn column_reverse( &self, row_index: Self::RowIndex ) -> Self::ColumnReverse {
        SimplexBoundaryDescend::new( row_index, self.ring_operator.clone() )
    }
}













#[cfg(test)]
mod tests {
    


}    