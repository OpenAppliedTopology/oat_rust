//!  Data structures for unfiltered simplices
//! 
//! # Advisory: consider ordinary Rust vectors
//! 





use ndarray::Order;
use pyo3::IntoPyObject;

use crate::algebra::rings::traits::RingOperations;
use crate::utilities::order::{is_sorted_strictly, OrderOperatorAuto};
use crate::utilities::sequences_and_ordinals::{SortedVec};

use std::cmp::Ordering;

use std::hash::Hash;



use std::fmt::Debug;



//  ===========================================================================
//  COMBINATORIAL SIMPLEX (UNWEIGHTED) -- DEFINITION
//  ===========================================================================


/// An unweighted simplex
/// 
/// The vertices should sorted in ascending order, however this condition is not
/// gauranteed.  Compare with [SimplexSafe], where strict ordering *is* gauranteed.
/// 
/// # Advisory
/// 
/// Plain Rust vectors can serve quite well to represent simplices, OAT provides
/// a substantial tool set to work with vector representations.  Using this more
/// primitive data format may simplifiy your code, and provide access to a 
/// broader tool set.
#[derive(Debug, PartialEq, Clone, Eq, Hash, IntoPyObject)]
pub struct Simplex< Vertex > 
{ 
    pub vertices: Vec< Vertex >     //  vertices should be sorted in ascending order
} 

impl    < Vertex > 
        Simplex
        < Vertex >   
        {
    
    pub fn number_of_vertices( &self ) -> usize { self.vertices.len() }
    pub fn dimension( &self ) -> usize { self.vertices.len() - 1 }
}        


impl    < Vertex >           
        PartialOrd for Simplex
        < Vertex >

    where   Vertex: Ord     {

    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl    < Vertex >   
        Ord for Simplex
        < Vertex >

    where Vertex: Ord   {

    fn cmp(&self, other: &Self) -> Ordering {

        // next compare simplex dimensions
        let comp = self.number_of_vertices().cmp( & other.vertices.len() );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        self.vertices.cmp( & other.vertices )
    }
}

impl    < Vertex >   
        IntoIterator for Simplex
        < Vertex >      {

    type Item = Vertex;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter { self.vertices.into_iter() }
}



//  ================================================================================
//  COMBINATORIAL SIMPLEX (UNWEIGHTED) -- DEFINITION
//  ================================================================================


/// An unweighted simplex; the vertices should sorted in ascending order.
/// 
/// This is simply a wrapper around a `SortedVec`.  However, this wrapper has a special
/// order (first by length, second by lexicographic order of elements) which is useful
/// in homology computations.
/// 
/// # Advisory
/// 
/// Plain Rust vectors can serve quite well to represent simplices, OAT provides
/// a substantial tool set to work with vector representations.  Using this more
/// primitive data format may simplifiy your code, and provide access to a 
/// broader tool set.  Moreover, the safety gaurantees associated with this data
/// struture can often be replicated by simpler means.
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct SimplexSafe< Vertex: Ord + PartialEq + Eq > 
{ 
    vertices: SortedVec< Vertex >     //  vertices should be sorted in ascending order
} 

impl    < Vertex: Ord > 
        SimplexSafe
        < Vertex >   
        {
    
    /// Construct a new simplex
    pub fn new( vertices: SortedVec< Vertex > ) -> Self { SimplexSafe { vertices } }
    /// Construct from a vector; throws an error if the vector is not sorted in striclty ascending order
    pub fn from_vec( vertices: Vec< Vertex > ) -> Result< Self, Vec<Vertex > > where Vertex: Ord {
        SortedVec::new( vertices ).map(|x| SimplexSafe{ vertices: x } )        
    }        
    /// Number of vertices
    pub fn number_of_vertices( &self ) -> usize { self.vertices.len() }
    /// Dimension of the simplex
    pub fn dimension( &self ) -> usize { self.vertices.len() - 1 }
    /// Immutable reference to the internally stored sorted vector
    pub fn sorted_vec( &self ) -> & SortedVec< Vertex > { & self.vertices }
    /// Immutable reference to the internally stored vector 
    pub fn vertices( &self ) -> & Vec< Vertex > { self.vertices.vec() }    
    /// Consumes the simplex, returning the sorted vector of vertices
    pub fn into_sorted_vec( self ) -> SortedVec< Vertex > { self.vertices }    
    /// Consumes the simplex, returning the vector of vertices
    pub fn into_vec( self ) -> Vec< Vertex > { self.vertices.into_vec() }        
}        


impl    < Vertex >           
        PartialOrd for SimplexSafe
        < Vertex >

    where   Vertex: Ord     {

    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl    < Vertex >   
        Ord for SimplexSafe
        < Vertex >

    where Vertex: Ord   {

    fn cmp(&self, other: &Self) -> Ordering {

        // next compare simplex dimensions
        let comp = self.number_of_vertices().cmp( & other.number_of_vertices() );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        return self.vertices().cmp( other.vertices() )
    }
}

impl    < Vertex: Ord >   
        IntoIterator for SimplexSafe
        < Vertex >      {

    type Item = Vertex;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter { self.into_vec().into_iter() }
}




/// Returns the index where a vertex can be inserted into a strictly sorted list to produce a new strictly sorted list.
/// 
/// # Arguments
/// 
/// * `facet` - a vector of vertices [v0 < .. < vn] representing the simplex, which should be sorted in strictly ascending order.
/// * `new_vertex` - the vertex v to be inserted into the simplex.
/// 
///  # Returns
/// 
/// * `Ok(k)`, where `k` is the index such that [v0 < .. < v < vk < .. < vn] is a strictly sorted list.
/// * `Err((Some(facet), None))` if the facet is not strictly sorted.
/// * `Err((None, Some(new_vertex)))` if the vertex already belongs to the facet.
///
/// # Example
/// 
///  ```
///  use oat_rust::topology::simplicial::simplices::unweighted::cofacet_vertex_insertion_locus;
/// 
///  let facet = vec![0, 2, 2];
///  let new_vertex = 2;
///  let result = cofacet_vertex_insertion_locus(&facet, new_vertex);
///  assert_eq!(result, Err((Some(&facet), Some(2)))); // facet is not strictly sorted, and vertex already exists
///  let new_vertex = 3;
///  let result = cofacet_vertex_insertion_locus(&facet, new_vertex);
///  assert_eq!(result, Err((Some(&facet), None))); // facet is not strictly sorted
/// 
///  let facet = vec![1, 2, 3];
///  let new_vertex = 2;
///  let result = cofacet_vertex_insertion_locus(&facet, new_vertex);
///  assert_eq!(result, Err((None, Some(2)))); // vertex already exists in the facet
///  let new_vertex = 4;
///  let result = cofacet_vertex_insertion_locus(&facet, new_vertex);
///  assert_eq!(result, Ok(3)); // vertex can be inserted at index 3
///  let new_vertex = 0;
///  let result = cofacet_vertex_insertion_locus(&facet, new_vertex);
///  assert_eq!(result, Ok(0)); // vertex can be inserted at index 0
///  ```
pub fn cofacet_vertex_insertion_locus< Vertex >(
        facet:              & Vec< Vertex >,
        new_vertex:         Vertex,
    )
    -> Result< usize, (Option< &Vec<Vertex> >, Option< Vertex > ) >

    where 
        Vertex:           Copy + Ord,
{
    // Return an error if the facet is not strictly sorted or the vertex already belongs to the facet.
    let facet_error = if is_sorted_strictly(facet, &OrderOperatorAuto) {
        None
    } else {
        Some(facet)
    };

    let vertex_error = if facet.contains(&new_vertex) {
        Some(new_vertex)
    } else {
        None
    };

    if facet_error.is_some() || vertex_error.is_some() {
        return Err( (facet_error, vertex_error) );
    }
    
    // otherwise compute the insertion locus
    let mut insertion_locus = facet.len();
    for p in 0 .. facet.len() {
        if facet[p] > new_vertex {
            insertion_locus = p;
            break;
        }
    }

    return Ok(insertion_locus)

}        



/// Returns the entry of the boundary matrix for the cofacet obtained by inserting a vertex into a facet.
/// 
/// Concretely, if the cofacet has form `[v0 < v1 < ... < vn]` and the new vertex is `vk`, then 
/// this function returns `([v0 < .. < vn], (-1)^(k))`
/// 
/// # Arguments
/// 
/// * `facet` - a vector of vertices representing the facet, which should be sorted in strictly ascending order.
/// * `new_vertex` - the vertex to be inserted into the facet.
/// * `ring_operator` - a ring operator that defines the coefficient field for the simplicial complex.
/// 
/// # Returns
/// 
/// * `Ok((cofacet,coefficient))` if the inputs are valid
/// * `Err((Some(facet), None))` if the facet is not strictly sorted.
/// * `Err((None, Some(new_vertex)))` if the vertex already belongs to the facet.
/// 
/// # Example
/// ```
/// use oat_rust::topology::simplicial::simplices::unweighted::coboundary_entry_for_facet_vertex_pair;
/// use oat_rust::algebra::rings::types::native::RingIsize;
/// 
/// let ring_operator = RingIsize::new();
/// 
/// 
/// let facet = vec![1, 2, 3];
/// let new_vertex = 4;
/// assert_eq!(
///     coboundary_entry_for_facet_vertex_pair(&facet, new_vertex, ring_operator),
///     Ok((vec![1, 2, 3, 4], -1))
/// );
/// 
/// let facet = vec![1, 2, 3];
/// let new_vertex = 0;
/// assert_eq!(
///     coboundary_entry_for_facet_vertex_pair(&facet, new_vertex, ring_operator),
///     Ok((vec![0, 1, 2, 3], 1))
/// );    
/// 
/// let facet = vec![1, 2, 3];
/// let new_vertex = 2;
/// assert_eq!(
///     coboundary_entry_for_facet_vertex_pair(&facet, new_vertex, ring_operator),
///     Err((None, Some(2))) // vertex already belongs to the facet
/// );    
/// 
/// let facet = vec![1, 1, 3];
/// let new_vertex = 2;
/// assert_eq!(
///     coboundary_entry_for_facet_vertex_pair(&facet, new_vertex, ring_operator),
///     Err((Some(&facet), None)) // facet is not strictly sorted
/// );        
/// 
/// let facet = vec![1, 1, 3];
/// let new_vertex = 1;
/// assert_eq!(
///     coboundary_entry_for_facet_vertex_pair(&facet, new_vertex, ring_operator),
///     Err((Some(&facet), Some(1))) // facet is not strictly sorted, and contains the vertex
/// );  
/// ```
pub fn coboundary_entry_for_facet_vertex_pair< Vertex, RingOperator >(
        facet:              & Vec< Vertex >,
        new_vertex:         Vertex,
        ring_operator:      RingOperator,
    ) -> Result< 
            ( Vec<Vertex>, RingOperator::Element  ),
            (Option< &Vec<Vertex> >, Option< Vertex > ) 
        >

    where 
        Vertex:           Copy + Ord,
        RingOperator:     RingOperations,
{        
    // get the insertion locus or return an error if the input is invalid
    let insertion_locus = cofacet_vertex_insertion_locus( facet, new_vertex )?;
    
    
    // costruct the cofacet
    let mut cofacet = Vec::with_capacity( facet.len() + 1 );
    cofacet.extend_from_slice( &facet[0..insertion_locus] );
    cofacet.push( new_vertex );
    cofacet.extend_from_slice( &facet[insertion_locus..] );


    // compute the coefficient
    let coefficient = ring_operator.minus_one_to_power( insertion_locus );
    
    // return
    return Ok(( cofacet, coefficient ))
}