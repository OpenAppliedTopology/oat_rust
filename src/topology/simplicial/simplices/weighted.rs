//! Tools to define and enumerate filtered and unfiltered simplices
//! 
//! In this library,
//! - filtered simplices (meaning simplices with assigned birth times) are represented by the [`WeightedSimplex`] struct
//! - unfiltered simplices (meaning simplices without birth times) are represented by vectors with *strictly ascending entries*

use derive_getters::Dissolve;
use derive_new::new;
use itertools::Itertools;
use num::Bounded;
use ordered_float::OrderedFloat;
use pyo3::types::{PyDict, PyDictMethods};
use pyo3::{pyclass, IntoPyObject, Py, PyAny};



use crate::algebra::matrices::query::MatrixOracle;
use crate::algebra::rings::traits::RingOperations;
use crate::algebra::vectors::entries::KeyValNew;
use crate::topology::simplicial::from::graph_weighted::{AgileBoundaryIteratorLexicographicOrder, AgileCoboundaryIterartorLexicographicOrder, VietorisRipsComplex};
use crate::utilities::order::{JudgeOrder, JudgePartialOrder};


use std::cmp::Ordering;

use std::hash::Hash;



use std::fmt::Debug;
use std::sync::Arc;


type Vertex = u16;




//  ================================================================================
//  TWIST ORDER OPERATOR - FILTERED SIMPLEX
//  ================================================================================


/// Orders simplices first by dimension (descending) then by filtration (ascending)
/// then by lexicographic order (ascending)
/// 
/// This order is specially designed to facilitate the `twist` optimization in 
/// persistent homology calculations.
/// 
/// This object implements the `JudgeOrder` and `JudgePartialOrder` traits.
#[derive(Clone,Debug,new)]
pub struct OrderOperatorTwistWeightedSimplex;

impl < FilVal: Ord + Clone + Debug >

    JudgePartialOrder
        < WeightedSimplex< FilVal > > for 

    OrderOperatorTwistWeightedSimplex 
{
    /// Orders simplices first by dimension (descending) then by filtration (ascending)
    /// then by lexicographic order (ascending)
    /// 
    /// This order is specially designed to facilitate the `twist` optimization in 
    /// persistent homology calculations.    
    fn judge_partial_cmp( & self, lhs: &WeightedSimplex< FilVal >, rhs: &WeightedSimplex< FilVal > ) -> Option<Ordering> {
        Some(self.judge_cmp(lhs,rhs))
    }
}

impl < FilVal: Ord + Clone + Debug >

    JudgeOrder
        < WeightedSimplex< FilVal > > for 

    OrderOperatorTwistWeightedSimplex 
{

    /// Orders simplices first by dimension (descending) then by filtration (ascending)
    /// then by lexicographic order (ascending)
    /// 
    /// This order is specially designed to facilitate the `twist` optimization in 
    /// persistent homology calculations.    
    fn judge_cmp( & self, lhs: &WeightedSimplex< FilVal >, rhs: &WeightedSimplex< FilVal > ) -> Ordering {
        // compare simplex dimensions
        let mut comp = rhs.vertices.len().cmp( & lhs.vertices.len() ); // note we put RHS on the LEFT here
        if comp != Ordering::Equal { return comp }

        // first compare filtration values
        comp = lhs.weight.cmp( & rhs.weight );
        if comp != Ordering::Equal { return comp }        

        // finally, compare simplices lexicographically
        lhs.vertices.cmp( & rhs.vertices )
    }
}







//  ================================================================================
//  COMBINATORIAL SIMPLEX (WEIGHTED) -- DEFINITION
//  ================================================================================



/// A weighted simplex
/// 
/// Concretely, this corresponds to a pair `([v0,..,vn],w)` where
/// - `[v0,..,vn]` is a list of vertices (which should be sorted in strictly ascending order)
/// - w is the weight
#[derive(Debug, PartialEq, Clone, Eq, Hash, Dissolve, new)]
pub struct WeightedSimplex<FilVal: Clone + Debug>{
    pub weight:     FilVal,     // the filtration value of simplex
    pub vertices:   Vec<Vertex> // a sorted vector in strictly descending order recording vertices
}

impl <FilVal: Clone + Debug>
    WeightedSimplex<FilVal> {

    /// Number of vertices.
    pub fn number_of_vertices(&self) -> usize { self.vertices.len() }

    /// Dimension of the simplex.
    pub fn dimension(&self) -> usize { self.vertices.len() - 1 }

    /// Weight of the simplex.
    pub fn weight(&self) -> FilVal { self.weight.clone() }    

    /// Filtraiton value of the simplex.
    /// 
    /// This is identical to the weight. We provide the same function with
    /// two different names for convenience, because in the overwhelming
    /// majority of cases, the filtration value is the same as the weight,
    /// and `x.filtraion()` will be more understandable than `x.weight()`
    /// for the average user.
    pub fn filtration(&self) -> FilVal { self.weight.clone() }

    /// The vertices of the simplex.
    pub fn vertices(&self) -> &Vec<Vertex> { &self.vertices }

    /// Returns entries of the boundary matrix column for this simplex, in ascending lexicographic order (filtration values are not sorted)
    /// 
    /// # Edge case 
    /// 
    /// If `self` is a 0-simplex, (i.e. contains only one vertex) then an emtpy iterator is returned.
    /// However, the object `AgileBoundaryIteratorLexicographicOrder` is capable of returning an iterator
    /// that runs over the -1 simplex (i.e. emtpy siplex) as the boundary of a 0-simplex; see the documentation
    /// for that object for more details.
    pub fn agile_boundary_iterator_lexicographic_order
                < DissimilarityMatrix, RingOperator >
                ( 
                    self, 
                    vietoris_rips_complex:  Arc< VietorisRipsComplex< DissimilarityMatrix, RingOperator > >,
                ) 
        -> 
        AgileBoundaryIteratorLexicographicOrder< DissimilarityMatrix, RingOperator >
        where
            DissimilarityMatrix:                        Clone + MatrixOracle<
                                                            ColumnIndex     =   usize,
                                                            RowIndex        =   usize,
                                                            Coefficient     =   FilVal, 
                                                            RowEntry:           KeyValNew,
                                                            Row:                Clone,
                                                        >,
            DissimilarityMatrix::Coefficient:           Ord + Bounded,            
            RingOperator:                               Clone + RingOperations, 
            FilVal:                                     Copy,
    {
        AgileBoundaryIteratorLexicographicOrder::new(
            vietoris_rips_complex,
            self,
            false, // don't allow the empty simplex to be returned as a boundary
        )
    }  

    /// Returns the coboundary of `self` as a vector **with entries sorted in ascending lexicograph order (filtration values are not sorted)**
    pub fn agile_coboundary_iterator_lexicographic_order
                < DissimilarityMatrix, RingOperator >
                ( 
                    self, 
                    vietoris_rips_complex:  Arc< VietorisRipsComplex< DissimilarityMatrix, RingOperator > >,
                ) 
        -> 
            Result< 
                AgileCoboundaryIterartorLexicographicOrder< DissimilarityMatrix, RingOperator >,
                Vec<Vertex>,
            >
        where
            DissimilarityMatrix:                        Clone + MatrixOracle<
                                                            ColumnIndex     =   usize,
                                                            RowIndex        =   usize,
                                                            Coefficient     =   FilVal, 
                                                            RowEntry:           KeyValNew,
                                                            Row:                Clone,
                                                        >,
            DissimilarityMatrix::Coefficient:           Ord + Bounded,            
            RingOperator:                               Clone + RingOperations, 
            FilVal:                                     Copy,
    {
        let iterator =  AgileCoboundaryIterartorLexicographicOrder::new(
            self.vertices().clone(),
            vietoris_rips_complex.dissimilarity_matrix.clone(),
            vietoris_rips_complex.ring_operator(),
        )?;

        if iterator.facet_filtration != self.filtration() {
            println!("Error: the weight of this simplex does not match its filtration value in the Vietoris-Rips complex.");
            return Err( self.vertices.clone() );
        }

        Ok(iterator)
    }

}

impl <FilVal>
    PartialOrd for WeightedSimplex<FilVal> 
    where FilVal: Eq + PartialEq + Ord + PartialOrd + Clone + Debug

{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl <FilVal>
    Ord for WeightedSimplex<FilVal> 
    where FilVal: Eq + PartialEq + Ord + PartialOrd + Clone + Debug
{
    /// Compares simplices first by filtration value, and second by lexicographic order.
    fn cmp(&self, other: &Self) -> Ordering {

        // first compare filtration values
        let comp = self.weight.cmp( & other.weight );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        self.vertices.cmp( & other.vertices )
    }
}