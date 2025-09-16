//! The filtered Clique complex (aka filtered Flag complex, filtered Vietoris-Rips complex); often constructed from a metric space or point cloud
//!
//! 
//! # Clique chain complex oracles
//! Mathematically, a filtered clique complex is determined by two pieces of data:
//! * a **dissimilarity matrix**
//!     * we represent this by a symmetric matrix *S* that has been flattened into a vector *v*
//! * a **threshold**  that determines where we stop growing the filtration
//!     * we represent this as a real number *t*
//! The boundary matrices of this chain complex oracle are indexed by `WeightedSimplex` objects which are
//! (essentially) pairs of form `(simplex, filtration_value)`.
//! 
//! Example:
//!     - construct a dissimilarity matrix
//!     - construct the associated clique complex
//!     - access a row + column of the boundary matrix
//!     - compute the dimension 1 barcode





use core::panic;
use std::collections::binary_heap::Iter;
use std::iter::Peekable;
use std::sync::Arc;

use std::hash::Hash;
use std::fmt::Debug;
use std::vec::IntoIter;

use derive_new::new;
use itertools::Itertools;
use num::iter::Range;
use num_traits::Bounded;
use ordered_float::OrderedFloat;
use sprs::vec;

use crate::algebra::chain_complexes::{ChainComplex, FilteredChainComplex};
use crate::algebra::vectors::operations::ChangeEntryType;
use crate::topology::simplicial::simplices::unweighted::{coboundary_entry_for_facet_vertex_pair, Simplex};
use crate::topology::simplicial::simplices::vector::insert_vertex;
use crate::topology::simplicial::simplices::weighted::{WeightedSimplex, OrderOperatorTwistWeightedSimplex};

use crate::algebra::vectors::entries::{KeyValGet, KeyValNew};


// use crate::matrix::{SmOracle, RingMetadata, MajorDimension};
// use crate::chx::{ChainComplex, ChxTransformKind};


use crate::algebra::matrices::operations::umatch::differential::DifferentialUmatch;
use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::matrices::query::{MatrixAlgebra, MatrixOracle};
use crate::algebra::rings::traits::DivisionRingOperations;
use crate::algebra::rings::traits::{RingOperations};

use crate::utilities::order::{ is_sorted_strictly, JudgeOrder, MakeNoneMaximum, OrderOperatorAuto, OrderOperatorByKey};
use crate::utilities::iterators::general::{minmax, symmetric_difference_of_ordered_iterators, TwoTypeIterator, IterWrappedArcVec, PeekUnqualified};
// use oat_rust::utilities::indexing_and_bijection::sparsevec_lookup;



// The specific type we use to store vertices
type Vertex = u16;
                                        










//  ===========================================================
//  FUNCTION TO EVALUATE DIAMETER
//  ===========================================================


/// Returns the filtration value of a simplex in the Vietoris-Rips complex, or `Err( vertices )`, if the input is invalid.
///
/// # Parameters
/// 
/// -`vertices`: vertices of the simplex
/// 
/// # Returns
/// 
/// - `Err( vertices.clone() )` if either
///   - `vertices` is not sorted in strictly ascending order, or
///   - the vietoris-rips complex does not contain the simplex, meaning that there exists a pair of vertices
///     `a, b` such that `dissimilarity_matrix.structural_nonzero_entry( &a, &b )` is `None`. 
///     **Note**  *we do not disallow dissimilarity matrices with missing diagonal entries*, so it is for `dissimilarity_matrix.structural_nonzero_entry( &a, &a )` to equal `None`.
///     In this case the function will return `Err( vertices.clone() )`.
/// 
/// - `Ok(DissimilarityMatrix::Coefficient::min_value())` if the vector `vertices` is empty; this is the minimum possible value that can be represented by an object of type  `DissimilarityMatrix::Coefficient`
/// - the maximum of the pairwise dissimilarities of the vertices in the simplex, otherwise.
///   If the simplex contains only one vertex `v`, then the value returned is `Ok(dissimilarity_matrix.structural_nonzero_entry( &v, &v ))`.
pub fn filtration_value_for_clique< DissimilarityMatrix >( 
        vertices:                       & Vec< Vertex >, 
        dissimilarity_matrix:           DissimilarityMatrix
    ) 
    -> Result< DissimilarityMatrix::Coefficient, Vec< Vertex >> 

    where
        DissimilarityMatrix:                MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
        DissimilarityMatrix::Coefficient:   Debug + Ord + num_traits::Bounded,
{
    if ! is_sorted_strictly(&vertices, & OrderOperatorAuto) {
        return Err( vertices.clone() ); // vertices must be sorted in ascending order
    }
    let mut diam = DissimilarityMatrix::Coefficient::min_value();
    let mut a;
    let mut b;
    for ii in 0..vertices.len() {
        a = usize::from( vertices[ii] ); 
        for jj in ii .. vertices.len() {
            b = usize::from( vertices[jj] ); 
            match dissimilarity_matrix.structural_nonzero_entry( &a, &b ) {
                None => { return Err( vertices.clone() ) } // the simplex never enters the filtration
                Some( diam_bound ) => {
                    if diam_bound > diam { diam = diam_bound.clone(); }
                }
            }
        }
    }
    Ok(diam)
}   



//  ===========================================================
//  SORTED ENUMERATION OF CLIQUES
//  ===========================================================




/// A vector of simplices sorted first by dimension (ascending `0..dimension_max`, including `dimension_max`) 
/// then by diameter (descending) then by lexicographic order (descending)
/// 
/// If `allow_empty_simplex = True`, then the empty simplex (which has dimension -1) will be included.
pub fn filtered_cliques_in_row_reduction_order< DissimilarityMatrix >( 
            dimension_max:                  isize,
            dissimilarity_matrix:           DissimilarityMatrix,
            dissimilarity_matrix_size:      usize,   
            allow_empty_simplex:            bool         
        ) 
        ->  Vec< WeightedSimplex< DissimilarityMatrix::Coefficient > >  

    where   DissimilarityMatrix:                    Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
            DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
            DissimilarityMatrix::RowEntry:          KeyValNew,   // required for DiagonalEntryIterator to be an iterator    
            DissimilarityMatrix::Row:               Clone,   
        
    {
        let order_operator = OrderOperatorTwistWeightedSimplex{};
        let dissimilarity_matrix_ref = & dissimilarity_matrix;
        let vec: Vec<_> = (-1..dimension_max+1).flat_map(|dimension|
            {
                let mut vec =   AgileSimplexIteratorLexicographicOrder::new(                          
                                        dissimilarity_matrix_ref, 
                                        dissimilarity_matrix_size,                                      
                                        dimension, // dimension is an isize, but we need a usize here   
                                        allow_empty_simplex, // whether or not to include the empty simplex                 
                                    )
                                    .collect_vec();
                vec.sort_by(|x,y| order_operator.judge_cmp(y,x) );
                vec
            })
        .collect();
        vec
}









//  ===========================================================
//  FILTERED VIETORIS RIPS BOUNDARY MATRIX 
//  ===========================================================



/// The chain complex of a filtered Vietoris-Rips complex
/// 
/// This struct represents a filtered Vietoris-Rips complex, as well as its boundary matrix.
/// The coefficient ring is encoded in the `RingOperator` field.
/// 
/// # Internal structure
/// 
/// The `VietorisRipsComplex` struct contains the following fields:
/// - `ring_operator`: a ring operator for the coefficient ring
/// - `dissimilarity_matrix`: symmetric matrix (oracle) representing pairwise dissimilarity values
/// - `dissimilarity_matrix_size`: the number of rows/columns in the dissimilarity matrix; this is necessary because the matrix oracle traits are passive, and do not give a means to directly count or generate the number of row/column indices
/// - `filtration_ordered_neighbor_lists`: a vector of vectors `V`. For each vertex `i`, the vector
///    `V[i] = [ (f0,v0) .. (fn,vn) ]` contains a list of `(filtration, neighbor)` pairs for vertex `i`,
///    such that `filtration` is the filtration value of the edge `(i,neighbor)` in the underlying weighted graph.
///    The list `[ (f0,v0) .. (fn,vn) ]` is sorted by filtration value, with ties broken by lexicographic order of the vertices.
/// 
/// # Examples 
/// 
/// This example shows how to construct a [VietorisRipsComplex], and use it to compute persistent homology.
/// 
/// 
/// ```
/// use oat_rust::algebra::rings::types::native::{FieldRational64};    
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use oat_rust::algebra::matrices::operations::umatch::differential::DifferentialUmatch;    
/// use oat_rust::algebra::matrices::types::third_party::IntoCSR;
/// use oat_rust::algebra::chain_complexes::barcode::get_barcode;
/// 
/// use oat_rust::topology::simplicial::from::graph_weighted::VietorisRipsComplex;        
/// use oat_rust::topology::point_cloud::unit_circle;            
/// 
/// use oat_rust::utilities::order::{ is_sorted_strictly, OrderOperatorByKeyCustom, OrderOperatorAuto, OrderOperatorAutoReverse}; 
/// use oat_rust::utilities::distances::{rowwise_distances};
/// 
/// 
/// use std::sync::Arc;
/// use itertools::Itertools;      
/// 
/// 
/// // Define a point cloud (points sampled from a unit circle)
/// //
/// // In this example we don't want any noise, but we can add noise
/// // to the circle by replacying None with Some(-1.0 .. 1.0); this will
/// // add noise uniformly sampled from [-1,1] to each x and y coordinate
/// 
/// let number_of_points = 5;
/// let pcloud = unit_circle( number_of_points, None);
/// 
/// // Get the distance matrix
/// 
/// let dissimilarity_matrix_data
///     = rowwise_distances(pcloud)
///         .into_iter()
///         .map(|x| x.into_iter().enumerate().collect_vec() )
///         .collect_vec()
///         .into_csr( number_of_points, number_of_points );
/// let dissimilarity_matrix = & dissimilarity_matrix_data;       
/// 
/// // Define the coefficient ring
/// 
/// let ring_operator = FieldRational64::new();
/// 
/// // Define the chain complex
/// 
/// let chain_complex_data = VietorisRipsComplex::new(
///     & dissimilarity_matrix,
///     number_of_points, ring_operator 
/// ).ok().unwrap();
/// let chain_complex = Arc::new(chain_complex_data);    
/// 
/// let dimension_max = 2;       
/// let row_index_vec = chain_complex.cliques_in_row_reduction_order(dimension_max); // runs over the row indices we want to use when computing the factorization
/// 
/// 
/// // Check that oracle has strictly sorted rows
/// 
/// for row_index in row_index_vec.iter() {        
///     assert!( is_sorted_strictly( 
///                                     & chain_complex.row(row_index).collect_vec() , // this line accesses a row 
///                                     & OrderOperatorByKeyCustom::new( OrderOperatorAuto ) 
///                                 ) );
/// }      
/// 
/// // Check that oracle has strictly sorted columns
/// 
/// 
/// for row_index in row_index_vec.iter() {
///     assert!( is_sorted_strictly(    & chain_complex.column_reverse(row_index).collect_vec() ,  // this line accesses a column
///                                     & OrderOperatorByKeyCustom::new( OrderOperatorAutoReverse::new() )  
///                                 ) );
/// }    
/// 
/// 
/// // Get a Umatch factorization
/// // --------------------------
/// 
/// let min_homology_dimension = 0;
/// let max_homology_dimension = 2;
/// let umatch = DifferentialUmatch::new(
///         chain_complex.clone(), 
///         min_homology_dimension,
///         max_homology_dimension,
///     );       
/// 
/// println!("Umatch factorization: {:#?}", &umatch);  
/// 
/// // Get the barcode
/// // ---------------
/// 
/// let return_cycle_representatives = true;
/// let return_bounding_chains = true;
/// let barcode =  get_barcode( 
///         &umatch, 
///         return_cycle_representatives,
///         return_bounding_chains,
///     );
///             
/// // Print the barcode in dimension 1
/// // --------------------------------
///             
/// let dimension = 1;
/// for (id,birth,death) in barcode.intervals_f64( dimension ) {
///     println!("Interval id {:?}:", id);
///     println!("    Dimension {:?}", dimension);
///     println!("    Birth:    {:?}", birth);
///     println!("    Death:    {:?}", death);
/// }
/// 
/// // This should print the following:
/// //
/// // Interval id 5:
/// //     Dimension 1
/// //     Birth:    1.1755705045849465
/// //     Death:    1.902113032590307
/// 
/// // Print the corresponding betti curve in dimension 1
/// // --------------------------------------------------
/// 
/// // The betti curve is formatted as a list of pairs (t0,b0) .. (tN, bN) such that
/// // the betti number for each half-closed interval [ti,ti+1) is bi.
/// 
/// println!(""); // new line
/// println!("Betti curve in dimension {:?}:", dimension);
/// for (t,b) in barcode.betti_curve(dimension) {
///     println!("    {:?} --> betti number {:?}", t.into_inner(), b);
/// }
/// 
/// // This should print the following:
/// //
/// // Betti curve in dimension 1:
/// //     1.1755705045849465 --> betti number 1
/// //     1.902113032590307 --> betti number 0        
/// 
/// 
/// // Print a cycle representative for the bar in dimension 1
/// // -------------------------------------------------------
/// //
/// // Note: calling `barcode.bar(5)` returns information abou the 5th
/// // bar in the barcode. We know that this is the bar we want, because
/// // the id number 5 appears in the barcode information printed above
/// 
/// let bar = barcode.bar(5);
/// println!(""); // new line
/// println!("Cycle representative for the bar in dimension 1:");
/// for entry in bar.cycle_representative().clone().unwrap() {
///     println!("    Simplex: {:?}", entry.0.vertices());
///     println!("       Filtration: {:?}", entry.0.filtration());            
///     println!("       Coefficient: {:?}", entry.1);            
/// }
/// 
/// // This should print the following
/// //
/// // Cycle representative for the bar in dimension 1:
/// //     Simplex: [3, 4]
/// //     Filtration: OrderedFloat(1.175570504584946)
/// //     Coefficient: Ratio { numer: -1, denom: 1 }
/// //     Simplex: [0, 1]
/// //     Filtration: OrderedFloat(1.1755705045849463)
/// //     Coefficient: Ratio { numer: -1, denom: 1 }
/// //     Simplex: [1, 2]
/// //     Filtration: OrderedFloat(1.1755705045849463)
/// //     Coefficient: Ratio { numer: -1, denom: 1 }
/// //     Simplex: [2, 3]
/// //     Filtration: OrderedFloat(1.1755705045849463)
/// //     Coefficient: Ratio { numer: -1, denom: 1 }
/// //     Simplex: [0, 4]
/// //     Filtration: OrderedFloat(1.1755705045849465)
/// //     Coefficient: Ratio { numer: 1, denom: 1 }
/// ```

#[derive(Clone,Debug,)]
pub struct VietorisRipsComplex< DissimilarityMatrix, RingOperator >
    where 
        DissimilarityMatrix:    MatrixOracle< ColumnIndex=usize, RowIndex=usize, >
{
    pub ring_operator:                              RingOperator, // ring meta data
	pub dissimilarity_matrix:                       DissimilarityMatrix,
    pub dissimilarity_matrix_size:                  usize,
    pub filtration_ordered_neighbor_lists:          Arc< Vec<  // one vec for every vertex in the whole Vietoris-Rips complex
                                                        Arc< Vec< (DissimilarityMatrix::Coefficient, Vertex) >
                                                    >>>,
    pub allow_empty_simplex:                        bool, // determines whether or not we include the emtpy simplex in the simplicial complex
}


/// Methods of VietorisRipsComplex struct
impl < DissimilarityMatrix, RingOperator  > 

    VietorisRipsComplex
        < DissimilarityMatrix, RingOperator > 

    where   DissimilarityMatrix:                MatrixOracle< ColumnIndex=usize, RowIndex=usize, >,
            DissimilarityMatrix::Coefficient:   Debug + Ord + num_traits::Bounded,
            RingOperator:                       RingOperations,
               
{

    /// Construct a (lazy) boundary matrix for the filtered VR complex of a symmetric matrix.
    /// - `dissimilarity_matrix` a symmetric matrix representing pairwise dissimilarity values
    /// - `dissimilarity_matrix_size` an integer representing the number of rows/columns of the dissimilarity matrix; this parameter is necessary because the matrix oracle traits are passive, and do not give a means to directly count or generate the number of row/column indices
    /// - `dissimilarity_value_max` maximum dissimilarity threshold
    /// - `dissimilarity_value_min` diameter of the empty simplex; this can be any number that is equal or less than every structural nonzero entry in the dissimilarity matrix
    /// - `ring_operator`: ring operator for the coefficient ring
    /// The filtration parameter of
    /// - the empty simplex `[]` is `dissimilarity_value_min`
    /// - a singletime `[i]` is `dissimilarity_matrix[i,i]`
    /// - an edge `[i,j]` is `dissimilarity_matrix[i,j]`
    /// - an arbitrary simplex `s` is the maximum of the filtration parameters of cardinality 0, 1, and 2 subsets of `s`
    /// subject to the condition that the filtration parameter of
    /// - a singleton `[i]` is undefined if `dissimilarity_matrix[i,i]` is structurally zero
    /// - an edge `[i,j]` is undefined if `dissimilarity_matrix[i,j]` is structurally zero
    /// - an arbitrary simplex `s` is undefined if the filtration parameter of any cardinality 0, 1, or 2 subset of `s` is undefined
    /// This means that some vertices may never enter the filtration.
    /// 
    /// # NB: Empty simplex
    /// 
    /// This constructor initializes the `allow_empty_simplex` field to `false`, meaning that the empty simplex is not included in the simplicial complex.
    /// However, for those who want to include the empty simplex, it is possible to set the `allow_empty_simplex` field to `true` after construction.
    /// 
    /// # Safety checks
    /// 
    /// The constructor checks that 
    /// - the input matrix is symmetric
    /// - the filtration parameter of every edge equals or exceeds the filtration parameter of its vertices.
    pub fn new(
            dissimilarity_matrix:               DissimilarityMatrix,
            dissimilarity_matrix_size:          usize,
            ring_operator:                      RingOperator,
        ) -> Result< Self, (usize,usize) >
    {

        let number_of_points = dissimilarity_matrix_size;

        // initialize the filtration ordered neighbor lists
        let mut filtration_ordered_neighbor_lists = Vec::with_capacity(number_of_points);

        // check that the dissimilarity matrix is symmetric
        for i in 0 .. number_of_points {

            // we will store the min filtration value of the row in this variable after iterating through the row
            let mut min_filtration = MakeNoneMaximum::from_opt(None); 
            
            // we will store the number of neighbors in this variable after iterating through the row
            let mut num_neighbors = 0;
            
            // iterate through the row
            for entry in dissimilarity_matrix.row( & i ) {   
                let n = entry.key(); 
                let f = entry.val(); 

                // update the min filtration
                min_filtration = min_filtration.min( MakeNoneMaximum::from_val(f.clone()) );   

                // increment the number of neighbors
                num_neighbors += 1;

                // check that the dissimilarity matrix has equal values for (i,j) and (j,i)         
                let entry_n_i = dissimilarity_matrix.structural_nonzero_entry( &n, &i );
                if Some( f.clone() ) != entry_n_i {
                    println!("\n\nError: The dissimilarity matrix passed to the `VietorisRipsComplex` constructor is not symmetric");
                    println!(" - Entry {:?} equals {:?}", (i,n) , f.clone() );
                    println!(" - However, entry {:?} equals {:?}.", (n,i), entry_n_i);
                    println!("\nThis message is generated by OAT.\n\n");                    
                    return Err( (i,n) );
                }
            }

            // check that the diagonal entry in row i is no greater than any other entry in the row
            let diagonal_entry = dissimilarity_matrix.structural_nonzero_entry( &i, &i );
            if min_filtration != MakeNoneMaximum::from_opt( diagonal_entry.clone() ) {
                println!("\n\nError during construction of the `VietorisRipsComplex`: Entry ({:?},{:?}) of the dissimilarity matrix is {:?} but the minimum structural nonzero entry in row {:?} is {:?}.  In this case `None` represents a value strictly greater than `Some(x)`, for every filtration value `x`. This indicates that input dissimilarity matrix is invalid.\nThis message is generated by OAT.\n\n",
                    i, 
                    i, 
                    diagonal_entry, 
                    i, 
                    min_filtration.into_inner() 
                );
                return Err( (i,i) );
            }

            // construct the filtration-ordered list of neighbors
            let mut filtration_ordered_neighbor_list = Vec::with_capacity(num_neighbors);
            for entry in dissimilarity_matrix.row( & i ) {   
                let neighbor = entry.key() as Vertex; 
                let filtration = entry.val(); 
                filtration_ordered_neighbor_list.push( (filtration, neighbor) );
            }
            // sort the list
            filtration_ordered_neighbor_list.sort();

            // push the list into our records, wrapping in Arc
            filtration_ordered_neighbor_lists.push(Arc::new(filtration_ordered_neighbor_list));
        }        

        Ok(
            VietorisRipsComplex { 
                ring_operator, 
                dissimilarity_matrix, 
                dissimilarity_matrix_size,
                filtration_ordered_neighbor_lists: Arc::new(filtration_ordered_neighbor_lists),
                allow_empty_simplex: false,
            }
        )

        
    }            




    /// Returns the diameter of a simplex
    ///
    /// # Parameters
    /// 
    /// -`vertices`: vertices of the simplex
    /// 
    /// # Returns
    /// 
    /// - `Err( vertices.clone() )` if either
    ///   - `vertices` is not sorted in strictly ascending order, or
    ///   - the vietoris-rips complex does not contain the simplex, meaning that there exists a pair of vertices
    ///     `a, b` such that `dissimilarity_matrix.structural_nonzero_entry( &a, &b )` is `None`. 
    ///     **Note**  *we do not disallow dissimilarity matrices with missing diagonal entries*, so it is for `dissimilarity_matrix.structural_nonzero_entry( &a, &a )` to equal `None`.
    ///     In this case the function will return `Err( vertices.clone() )`.
    /// 
    /// - `Ok(DissimilarityMatrix::Coefficient::min_value())` if the `vertices` is empty; this is the minimum possible value that can be represented by an object of type  `DissimilarityMatrix::Coefficient`
    /// - the maximum of the pairwise dissimilarities of the vertices in the simplex, otherwise.
    ///   If the simplex contains only one vertex `v`, then the value returned is `Ok(dissimilarity_matrix.structural_nonzero_entry( &v, &v ))`.
    pub fn filtration_value_for_clique( 
            & self, 
            vertices: & Vec< Vertex > 
        ) -> Result< DissimilarityMatrix::Coefficient, Vec< Vertex >> 
    {
        filtration_value_for_clique( vertices, self.dissimilarity_matrix_ref() )
    }    


    /// Add a filtration value to a simplex, returning a `WeightedSimplex` object
    /// 
    /// # Parameters
    /// 
    /// - `vertices`: vertices of the simplex
    /// 
    /// # Returns
    /// 
    /// - `Ok(WeightedSimplex)` if the filtration value is defined for the simplex
    /// - `Err( vertices.clone() )` if the filtration value is not defined for the simplex,
    ///    meanint that either the vertices are not sorted in strictly ascending order, or
    ///    the dissimilarity matrix lacks a value for one or more pairs of vertices
    pub fn add_filtration_value_to_simplex( &self, vertices: &Vec<Vertex> ) 
        -> Result<
                WeightedSimplex< DissimilarityMatrix::Coefficient >,
                Vec<Vertex>
            >
    {
        let filtration_value = self.filtration_value_for_clique( vertices )?;
        Ok(WeightedSimplex { vertices: vertices.clone(), weight: filtration_value })
    }


    /// Returns a reference to the internally stored symmetric dissimilarity matrix, wrapped in a struct that allows query operations to be performed.
    pub fn dissimilarity_matrix_ref(&self) -> & DissimilarityMatrix { & self.dissimilarity_matrix }

    /// Returns a reference to the internally stored symmetric dissimilarity matrix, wrapped in a struct that allows query operations to be performed.
    pub fn dissimilarity_matrix_size(&self) -> usize { self.dissimilarity_matrix_size }    

    /// Returns the dissimilarity value of two vertices, wrapped in a `MakeNoneMaximum` struct
    /// 
    /// The purpose of the `MakeNoneMaximum` struct is to "switch" the order of None values,
    /// so that they are regarded as strictly greater (rather than strictly less) than all `Some(x)` values.
    pub fn dissimilarity_as_nonegreater( &self, a: Vertex, b: Vertex ) -> MakeNoneMaximum<DissimilarityMatrix::Coefficient> 
        { MakeNoneMaximum { opt: self.dissimilarity_matrix.structural_nonzero_entry(& (a as usize), & (b as usize) ) } }

    /// The ring operator for the coefficient ring.
    pub fn ring_operator( &self ) -> RingOperator where RingOperator: Clone { self.ring_operator.clone() }

    /// A vector of simplices appearing in the order that rows are visited during row reduction of the boundary matrix for U-match factorization.
    /// 
    /// Simplices sorted first by dimension (ascending `0..dimension_max`, including `dimension_max`) 
    /// then by diameter (descending) then by lexicographic order (descending)
    pub fn cliques_in_row_reduction_order( &self, dimension_max: isize ) 
        -> 
        Vec< WeightedSimplex< DissimilarityMatrix::Coefficient > > 
        where 
            DissimilarityMatrix::RowEntry:      KeyValNew,   // required for DiagonalEntryIterator to be an iterator    
            DissimilarityMatrix::Row:           Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator   
    { 
        return filtered_cliques_in_row_reduction_order( 
            dimension_max,  
            self.dissimilarity_matrix_ref(), 
            self.dissimilarity_matrix_size,
            self.allow_empty_simplex, // determines whether or not we include the empty simplex in the simplicial complex
        )
    }    

    /// A vector of simplices sorted in the order imposed on the boundary matrix.
    /// 
    /// Simplices are sorted first by dimension (ascending `0..dimension_max`, including `dimension_max`),
    /// then by filtration value (ascending), then by lexicographic order (ascending).
    pub fn cliques_in_boundary_matrix_order_fixed_dimension( &self, dimension_max: isize ) 
        -> Vec< WeightedSimplex< DissimilarityMatrix::Coefficient > > 
        where 
            DissimilarityMatrix::RowEntry:      KeyValNew,   // required for DiagonalEntryIterator to be an iterator    
            DissimilarityMatrix::Row:           Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator    
    { 
        let mut simplices = self.cliques_in_lexicographic_order_fixed_dimension(dimension_max).collect_vec();
        simplices.sort_by( |x,y| OrderOperatorTwistWeightedSimplex{}.judge_cmp( x, y ) );
        simplices
    }

    /// An iterator that runs over simplices in ascending lexicographic order (not by filtration, etc.).
    pub fn cliques_in_lexicographic_order_fixed_dimension( &self, dimension: isize ) 
        -> AgileSimplexIteratorLexicographicOrder< & DissimilarityMatrix > 
        where 
            DissimilarityMatrix::RowEntry:      KeyValNew,   // required for DiagonalEntryIterator to be an iterator    
            DissimilarityMatrix::Row:           Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator    
    { 
        AgileSimplexIteratorLexicographicOrder::new(                          
            self.dissimilarity_matrix_ref(), // a simple trick to make a copy of the dissimilarity matrix that implements the Copy trait
            self.dissimilarity_matrix_size,                
            dimension,   
            self.allow_empty_simplex,                            
        )
    }        

    /// Returns a factorization of the boundary matrix restricted to simplices
    /// of dimension 0 .. dimension_max+1 (inclusive)
    /// 
    /// This factorization can be used to compute betti numbers, cycle representatives,
    /// and more.  See the documenation for `FactoredChainComplex`.
    pub fn factor_from_arc( self, max_homology_dimension: isize ) 
            ->
            DifferentialUmatch< 
                    Arc< VietorisRipsComplex< DissimilarityMatrix, RingOperator > >,                                     
                >    

        where
            RingOperator::Element:                  Debug + PartialEq, 
            RingOperator:                           Clone + DivisionRingOperations,
            DissimilarityMatrix:                    Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize, >,    
            DissimilarityMatrix::RowEntry:          KeyValNew,   // required for DiagonalEntryIterator to be an iterator               
            DissimilarityMatrix::Row:               Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator
            DissimilarityMatrix::Coefficient:       Copy + Clone + Debug + Ord + Hash + Bounded,             
    {
        let min_homology_dimension = 0;
        let arc = Arc::new(self);
        DifferentialUmatch::new(
            arc,
            min_homology_dimension,
            max_homology_dimension,
        )
    }



}




















//  ===========================================================
//  CHAIN COMPLEX TRAIT
//  ===========================================================




impl  < DissimilarityMatrix, RingOperator >

    ChainComplex for
    
    Arc< 
            VietorisRipsComplex
                < DissimilarityMatrix, RingOperator > 
        >    

    where
        DissimilarityMatrix:                        Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,     
        DissimilarityMatrix::Coefficient:           Copy + Debug + Ord + Bounded, 
        DissimilarityMatrix::RowEntry:              KeyValNew,   // required for DiagonalEntryIterator to be an iterator
        DissimilarityMatrix::Row:                   Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator
        RingOperator:                               Clone + RingOperations,
        RingOperator::Element:                      Debug + PartialEq,
{

    type BasisVectorIndicesIterable = Vec< WeightedSimplex< DissimilarityMatrix::Coefficient > >; 


    fn basis_vector_indices_for_dimension( &self, dimension: isize ) -> Self::BasisVectorIndicesIterable {
        self.cliques_in_boundary_matrix_order_fixed_dimension( dimension )
    }


    fn dimension_for_basis_vector_with_index( &self, index: & Self::RowIndex ) -> Result<isize, Self::RowIndex> {        

        match self.has_row_for_index( index ) {
            false => Err( index.clone() ), // if input is invalid return an error
            true => Ok( index.dimension() as isize ), // otherwise return the dimension of the simplex
        }
    }
}







//  ===========================================================
//  *FILTERED* CHAIN COMPLEX TRAIT
//  ===========================================================




impl  < DissimilarityMatrix, RingOperator >

    FilteredChainComplex for
    
    Arc< 
            VietorisRipsComplex
                < DissimilarityMatrix, RingOperator > 
        >    

    where
        DissimilarityMatrix:                        Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,     
        DissimilarityMatrix::Coefficient:           Copy + Debug + Ord + Bounded, 
        DissimilarityMatrix::RowEntry:              KeyValNew,   // required for DiagonalEntryIterator to be an iterator
        DissimilarityMatrix::Row:                   Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator
        RingOperator:                               Clone + RingOperations,
        RingOperator::Element:                      Debug + PartialEq,
{

    type FiltrationValue = DissimilarityMatrix::Coefficient; // the filtration value is the same as the dissimilarity value

    /// Returns the filtration value of the basis vector (i.e. simplex) with the given index.
    /// 
    /// This function returns `Ok(filtration_value)` if the input is a correctly formatted simplex in the Vietoris-Rips complex
    /// (with correctly calculated weight), and `Err(index.clone())` otherwise.
    fn filtration_value_for_basis_vector_with_index( 
            &self, 
            index: & Self::RowIndex 
        ) ->    
            Result<
                Self::FiltrationValue, 
                Self::RowIndex 
            > {
        
        match self.has_row_for_index( index ) {
            false => Err( index.clone() ), // if input is invalid return an error
            true => {
                Ok( index.filtration() )
            }            
        }
    }
}


















//  ===========================================================
//  FACET ITERATORS
//  ===========================================================



/// A iterator that runs over all entries in the boundary of a simplex, in ascending lexicographic order.
/// 
/// If the original simplex is `[0,1,2]`, it returns entries in the following order: `([0,1],1), ([0,2],-1), ([1,2],1)`.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::topology::simplicial::from::graph_weighted::{AgileBoundaryIteratorLexicographicOrder, VietorisRipsComplex};
/// use oat_rust::topology::simplicial::simplices::weighted::WeightedSimplex;   
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::rings::types::native::FieldFloat64;
/// use ordered_float::OrderedFloat;
/// use std::sync::Arc;
///
/// // define the ring operator
/// let ring_operator = FieldFloat64::new();        
/// 
/// // define the dissimilarity matrix
/// // (this is a convenient way to build a VecOfVec from a "dense matrix")
/// let dissimilarity_matrix = VecOfVec::from_ragged(
///     & vec![
///         vec![OrderedFloat(0.0), OrderedFloat(1.0), OrderedFloat(2.0)],
///         vec![OrderedFloat(1.0), OrderedFloat(0.0), OrderedFloat(3.0)],
///         vec![OrderedFloat(2.0), OrderedFloat(3.0), OrderedFloat(0.0)],
///     ]
/// );
/// 
/// // define the Vietoris-Rips complex
/// let vietoris_rips_complex = Arc::new( 
///     VietorisRipsComplex::new(
///         & dissimilarity_matrix,
///         3, // dissimilarity_matrix_size
///         ring_operator,
///     ).ok().unwrap()
/// );
/// 
/// // define the simplex
/// let simplex = WeightedSimplex {
///     vertices: vec![0, 1, 2],
///     weight: OrderedFloat(2.0),
/// }; 
/// 
/// // create the iterator
/// let mut boundary = AgileBoundaryIteratorLexicographicOrder::new(
///     vietoris_rips_complex,
///     simplex,
///     false, // this indicates we don't consider the emtpy simplex to be a "real" simplex
/// );   
/// 
/// // check the results
/// assert!( 
///     boundary.eq(
///         vec![
///             ( WeightedSimplex { vertices: vec![0, 1], weight: OrderedFloat(1.0) },  1.0 ),
///             ( WeightedSimplex { vertices: vec![0, 2], weight: OrderedFloat(2.0) }, -1.0 ),
///             ( WeightedSimplex { vertices: vec![1, 2], weight: OrderedFloat(3.0) },  1.0 ),
///         ]
///     )         
/// );
/// ```
/// A iterator that runs over all entries in the boundary of a simplex, in **descending lexicographic order**
/// 
/// If the original simplex is `[0,1,2]`, it return entries in the following order: `([1,2],1), ([0,2],-1), ([0,1],1)`.
#[derive(Clone, Debug)]
pub struct AgileBoundaryIteratorLexicographicOrder
        < DissimilarityMatrix, RingOperator >
    where
        DissimilarityMatrix:                            MatrixOracle< ColumnIndex=usize, RowIndex=usize >,       
        RingOperator:                                   RingOperations, 
{
    pub vietoris_rips_complex: Arc< VietorisRipsComplex< DissimilarityMatrix, RingOperator > >,
	pub simplex: WeightedSimplex<DissimilarityMatrix::Coefficient>,
	pub removal_location_for_next_facet: Option<usize>,
    pub ring_operator: RingOperator,
}


impl< DissimilarityMatrix, RingOperator > 
    
    AgileBoundaryIteratorLexicographicOrder
        < DissimilarityMatrix, RingOperator >
    where // we impose all these conditions so that the VietorisRipsComplex will implement MatrixAlgebra, which lets us extract the coefficient ring
        DissimilarityMatrix::Coefficient:           Copy + Debug + Ord + Bounded,
        DissimilarityMatrix::RowEntry:              KeyValNew,   // required for DiagonalEntryIterator to be an iterator
        DissimilarityMatrix::Row:                   Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator 
        DissimilarityMatrix:                        Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,         
        RingOperator:                               Clone + RingOperations,
        RingOperator::Element:                      Debug + PartialEq,
{

    /// Create a new [AgileBoundaryIteratorLexicographicOrder]
    /// 
    ///  # Parameters
    /// 
    ///  - `vietoris_rips_complex`: an Arc-wrapped Vietoris-Rips complex
    ///  - `simplex`: a `WeightedSimplex` object representing the simplex for which the boundary is computed
    ///  - `allow_empty_simplex`: a boolean flag that indicates whether the iterator should include the 
    ///     empty simplex (of dimension -1) in the boundary of a singleton simplex.
    /// 
    /// # See also
    /// 
    /// The documentation for [AgileBoundaryIteratorLexicographicOrder] for more details on the iterator.
    pub fn new(
            vietoris_rips_complex: Arc< VietorisRipsComplex< DissimilarityMatrix, RingOperator > >,
            simplex: WeightedSimplex<DissimilarityMatrix::Coefficient>,
            allow_empty_simplex: bool,
        ) -> Self
    {
        // get the initial removal location for the next facet; this should be the last vertex in the simplex
        let removal_location_for_next_facet = match simplex.vertices.len() {
            0 => {
                // if the simplex is empty, there is no next vertex to remove
                None
            } 
            1 => {
                // if the simplex is a singleton, we can remove it to get the empty simplex, if that's allowed
                if allow_empty_simplex {
                    Some(0) 
                } else {
                    None 
                }
            } 
            _ => {
                Some(simplex.vertices.len() - 1)
            }
        };

        let ring_operator = vietoris_rips_complex.ring_operator();
        AgileBoundaryIteratorLexicographicOrder {
            vietoris_rips_complex,
            simplex,
            removal_location_for_next_facet,
            ring_operator,
        }
    }
}        




/// implement standard methods of Iterator for AgileBoundaryIteratorLexicographicOrder struct
impl< DissimilarityMatrix, RingOperator >

    Iterator for 
    
    AgileBoundaryIteratorLexicographicOrder
        < DissimilarityMatrix, RingOperator >
    where
        DissimilarityMatrix:                            MatrixOracle< 
                                                            ColumnIndex     =   usize, 
                                                            RowIndex        =   usize,
                                                            Coefficient:        Ord + Bounded,
                                                        >,   
        RingOperator:                                   RingOperations, 
{
	type Item =     (
                        WeightedSimplex
                            <DissimilarityMatrix::Coefficient>, 
                        RingOperator
                            ::Element
                    );

	fn next( &mut self ) -> Option<Self::Item> {

        // return None if there is no "next" vertex to delete; otherwise get a pointer
        // to the next vertex to delete
        let removal_location_for_next_facet = self.removal_location_for_next_facet?;
        
        // remove the vertex and compute the filtraiton and coefficient
        let mut simplex = self.simplex.clone();
        simplex.vertices.remove(removal_location_for_next_facet );
        simplex.weight = self.vietoris_rips_complex.filtration_value_for_clique(&simplex.vertices).unwrap();
        let coeff = self.ring_operator.minus_one_to_power( removal_location_for_next_facet );

        // update the removal location for the next facet
        match removal_location_for_next_facet {
            0 => { self.removal_location_for_next_facet = None; }
            _ => { self.removal_location_for_next_facet = Some(removal_location_for_next_facet - 1); }
        }

        // return the facet and coefficient
        Some((simplex, coeff))        
    }
}


//  ---------------------------------------------------------------------------
//  IMPLEMENT MATRIX ORACLE CLIQUE BOUNDARY MATRIX + DEFINE SOME SORTED (CO)BOUNDARY ITERATORS
//  ---------------------------------------------------------------------------



//  MATRIX ORACLE
//  ------------------------------------------






// Matrix Oracle
impl  < DissimilarityMatrix, RingOperator >

    MatrixOracle for
    
    Arc< 
            VietorisRipsComplex
                < DissimilarityMatrix, RingOperator > 
        >    

    where
        DissimilarityMatrix::Coefficient:           Copy + Debug + Ord + Bounded, 
        DissimilarityMatrix::RowEntry:              KeyValNew,   // required for DiagonalEntryIterator to be an iterator
        DissimilarityMatrix::Row:                   Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator
        RingOperator:                               Clone + RingOperations,
        RingOperator::Element:                      Debug + PartialEq,
        DissimilarityMatrix:                        Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >, 
{   

    type Coefficient            =   RingOperator::Element;// The type of coefficient stored in each entry of the matrix    
    type RowIndex               =   WeightedSimplex<DissimilarityMatrix::Coefficient>;// The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type RowEntry               =   ( Self::RowIndex, Self::Coefficient );// The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnIndex            =   WeightedSimplex<DissimilarityMatrix::Coefficient>;// The type of column indices    
    type ColumnEntry            =   ( Self::RowIndex, Self::Coefficient );// The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =   AgileCoboundaryIteratorFiltrationOrder< DissimilarityMatrix, RingOperator >; // An iterator that runs over the entries in a row, sorted by filtration value (ascending) and then by lexicographic order (descending)
    type RowReverse             =   std::vec::IntoIter<(WeightedSimplex<DissimilarityMatrix::Coefficient>, RingOperator::Element)>; // What you get when you ask for a row with the order of entries reversed
    type Column                 =   std::vec::IntoIter<(WeightedSimplex<DissimilarityMatrix::Coefficient>, RingOperator::Element)>; // What you get when you ask for a column
    type ColumnReverse          =   std::vec::IntoIter<(WeightedSimplex<DissimilarityMatrix::Coefficient>, RingOperator::Element)>;
    
    fn row(                     &   self, row_index: & Self::RowIndex   )   -> Self::Row {
        self.row_result(row_index).ok().unwrap()
    }
    
    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
        let mut row: Vec<_>     =   self.row( index ).collect();
        (&mut row).reverse();
        row.shrink_to_fit();
        return row.into_iter()
    }

    fn row_result(                 &   self, index: & Self::RowIndex   )   -> Result< Self::Row, Self::RowIndex > {
        let row_result =         AgileCoboundaryIteratorFiltrationOrder::new(
            index.vertices().clone(),
            self.dissimilarity_matrix.clone(),
            self.dissimilarity_matrix_size,
            self.filtration_ordered_neighbor_lists.clone(),
            self.ring_operator(),
        );
        match row_result {
            Ok(row) => Ok(row),
            Err(_) => {
                // if the row is not defined, we return an error
                Err(index.clone())
            }
        }
    }    
    
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {

        let cardinality = index.number_of_vertices();

        let iter = AgileBoundaryIteratorLexicographicOrder::new(
            self.clone(), // vietoris_rips_complex: 
            index.clone(), // simplex:  
            self.allow_empty_simplex, // the AgileBoundaryIteratorLexicographicOrder will be empty is this value is set to false (and self.allow_empty_simplex *IS* set to false by default)
        );        
        // initialize empty list
        let mut vec = Vec::with_capacity( cardinality );
        // fill with entries from the iterator
        vec.extend( iter );
        // sort
        vec.sort_by(|x,y| x.0.cmp(&y.0)); 
        // return as iterator
        vec.into_iter()
    }
    
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
        let cardinality = index.number_of_vertices();
        let iter = AgileBoundaryIteratorLexicographicOrder::new(
            self.clone(), // vietoris_rips_complex: 
            index.clone(), // simplex: 
            self.allow_empty_simplex, // the AgileBoundaryIteratorLexicographicOrder will be empty is this value is set to false (and self.allow_empty_simplex *IS* set to false by default)
        );        
        // initialize empty list
        let mut vec = Vec::with_capacity( cardinality );
        // fill with entries from the iterator
        vec.extend( iter );
        // sort
        vec.sort_by(|x,y| y.0.cmp(&x.0)); // sort in REVERSE order
        // return as iterator
        vec.into_iter()
    }
    
    /// Returns `false` if any of the following conditions fails:
    /// - the list of vertices is sorted in strictly ascending order (in particular, it does not contain repeated entries)
    /// - every pair of vertices is adjacent in the underlying graph
    /// - the filtration value of `index` is the maximum of the weights of the edges (and the vertices) in `index`. 
    ///   In practice, the weights of the vertices come into play only in the special case where there are no
    ///   edges, i.e. when the simplex is a singleton.
    /// If all these conditions are satisfied, then the function returns `true`.
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {

        match self.filtration_value_for_clique( index.vertices() ) {
            Err(_) => false, // if the filtration value is not defined, then the row does not exist
            Ok(filtration_value) => {
                filtration_value == index.weight
            }            
        }
    }
    
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.has_row_for_index( index )
    }
    
    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {

        // both indices must be valid
        if !  self.has_row_for_index(row)  {
            panic!("Attempted to look up the entry in row {:?}, column {:?} of a Vietoris Rips boundary matrix, but {:?} is either an invalid simplex (not sorted or has a repeat entry) or it is not contained in the complex", row, column, row );
        }
        if !  self.has_column_for_index(row)  {
            panic!("Attempted to look up the entry in row {:?}, column {:?} of a Vietoris Rips boundary matrix, but {:?} is either an invalid simplex (not sorted or has a repeat entry) or it is not contained in the complex", row, column, column );
        }        

        // this criterion has to be satisfied to have a facet-cofacet pair
        if column.number_of_vertices() != row.number_of_vertices() + 1 { 
            return None 
        }

        // get the symmetric difference
        let mut symmetric_difference                        =   symmetric_difference_of_ordered_iterators( row.vertices.iter(), column.vertices.iter() );
        let first_different_vertex                      =   symmetric_difference.next(); 
        
        if symmetric_difference.next().is_some() {
            // in this case the row and column indices differ by more than one vertex, so they cannot be a facet-cofacet pair
            return None
        } 
                
        // if we make it to here then we do have a facet-cofacet pair, so we have to calculate the corresponding coefficient by finding 
        // which vertex in column is not contained in row, and then calculating the corresponding coefficient
        for (counter, vertex) in column.vertices.iter().enumerate() {
            if Some(vertex) == first_different_vertex {
                let coefficient         =   self.ring_operator.minus_one_to_power( counter );
                return Some( coefficient )
            }
        }    

        panic!("Something went wrong in the calculation of a structural nonzero");
    } 

} 









// Matrix Algebra
impl  < DissimilarityMatrix, RingOperator >

    MatrixAlgebra for
    
    Arc< 
            VietorisRipsComplex
                < DissimilarityMatrix, RingOperator > 
        >    

    where
        DissimilarityMatrix:                        Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,     
        DissimilarityMatrix::Coefficient:           Copy + Debug + Ord + Bounded, 
        DissimilarityMatrix::RowEntry:              KeyValNew,   // required for DiagonalEntryIterator to be an iterator
        DissimilarityMatrix::Row:                   Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator
        RingOperator:                               Clone + RingOperations,
        RingOperator::Element:                      Debug + PartialEq,
{
    type RingOperator                                   =   RingOperator;

    type OrderOperatorForRowEntries                     =   OrderOperatorByKey;

    type OrderOperatorForRowIndices                     =   OrderOperatorAuto;

    type OrderOperatorForColumnEntries                  =   OrderOperatorByKey;

    type OrderOperatorForColumnIndices                  =   OrderOperatorAuto;

    fn ring_operator( &self ) -> Self::RingOperator {
        self.ring_operator.clone()
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        OrderOperatorByKey::new()
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        OrderOperatorAuto::new()
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        OrderOperatorByKey::new()
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        OrderOperatorAuto::new()
    }
}









//  Matrix Oracle Operations
//  ------------------------------------------


impl < DissimilarityMatrix, RingOperator >

    MatrixOracleOperations for
    
    Arc< 
        VietorisRipsComplex
            < DissimilarityMatrix, RingOperator > 
    >    

    where
        DissimilarityMatrix:                        Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,     
        DissimilarityMatrix::Coefficient:           Copy + Debug + Ord + Bounded, 
        DissimilarityMatrix::RowEntry:              KeyValNew,   // required for DiagonalEntryIterator to be an iterator
        DissimilarityMatrix::Row:                   Clone,       // required for AgileSimplexIteratorLexicographicOrder to be an iterator
        RingOperator:                               Clone + RingOperations,
        RingOperator::Element:                      Debug + PartialEq,    

{
    fn bottom_entry_for_column( &self, column_index: &<Arc<VietorisRipsComplex<DissimilarityMatrix, RingOperator>> as crate::algebra::matrices::query::MatrixOracle>::ColumnIndex ) 
        -> Option< 
                <Arc<VietorisRipsComplex<DissimilarityMatrix, RingOperator>> as MatrixOracle>::ColumnEntry
        > {

        let order_operator = self.order_operator_for_column_entries();
        let iter = AgileBoundaryIteratorLexicographicOrder::new(
            self.clone(),
            column_index.clone(),
            false, // don't allow the empty simplex to be returned as a boundary
        );
        return iter.max_by(|x, y| order_operator.judge_cmp( &x, &y ) );
    }

}






//  ===========================================================================
//  SMALL COFACET ITERATOR
//  ===========================================================================





/// Iterates over all vertices which can be added to a facet to form a cofacet of equal filtration value
/// 
/// Internally, this iterator essentially works by intersecting the set of all neighborhoods of the vertices
/// in the facet.
#[derive(Clone,PartialEq)]
pub struct SmallCofacetVertexIterator< DissimilarityMatrix >
    where
        DissimilarityMatrix:                    MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
        DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
        DissimilarityMatrix::RowEntry:          KeyValNew,   // required for DiagonalEntryIterator to be an iterator
{
    pub proveably_exhausted:                    bool,
    pub neighbor_iterators:                     Vec< DissimilarityMatrix::Row >, 
    pub facet:                                  Arc< Vec< Vertex > >,
    pub facet_filtration:                       DissimilarityMatrix::Coefficient,
}


impl < DissimilarityMatrix > 

    SmallCofacetVertexIterator
        < DissimilarityMatrix >

    where
        DissimilarityMatrix:                    MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
        DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
        DissimilarityMatrix::RowEntry:          KeyValNew,   // required for DiagonalEntryIterator to be an iterator        
{


    /// Constructs a new `SmallCofacetVertexIterator` iterator.
    /// 
    /// # Parameters
    /// 
    /// - `dissimilarity_matrix`: a symmetric matrix representing pairwise dissimilarity values
    /// - `dissimilarity_matrix_size`: an integer representing the number of rows/columns of the dissimilarity matrix; this parameter is necessary because the matrix oracle traits are passive, and do not give a means to directly count or generate the number of row/column indices
    /// - `facet`: a vector of vertices that form the facet
    /// - `facet_filtration`: the filtration value of the facet, which is the maximum of the filtration values of all vertices in the facet
    /// 
    /// # Returns
    /// 
    /// A new `SmallCofacetVertexIterator` iterator that can be used to iterate over all vertices which can be added to the facet to form a cofacet of equal filtration value (all other vertices are skipped).
    /// 
    /// # Note
    /// 
    /// We do not that the inputs are valid, i.e. that
    /// - the facet is sorted in strictly ascending order
    /// - the user-provided filtration is correct
    /// - all vertices in the facet are adjacent
    pub fn new(
            dissimilarity_matrix: DissimilarityMatrix,
            facet: Arc< Vec< Vertex > >,
            facet_filtration: DissimilarityMatrix::Coefficient,
        ) -> Self 
    {

        // we haven't checked any vertices yet, so we ONLY mark the iterator as proveably exhausted if the facet is empty
        // (if it is, then by convention SmallCofacetVertexIterator should act as an empty iterator)
        let proveably_exhausted =  facet.is_empty();


        // collect the neighbor iterartors
        let mut neighbor_iterators = Vec::with_capacity(facet.len());
        for vertex in facet.iter() {
            let row = dissimilarity_matrix.row(& (*vertex as usize ) );
            neighbor_iterators.push(row);
        }

        SmallCofacetVertexIterator {
            proveably_exhausted,
            neighbor_iterators,
            facet,
            facet_filtration,
        }
    }


    /// Reurns the lexicographically smallest vertex that lies within distance `self.facet_filtration` of vertex `self.facet[k]`
    pub fn next_plausible_candidate_for_kth_facet_vertex( &mut self, k: usize )
        -> Option< Vertex >
    {   
        let neighbor_iterator = &mut self.neighbor_iterators[k];
        while let Some(entry) = neighbor_iterator.next() {                    

            let neighbor = entry.key() as Vertex;
            let filtration = entry.val();

            // check that the neighbor is not already in the facet
            if self.facet.contains(&neighbor) {
                continue; // if the neighbor is already in the facet, we skip it
            }

            // if the neighbor equals candidate but filtration is greater than the facet, we skip this neighbor and seek a new candidate
            if filtration > self.facet_filtration { 
                continue
            } else {
                return Some( neighbor );
            }

        }     
        return None;   
    }
}




impl < DissimilarityMatrix >

    Iterator for 

    SmallCofacetVertexIterator
        < DissimilarityMatrix >

    where
        DissimilarityMatrix:                    MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
        DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
        DissimilarityMatrix::RowEntry:          KeyValNew,   // required for DiagonalEntryIterator to be an iterator        
{
    type Item = Vertex;

    fn next( &mut self ) -> Option<Self::Item> {

        if self.proveably_exhausted {
            return None; // if we are exhausted, we return None
        }

        // NB: 

        let mut candidate_neighbor: u16                       = self.next_plausible_candidate_for_kth_facet_vertex( 0 )?;
        let mut number_of_facet_vertices_checked            = 1;
        let mut next_facet_vertex_to_check: usize           = 0 % self.facet.len(); // we will increment this variable in the first phase of the outer loop
        
        'outer_loop: loop {               

            // if we have checked all vertices then the candidate is valid; return it
            if number_of_facet_vertices_checked == self.facet.len() {
                return Some( candidate_neighbor );
            }

            // otherwise get a pointer to the next vertex to check
            // we do this by incrementing the vertex to check by 1, modulo the number of vertices in the facet
            next_facet_vertex_to_check += 1;
            next_facet_vertex_to_check %= self.facet.len();    

            while let Some( neighbor ) = self.next_plausible_candidate_for_kth_facet_vertex( next_facet_vertex_to_check ) {
                if neighbor < candidate_neighbor {
                    // in this case neighbor can't be a cofacet vertex, since candidate_neighbor is the lowest possible number, so move on
                    continue 
                }
                if neighbor > candidate_neighbor {
                    // in this case candidate_neighbor is not viable; neighbor is the new lowest-possible value candidate value
                    candidate_neighbor = neighbor;
                    number_of_facet_vertices_checked = 1; // reset the number of facet vertices checked
                    continue 'outer_loop;
                } 
                // if we make it here, then neighbor == candidate_neighbor, candidate_neighbor is sufficiently close
                // to vertex `next_vertex_to_check`
                // therefore we can increment the number of facet vertices checked, and move on to the next facet vertex
                number_of_facet_vertices_checked += 1;
                continue 'outer_loop;
            }

            // if we make it here then next_facet_vertex_to_check has no more viable neighbors, so no more cofacet vertices remain to be found
            self.proveably_exhausted = true;
            return None;
        }
    }
}












//  ===========================================================================
//  DIAGONAL ENTRY ITERATOR  -- ordered lexicographically
//  ===========================================================================







/// Runs over all structural nonzero diagonal entries in a dissimilarity matrix
/// 
/// Items returned have type `Option< DissimilarityMatrix::RowEntry >`
#[derive(Clone,Debug,PartialEq, Eq)]
pub struct DiagonalEntryIterator< DissimilarityMatrix >
    where   DissimilarityMatrix:                    MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
            DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
{
    pub row_number_iterator:        std::ops::Range<usize>,
    pub dissimilarity_matrix:       DissimilarityMatrix,    
}



impl < DissimilarityMatrix >

    DiagonalEntryIterator
        < DissimilarityMatrix > 

    where   DissimilarityMatrix:                    MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
            DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
{
    pub fn new( dissimilarity_matrix: DissimilarityMatrix, size: usize ) -> Self {
        DiagonalEntryIterator {
            row_number_iterator: 0..size,
            dissimilarity_matrix,
        }
    }
}




impl < DissimilarityMatrix >

    Iterator for

    DiagonalEntryIterator
        < DissimilarityMatrix > 

    where   DissimilarityMatrix:                    MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
            DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
            DissimilarityMatrix::RowEntry:          KeyValNew,
{
    type Item =  DissimilarityMatrix::RowEntry;

    fn next( &mut self ) -> Option<Self::Item> {        
        while let Some(row_number) = self.row_number_iterator.next() {
            match self.dissimilarity_matrix.structural_nonzero_entry( & row_number, & row_number ) {
                None => continue, // skip to the next row if this one has no entry
                Some(coefficient) => {
                    return Some( 
                        DissimilarityMatrix::RowEntry::new( row_number, coefficient ) 
                    );
                }
            }
        }
        return None;        
    }
}












//  ===========================================================================
//  SIMPLEX ITERATOR  -- ordered lexicographically
//  ===========================================================================










/// An iterator that runs over all simplices of given dimension in a clique complex, in **ascending lexicographic order**.
/// 
/// # How it works
/// 
/// Internally, this object stores a list of iterators, a list of filtration values, and a list of vertices.
/// The structure is as shown in the following example:
/// 
/// ```text
/// 
/// Iterators       =   [ D   I1  I2  I3  I4  I5 ]
/// Vertices        =   [ V0  V1  V2  V3  V4     ]
/// Filtrations     =   [ F0  F1  F2  F3  F4     ]
/// ```
/// 
/// Here `D` is an iterator that runs over all pairs `(v,f)` such that `dissimilarity_matrix[v,v]=f` is
/// a structural nonzero entry in the dissimilarity matrix. The iterator runs over vertices in
/// lexicographic order. The value `f` is the filtration value of the vertex.
/// 
/// The each iterator `I1, I2, ..` is a (partially expended) copy of row `V0` of the dissimilarity matrix.
/// 
/// For each `k`
/// - the sequence `[V0, V1, ..., Vk]` is a valid simplex of dimension `k`, in the Vietoris-Rips complex
/// - the filtration value of this simplex is `fk`
/// 
/// **The number of iterators always exceeds the number of vertices (respectively, filtrations) by one.**
/// 
/// When the function `self.next()` is called, this object will pass through the following stages.
/// 
/// **(A)** Begin drawing elements `(V,F)` from the last iterator in `Iterators`.
/// Here `V` is a neighbor of `V0` and `F` is the filtration value of the simplex `[V0, V]`. 
/// 
/// **(A1)** If `[V0 .. V4 V]` is a valid simplex with filtration value `K`, then ..
/// 
/// - If `5 = self.dimension`, then the simplex `[V0, V1, V2, V3, V4, V]` is returned, and the process
///   terminates until the next call to `self.next()`.
/// 
/// - If `5 < self.dimension`, then we push V to the end of the list of vertice and F to the end of the list of filtrations.
///   We then make an identical copy of I5, and push it to the end of the list of iterators. Thus the structure becomes
///  as follows:
/// 
/// ```text
/// Iterators       =   [ D   I1  I2  I3  I4  I5  I6 ]
/// Vertices        =   [ V0  V1  V2  V3  V4  V5     ]
/// Filtrations     =   [ F0  F1  F2  F3  F4  F5     ]
/// ``` 
/// Then we return to step (A) and continue drawing elements from the last iterator `I6`.
/// 
/// **(A2)** If `[V0 .. V4 V]` is **not** a valid simplex, then we continue drawing elements from the last iterator `I6`
/// until we find a valid simplex (in which case we proceed with (A1)), or until the iterator is exhausted. 
/// 
/// **(A3)** If iterator `I5` is exhausted then we remove `I5`, `V4`, and `F4` from their respective lists.
/// Thus the structure becomes:   
///   
/// ```text
/// Iterators       =   [ D   I1  I2  I3  I4 ]
/// Vertices        =   [ V0  V1  V2  V3     ]
/// Filtrations     =   [ F0  F1  F2  F3     ]  
/// ```
/// 
/// We then begin fresh with step (A).
/// 
/// **(B)** Once we exhaust iterator `D`, the process terminates. At this point we have generated all simplices of
/// dimension `self.dimension`
/// 
/// 
/// 
/// 
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::third_party::IntoCSR;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;   
/// use oat_rust::algebra::rings::types::native::FieldRationalSize;     
/// use oat_rust::topology::simplicial::simplices::weighted::WeightedSimplex;
/// use oat_rust::topology::simplicial::from::graph_weighted::VietorisRipsComplex;
/// use std::sync::Arc;   
/// use itertools::Itertools;                   
/// 
/// let number_of_points = 10;
/// let maxdim = 2;       
/// 
/// // initialize a random dissimilarity matrix
/// let dissimilarity_matrix_vecvec = VecOfVec::random_symmetric_zero_diagonal_with_enclosing_radius_threshold(number_of_points);
/// let dissimilarity_matrix =  Arc::new( dissimilarity_matrix_vecvec.into_csr(number_of_points, number_of_points) );                
/// 
/// // create the Vietoris-Rips complex
/// let ring_operator = FieldRationalSize::new();
/// let vietoris_rips_complex_data = VietorisRipsComplex::new( dissimilarity_matrix, number_of_points, ring_operator ).ok().unwrap();
/// let vietoris_rips_complex = Arc::new(vietoris_rips_complex_data);        
/// 
/// // generate simplices using the agile iterator
/// let mut simplices_agile = vietoris_rips_complex.cliques_in_row_reduction_order(maxdim);
/// simplices_agile.sort();
/// 
/// // generate simplices using a brute-force approach
/// let mut simplices_brute_force = Vec::new();
/// for dimension in 0 .. maxdim as usize + 1 {
///     for combination in (0..number_of_points as u16).combinations(dimension+1) {
///         if let Ok( filtration ) = vietoris_rips_complex.filtration_value_for_clique(&combination) {
///             simplices_brute_force.push(WeightedSimplex {
///                 vertices: combination,
///                 weight: filtration,
///             });
///         }
///     }
/// }
/// simplices_brute_force.sort();        
/// 
/// // verify that both generation methods align
/// assert!(
///     simplices_agile.eq( & simplices_brute_force ),
///     "The agile simplex enumeration does not match the brute-force enumeration"
/// );
/// ```
/// 
/// An empty iterator:
/// 
/// ```
/// use sprs::CsMatBase;
/// use ordered_float::OrderedFloat;
/// use oat_rust::topology::simplicial::from::graph_weighted::AgileSimplexIteratorLexicographicOrder;
/// use oat_rust::topology::simplicial::simplices::weighted::WeightedSimplex;
/// 
/// 
/// let dissimilarity_matrix = CsMatBase::new( (2,2), vec![0,1,2], vec![0,1], vec![OrderedFloat(0.0),OrderedFloat(0.0)] );
///     
/// let simplices: Vec<_>    =   AgileSimplexIteratorLexicographicOrder::new( 
///                                  & dissimilarity_matrix, 
///                                  2,          //  matrix size
///                                  1,          //  dimension
///                                  false,      //  do not include the empty simplex (this only applies when dimension is -1)
///                              ).collect();
/// let empty_list           =   Vec::< WeightedSimplex<OrderedFloat<f64>> >::new();
///             
/// assert!( empty_list.eq( &simplices ) ) 
/// ```
#[derive(Clone,Debug)]
pub struct AgileSimplexIteratorLexicographicOrder< DissimilarityMatrix > 
    
    where   DissimilarityMatrix:                    MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
            DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
            DissimilarityMatrix::RowEntry:          KeyValNew,            // apparently to be inside an itertwotype we require an object to be an iterator; and DiagonalEntryIterator is only an iterator if DissimilarityMatrix::RowEntry implements KeyValNew              
{
    pub dissimilarity_matrix:       DissimilarityMatrix,
    pub neighbor_iterators:         Vec< 
                                        TwoTypeIterator<
                                            DiagonalEntryIterator< DissimilarityMatrix >,
                                            DissimilarityMatrix::Row,
                                        >
                                    >,
    pub face_filtrations:           Vec< DissimilarityMatrix::Coefficient >, // the filtration values of the simplices we are building
    pub vertices:                   Vec<Vertex>, // the vertices of the simplex we are building
    pub dimension:                  isize,
}

impl < DissimilarityMatrix >

    AgileSimplexIteratorLexicographicOrder
        < DissimilarityMatrix > 

    where   DissimilarityMatrix:                    Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
            DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
            DissimilarityMatrix::RowEntry:          KeyValNew,                        
    {

    /// Creates a new AgileSimplexIteratorLexicographicOrder object that runs over all simplices of given dimension in a clique complex
    /// 
    /// # Arguments
    /// 
    /// * `dissimilarity_matrix` - a dissimilarity matrix that defines the clique complex; it must be indexed by `usize` and have a `Coefficient` type that implements `Ord + Debug + Bounded`
    /// * `dissimilarity_matrix_size` - the number of rows and columns of the dissimilarity matrix
    /// * `dimension` - the dimension of simplices to iterate over; must be a non-negative integer
    /// 
    /// # Returns
    /// 
    /// A new AgileSimplexIteratorLexicographicOrder object that runs over all simplices of given dimension in a clique complex
    /// 
    /// # Initialization
    /// 
    /// See the documentation of the AgileSimplexIteratorLexicographicOrder struct on internal structure.  We initialize the relevant lists as follows:
    /// ```text
    /// Iterators       =   [ D   ]
    /// Vertices        =   [     ]
    /// Filtrations     =   [     ]  
    /// ```
    pub fn new( 
            dissimilarity_matrix:           DissimilarityMatrix,
            dissimilarity_matrix_size:      usize, // the size of the dissimilarity matrix
            dimension:                      isize,
            allow_empty_simplex:            bool,
        ) 
        -> 
        Self 
        where
            DissimilarityMatrix:                            MatrixOracle< ColumnIndex=usize, RowIndex=usize >
    {      

        // if the dimension is negative, initialize all internal data as empty, with one caveat
        if dimension < 0 {
            let neighborhood_iterators = Vec::with_capacity(0);
            let face_filtrations = Vec::with_capacity(0);
            // if dimension is -1 (corresponding to the empty simplex), then we add one unit of capacity to vertices as a flag;
            // the first time we call `next()` on the iterator, we will return an empty simplex with minimum filtration value, 
            // then reduce capacity to zero
            let mut vertices= match dimension {
                -1 => Vec::with_capacity(1), 
                _ => Vec::with_capacity(0), 
            };
            // if directed to exclude the empty clique, then we reset the flag from the previous step
            if ! allow_empty_simplex {
                vertices = Vec::with_capacity(0); 
            }
            return AgileSimplexIteratorLexicographicOrder {
                dimension,
                dissimilarity_matrix,
                neighbor_iterators: neighborhood_iterators,
                face_filtrations,
                vertices,
            };
        }

        // in this case dimension is non-negative, so we can convert to usize for convenience
        let dimension = dimension as usize;

        // runs over vertices and their filtration values
        let diagonal_entry_iterator         =   TwoTypeIterator::Version1(
                                                    DiagonalEntryIterator::new( dissimilarity_matrix.clone(), dissimilarity_matrix_size )
                                                );
        let mut neighbor_iterators          =   Vec::with_capacity( dimension + 1 );
        neighbor_iterators.push( diagonal_entry_iterator ); // the first neighbor iterator is the diagonal entry iterator, which runs over all vertices and their filtration values

        let face_filtrations            =   Vec::with_capacity(dimension+1); 

        let vertices                    =   Vec::with_capacity(dimension+1);


        AgileSimplexIteratorLexicographicOrder {
            dimension: dimension as isize, 
            dissimilarity_matrix,
            neighbor_iterators,
            face_filtrations,
            vertices,
        }
    }




    // Returns `Some(n-1)` if `self.neighbor_iterators` has length `n > 0`, otherwise returns `None`.
    pub fn last_neighborhood_ordinal( &self ) -> Option<usize> {
        match self.neighbor_iterators.len() {
            0 => None, // if there are no neighbor iterators, we are done
            n => Some( n - 1 ), // otherwise return the ordinal of the last neighbor iterator
        }
    }

}







impl< DissimilarityMatrix > 

    Iterator for 
    
    AgileSimplexIteratorLexicographicOrder
        < DissimilarityMatrix > 
        
    where   DissimilarityMatrix:                    Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
            DissimilarityMatrix::Coefficient:       Ord + Debug + Bounded,
            DissimilarityMatrix::RowEntry:          KeyValNew,   // required for DiagonalEntryIterator to be an iterator    
            DissimilarityMatrix::Row:               Clone,   
{
    type Item = WeightedSimplex<DissimilarityMatrix::Coefficient>;


    /// This implementation requires `DissimilarityMatrix` to implement Clone, because we want to be able to clone
    /// an `TwoTypeIterator< DiagonalEntryIterator, DissimilarityMatrix::Row >` object, and the `DiagonalEntryIterator` requires
    /// `DissimilarityMatrix` to implement `Clone` in order to be able to clone itself.
    /// 
    /// In practice this hasn't been a problem, since we usually use an `Arc<T>` or `& T` for our dissimilarity matrix.
    fn next(&mut self) -> Option<Self::Item> {



        // first handle some the edge cases where dimension < 0
        if self.dimension < 0 {
            if self.vertices.capacity() == 0 {
                // we use self.vertices.capacity() as a flag to indicate that, even if we do count the empty simplex as a simplex of dimension -1, we have already returned it
                return None;
            } else {
                if self.vertices.len() > 0 {
                    panic!("AgileSimplexIteratorLexicographicOrder::next() called on a simplex iterator with dimension < 0, but the vertices list is not empty. This indicates that the simplex iterator was not properly initialized. This error message originates in the oat_rust::topology::simplicial::from::graph_weighted module, in the AgileSimplexIteratorLexicographicOrder struct.");
                }
                self.vertices.shrink_to(0); 
                return  Some( WeightedSimplex { 
                            vertices: self.vertices.clone(), 
                            weight: DissimilarityMatrix::Coefficient::min_value(), 
                        } );
            }
        }



        'outer_loop: loop {

            
            // if ZERO neighborhoods are stored in memory, then there are no more neighbors to check, so there are no more simplices to return, and the iterator is exhausted
            let last_neighborhood_ordinal                =   self.last_neighborhood_ordinal()?; 


            // this while-loop scans through elements of last_neighbor_iterator
            'search_within_neighborhood: while let Some(entry) = ( &mut self.neighbor_iterators[last_neighborhood_ordinal] ).next() { // get the next entry from the last neighbor iterator
                let v = entry.key(); // get the vertex of the next neighbor
                let mut f = entry.val(); // get the filtration value of the next neighbor

                // update f to be at least as large as the filtration value of the simplex we are building
                // if the last_ordinal is 0 there's no prior filtration to compare
                // if the last_ordinal is 1 then the only prior filtration to compare is the filtration of vertex 0, which is
                // is no greater than f since f is the filtration parameter of an edge incident to vertex 0
                // thus we only update f if the last_ordinal is greater than 1
                if last_neighborhood_ordinal > 1 { 
                    // update the filtration value to be at least as large as the previous filtration value
                    // recall that the previous filtration value is stored int he slot ONE TO THE LEFT of the last neighborhood ordinal
                    f = f.max(self.face_filtrations[last_neighborhood_ordinal-1].clone()); 
                }


                // update the filtration to equal the maximum of the distances to the other vertices in the facet
                // !! NB: we can skip the first vertex, because the original value of f was the distance from that vertex to the new vertex v
                if last_neighborhood_ordinal > 0 {
                    for vertex in self.vertices[1..last_neighborhood_ordinal].iter().cloned().map(|x| x as usize) {
                        match self.dissimilarity_matrix.structural_nonzero_entry( &vertex, &v) {
                            Some(dissimilarity) => {
                                f = f.max(dissimilarity); // update the filtration value to be at least as large as the dissimilarity value
                            },
                            None => { 
                                continue 'search_within_neighborhood; // this vertex can't be added; continue scanning for the next neighbor in the neighborhood
                            }
                        }
                    }
                }


                // if we make it hear then we have checked that v is adjacent to every vertex in self.vertices, so we can add it to create a valid simplex in the VR complex
                // NP: no matter what happens in the following match statement, we will exit the ambient while-loop
                match last_neighborhood_ordinal == self.dimension as usize  {
                    true => { 
                        // the new simplex has dimension self.dimension
                        let mut vertices = self.vertices.clone();
                        vertices.push(v as u16);

                        return Some( WeightedSimplex { vertices, weight: f, } );

                    } false => { 
                        // the new siplex has sub-full dimension
                        //
                        // in this case we record the new vertex and filtration
                        self.vertices.push(v as u16); // add the new vertex to the simplex                        
                        self.face_filtrations.push(f); // update the filtration value of the simplex

                        // we also push a new neighbor iterator to our list
                        match last_neighborhood_ordinal == 0 {
                            // if the current list of iterators only includes the diagonal entry iterator, grab the neighors of the first vertex
                            true => {
                                let mut new_neighbor_iterator = self.dissimilarity_matrix.row( &v ); // get the iterator for the new vertex
                                // remove all entries for vertices that come before v
                                while let Some(entry) = new_neighbor_iterator.next() { 
                                    if entry.key() == v { break }
                                }
                                // now we can push the new neighbor iterator to the list of iterators
                                self.neighbor_iterators.push( 
                                    TwoTypeIterator::Version2(
                                        new_neighbor_iterator
                                    )
                                ); 
                            } false => { // otherwise just clone the last neighbor iterator; it will start in the right place, since it just returned vertex v
                                self.neighbor_iterators.push( 
                                    self.neighbor_iterators[last_neighborhood_ordinal].clone()
                                ); 
                            }
                        }

                        // now we can keep scanning for the next vertex to add
                        continue 'outer_loop; 

                    }                    
                }
            }
            // at this point our neighbor iterator is exhausted and we haven't found a new vertex to add
            // so first we have to step back in the search tree:
            self.neighbor_iterators.pop(); 
            self.face_filtrations.pop();
            self.vertices.pop(); 
            // now we can re-start the outer loop           

        }

    }

}












//  ===========================================================================
//  EDGE IDENTIFIER ITERATOR  -- for cofacets ordered by filtration
//  ===========================================================================

/// Iterates over the "big cofacet edge identifiers" of a facet in a Vietoris-Rips complex.
/// 
/// 
/// A **big cofacet** is a cofacet with filtration strictly greater than the given facet.
/// The **edge identifier** for a big cofacet is the unique edge whose length equals
/// the filtration of the cofacet, and which is lexicographically smallest
/// among all edges with that property.
/// There is a 1-1 correspondence between big cofacets and edge identifiers.
/// 
/// We can impose a total order on edge identifiers by
/// 
/// (1) mapping each edge identifier to the ordered pair `(filtration, vertex)` where `vertex` is the vertex that is not in the facet,
/// and `filtration` is the filtration value of the edge identifier (equivalently, of the associated cofacet), 
/// (2) ordering these pairs lexcicographically
/// 
/// This order corresponds exactly to the order we impose on cofacets (first by filtraiton, with ties broken by lexicographic order of vertices).
/// 
/// **This iterator returns edge identifiers in for big cofacets, in the order described.**
#[derive(Clone,Debug,PartialEq, Eq)]
pub struct BigCofacetEdgeIterator
        < DissimilarityMatrix >

    where
        DissimilarityMatrix:                MatrixOracle< ColumnIndex=usize, RowIndex=usize >,        
{
    facet_filtration:                       DissimilarityMatrix::Coefficient,
    facet:                                  Arc<Vec< Vertex >>,
    filtration_ordered_neighborhood_lists:  Vec<
                                                    IterWrappedArcVec<
                                                    (
                                                        DissimilarityMatrix::Coefficient, 
                                                        Vertex,
                                                    )
                                                >
                                            >,
    dissimilarity_matrix:                   DissimilarityMatrix,
    cofacet_edge_identifiers:               Vec< MakeNoneMaximum<(DissimilarityMatrix::Coefficient, Vertex)> >,
}


impl < DissimilarityMatrix >   

    BigCofacetEdgeIterator
        < DissimilarityMatrix >
    where
        DissimilarityMatrix:                MatrixOracle< ColumnIndex=usize, RowIndex=usize >,        
        DissimilarityMatrix::Coefficient:   Ord + Copy,
{


    fn new( 
        facet:                                              Arc<Vec<Vertex>>,
        dissimilarity_matrix:                               DissimilarityMatrix,
        facet_filtration:                                   DissimilarityMatrix::Coefficient,
        all_ambient_filtration_ordered_neighborhood_lists:  Arc< Vec<  // one vec for every vertex in the whole Vietoris-Rips complex
                                                                Arc< Vec< (DissimilarityMatrix::Coefficient, Vertex) >
                                                            >>>, 
    ) -> Self {

        let facet_cardinality = facet.len();
        
        // look up the neighborhood iterators for each vertex in the facet, and place them all in a vector
        // (recall these iterators run over (filtration, neighbor) pairs in filtration order, which in this case coincides with lexicographic order)
        let mut neighbor_lists_filtration_order = Vec::with_capacity( facet_cardinality );
        for v in facet.iter() {
            neighbor_lists_filtration_order.push(
                IterWrappedArcVec::new(
                    all_ambient_filtration_ordered_neighborhood_lists[ *v as usize ].clone()
                )                
            );
        }

        // initialize the iterator
        let mut initialized_cofacet_iterator = BigCofacetEdgeIterator {
            facet_filtration: facet_filtration,
            facet,
            filtration_ordered_neighborhood_lists: neighbor_lists_filtration_order,
            dissimilarity_matrix,
            cofacet_edge_identifiers: Vec::new(),
        };


        // find the first valid cofacet edge identifier for each vertex in the facet
        let mut cofacet_edge_identifiers: Vec<_> = Vec::with_capacity( facet_cardinality );
        for k in 0 .. facet_cardinality {
            cofacet_edge_identifiers.push(
                MakeNoneMaximum::from_opt(
                    initialized_cofacet_iterator
                        .find_next_identifier_for_the_kth_vertex( k )
                )
            );
        }

        // insert these identifiers into the initialized object
        initialized_cofacet_iterator.cofacet_edge_identifiers = cofacet_edge_identifiers;

        return initialized_cofacet_iterator;
    }


    /// The number of vertices in the facet
    fn facet_cardinality( &self ) -> usize {
        self.facet.len()
    }

    /// The dimension of the facet
    fn facet_dimension( &self ) -> isize {
        ( self.facet.len() as isize ) - 1
    }    



    /// Checks if a given edge is a valid cofacet identifier for the facet
    /// 
    /// Specifically, the edge must be longer than the facet filtration, 
    /// and it must be farther from the vertex vk than any other vertex in the facet.
    /// If it is equidistant to vk and any other vertex vj in the facet,
    /// then it is valid only if k < j.
    fn edge_is_a_valid_cofacet_identifier( 
        &self, 
        k: usize, 
        neighbor_of_vk: Vertex, 
        distance_vk_to_neighbor_of_vk: DissimilarityMatrix::Coefficient,
    ) -> bool {

        // only edges strictly longer than the facet filtration are valid cofacet identifiers
        if distance_vk_to_neighbor_of_vk <= self.facet_filtration { 
            return false; 
        }        

        // check that the edge is farther from vk than any other vertex in the facet
        for ( facet_vertex_ordinal, facet_vertex ) in self.facet.iter().enumerate() {

            // if neighbor_of_vk is already in the facet, then we can't add it to form a cofacet
            if *facet_vertex == neighbor_of_vk { 
                continue; 
            }

            // check the filtration value
            if facet_vertex_ordinal == k { continue }   // in this case we've already done the relevant check for vk, so we skip it
            match self.dissimilarity_matrix.structural_nonzero_entry( 
                & ( neighbor_of_vk as usize), 
                &( * facet_vertex as usize ) 
            ) {
                // case 1: facet_vertex is not adjacent to neighbor_of_vk, so we exclude neighbor_of_vk
                None => {
                    return false;                                        
                }
                // case 2: facet_vertex is adjacent to neighbor_of_vk
                Some(f) => {
                    match f.cmp(&distance_vk_to_neighbor_of_vk) {
                        // if neighbor_of_vk is strictly closer to facet_vertex than to vk, then we don't discard neighbor_of_vk
                        std::cmp::Ordering::Less => {
                            continue
                        }
                        // if neighbor_of_vk is strictly closer to vk than to facet_vertex, then we discard neighbor_of_vk
                        std::cmp::Ordering::Greater => {
                            return false;                                        
                        }          
                        // if there is a tie, we give ownership of neighbor_of_vk to the facet vertex with the lowest index                          
                        std::cmp::Ordering::Equal => {
                            if facet_vertex_ordinal < k {
                                return false;
                            } else {
                                continue
                            }
                        }
                    }
                }
            }
        }
        // if we make it to here, then the edge is a valid cofacet identifier
        return true
    }




    fn find_next_identifier_for_the_kth_vertex( &mut self, facet_vertex_ordinal:usize ) -> Option< (DissimilarityMatrix::Coefficient, Vertex) >{

        // ordinarily we could do a standard find operation to pull out the next valid edge,
        // but rust doesn't like certain borrow conflicts, so we do it by hand
        let mut replacement_identifier: Option< (DissimilarityMatrix::Coefficient, Vertex) >;
        loop {
            replacement_identifier  =   self.filtration_ordered_neighborhood_lists[facet_vertex_ordinal].next();
            // if the next candidate is Some, then we check if it is a valid replacement identifier
            if let Some(  (f, v) ) = replacement_identifier {
                if  self.edge_is_a_valid_cofacet_identifier(
                        facet_vertex_ordinal, 
                        v, 
                        f
                    )                 
                {
                    break;
                // otherwise we keep looking for a valid replacement identifier                    
                } else {
                    continue; 
                }
            // if there are replacement_identifier is None, then we are out of candidates!
            } else {
                break;
            }
        }       
        replacement_identifier 
    }

}





impl < DissimilarityMatrix >  

    Iterator for

    BigCofacetEdgeIterator
        < DissimilarityMatrix >
    where
        DissimilarityMatrix:                MatrixOracle< ColumnIndex=usize, RowIndex=usize >,        
        DissimilarityMatrix::Coefficient:   Ord + Copy,
{
    type Item = (DissimilarityMatrix::Coefficient, Vertex);

    fn next(&mut self) -> Option<Self::Item> {

        // find the next valid cofacet identifier
        let (facet_vertex_ordinal, _next_identifier_ref)    =   self.cofacet_edge_identifiers
                                                        .iter()
                                                        .enumerate()
                                                        .min_by_key( 
                                                            |(_facet_vertex_ordinal, identifier)| {
                                                                // we use the filtration as the primary key, and the vertex number as the secondary key
                                                                **identifier
                                                            }
                                                        )?;
        let replacement_identifier =    MakeNoneMaximum::from_opt(
            self.find_next_identifier_for_the_kth_vertex( 
                facet_vertex_ordinal 
            )
        );
        
        // pull out the selected identifier, and replace it with its successor
        let selected_identifier =   std::mem::replace(
                                        &mut self.cofacet_edge_identifiers[facet_vertex_ordinal], // element we pull from the vector
                                        replacement_identifier                                     // element we put in the vector
                                    );
        
        return selected_identifier.into_inner()
        
    }
}        











//  ===========================================================================
//  COBOUNDARY ITERATOR ----   FILTRATION ORDER
//  ===========================================================================









/// Iterates over all `(WeightedSimplex, RingOperator::Element)` pairs in the coboundary
/// of a simplex, in **filtration order** (with ties broken by lexicographic order of vertices).
/// 
/// # Internal structure
/// 
/// Stores Option<S> and Option<T>, where S and T are iterators that run over vertices which can be added to
/// form valid cofacets.  S runs over vertices for "small cofacet vertices" (those with filtration equal to the facet filtration),
/// and T runs over vertices for "big cofacet edge identifiers" (those with filtration strictly greater than the facet filtration).
///
/// Both Option<S> and Option<T> are initialized as `None`, since there is a small startup cost to initializing them,
/// and in many applications not all entries are needed.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::topology::simplicial::from::graph_weighted::{AgileBoundaryIteratorLexicographicOrder, AgileCoboundaryIteratorFiltrationOrder, VietorisRipsComplex};
/// use oat_rust::topology::simplicial::simplices::weighted::WeightedSimplex;   
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use oat_rust::algebra::rings::types::native::FieldFloat64;
/// use ordered_float::OrderedFloat;
/// use std::sync::Arc;
/// 
/// // define the ring operator
/// let ring_operator = FieldFloat64::new();        
/// 
/// // define the dissimilarity matrix
/// // (this is a convenient way to build a VecOfVec from a "dense matrix")
/// let dissimilarity_matrix = VecOfVec::from_ragged(
///     & vec![
///         vec![OrderedFloat(0.0), OrderedFloat(1.0), OrderedFloat(2.0), OrderedFloat(3.0)],
///         vec![OrderedFloat(1.0), OrderedFloat(0.0), OrderedFloat(1.0), OrderedFloat(2.0)],
///         vec![OrderedFloat(2.0), OrderedFloat(1.0), OrderedFloat(0.0), OrderedFloat(1.0)],
///         vec![OrderedFloat(3.0), OrderedFloat(2.0), OrderedFloat(1.0), OrderedFloat(1.0)],                
///     ]
/// );
/// 
/// // define the Vietoris-Rips complex
/// let vietoris_rips_complex = Arc::new( 
///     VietorisRipsComplex::new(
///         & dissimilarity_matrix,
///         3, // dissimilarity_matrix_size
///         ring_operator,
///     ).ok().unwrap()
/// );
/// 
/// // Check a boundary iterator
/// // ---------------------------------------------------------------
/// 
/// // define the simplex
/// let simplex = WeightedSimplex {
///     vertices: vec![0, 1, 2],
///     weight: OrderedFloat(2.0),
/// };         
/// 
/// // create the iterator
/// let boundary = AgileBoundaryIteratorLexicographicOrder::new(
///     vietoris_rips_complex.clone(),
///     simplex,
///     false, // this indicates we don't consider the emtpy simplex to be a "real" simplex
/// );   
/// 
/// // check the results
/// assert!( 
///     boundary.eq(
///         vec![
///             ( WeightedSimplex { vertices: vec![0, 1], weight: OrderedFloat(1.0) },  1.0 ),
///             ( WeightedSimplex { vertices: vec![0, 2], weight: OrderedFloat(2.0) }, -1.0 ),
///             ( WeightedSimplex { vertices: vec![1, 2], weight: OrderedFloat(1.0) },  1.0 ),
///         ]
///     )         
/// );
/// 
/// // Check a coboundary iterator
/// // ---------------------------------------------------------------   
/// 
/// // Here we use the constructor from the Vietoris-Rips complex, which is more convenient  
/// 
/// // define the simplex
/// let simplex = WeightedSimplex {
///     vertices: vec![0, 2],
///     weight: OrderedFloat(2.0),
/// };  
/// 
/// // create the iterator
/// let coboundary = vietoris_rips_complex.row( & simplex );
/// 
/// // check the results
/// assert!( 
///     coboundary.eq(
///         vec![
///             ( WeightedSimplex { vertices: vec![0, 1, 2], weight: OrderedFloat(2.0) },  -1.0 ),
///             ( WeightedSimplex { vertices: vec![0, 2, 3], weight: OrderedFloat(3.0) },   1.0 ),
///         ]
///     )         
/// );
/// ```
#[derive(Clone)]
pub struct AgileCoboundaryIteratorFiltrationOrder
        < DissimilarityMatrix, RingOperator >
    where
        DissimilarityMatrix:                                MatrixOracle< ColumnIndex=usize, RowIndex=usize >,  
        DissimilarityMatrix::Row:                           Clone, 
        DissimilarityMatrix::RowEntry:                      KeyValNew,     
        DissimilarityMatrix::Coefficient:                   Ord + Copy + Debug + Bounded,      
        RingOperator:                                       RingOperations,  
{
    facet:                                                  Arc<Vec<Vertex>>,
    facet_filtration:                                       DissimilarityMatrix::Coefficient,
    small_cofacet_vertices:                                 Option< SmallCofacetVertexIterator< DissimilarityMatrix > >,
    big_cofacet_edge_identifiers:                           Option< BigCofacetEdgeIterator< DissimilarityMatrix > >,
    dissimilarity_matrix:                                   DissimilarityMatrix,
    dissimilarity_matrix_size:                              usize, // the size of the dissimilarity matrix
    zero_dimensional_cofacets_opt:                          Option< 
                                                                IntoIter<
                                                                    WeightedSimplex< 
                                                                        DissimilarityMatrix::Coefficient
                                                                    >
                                                                >
                                                            >,  
    all_ambient_filtration_ordered_neighborhood_lists:      Arc< Vec<  // one vec for every vertex in the whole Vietoris-Rips complex
                                                                     Arc< Vec< 
                                                                        (DissimilarityMatrix::Coefficient, Vertex) 
                                                                    >>
                                                            >>,
    ring_operator:                                          RingOperator,
}



impl < 'a, DissimilarityMatrix, RingOperator >

    AgileCoboundaryIteratorFiltrationOrder
        < DissimilarityMatrix, RingOperator >
    where
        DissimilarityMatrix:                Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,        
        DissimilarityMatrix::Row:           Clone,   
        DissimilarityMatrix::RowEntry:      KeyValNew,             
        DissimilarityMatrix::Coefficient:   Ord + Copy + Debug + Bounded,   
        RingOperator:                       RingOperations,
{
    /// Creates a new AgileCoboundaryIteratorFiltrationOrder object
    /// 
    /// # Arguments
    /// 
    /// * `facet` - a vector of vertices that form the facet; must be sorted in strictly ascending order
    /// * `dissimilarity_matrix` - a dissimilarity matrix that defines the clique complex; it must be indexed by `usize` and have a `Coefficient` type that implements `Ord + Debug + Bounded`
    /// * `all_ambient_filtration_ordered_neighborhood_lists` - a vector of vectors, where each inner vector is an `Arc<Vec<(DissimilarityMatrix::Coefficient, Vertex)>>` that contains the filtration-ordered neighborhood lists for each vertex in the whole Vietoris-Rips complex
    /// * `ring_operator` - a ring operator that defines the coefficients of the chain complex
    /// 
    /// # Returns
    /// 
    /// A new AgileCoboundaryIteratorFiltrationOrder object that runs over all nonzero entries in the boundary of the facet in filtraton order (with ties split by lexicographic order of vertices)
    /// 
    /// # Errors
    /// 
    /// Returns an error if the facet is invalid, either because the vertices are not sorted in strictly ascending order,
    /// or because it does not represent a clique in the underlying weighted graph.
    pub fn new( 
        facet:                                              Vec<Vertex>,
        dissimilarity_matrix:                               DissimilarityMatrix,
        dissimilarity_matrix_size:                          usize,
        all_ambient_filtration_ordered_neighbor_lists:      Arc< Vec<  // one vec for every vertex in the whole Vietoris-Rips complex
                                                                Arc< Vec< (DissimilarityMatrix::Coefficient, Vertex) >
                                                            >>>, 
        ring_operator:                                      RingOperator,
    ) -> Result< Self, Vec<Vertex> > {


        // first place the facet inside a vec that has capacity for exactly one additional vertex
        // (this is just to ensure we don't use unnecessary memory)
        let mut exact_facet = Vec::with_capacity( facet.len() + 1 );
        exact_facet.extend_from_slice( &facet );
        let facet = Arc::new( exact_facet );



        // returns an error if the facet is invalid, either because the vertices are not sorted in strictly ascending order, or because it does not represent a clique in the underlying graph
        let facet_filtration = filtration_value_for_clique(&(*facet), &dissimilarity_matrix)?;


        // we initialize the dissimilarity matrix diagonal iterator iff the facet is empty
        //
        // if it *is* the empty facet ([],-MAX), then the cofacets are essentially ([vertex], filtration) pairs such that dissimilarity_matrix[vertex, vertex] is structurally nonzero
        // moreover, the coefficient of every entry in the coboundary is 1
        let zero_dimensional_cofacets_opt = match facet.is_empty() {
            true => {
                let mut cofacets = Vec::new();
                for entry in  DiagonalEntryIterator::new( 
                        dissimilarity_matrix.clone(), 
                        dissimilarity_matrix_size
                ){
                    let vertex = entry.key() as Vertex; // get the vertex from the entry
                    let filtration = entry.val(); // get the filtration value from the entry
                    let simplex = vec![vertex];
                    let cofacet = WeightedSimplex::new( filtration, simplex );
                    cofacets.push( cofacet );
                }
                cofacets.sort(); // sort the cofacets
                Some( cofacets.into_iter() )
            }
            false => None, 
        };
        
        Ok( 
            AgileCoboundaryIteratorFiltrationOrder {
                facet,
                facet_filtration,
                dissimilarity_matrix,
                dissimilarity_matrix_size,
                ring_operator,                
                small_cofacet_vertices:                                 None,
                big_cofacet_edge_identifiers:                           None,
                zero_dimensional_cofacets_opt:                          zero_dimensional_cofacets_opt,
                all_ambient_filtration_ordered_neighborhood_lists:      all_ambient_filtration_ordered_neighbor_lists,
            } 
        )
    }


    /// Generates the entry of the coboundary of `self.facet` corresponding to the given vertex.
    /// 
    /// # Arguments
    /// 
    /// * `vertex` - the vertex to add to the facet
    /// * `filtration` - the filtration value of the cofacet that is generated by adding the vertex to the facet
    /// 
    /// # Returns
    /// 
    /// A tuple containing the weighted simplex and the coefficient of the entry in the coboundary
    pub fn coboundary_entry_for_exogenous_vertex(
            & self,
            vertex:             Vertex,
            filtration:         DissimilarityMatrix::Coefficient,
        ) -> ( WeightedSimplex<DissimilarityMatrix::Coefficient>, RingOperator::Element ) 

        where
            RingOperator:               Clone,    
    {        

        let (cofacet_vertices, coefficient) = coboundary_entry_for_facet_vertex_pair(
            & self.facet,
            vertex,
            self.ring_operator.clone(),
        ).ok().unwrap();   

        // create the weighted simplex
        let weighted_simplex = WeightedSimplex {
            vertices: cofacet_vertices,
            weight: filtration,
        };

        return ( weighted_simplex, coefficient );
    }


}





impl < DissimilarityMatrix, RingOperator >

    Iterator for

    AgileCoboundaryIteratorFiltrationOrder
        < DissimilarityMatrix, RingOperator >
    where
        DissimilarityMatrix::Coefficient:       Clone + Copy + Debug + Ord + Bounded,
        DissimilarityMatrix::RowEntry:          KeyValNew,   // required for DiagonalEntryIterator to be an iterator
        DissimilarityMatrix::Row:               Clone,
        DissimilarityMatrix:                    Clone + MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
        RingOperator:                           Clone + RingOperations,
{
    type Item = (WeightedSimplex<DissimilarityMatrix::Coefficient>, RingOperator::Element);

    fn next(&mut self) -> Option<Self::Item> {


        // check for 0-dimensional cofacets first
        // -----------------------------------------------------    

        // The field `self.zero_dimensional_cofacets_opt` is Some iff the facet is empty (ie has dimension -1).

        if let Some( zero_dimensional_cofacets ) = &mut self.zero_dimensional_cofacets_opt {

            return zero_dimensional_cofacets
                .next()
                .map( |weighted_simplex| {
                    ( weighted_simplex, RingOperator::one() )
                })
        }


        // check for small cofacets next
        // -----------------------------------------------------
        //
        // we define "small cofacets" as cofacets with filtration equal to the facet filtration,

        // if we haven't yet initialized the small cofacet vertices iterator, do so now
        if self.small_cofacet_vertices.is_none() {
            self.small_cofacet_vertices = Some(
                SmallCofacetVertexIterator::new(
                    self.dissimilarity_matrix.clone(),
                    self.facet.clone(),
                    self.facet_filtration, // the filtration value of the facet
                )
            );
        }

        // now self.small_cofacet_vertices is initialized
        // check if it contains any remaining small cofacets, then return the next one
        if let Some( small_cofacet_vertices ) = &mut self.small_cofacet_vertices {
            if let Some( neighbor ) = small_cofacet_vertices.next() {

                // let filtration = match self.facet.is_empty() {
                //     true => { self.dissimilarity_matrix.structural_nonzero_entry(& (neighbor as usize),& (neighbor as usize)).unwrap() }
                //     false => { self.facet_filtration.clone() } 
                // };

                return Some( self.coboundary_entry_for_exogenous_vertex( 
                    neighbor, 
                    self.facet_filtration.clone(), // the filtration value of small cofacets equals that of the facet by definition
                ) );
            }
        }


        // check for big cofacets next
        // -----------------------------------------------------    
        //
        // we define "big cofacets" as cofacets with filtration strictly greater than the facet filtration,
     

        // if we haven't yet initialized the big cofacet edge identifiers iterator, do so now                   
        if self.big_cofacet_edge_identifiers.is_none() {
            self.big_cofacet_edge_identifiers = Some(
                BigCofacetEdgeIterator::new(
                    self.facet.clone(),
                    self.dissimilarity_matrix.clone(),
                    self.facet_filtration, // the filtration value of the facet
                    self.all_ambient_filtration_ordered_neighborhood_lists.clone()
                )
            );
        }

        // now self.big_cofacet_vertices is initialized        
        // if there are remaining big cofacets, then return the next one
        if let Some( big_cofacet_edge_identifiers ) = &mut self.big_cofacet_edge_identifiers {
            if let Some( (cofacet_filtration, vertex) ) = big_cofacet_edge_identifiers.next() {
                return Some( self.coboundary_entry_for_exogenous_vertex( 
                    vertex, 
                    cofacet_filtration // the new filtration will be the length of the edge identifier
                ) );
            }
        }


        // otherwise all entries in the coboundary are exhausted
        return None
    }
}








//  ===========================================================================
//  COBOUNDARY ITERATOR ----   LEXICOGRAPHIC ORDER
//  ===========================================================================




/// An iterator that runs over all cofacets of a given simplex, in ascending lexicographic order, regardless of filtration value
pub struct AgileCoboundaryIterartorLexicographicOrder
        < DissimilarityMatrix, RingOperator >
    where
        DissimilarityMatrix:                        MatrixOracle< ColumnIndex=usize, RowIndex=usize, >,   
{
	pub facet:                                      Vec<Vertex>,
    pub facet_filtration:                           DissimilarityMatrix::Coefficient,
	pub dissimilarity_matrix:                       DissimilarityMatrix,
    pub neighbors_of_vertex_0:                      DissimilarityMatrix::Row,    
    pub ring_operator:                              RingOperator,
}



impl < DissimilarityMatrix, RingOperator >

    AgileCoboundaryIterartorLexicographicOrder
        < DissimilarityMatrix, RingOperator >

    where
        DissimilarityMatrix:                        MatrixOracle< 
                                                        ColumnIndex         =   usize, 
                                                        RowIndex            =   usize, 
                                                        Coefficient:            Ord + Bounded, 
                                                    >,   
{
    pub fn new(
        facet:                                              Vec<Vertex>,
        dissimilarity_matrix:                               DissimilarityMatrix,
        ring_operator:                                      RingOperator,
    ) -> Result<Self, Vec<Vertex>> {

        let facet_filtration = filtration_value_for_clique(&facet, &dissimilarity_matrix)?;
        let neighbors_of_vertex_0 = dissimilarity_matrix.row( & (facet[0] as usize) );

        return Ok(
            AgileCoboundaryIterartorLexicographicOrder {
                facet,
                facet_filtration,
                dissimilarity_matrix,
                neighbors_of_vertex_0,
                ring_operator,
            }
        );
    }
}        




/// implement standard methods of Iterator for AgileCoboundaryIterartorLexicographicOrder struct
impl< 'a, DissimilarityMatrix, RingOperator >

    Iterator for 
    
    AgileCoboundaryIterartorLexicographicOrder
        < DissimilarityMatrix, RingOperator >
        
    where
        DissimilarityMatrix::Coefficient:       Debug + Ord,
        DissimilarityMatrix:                    MatrixOracle< ColumnIndex=usize, RowIndex=usize >,
        RingOperator:                           Clone + RingOperations,    
{
	type Item = (
                    WeightedSimplex
                        < DissimilarityMatrix::Coefficient >, 
                    RingOperator
                        ::Element
                );

    fn next( &mut self ) -> Option<Self::Item> {
        
        'outer_loop: while let Some( neighbor_entry ) = self.neighbors_of_vertex_0.next() {
            
            // get the next neighbor of vertex 0
            let neighbor = neighbor_entry.key();
            let mut filtration = neighbor_entry.val();
            filtration = filtration.max(self.facet_filtration.clone()); // update the filtration value to be at least as large as the facet filtration

            // skip if the vertex is already in the facet
            if self.facet.contains(& (neighbor as Vertex )) {
                continue 'outer_loop; // skip this vertex, since it is already in the facet
            }

            // exclude the vertex if it is too far from any facet vertices
            // otherwise update the filtration value to be the maximum of the distances to the facet vertices            
            for facet_vertex in self.facet[1..].iter().cloned() {                
                match self.dissimilarity_matrix.structural_nonzero_entry( & (facet_vertex as usize), & neighbor ) {
                    None => continue 'outer_loop, // skip this vertex, since it is not adjacent to the facet vertex
                    Some(new_filtration) => filtration = filtration.max(new_filtration), // the vertex is adjacent to the facet vertex

                }
            }

            // if we make it here, then the vertex is neighbor is adjacent to all facet vertices

            // create the cofacet and get its coefficient
            // !! NB: this will throw an error if the facet is not sorted or the vertex already belongs to the facet
            let (cofacet, coefficient) = coboundary_entry_for_facet_vertex_pair( 
                &self.facet, 
                neighbor as  Vertex, 
                self.ring_operator.clone(),
            ).ok().unwrap();

            // wrap the cofacet in a WeightedSimplex
            let weighted_cofacet = WeightedSimplex {
                vertices: cofacet,
                weight: filtration,
            };

            return Some((weighted_cofacet, coefficient));
        }

        return None; // no more cofacets to return
	}
}

























//  =========================================================================
//  UNIT TESTS
//  =========================================================================

#[cfg(test)]
mod tests {
    use crate::{algebra::matrices::{debug::matrix_oracle_is_internally_consistent, operations::umatch::differential::DifferentialUmatch}, topology::{point_cloud, simplicial::simplices::{self, weighted::WeightedSimplex}}, utilities::order::OrderOperatorAuto};

    

    





    #[test]
    fn check_that_some_basic_functions_run_without_error_small() {


        use crate::topology::simplicial::{simplices::weighted::WeightedSimplex, from::graph_weighted::{VietorisRipsComplex}};        
        use crate::topology::point_cloud::unit_circle;    

        use crate::algebra::vectors::entries::KeyValGet;
        use crate::algebra::rings::types::native::{FieldRational64};    
        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};    
        use crate::algebra::matrices::types::third_party::IntoCSR;
        
        use crate::utilities::order::{ is_sorted_strictly, OrderOperatorByKeyCustom, OrderOperatorAutoReverse};
        use crate::utilities::iterators::general::minmax;    
        use crate::utilities::distances::{rowwise_distances};

        use std::sync::Arc;
        use ordered_float::OrderedFloat;  
        use itertools::Itertools;      

        let number_of_points = 2;
        let min_homology_dimension = 0;
        let max_homology_dimension = 2;


        let pcloud: Vec<Vec<f64>> = vec![vec![0.], vec![1.]];

        let dissimilarity_matrix_data
            = rowwise_distances(pcloud.clone())
                .into_iter()
                .map(|x| x.into_iter().enumerate().collect_vec() )
                .collect_vec()
                .into_csr( number_of_points, number_of_points );
        let dissimilarity_matrix = & dissimilarity_matrix_data;

        let dissimilarity_value_min = OrderedFloat(0.0);        
        let dissimilarity_value_max = 
        minmax( 
                (0..number_of_points).map(
                        |x| 
                        dissimilarity_matrix.row(&x).into_iter().map(
                                |x| 
                                x.val()
                            ) 
                    ) 
            ).unwrap_or( dissimilarity_value_min.clone() );         

    
        let ring_operator = FieldRational64::new();
        let chain_complex_data = VietorisRipsComplex::new( & dissimilarity_matrix, number_of_points, ring_operator ).ok().unwrap();
        // let chain_complex_ref = & chain_complex; 
        // let chain_complex = VietorisRipsComplexArc{ arc: Arc::new(chain_complex_data) };   
        let chain_complex = Arc::new(chain_complex_data);           
        let mut row_index_vec = chain_complex.cliques_in_row_reduction_order(max_homology_dimension);
        let mut column_index_vec = chain_complex.cliques_in_row_reduction_order(max_homology_dimension+1);
        
        row_index_vec.reverse();
        column_index_vec.reverse();

   
        // CHECK SIMPLICES ARE FORMATTED CORRECTLY
        // ---------------------------------------------------------------

        // the chain_complex.filtration_value check that the simplex is valid (strictly sorted and a true clique)
        for simplex in row_index_vec.iter() {
            let filtration = chain_complex.filtration_value_for_clique(simplex.vertices()).ok().unwrap();
            assert!( filtration == simplex.filtration() );
        }
        for simplex in column_index_vec.iter() {
            let filtration = chain_complex.filtration_value_for_clique(simplex.vertices()).ok().unwrap();
            assert!( filtration == simplex.filtration() );
        }        


        // CHECK MATRIX IS INTERNALLY VALID
        // ---------------------------------------------------------------
        assert!(
            crate::algebra::matrices::debug::matrix_oracle_is_internally_consistent(
                chain_complex.clone(),
                row_index_vec.iter().cloned(),
                column_index_vec.iter().cloned(),                
            )
        );

        // CHECK ORDER OPERATORS ARE INTERNALLY VALID
        // ---------------------------------------------------------------        
        let mut row_index_vec_sorted = row_index_vec.clone();
        let mut column_index_vec_sorted = column_index_vec.clone();
        row_index_vec_sorted.sort();
        column_index_vec_sorted.sort();

        assert!(
            crate::algebra::matrices::debug::matrix_order_operators_are_internally_consistent(
                chain_complex.clone(), 
                row_index_vec_sorted.iter().cloned(),
                column_index_vec_sorted.iter().cloned(),                
            ).is_ok()
        );
        
      
    

    

    
        println!("starting umatch");
        let umatch = DifferentialUmatch::new(
                chain_complex.clone(), 
                min_homology_dimension,
                max_homology_dimension,
            );      
    
    
        println!("setting up to unpack");  
        let return_cycle_representatives = true;
        let return_bounding_chains = true;
        let _barcode = crate::algebra::chain_complexes::barcode::get_barcode( &umatch, return_cycle_representatives, return_bounding_chains );
    
        println!("getting intervals, reps, bounding chains");
        // let (intervals, representatives, bounding_chains ) = barcode.unwrap();
    
        // let convert_coefficients = |x: Ratio<i64> | (x.numer().clone(), x.denom().clone());
        // let convert_simplex = |x: WeightedSimplex< OrderedFloat<f64> > | x.vertices.iter().map(|y| y.clone() as usize ).collect_vec();
        
        // println!("start reformat reps");     
        // let represntatives_new = representatives.unwrap().iter().map(
        //             |x| // each x is a list of cycle reps
        //             x.iter().map(
        //                 |y| // each y is a cycle rep
        //                 y.iter().map(
        //                     |z| // each z is an entry in a cycle rep
        //                     ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
        //                 ).collect_vec()
        //             ).collect_vec()
        //         ).collect_vec();
    
        // println!("start reformat bounding chains");                 
        // let bounding_chains_new = bounding_chains.unwrap().iter().map(
        //             |x| // each x is a list of cycle reps
        //             x.iter().map(
        //                 |y| // each y is a cycle rep
        //                 y.iter().map(
        //                     |z| // each z is an entry in a cycle rep
        //                     ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
        //                 ).collect_vec()
        //             ).collect_vec()
        //         ).collect_vec();              
    }






    









  
    
    #[test]
    fn check_that_some_basic_functions_run_without_error() {


        use crate::topology::simplicial::{simplices::weighted::WeightedSimplex, from::graph_weighted::{VietorisRipsComplex}};        
        use crate::topology::point_cloud::unit_circle;    

        use crate::algebra::vectors::entries::KeyValGet;
        use crate::algebra::rings::types::native::{FieldRational64};    
        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};    
        use crate::algebra::matrices::types::third_party::IntoCSR;
        
        use crate::utilities::order::{ is_sorted_strictly, OrderOperatorByKeyCustom, OrderOperatorAutoReverse};
        use crate::utilities::iterators::general::minmax;    
        use crate::utilities::distances::{rowwise_distances};

        use std::sync::Arc;
        use ordered_float::OrderedFloat;  
        use itertools::Itertools;      

        let number_of_points = 8;
        let min_homology_dimension = 0;
        let max_homology_dimension = 2;


        let pcloud = unit_circle( number_of_points, Some(-1.0 .. 1.0));

        let dissimilarity_matrix_data
            = rowwise_distances(pcloud.clone())
                .into_iter()
                .map(|x| x.into_iter().enumerate().collect_vec() )
                .collect_vec()
                .into_csr( number_of_points, number_of_points );
        let dissimilarity_matrix = & dissimilarity_matrix_data;    
    
        let ring_operator = FieldRational64::new();
        let chain_complex_data = VietorisRipsComplex::new( & dissimilarity_matrix, number_of_points, ring_operator ).ok().unwrap();
        // let chain_complex_ref = & chain_complex; 
        // let chain_complex = VietorisRipsComplexArc{ arc: Arc::new(chain_complex_data) };   
        let chain_complex = Arc::new(chain_complex_data);           
        let mut row_index_vec = chain_complex.cliques_in_row_reduction_order(max_homology_dimension);
        let mut column_index_vec = chain_complex.cliques_in_row_reduction_order(max_homology_dimension+1);
        
        row_index_vec.reverse();
        column_index_vec.reverse();


        // CHECK SIMPLICES ARE FORMATTED CORRECTLY
        // ---------------------------------------------------------------

        // the chain_complex.filtration_value check that the simplex is valid (strictly sorted and a true clique)
        for simplex in row_index_vec.iter() {
            let filtration = chain_complex.filtration_value_for_clique(simplex.vertices()).ok().unwrap();
            assert!( filtration == simplex.filtration() );
        }
        for simplex in column_index_vec.iter() {
            let filtration = chain_complex.filtration_value_for_clique(simplex.vertices()).ok().unwrap();
            assert!( filtration == simplex.filtration() );
        }       

        // CHECK MATRIX INTERNALLY VALID
        // ---------------------------------------------------------------        

        assert!(
            crate::algebra::matrices::debug::matrix_oracle_is_internally_consistent(
                chain_complex.clone(),
                row_index_vec.iter().cloned(),
                column_index_vec.iter().cloned(),                
            )
        );
        

        // CHECK ORDER OPERATORS ARE INTERNALLY VALID
        // ---------------------------------------------------------------        
        let mut row_index_vec_sorted = row_index_vec.clone();
        let mut column_index_vec_sorted = column_index_vec.clone();
        row_index_vec_sorted.sort();
        column_index_vec_sorted.sort();

        assert!(
            crate::algebra::matrices::debug::matrix_order_operators_are_internally_consistent(
                chain_complex.clone(), 
                row_index_vec_sorted.iter().cloned(),
                column_index_vec_sorted.iter().cloned(),                
            ).is_ok()
        ); 
    

        println!("starting umatch");
        let iter_row_index = chain_complex.cliques_in_row_reduction_order(max_homology_dimension);  
        let differential_umatch = DifferentialUmatch::new(
            chain_complex.clone(),
            min_homology_dimension,
            max_homology_dimension,
        );   
 
    
    
        println!("setting up to unpack");  
        let return_cycle_representatives = true;
        let return_bounding_chains = true;
        let _barcode = crate::algebra::chain_complexes::barcode::get_barcode( 
            &differential_umatch,
            return_cycle_representatives,
            return_bounding_chains,
        );
    
        println!("getting intervals, reps, bounding chains");
        // let (intervals, representatives, bounding_chains ) = barcode.unwrap();
    
        // let convert_coefficients = |x: Ratio<i64> | (x.numer().clone(), x.denom().clone());
        // let convert_simplex = |x: WeightedSimplex< OrderedFloat<f64> > | x.vertices.iter().map(|y| y.clone() as usize ).collect_vec();
        
        // println!("start reformat reps");     
        // let represntatives_new = representatives.unwrap().iter().map(
        //             |x| // each x is a list of cycle reps
        //             x.iter().map(
        //                 |y| // each y is a cycle rep
        //                 y.iter().map(
        //                     |z| // each z is an entry in a cycle rep
        //                     ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
        //                 ).collect_vec()
        //             ).collect_vec()
        //         ).collect_vec();
    
        // println!("start reformat bounding chains");                 
        // let bounding_chains_new = bounding_chains.unwrap().iter().map(
        //             |x| // each x is a list of cycle reps
        //             x.iter().map(
        //                 |y| // each y is a cycle rep
        //                 y.iter().map(
        //                     |z| // each z is an entry in a cycle rep
        //                     ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
        //                 ).collect_vec()
        //             ).collect_vec()
        //         ).collect_vec();              
    }




    #[test]
    fn test_empty_simplex_iter() {
        use sprs::CsMatBase;
        use ordered_float::OrderedFloat;
        use crate::topology::simplicial::from::graph_weighted::AgileSimplexIteratorLexicographicOrder;
        use crate::algebra::matrices::query::MatrixOracle;

        
        let dissimilarity_matrix = CsMatBase::new( (2,2), vec![0,1,2], vec![0,1], vec![OrderedFloat(0.0),OrderedFloat(0.0)] );

        println!("entry (0,1) = {:?}", (& dissimilarity_matrix).structural_nonzero_entry(&0, &1));
            
        let simplices: Vec<_>    =   AgileSimplexIteratorLexicographicOrder::new( 
                            & dissimilarity_matrix, 
                            2,          //  matrix size
                            1,          //  dimension
                            false,  // include the empty simplex
                        ).collect();
        println!("generated simplices = {:?}", &simplices );
        let empty =  Vec::< WeightedSimplex<OrderedFloat<f64>> >::new();
                    
        assert!( empty.eq( &simplices ) )        
    }



    
    #[test]
    fn test_simplex_iterators () {
        use crate::topology::simplicial::from::graph_weighted::{AgileBoundaryIteratorLexicographicOrder, AgileCoboundaryIteratorFiltrationOrder, VietorisRipsComplex};
        use crate::topology::simplicial::simplices::weighted::WeightedSimplex;   
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::rings::types::native::FieldFloat64;
        use ordered_float::OrderedFloat;
        use std::sync::Arc;

        // define the ring operator
        let ring_operator = FieldFloat64::new();        

        // define the dissimilarity matrix
        // (this is a convenient way to build a VecOfVec from a "dense matrix")
        let dissimilarity_matrix = VecOfVec::from_ragged(
            & vec![
                vec![OrderedFloat(0.0), OrderedFloat(1.0), OrderedFloat(2.0), OrderedFloat(3.0)],
                vec![OrderedFloat(1.0), OrderedFloat(0.0), OrderedFloat(1.0), OrderedFloat(2.0)],
                vec![OrderedFloat(2.0), OrderedFloat(1.0), OrderedFloat(0.0), OrderedFloat(1.0)],
                vec![OrderedFloat(3.0), OrderedFloat(2.0), OrderedFloat(1.0), OrderedFloat(1.0)],                
            ]
        );

        // define the Vietoris-Rips complex
        let vietoris_rips_complex = Arc::new( 
            VietorisRipsComplex::new(
                & dissimilarity_matrix,
                3, // dissimilarity_matrix_size
                ring_operator,
            ).ok().unwrap()
        );

        // Check a boundary iterator
        // ---------------------------------------------------------------

        // define the simplex
        let simplex = WeightedSimplex {
            vertices: vec![0, 1, 2],
            weight: OrderedFloat(2.0),
        };         

        // create the iterator
        let boundary = AgileBoundaryIteratorLexicographicOrder::new(
            vietoris_rips_complex.clone(),
            simplex,
            false, // this indicates we don't consider the emtpy simplex to be a "real" simplex
        );   

        // check the results
        assert!( 
            boundary.eq(
                vec![
                    ( WeightedSimplex { vertices: vec![0, 1], weight: OrderedFloat(1.0) },  1.0 ),
                    ( WeightedSimplex { vertices: vec![0, 2], weight: OrderedFloat(2.0) }, -1.0 ),
                    ( WeightedSimplex { vertices: vec![1, 2], weight: OrderedFloat(1.0) },  1.0 ),
                ]
            )         
        );

        // Check a coboundary iterator
        // ---------------------------------------------------------------   

        // Here we use the constructor from the Vietoris-Rips complex, which is more convenient  

        // define the simplex
        let simplex = WeightedSimplex {
            vertices: vec![0, 2],
            weight: OrderedFloat(2.0),
        };  

        // create the iterator
        let coboundary = vietoris_rips_complex.row( & simplex );

        // check the results
        assert!( 
            coboundary.eq(
                vec![
                    ( WeightedSimplex { vertices: vec![0, 1, 2], weight: OrderedFloat(2.0) },  -1.0 ),
                    ( WeightedSimplex { vertices: vec![0, 2, 3], weight: OrderedFloat(3.0) },   1.0 ),
                ]
            )         
        );
    } 





}    