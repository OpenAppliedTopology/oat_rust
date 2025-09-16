//! Bars and barcodes for persistence modules

use itertools::Itertools;
use ordered_float::OrderedFloat;
use rand::seq::index;

use crate::algebra::chain_complexes::{ChainComplex, FilteredChainComplex};
use crate::algebra::matrices::operations::umatch::differential::{DifferentialComb, DifferentialUmatch};
use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::vectors::entries::KeyValPair;
use crate::algebra::matrices::operations::multiply::multiply_column_vector_with_matrix_and_return_reversed;
use crate::algebra::matrices::operations::umatch::row_major::Umatch;
use crate::algebra::matrices::query::{MatrixAlgebra, MatrixOracle};
use crate::algebra::rings::traits::DivisionRingOperations;

use crate::algebra::vectors::operations::VectorOperations;
use crate::utilities::order::{JudgeOrder, JudgePartialOrder};

// use polars::prelude::*;

use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::hash::Hash;

use derive_getters::{Getters, Dissolve};




/// A single bar in a persistent homology barcode (optionally including generators)
/// 
/// This object stores information about
/// 
/// - `id_number`: an integer used to distinguish this bar from all other bars in the barcode
/// - `birth`: the left endpoint of the bar
/// - `death`: the right endpoint of the bar (if it exists)
/// - `dimension`: the dimension of the homology group
/// - `birth_column`: in the special case of simplicial complexes, this is called the "birth simplex"
/// - `death_column`: in the special case of simplicial complexes, this is called the "death simplex" (if it exists)
/// - `cycle_representative`: a cycle representative in homology, which is born when the bar is born
/// - `bounding_chain`: a bounding chain, which is born when the bar dies
/// 
/// The primary motivation for this object is convenience -- reducing the number and complexity of
/// method calls needed to access commonly used properties of the barcode.
#[derive(Getters,Dissolve)]
#[derive(Debug,Clone,PartialEq, Eq, Hash)]
pub struct Bar< Index, Entry > {
    id_number:              usize,
    #[getter(skip)]
    dimension:              isize,
    birth:                  OrderedFloat< f64 >,    
    birth_column:           Index,
    death:                  Option< OrderedFloat< f64 > >,    
    death_column:           Option< Index >,
    cycle_representative:   Option< Vec< Entry > >,
    bounding_chain:         Option< Vec< Entry > >,
}

impl < Index, Entry >

    Bar
        < Index, Entry > {            
    pub fn dimension( &self ) -> isize { self.dimension }
    
    
    /// The birth fitration value, as an `OrderedFloat<f64>`
    pub fn birth_ordf64(&self)  -> OrderedFloat<f64> { self.birth }
    /// The death fitration value, as an `OrderedFloat<f64>`
    pub fn death_ordf64(&self)  -> OrderedFloat<f64> { self.death.unwrap_or( OrderedFloat( f64::INFINITY ) ) }    
    /// The length of the interval, as an `OrderedFloat<f64>`
    pub fn length_ordf64(&self)  -> OrderedFloat<f64> { self.death_ordf64() - self.birth_ordf64() }        
    /// The birth-death pair as `(OrderedFloat<f64>,OrderedFloat<f64>)`
    pub fn interval_ordf64( &self ) -> (OrderedFloat<f64>,OrderedFloat<f64>) {  (self.birth_ordf64(), self.death_ordf64()) }    

    
    /// The birth fitration value, as an `f64`
    pub fn birth_f64( &self ) -> f64 { self.birth_ordf64().into_inner() }
    /// The death fitration value, as an `f64`
    pub fn death_f64( &self ) -> f64 { self.death_ordf64().into_inner() }
    /// The length of the interval, as an `f64`
    pub fn length_f64(&self)  -> f64 { self.length_ordf64().into_inner() }            
    /// The birth-death pair as `(f64,f64)`
    pub fn interval_f64( &self ) -> (f64,f64) {  (self.birth_f64(), self.death_f64()) }    
    
}

/// Wrapper for a vector of [`Bar`]s, with handy utilities.
#[derive(Debug,Clone,PartialEq, Eq, Hash)]
pub struct Barcode
            < Index, Entry >
    { bars: Vec< Bar< Index, Entry > > }

impl < Index, Entry > 
    
    Barcode 
        < Index, Entry > {

    /// Construct a barcode from a list of bars
    pub fn new< I >( bars: I ) -> Self 
        where 
            I: IntoIterator< Item=Bar<Index,Entry> >
    { Barcode{ bars: bars.into_iter().collect() } }

    /// The number of bars in the barcode
    pub fn len( &self ) -> usize { self.bars.len() }

    /// An iterator that runs over all the bars
    pub fn iter(&self) -> std::slice::Iter<'_, Bar<Index, Entry>> {
        self.bars.iter()
    }

    /// Immutable reference to the internally stored vector of [`Bar`]s
    pub fn bars( &self ) -> & Vec< Bar< Index, Entry > > { &self.bars }

    /// Vector containing references to all bars of dimension `dim`
    pub fn bars_in_dim( &self, dim:isize ) -> Vec< & Bar< Index, Entry > > { 
        self.bars.iter().filter(|x| x.dimension()==dim ).collect() 
    }

    /// Immutable reference to the `bafr_id_number`th internally stored bar.
    pub fn bar( &self, bar_id_number: usize ) -> & Bar< Index, Entry > { &self.bars[bar_id_number] }

    /// Return a sorted list of all endpoints of intervals in the barcode.
    pub fn endpoints_ordf64_within_dim( &self, dim: isize ) -> Vec< OrderedFloat<f64> > {
        let iter_a = self.bars.iter().filter(|x|x.dimension()==dim).map(|x| x.birth_ordf64() );
        let iter_b = self.bars.iter().filter(|x|x.dimension()==dim).map(|x| x.death_ordf64() );        
        let mut finite_endpoints = HashSet::new();
        for val in iter_a { finite_endpoints.insert( val ); }
        for val in iter_b { finite_endpoints.insert( val ); }
        let mut finite_endpoints = finite_endpoints.drain().collect_vec();
        finite_endpoints.sort();
        finite_endpoints
    }

    /// Return a sorted list of all endpoints of intervals in the barcode
    /// 
    /// The search runs over all bars in the barcode
    pub fn endpoints_ordf64( &self ) -> Vec< OrderedFloat<f64> > {
        let iter_a = self.bars.iter().map(|x| x.birth_ordf64() );
        let iter_b = self.bars.iter().map(|x| x.death_ordf64() );        
        let mut finite_endpoints = HashSet::new();
        for val in iter_a { finite_endpoints.insert( val ); }
        for val in iter_b { finite_endpoints.insert( val ); }
        let mut finite_endpoints = finite_endpoints.drain().collect_vec();
        finite_endpoints.sort();
        finite_endpoints
    }    

    /// Either the `Some( maximum of the set { finite endpoints in intervals in the barcode })`,
    /// or `None`, if there exists no finite endpoint.
    /// 
    /// The search runs over bars of every dimension.
    pub fn max_finite_endpoint( &self ) -> Option< OrderedFloat<f64> > {
        let endpoints = self.endpoints_ordf64();
        let limit = &OrderedFloat(f64::INFINITY);
        return endpoints.iter().cloned().filter(|x| x < limit ).max()
    }

    /// Returns a vector of triples `(birth, death)`, where `id` is the uniue id of the bar.
    pub fn intervals_f64( &self, dim: isize ) -> Vec< ( usize, f64, f64 ) > { 
        self.bars
            .iter()
            .filter_map(|x|
                match x.dimension == dim {
                    true      => Some( ( *x.id_number(), x.birth_f64(),  x.death_f64(), ) ) ,
                    false     => None
                }                
            ) 
            .collect_vec()
    }

    /// The betti curve for homology in a given dimension
    /// 
    /// The output is formatted as a list of pairs (t0,b0) .. (tN, bN). 
    /// - The betti number for each half-closed interval [ti,ti+1) is bi
    /// - The betti number of [-infinity, t0) is 0
    /// - The betti number of [tN, +infinity) is 0
    pub fn betti_curve( &self, dim: isize ) -> Vec< ( OrderedFloat<f64>, usize ) > {

        let endpoints = self.endpoints_ordf64_within_dim(dim); // a list of all interval endpoints
        let mut bettis  =  vec![0; endpoints.len()]; // initialize with zeros
        for bar in self.bars() {
            if dim != bar.dimension() { continue }                      // exclude if dimension is wrong
            let birth = bar.birth_ordf64();
            for (ordinal, endpoint) in endpoints.iter().enumerate() {
                if endpoint < &birth { continue }                       // exclude if birth too late
                if let & Some(death) = bar.death() {
                    if death <= *endpoint { break }                     // exclude if death too early
                }
                bettis[ordinal] += 1;                                   // otherwise increment betti number
            }
        }
        return endpoints.into_iter()
                .zip(bettis)
                .collect()
    }
}

/// Extract the barcode of a filtered chain complex from a U-match factorization of its boundary matrix.
/// 
/// The output is formatted as a vector of ['Bar']s.
/// 
/// - `simplices_in_row_reduction_order` is an iterator that runs over every row-index `p` of the boundary matrix
/// such that dim(`p`) < max { dim(`q`) : `q` is a row-index }.
/// - `get_dimension` returns the dimension of each index
/// - `get_filtration_value` returns the birth time of each index
/// - `return_cycle_representatives` determines wheter cycle repersentatives are calculated
/// - `return_bounding_chains` determines whether bounding chains are returned
/// 
/// **Note** if both bounding chains and cycle representatives are returned, then this function
/// internally checks that the boundary of the bounding chain is indeed the corresponding
/// cycle representative.
pub fn get_barcode< 
            BoundaryMatrix, 
            OrderOperatorForRowAndColumnEntries, 
            OrderOperatorForRowAndColumnIndices, 
            IndexForRowsAndColumns,
            EntryForRowsAndColumns,           
        >
    ( 
        differential_umatch:            & DifferentialUmatch<  BoundaryMatrix >, 
        return_cycle_representatives:   bool,
        return_bounding_chains:         bool,
    ) 
    -> Barcode<
                BoundaryMatrix::RowIndex, 
                BoundaryMatrix::ColumnEntry 
            > 
    where
        BoundaryMatrix:                             MatrixAlgebra< 
                                                        RowEntry                        =   EntryForRowsAndColumns, 
                                                        RowIndex                        =   IndexForRowsAndColumns, 
                                                        OrderOperatorForRowEntries      =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForRowIndices      =   OrderOperatorForRowAndColumnIndices,
                                                        ColumnEntry                     =   EntryForRowsAndColumns, 
                                                        ColumnIndex                     =   IndexForRowsAndColumns,                                                        
                                                        OrderOperatorForColumnEntries   =   OrderOperatorForRowAndColumnEntries,
                                                        OrderOperatorForColumnIndices   =   OrderOperatorForRowAndColumnIndices,                                                        
                                                    >
                                                    + MatrixOracleOperations
                                                    + ChainComplex 
                                                    + FilteredChainComplex<
                                                            FiltrationValue = OrderedFloat<f64>,
                                                        >,
        OrderOperatorForRowAndColumnEntries:        Clone + JudgePartialOrder< EntryForRowsAndColumns >,
        OrderOperatorForRowAndColumnIndices:        Clone + JudgeOrder< IndexForRowsAndColumns >,
        IndexForRowsAndColumns:                     Clone + Debug + Hash + Eq, // required for the hashing performed by the generalized matching array   
        BoundaryMatrix::Coefficient:            Debug,       
        BoundaryMatrix::RowEntry:               KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >, 
        BoundaryMatrix::ColumnEntry:            KeyValPair< Key = IndexForRowsAndColumns, Val = BoundaryMatrix::Coefficient >,  
        BoundaryMatrix::RingOperator:           DivisionRingOperations,   
        EntryForRowsAndColumns:                     Clone + Debug + PartialEq,                
       
{    
    let differential_comb                 =  differential_umatch.differential_comb(); // DifferentialComb::new( umatch );    
    let asymmetric_umatch = differential_umatch.asymmetric_umatch(); // get a copy of the asymmetric U-match
    let matching_matrix = asymmetric_umatch.generalized_matching_matrix_ref(); 
    let boundary_matrix = differential_umatch.boundary_matrix();
    let mut get_filtration_value = | index: & IndexForRowsAndColumns | {
        boundary_matrix.filtration_value_for_basis_vector_with_index( index )
            .expect("barcode: filtration value for index not found")
    };

    let mut barcode                 =   Vec::new();
    for row_index in differential_umatch.indices_in_homologically_valid_dimensions() {
        if matching_matrix.has_a_match_for_column_index( & row_index ) { continue } // in this case we ignore the key, since it doesn't correspond to a cycle

        let death_column = matching_matrix.column_index_for_row_index( & row_index );
        let death = death_column.as_ref().map( &mut get_filtration_value );
        let birth = get_filtration_value( & row_index );        
        
        if death == Some( birth ) { continue } // if birth = death, don't include the bar

        let cycle_representative   =   match return_cycle_representatives {
            false   =>  { None },
            true    =>  { Some( differential_comb.column( & row_index ).collect_vec() ) }
        };

        // compute bounding chain and verify that it does bound the 
        let bounding_chain  =   match return_bounding_chains {
            false   =>  { None },
            true    =>  {
                            if let Some( column_index ) = matching_matrix.column_index_for_row_index( & row_index ) {
                                // the chain in question is a boundary
                                // the differential comb gives us a pair of columns x, y, and a scalar alpha, such that
                                // alpha * y = Dx
                                // here y is the cycle representative, and x is (almost) the bounding chain
                                // we just have to scale x by alpha_inverse
                                let ring_operator = boundary_matrix.ring_operator();
                                let alpha = matching_matrix
                                                                .coefficient_for_row_index( & row_index );
                                let alpha_inverse = ring_operator.invert( alpha );
                                let bounding_chain  =   differential_comb
                                                            .column( &column_index ) // get column x
                                                            .scale_by(  // scale by alpha_inverse
                                                                alpha_inverse,
                                                                ring_operator
                                                            )
                                                            .collect_vec();                           
                                let boundary_vec    =   boundary_matrix
                                                            .multiply_with_column_vector( & bounding_chain )
                                                            .collect_vec();
                                
                                if return_bounding_chains && return_cycle_representatives {
                                    assert_eq!( cycle_representative, Some(boundary_vec) ); // check that the bounding chain bounds the right cycle
                                }
                                Some( bounding_chain )
                            } else {
                                None
                            }
                        }
        };

        let dimension = boundary_matrix.dimension_for_basis_vector_with_index( & row_index ).ok().unwrap();
        let id_number = barcode.len();
        let birth_column    =   row_index;
        let bar     =   Bar{ 
                                                            id_number, 
                                                            dimension,
                                                            birth,
                                                            birth_column,
                                                            death,
                                                            death_column,
                                                            cycle_representative,
                                                            bounding_chain 
                                                        };

        barcode.push(bar);
    }
    Barcode{ bars: barcode }
}









//  ===========================================================================
//  ===========================================================================
//  TESTS
//  ===========================================================================
//  ===========================================================================





#[cfg(test)]
mod tests {

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    use crate::algebra::vectors::entries::KeyValGet;    
    use crate::algebra::matrices::query::MatrixOracle;
    
    use crate::algebra::matrices::operations::umatch::differential::DifferentialUmatch;

    use crate::topology::point_cloud::unit_circle;
    use crate::topology::simplicial::simplices::weighted::WeightedSimplex;
    use crate::topology::simplicial::from::graph_weighted::VietorisRipsComplex;

    use crate::utilities::distances::rowwise_distances;
    use crate::utilities::iterators::general::minmax;
    

    use std::sync::Arc;
    use sprs::CsMatBase;




    #[test]
    fn test_barcode_random_symmetric_matrix() {

    
        use crate::algebra::matrices::types::{third_party::IntoCSR, vec_of_vec::sorted::VecOfVec};
        
        let number_of_points = 20;
        let min_homology_dimension =  0;                
        let max_homology_dimension =  1;

        let dissimilarity_matrix_vecvec = VecOfVec::random_symmetric_zero_diagonal_with_enclosing_radius_threshold(number_of_points);
        let dissimilarity_matrix =  Arc::new( dissimilarity_matrix_vecvec.into_csr(number_of_points, number_of_points) );                        


        let ring_operator = crate::algebra::rings::types::native::FieldRationalSize::new();
        let boundary_matrix_data = VietorisRipsComplex::new( dissimilarity_matrix, number_of_points, ring_operator ).ok().unwrap();
        let boundary_matrix = Arc::new(boundary_matrix_data);

      
            
        // println!("starting umatch");
        let factored = DifferentialUmatch::new(
                    boundary_matrix, 
                    min_homology_dimension,
                    max_homology_dimension,
                );
    
        // println!("setting up to unpack");  
        let return_cycle_representatives = true;
        let return_bounding_chains = true;         
        let _barcode = get_barcode( &factored, return_cycle_representatives , return_bounding_chains, );
        // for bar in barcode.bars() {  println!("{:#?}", bar)}
    }    


    #[test]
    fn test_barcode_random_symmetric_matrix_with_theshold() {

    
        use crate::algebra::matrices::types::third_party::IntoCSR;
        
        let number_of_points = 20;    
        let min_homology_dimension =  0;          
        let max_homology_dimension =  1;
        let dissimilarity_value_max = OrderedFloat(1.1);


        let dissimilarity_matrix_data = crate::utilities::random::random_symmetric_matrix_zero_diag(number_of_points).into_csr( number_of_points, number_of_points );
        let mut tri = sprs::TriMat::new( (number_of_points, number_of_points) );
        for (v, (i,j) ) in dissimilarity_matrix_data.iter() { 
            if v.clone() <= dissimilarity_value_max{
                tri.add_triplet(i, j, v.clone()); 
            }            
        }
        let dissimilarity_matrix_data: CsMatBase<_,_,_,_,_> = tri.to_csr();
        let dissimilarity_matrix = & dissimilarity_matrix_data;

        for i in 0 .. number_of_points {
            for j in i .. number_of_points {
                assert_eq!( dissimilarity_matrix.structural_nonzero_entry( & i, & j ), dissimilarity_matrix.structural_nonzero_entry( & j, & i ) );
            }
        }
     
        let ring_operator = crate::algebra::rings::types::native::FieldRationalSize::new();
        let boundary_matrix_data = VietorisRipsComplex::new( dissimilarity_matrix, number_of_points, ring_operator ).ok().unwrap();
        let boundary_matrix = Arc::new(boundary_matrix_data);
            
        // println!("starting umatch");
        let factored = DifferentialUmatch::new(
                    boundary_matrix, 
                    min_homology_dimension,
                    max_homology_dimension,
                );
    
        // println!("setting up to unpack");  
        let return_cycle_representatives = true;
        let return_bounding_chains = true;         
        let _barcode = get_barcode( &factored, return_cycle_representatives , return_bounding_chains, );
        // for bar in barcode.bars() {  println!("{:#?}", bar)}
    }    


    #[test]
    fn test_barcode_circle() {

    
        use crate::algebra::matrices::types::third_party::IntoCSR;
        
        let number_of_points = 50;
        let min_homology_dimension =  0;                
        let max_homology_dimension =  1;

        let cloud = unit_circle(number_of_points, Some(0.0 .. 0.1) );
        let dissimilarity_matrix_data = rowwise_distances(cloud).into_csr( number_of_points, number_of_points );
        let dissimilarity_matrix = & dissimilarity_matrix_data;

        for i in 0 .. number_of_points {
            for j in i .. number_of_points {
                assert_eq!( dissimilarity_matrix.structural_nonzero_entry( & i, & j), dissimilarity_matrix.structural_nonzero_entry( & j, & i ) );
            }
        }
     
        let ring_operator = crate::algebra::rings::types::native::FieldRationalSize::new();

        let boundary_matrix_data = VietorisRipsComplex::new( dissimilarity_matrix, number_of_points, ring_operator ).ok().unwrap();
        let boundary_matrix = Arc::new(boundary_matrix_data);

            
        // println!("starting umatch");
        let factored = DifferentialUmatch::new(
                    boundary_matrix, 
                    min_homology_dimension,
                    max_homology_dimension,
                );
    
        // println!("setting up to unpack");  
        let return_cycle_representatives = true;
        let return_bounding_chains = true;         
        let _barcode = get_barcode( &factored, return_cycle_representatives , return_bounding_chains, );
        // for bar in barcode.bars() {  println!("{:#?}", bar)}
        // println!("number of pairs: {:?}", factored.umatch().generalized_matching_matrix_ref().number_of_structural_nonzeros() );
    }        




}









    // /// The betti numbers of a chain complex.
    // /// 
    // /// The betti numbers of a chain complex with boundary matrix `D` can be computed as 
    // /// follows, c.f. [U-match Factorization](https://arxiv.org/abs/2108.08831) (recall
    // /// that `D` has rows and columns for chains in every dimension): (i) obtain
    // /// a U-match factorization of `D` with matching matrix `M`, (ii) the `k`th betti number
    // /// is the number of indices `i` such that `M[i,:]=0` and `M[:,i]=0`, and index `i`
    // /// corrsponds to a chain of dimension `k`.
    // /// 
    // /// This function computes betti numbers according to this formula, assuming that `Self`
    // /// is the matching matrix of some boundary matrix `D`.  Argument `I` is an iterator
    // /// that runs over all the row (equivalently column) indices of `D`.  If you only need
    // /// homology up through dimension `d`, you only need to include indices for chains
    // /// of dimension `d` and below.  Argument `get_dimension`
    // /// is a function that returns the dimension of the chain associated with each index.
    // /// 
    // /// *Remark* under the hood this method is identical to [histogram_of_unmatched_indices].  It is 
    // /// included as a separate function primarily for the purpose of readability.
    // pub fn homology_dimensions< I, F >( &self, iter_IndexForRowsAndColumns: I, get_dimension: F ) 
    //         -> 
    //         Vec< usize > 
    //     where
    //         I: Iterator<Item=IndexForRowsAndColumns>, 
    //         F: FnMut(IndexForRowsAndColumns)->usize  
    //     {
    //     return self.generalized_matching_matrix_ref().histogram_of_unmatched_indices(iter_IndexForRowsAndColumns, get_dimension)
    // }   