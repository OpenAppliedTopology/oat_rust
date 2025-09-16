//! Random object generators


use itertools::Itertools;
use ordered_float::OrderedFloat;
use rand::Rng;

use crate::utilities::sequences_and_ordinals::SortedVec;



/// A sequence of combinations sampled with given probability.
/// 
/// We initialize an empty vector `v`.
/// For each `k` in `cardinalities` we generate the sequence of all cardinality-`k` subsets of  `{0.. size_of_ambient_set}`, in lexicographic order.
/// Subsets are represented as sorted vectors.
/// For each subset `s`, we push `s` into `v` with probability `probability`.
/// The vector `v` is returned to the user when this process terminates.
pub fn random_sequences< T: Iterator<Item=usize> >( 
        size_of_ambient_set: usize, 
        cardinalities: T, 
        probability: f64 
    ) -> Vec< SortedVec< usize > > 
{
    use rand::distributions::{Bernoulli, Distribution};

    let d = Bernoulli::new(  probability ).unwrap();
    let mut rng = rand::thread_rng();    
    
    let mut combinations    =   vec![];
    for cardinality in cardinalities {
        for combo in (0 .. size_of_ambient_set ).combinations( cardinality ) {
            if d.sample(&mut rng) {
                combinations.push( SortedVec::new(combo).ok().unwrap() )
            }
        }
    }  

    combinations

}



/// Generates an m x n matrix with entries sampled from the uniform 
/// distribution on [0,1]
pub fn random_matrix( m: usize, n:usize ) -> Vec< Vec< OrderedFloat<f64> > > {
    let mut rng = rand::thread_rng();

    let mut matrix = vec![ Vec::with_capacity(n); m ];
    for rowind in 0..m {
        for _ in 0..n {
            matrix[rowind].push( OrderedFloat(rng.gen_range(0.0 .. 1.0)) )
        }
    }
    matrix
}

/// Generates an m x n symmetric matrix of form (D+D')/2, where D is
/// a random matrix with entries sampled from the uniform  distribution on [0,1]
pub fn random_symmetric_matrix( m: usize ) -> Vec< Vec< OrderedFloat<f64> > > {

    let mut matrix = random_matrix(m,m);
    for p in 0..m {
        for q in p..m {
            let t = (matrix[p][q] + matrix[q][p]) / 2.0 ;
            matrix[p][q] = t;
            matrix[q][p] = t;
        }        
    }
    matrix
}

/// Generates an m x n symmetric matrix of form (D+D')/2, where D is
/// a random matrix with entries sampled from the uniform  distribution on [0,1]
pub fn random_symmetric_matrix_zero_diag( m: usize ) -> Vec< Vec< OrderedFloat<f64> > > {

    let mut matrix = random_matrix(m,m);
    for p in 0..m {
        for q in p..m {
            let t = (matrix[p][q] + matrix[q][p]) / 2.0 ;
            matrix[p][q] = t;
            matrix[q][p] = t;
        }        
    }
    for p in 0 .. m { matrix[p][p] = OrderedFloat(0.0) }
    matrix
}