//! Some functions that return matrices, e.g. `f(m) = random_m_by_m_sparse_matrix`.

use crate::matrices::types::vec_of_vec::sorted::VecOfVec;
use rand::Rng;                                          // we use this module to generate random elements
use rand::distributions::{Bernoulli, Distribution};     // we use this module to generate random elements



/// Generate a random invertible upper triangular matrix with integer coefficients.
/// 
/// - Entries on the diagonal are drawn from iid uniform distributions on [1, .., `modulus`)
/// - Each entry strictly above the diagonal is either (i) strucutrually zero, 
/// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
/// Algebraically nonzero entries are drawn iid from the uniform distribution on [1, .., `modulus`).
/// 
/// Nonzero entries are drawn in an iid fasion from [1, .., `modulus`).
pub fn random_upper_triangular( matrix_size: usize, modulus: usize) -> VecOfVec< usize, usize > {

    let mut rng = rand::thread_rng(); // this generates random integers
    let mut vec_of_vec = vec![];
    for row_index in 0 .. matrix_size {
        let coefficient_leading         =   rng.gen_range( 1 .. modulus );
        let mut new_vec     =   vec![ (row_index, coefficient_leading) ]; // start with a nonzero diagonal element            
        for q in row_index+1 .. matrix_size { // fill out the rest of this row of the matrix
            let coefficient   = rng.gen_range( 0 .. modulus );
            let flag = rng.gen_range(0usize .. 3usize);
            if      flag == 0 { new_vec.push( ( q, 0 )           ) }
            else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
            else              { continue }
        }
        vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
    }

    VecOfVec::new( vec_of_vec ) // formally wrap the matrix in a VecOfVec struct
}

/// Generate a random invertible upper triangular matrix with integer coefficients.
/// 
/// - Entries on the diagonal are equal to 1.
/// - Each entry strictly above the diagonal is either (i) strucutrually zero, 
/// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
/// Algebraically nonzero entries are drawn iid from the uniform distribution on [1, .., `modulus`).
/// 
/// Nonzero entries are drawn in an iid fasion from [1, .., `modulus`).
pub fn random_upper_unitriangular( matrix_size: usize, modulus: usize) -> VecOfVec< usize, usize > {

    let mut rng = rand::thread_rng(); // this generates random integers
    let mut vec_of_vec = vec![];
    for row_index in 0 .. matrix_size {
        let coefficient_leading         =   1;
        let mut new_vec     =   vec![ (row_index, coefficient_leading) ]; // start with a nonzero diagonal element            
        for q in row_index+1 .. matrix_size { // fill out the rest of this row of the matrix
            let coefficient   = rng.gen_range( 0 .. modulus );
            let flag = rng.gen_range(0usize .. 3usize);
            if      flag == 0 { new_vec.push( ( q, 0 )           ) }
            else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
            else              { continue }
        }
        vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
    }

    VecOfVec::new( vec_of_vec ) // formally wrap the matrix in a VecOfVec struct
}



/// Generate a random `VecOfVec< usize, usize >`
/// 
/// Nonzero entries are drawn in an iid fasion from [0, .., `modulus`) if `allow_nonstructural_zero == true` and from [1, .., `modulus`)
/// if `allow_nonstructural_zero == false`.
pub fn random_vec_of_vec_simple( 
            num_indices_major:          usize, 
            num_indices_minor:          usize, 
            approximate_density:        f64, 
            modulus:                    usize,
            allow_nonstructural_zero:   bool,
        ) 
        -> 
        VecOfVec< usize, usize > 
{

    let mut rng = rand::thread_rng(); // this generates random integers

    let d = Bernoulli::new( approximate_density ).unwrap(); // bernoulli random variable returning `true` with probability `approximate_density`
    let v = d.sample(&mut rand::thread_rng());
    println!("{} is from a Bernoulli distribution", v);

    let mut vecvec = Vec::new(); // initialize empty vector of vectors
    for row_index in 0 .. num_indices_major {
        vecvec.push( Vec::new() ); // push a new vector representing a row
        for column_index in 0 .. num_indices_minor { // for each column index
            if d.sample( &mut rand::thread_rng() ) { // add a structural nonzero entry with probability `approximate_density`
                let coefficient   = match allow_nonstructural_zero{ 
                    true => { rng.gen_range( 0 .. modulus ) },
                    false => { rng.gen_range( 1 .. modulus ) }
                };
                vecvec[ row_index ].push( (column_index, coefficient) );
            }
        }
        vecvec[ row_index ].shrink_to_fit();
    }

    return VecOfVec::new( vecvec )

}





//  RANDOM MATRICES
//  ===========================================================================


/// Generate a random upper-triangular row-major matrix of size `matrix_size` x `matrix_size`, with coefficients
/// in `{ 0, .., p-1 }`, having no zeros on the diagonal.
/// 
/// - Diagonal entries are drawn independently from the uniform distribution on `{ 1, .., p-1 }`.
/// - Entries strictly above the diagonal are equally likely to be (i) structurally zero (and therefore algebraically) zero, 
/// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
/// In case (iii), entries are drawn independently from the uniform distribution on `{ 1, .., p-1 }`.
pub fn random_upper_triangular_matrix_mod_p
        ( matrix_size: usize, modulus: usize ) 
        ->
        VecOfVec< usize, usize >
{
    use rand::Rng;        // we use this module to generate random elements

    let mut rng = rand::thread_rng(); // this generates random integers
    let mut vec_of_vec = vec![];
    for row_index in 0 .. matrix_size {
        let coefficient_leading         =   rng.gen_range( 1 .. modulus );
        let mut new_vec     =   vec![ (row_index, coefficient_leading) ]; // start with a nonzero diagonal element            
        for q in row_index+1 .. matrix_size { // fill out the rest of this row of the matrix
            let coefficient   = rng.gen_range( 1 .. modulus );
            let flag = rng.gen_range(0usize .. 3usize);
            if      flag == 0 { new_vec.push( ( q, 0 )           ) }
            else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
            else              { continue }
        }
        vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
    }

    return VecOfVec::new(vec_of_vec); // formally wrap the matrix in a VecOfVec struct
}


/// Generate a random row-major matrix of size `num_rows` x `num_cols` with coefficients
/// in `{ 0, .., p-1 }`.
/// 
/// Each entry is equally likely to be (i) structurally zero, 
/// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
/// In case (iii), entries are drawn independently from the uniform distribution on `{ 1, .., p-1 }`.
pub fn random_m_by_n_matrix
            ( num_rows: usize, num_cols: usize, modulus: usize ) 
            ->
            VecOfVec< usize, usize >
{
    use rand::Rng;        // we use this module to generate random elements

    let mut rng = rand::thread_rng(); // this generates random integers
    let mut vec_of_vec = vec![];
    for row_index in 0 .. num_rows {
        let mut new_vec     =   vec![]; // start with an empty vector
        for q in row_index+1 .. num_cols { // fill it in
            let coefficient   = rng.gen_range( 1 .. modulus );
            let flag = rng.gen_range(0usize .. 3usize);
            if      flag == 0 { new_vec.push( ( q, 0 )           ) }
            else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
            else              { continue }
        }
        vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
    }

    return VecOfVec::new(vec_of_vec); // formally wrap the matrix in a VecOfVec struct
}