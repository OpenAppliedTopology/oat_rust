//! Variant of vec-of-vec where rows iterate over entries of form `& (column_index,coefficient)`, not `(column_index,coefficient)`
//! 
//! See the [parent module](super) for details on what vector-of-vectors format means.
//! 
//! This module is nearly identical to [sorted](super::sorted), but the oracle returns references to entries, not clones.

use crate::algebra::vectors::entries::KeyValGet;
use crate::algebra::matrices::query::MatrixOracle;

                                    use crate::utilities::binary_search::{find_sorted_binary_oracle};
                                    use crate::utilities::order::{JudgePartialOrder, is_sorted_strictly, OrderOperatorByKey, };
                                    


use rand::Rng;                                          // we use this module to generate random elements
use rand::distributions::{Bernoulli, Distribution};     // we use this module to generate random elements

use std::iter::Rev;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::slice::Iter;
use itertools::Itertools;

use super::super::bimajor::MatrixBimajorData;



/// Vector-of-vectors sparse matrix format, with internal vectors stored in sorted order
/// 
/// See the [parent module](super) for details on what vector-of-vectors format means.
/// 
/// - This VecOfVec structure is nearly identical to [sorted::VecOfVec](super::sorted::VecOfVec), but the oracle returns references to entries, not clones.
/// - Entries in each internal vector are stored in *strictly* ascending order of index (no repeat indices).
/// - Order of incides is determined by Rust's `PartialOrd` trait.
/// 
/// 
/// # See also
/// 
/// If you need a customized order for entries, see the [variant for custom orders](super::sorted_custom).
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted_ref::VecOfVec;
/// use oat_rust::algebra::matrices::query::*;
/// use std::marker::PhantomData;
/// 
/// // Create a new matrix.
/// let matrix  =   VecOfVec::new(
///                     vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                 );
/// 
/// 
/// 
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VecOfVec
            < ColumnIndex, Coefficient >

{
    vec_of_vec:         Vec< Vec< ( ColumnIndex, Coefficient ) > >,
    order_operator:     OrderOperatorByKey,
    // phantom_lifetime:   PhantomData< &'a ColumnIndex >,
}

impl < ColumnIndex, Coefficient >
        
        VecOfVec
            < ColumnIndex, Coefficient > 

{
    /// Make a new `VecOfVec`.  Throws an error if one of the internal vectors is not sorted
    /// in strictly asecending order, by index.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// // This line should not panic.
    /// let matrix = VecOfVec::new( vec![ vec![ (0,5), (1,6) ] ] ); // no panic
    /// 
    /// // This line should panic.
    /// // let matrix = VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ] );
    /// // And here is the test that confirms it.
    /// let result = std::panic::catch_unwind(|| {VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ] )} );
    /// assert!(result.is_err()); 
    /// ```
    pub fn new( vecvec: Vec < Vec < ( ColumnIndex, Coefficient ) > > ) -> Self  
        where   ColumnIndex: PartialOrd,
                ( ColumnIndex, Coefficient ):     KeyValGet < Key = ColumnIndex, Val = Coefficient >, // if we comment this out then  we get an error sayint that `OrderOperatorByKey` doesn't implement the `JudgePartialOrder` trait; probably this has something to do with the fact that `OrderOperatorByKey` autoimplements `OrderOperatorByKey` for structs that implement `KeyValGet`
    {
        for vec in vecvec.iter() {
            if ! is_sorted_strictly( vec, &mut OrderOperatorByKey::new() ) {
                panic!("Attempt to construct `VecOfVec` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
            }
        }        

        VecOfVec{   
                    vec_of_vec:         vecvec,   
                    order_operator:   OrderOperatorByKey::new(),                 
                    // phantom_lifetime:   PhantomData,                  
                }
    }

    /// Returns a clone of the internally stored order comparator.
    pub fn clone_order_operator( &self ) -> OrderOperatorByKey 
        where OrderOperatorByKey :   Clone
    { self.order_operator.clone() }

    /// Returns an immutable reference to the `Vec< Vec< (usize, Coefficient) > >` that stores the entries of of the matrix, internally.
    pub fn vec_of_vec_borrowed( &self ) -> & Vec< Vec< (ColumnIndex, Coefficient) > > { & self.vec_of_vec }


    /// Returns the internal `Vec< Vec< (usize, Coefficient) > >` (consuming `self`)
    pub fn vec_of_vec_unwrap( self ) -> Vec< Vec< (ColumnIndex, Coefficient) > > { self.vec_of_vec }    


    /// Append a vector to the internally stored `Vec<Vec<EntryType>>` struct.  The new vector must be sorted in strictly ascending order; the function panics, otherwise.
    pub fn push_vec( &mut self, new_sorted_vec: Vec < (ColumnIndex, Coefficient) > ) 
        where 
            OrderOperatorByKey: JudgePartialOrder< (ColumnIndex, Coefficient)>
    {
        if ! is_sorted_strictly( & new_sorted_vec, &mut self.order_operator ) {
            panic!("Attempt to append a non-strictly-sorted vector to `VecOfVec`.  To proceed, the new vector must be sorted in *strictly* ascending order, by index.") 
        }        
        self.vec_of_vec.push( new_sorted_vec );
    }    

    /// Ensures that `self` has at least `m` rows.
    /// 
    /// If `self` has fewer than `m` major slices, then
    /// - (if necessary) extend the capacity of `self.vec_of_vec` to `m` using `reserve_exact`
    /// - push as many empty vectors to `self.vec_of_vec` as necessary, to ensure there are `m` rows
    pub fn ensure_number_of_rows( &mut self, m: usize ) {
        if m > self.vec_of_vec.len() {
            if m > self.vec_of_vec.capacity() {
                self.vec_of_vec.reserve_exact( m-self.vec_of_vec.capacity() ) // reserve the exact right amount of extra capacity, if needed
            }
        }
        for _ in 0 .. (m-self.vec_of_vec.len()) {
            self.vec_of_vec.push( Vec::with_capacity(0) );
        }

    }

    /// Creates a [`VecOfVec`] from an iterable that runs over iterables that run over tuples.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted_ref::VecOfVec;
    ///     
    /// let iter = (0..2).map( |x| vec![(x,x)] );
    /// let vec_of_vec = VecOfVec::from_iterable_of_iterables( iter );
    /// let ground_truth =  VecOfVec::new(
    ///                         vec![ 
    ///                             vec![ (0,0)        ],
    ///                             vec![        (1,1) ],
    ///                         ]
    ///                     );
    /// assert_eq!( vec_of_vec, ground_truth  );
    /// ```
    pub fn from_iterable_of_iterables< I >( iter: I ) -> VecOfVec< ColumnIndex, Coefficient > 
        where
            I:          IntoIterator,
            I::Item:    IntoIterator< Item = (ColumnIndex, Coefficient) >,
            ColumnIndex:     Clone + PartialOrd,
            Coefficient:     Clone,
    {
        let vecvec =    iter.into_iter().map( |x| x.into_iter().collect_vec() ).collect_vec();
        VecOfVec::new( vecvec )
    }

    pub fn vec_of_vec( &self ) -> &Vec< Vec< ( ColumnIndex, Coefficient ) > > { &self.vec_of_vec }


}



impl < Coefficient >

    VecOfVec< usize, Coefficient > 
{

    /// Maximum column index of an explicitly stored entry (or None, if no entries are stored)
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// let a: Vec<Vec<(usize,f64)>>    =   vec![vec![ (0,1.) ]];   // matrix with 1 entry
    /// let b: Vec<Vec<(usize,f64)>>    =   Vec::new();             // empty matrix
    /// 
    /// assert_eq!( VecOfVec::new(a).ok().unwrap().max_column_index(), Some(0) );
    /// assert_eq!( VecOfVec::new(b).ok().unwrap().max_column_index(), None    );
    /// ```
    pub fn max_column_index( &self ) -> Option< usize > {
        let mut max_col = None;
        for row in self.vec_of_vec.iter() {
            for col_index in row.iter().map(|x| x.0) {
                max_col = max_col.max( Some(col_index) )
            }
        }     
        max_col   
    }


    /// Prints a dense form of the matrix
    /// 
    /// - Entries that are not explicitly stored will appear as `filler_coefficient`.
    /// - The number of columns will equal the maximum column index of any explicitly stored entry.
    /// - If no entry is explicitly stored, then this function prints `Matrix contains no explicitly stored entries`
    pub fn print_dense( &self, filler_coefficient: Coefficient )
            where
                Coefficient:    std::fmt::Debug + Clone,
    {
        match self.max_column_index() {
            None => { println!("Matrix contains no explicitly stored entries") }
            Some( max_column_index ) => {
                let mut to_print = vec![ filler_coefficient.clone(); max_column_index + 1 ]; // we need +1 becuase of zero indexing
                for row_data in self.vec_of_vec.iter() {
                    for entry in to_print.iter_mut() { *entry = filler_coefficient.clone() }; // clear the entries
                    for (j,v) in row_data.iter().cloned() { to_print[j]=v } // fill with the nonzero entries
                    println!("{:?}", &to_print);
                }
            }
        }
    }

    /// Generates a new `VecOrVecSimple` representing the transpose, with copied (not borrowed) data.
    /// 
    /// Compare this method with [transpose](crate::algebra::matrices::operations::MatrixOracleOperations::transpose), which
    /// generates a *lazy* transpose.
    /// 
    /// # Arguments
    /// 
    /// If `num_rows_of_transpose` is less than 1 + the maximum column index of a
    /// any explicitly stored entry in `self`, then returns `None`.  You can check the maximum column index with the method
    /// [VecOfVec::max_column_index]. Otherwise returns the transpose of a
    /// matrix of `self`, where `self` is regarded as a matrix of size `(self.vec_of_vec.len()) x num_rows_of_transpose`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::operations::MatrixOracleOperations;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted_ref::VecOfVec;
    /// use oat_rust::algebra::matrices::query::MatrixOracle;
    /// 
    /// use itertools::Itertools;
    ///  
    /// let matrix     
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"0,0"), (1,"0,1"), ],
    ///                 vec![ (0,"1,0"), (1,"1,1"), ],
    ///             ]
    ///         );
    /// 
    /// let transpose   
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"0,0"), (1,"1,0"), ],
    ///                 vec![ (0,"0,1"), (1,"1,1"), ],
    ///             ]
    ///         );
    /// 
    /// let transpose_lazy  =   matrix.transpose();
    /// let transpose_deep  =   & matrix.transpose_deep(2).unwrap();
    ///
    ///  
    /// for p in 0..2 {
    ///     let a   =   transpose_deep.row(&p).cloned().collect_vec();
    ///     let b   =   transpose_lazy.row(&p).collect_vec();
    ///     let c   =   transpose.row(&p).cloned().collect_vec();
    ///     assert_eq!(a, b);
    ///     assert_eq!(a, c);
    /// }
    /// 
    /// // check that matrix returns None when too few rows are specified
    /// assert_eq!( matrix.antitranspose_deep(1), None );
    /// ```
    pub fn transpose_deep( &self, num_rows_of_transpose: usize ) -> Option< Self >
    
        where
            Coefficient:  Clone
    {
        // count the explicit entries in each column
        let mut hist = vec![0; num_rows_of_transpose]; // a histogram that counts the number of nonzero entries in each column
        for row in self.vec_of_vec.iter() {
            for col_index in row.iter().map(|x| x.0) {
                if col_index + 1 > num_rows_of_transpose {
                    return None // in this case we haven't reserved enough space
                } else {
                    hist[ col_index ] +=1; // otherwise count the column in the historgram
                }
            }
        }

        // allocate exactly the right amount of space in the new matrix
        let mut transpose_data = Vec::with_capacity(num_rows_of_transpose);
        for num_indices in hist.into_iter() {
            transpose_data.push( Vec::with_capacity(num_indices) ); 
        }

        // fill the new matrix
        for row in 0 .. self.vec_of_vec.len() {
            for (col,val) in self.vec_of_vec[row].iter() {
                transpose_data[*col].push( (row,val.clone()) )
            }
        }

        return Some( VecOfVec { vec_of_vec: transpose_data, order_operator: OrderOperatorByKey::new() } )
    }


    /// Generates a new [VecOrVec] representing the anti-transpose, with copied (not borrowed) data.
    /// 
    /// # Arguments
    /// 
    /// If `num_rows_of_antitranspose` is less than 1 + the maximum column index of a
    /// any explicitly stored entry in `self`, then `None` is returned.
    /// You can check the maximum column index of `self` with the method [VecOfVec::max_column_index].
    /// Otherwise returns the antitranspose of a
    /// matrix of `self`, where `self` is regarded as a matrix of size `(self.vec_of_vec.len()) x num_rows_of_antitranspose`.
    /// 
    /// # Caution
    /// 
    /// This method differs from [order_antitranspose](crate::algebra::matrices::operations::MatrixOracleOperations::order_antitranspose) in several ways:
    /// - it makes a copy of the underlying data
    /// - it changes the indices of the nonzero entires; by contrast, the order antitranspose only changes the order in which entries are returned by iterators
    ///   without changing the underlying indices.
    ///   For details see  the documentation for [OrderAntiTranspose](crate::algebra::matrices::types::transpose::OrderAntiTranspose). 
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::operations::MatrixOracleOperations;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted_ref::VecOfVec;
    /// use oat_rust::algebra::matrices::query::MatrixOracle;
    /// 
    /// use itertools::Itertools;
    /// 
    /// let matrix     
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"0,0"), (1,"0,1"), (2,"0,2"), ],
    ///                 vec![ (0,"1,0"), (1,"1,1"), (2,"1,2"), ],
    ///             ]
    ///         );
    ///                 
    /// let antitranspose   
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"1,2"), (1,"0,2"), ],
    ///                 vec![ (0,"1,1"), (1,"0,1"), ],
    ///                 vec![ (0,"1,0"), (1,"0,0"), ],                        
    ///             ]
    ///         );
    ///                     
    /// let antitranspose_deep  =  matrix.antitranspose_deep(3).unwrap();
    /// let antitranspose_deep  =  & antitranspose_deep;
    ///                     
    /// // check that calculation is correct
    /// for p in 0..3 {
    ///     let a   =   antitranspose.row(&p).collect_vec();
    ///     let b   =   antitranspose_deep.row(&p).collect_vec();
    ///     assert_eq!(a, b)
    /// }
    /// 
    /// // check that matrix returns None when too few rows are specified
    /// assert_eq!( matrix.antitranspose_deep(2), None );
    /// ```
    pub fn antitranspose_deep( &self, num_rows_of_antitranspose: usize ) -> Option< Self >
    
        where
            Coefficient:  Clone + std::fmt::Debug,
    {
        let num_rows_original = self.vec_of_vec.len();
        
        let transpose   =   self.transpose_deep(num_rows_of_antitranspose)?;
        let mut transpose_vecvec    =   transpose.vec_of_vec_unwrap();
        ( &mut transpose_vecvec ).reverse();
        let num_rows_in_transpose = transpose_vecvec.len();
        for p in 0 .. num_rows_in_transpose {
            transpose_vecvec[ p ].reverse();
        }        
        for row in transpose_vecvec.iter_mut() {  
            row.iter_mut()
                .for_each(
                    |(j,_v)|
                    *j = num_rows_original - *j -1
                );      
            println!("row after subtracting: {:?}", row);                 
        }
        for p in 0..transpose_vecvec.len() {
            println!("raw row before: {:?}", &transpose_vecvec[p] );
        }                
        return Some( VecOfVec::new( transpose_vecvec ) )
    }    

    /// Returns a [MatrixBiajorData](crate::algebra::matrices::types::bimajor::MatrixBimajorData) that contains both row-major and colum-major copies of `self`
    /// 
    /// This doubles memory usage (it is not lazy); however it allows efficient access to **both rows and columns**.
    /// 
    /// # Arguments
    /// 
    /// If `num_rows_of_transpose` is less than 1 + the maximum column index of a
    /// any explicitly stored entry in `self`, then `None` is returned.  Otherwise returns the antitranspose of a
    /// matrix of `self`, where `self` is regarded as a matrix of size `(self.vec_of_vec.len()) x num_rows_of_transpose`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::matrices::query::MatrixOracle;      
    /// 
    /// use itertools::Itertools;
    /// 
    /// let matrix   
    ///     =   VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"a"), (1,"b"), (2,"c"), ],
    ///                 vec![ (0,"d"), (1,"e"), (2,"f"), ],
    ///             ]
    ///         ).ok().unwrap();
    ///     
    /// let bimajor =   matrix.clone().bimajor(3).unwrap();
    ///     
    /// for p in 0..2 {
    ///     assert_eq!(     ( & matrix  ).row(&p).collect_vec(),
    ///                     ( & bimajor ).row(&p).collect_vec(),            );
    ///     assert_eq!(     ( & matrix  ).row(&p).collect_vec(),
    ///                     ( & bimajor ).row(&p).collect_vec(),     );
    ///     assert_eq!(     ( & matrix  ).row_reverse(&p).collect_vec(),
    ///                     ( & bimajor ).row_reverse(&p).collect_vec(),    );
    /// }
    /// 
    /// for p in 0..3 {
    ///     assert_eq!(     ( & matrix  ).column(&p).collect_vec(),
    ///                     ( & bimajor ).column(&p).collect_vec(),            );
    ///     assert_eq!(     ( & matrix  ).column(&p).collect_vec(),
    ///                     ( & bimajor ).column(&p).collect_vec(),     );
    ///     assert_eq!(     ( & matrix  ).column_reverse(&p).collect_vec(),
    ///                     ( & bimajor ).column_reverse(&p).collect_vec(),    );
    /// }     
    /// ```
    pub fn bimajor( self, num_rows_of_transpose: usize ) -> Option< MatrixBimajorData< Self, Self > >

        where
            Coefficient:  Clone + std::fmt::Debug,
    {        
        let transpose = self.transpose_deep(num_rows_of_transpose);
        match transpose {
            None => return None,
            Some(transpose) => return Some( MatrixBimajorData{ matrix_rows_data: self, matrix_columns_data: transpose } )
        }     
    }        




}


impl VecOfVec< usize, usize > {

    //  RANDOM MATRICES
    //  ===========================================================================    

    //! Some functions that return matrices, e.g. `f(m) = random_m_by_m_sparse_matrix`.



    /// Random invertible upper-triangular row-major matrix with coefficients mod p
    /// 
    /// - The returned matrix has size `matrix_size` x `matrix_size`
    /// - Diagonal entries are drawn independently from the uniform distribution on `{ 1, .., p-1 }`.
    /// - Entries strictly above the diagonal are equally likely to be (i) structurally zero (and therefore algebraically) zero, 
    /// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
    /// In case (iii), entries are drawn independently from the uniform distribution on `{ 1, .., p-1 }`.
    pub fn random_mod_p_upper_triangular_invertible( 
            matrix_size: usize, 
            modulus: usize 
        ) ->
            VecOfVec< usize, usize >
    {
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

        VecOfVec::new(vec_of_vec)// formally wrap the matrix in a VecOfVec struct
    }

    /// Generate a random invertible upper triangular matrix with integer coefficients.
    /// 
    /// - Entries on the diagonal are equal to 1.
    /// - Each entry strictly above the diagonal is either (i) strucutrually zero, 
    /// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
    /// Algebraically nonzero entries are drawn iid from the uniform distribution on [1, .., `modulus`).
    /// 
    /// Nonzero entries are drawn in an iid fasion from [1, .., `modulus`).
    pub fn random_mod_p_upper_unitriangular( 
            matrix_size: usize, 
            modulus: usize
        ) -> 
            VecOfVec< usize, usize > {

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


    /// Generate a random row-major matrix of size `num_rows` x `num_cols` with coefficients
    /// in `{ 0, .., p-1 }`.
    /// 
    /// Each entry is equally likely to be (i) structurally zero, 
    /// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
    /// In case (iii), entries are drawn independently from the uniform distribution on `{ 1, .., p-1 }`.
    pub fn random_mod_p( 
            num_rows: usize, 
            num_cols: usize, 
            modulus: usize 
        ) ->
            VecOfVec< usize, usize >
    {
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

        VecOfVec::new(vec_of_vec)// formally wrap the matrix in a VecOfVec struct
    }

    /// Generate a random `VecOfVec< usize, usize >`
    /// 
    /// Nonzero entries are drawn in an iid fasion from [0, .., `modulus`) if `allow_nonstructural_zero == true` and from [1, .., `modulus`)
    /// if `allow_nonstructural_zero == false`.
    pub fn random_mod_p_with_density( 
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

        VecOfVec::new( vecvec )

    }    
}


//  ORACLE IMPLEMENTATIONS FOR &'a VecOfVec
//  ---------------------------------------------

impl < 'a, ColumnIndex, Coefficient > 

    MatrixOracle for 
    &'a VecOfVec < ColumnIndex, Coefficient >

    where   ColumnIndex:        Clone + Debug + Ord + Eq,
            Coefficient:        Clone + Debug + Ord + Eq,

{ 
    type Coefficient        =   Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   usize;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   ColumnIndex;       // The type of column indices
    
    type RowEntry           =   &'a ( ColumnIndex, Coefficient );          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   ( usize, Coefficient );       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                =   Iter< 'a, (ColumnIndex, Coefficient) >;  // What you get when you ask for a row.
    type RowReverse         =   Rev<std::slice::Iter<'a, (ColumnIndex, Coefficient)>>;  // What you get when you ask for a row with the order of entries reversed
    
    type Column             =   VecOfVecMatrixColumn< 'a, ColumnIndex, Coefficient >; // What you get when you ask for a column
    type ColumnReverse      =   VecOfVecMatrixColumnReverse< 'a, ColumnIndex, Coefficient >; // What you get when you ask for a column with the order of entries reversed 

    // entry lookup
    fn structural_nonzero_entry( & self, row: &Self::RowIndex, column: &Self::ColumnIndex )   ->  Option< Self::Coefficient > {
        let row = & self.vec_of_vec[ * row ];
        find_sorted_binary_oracle( 
                    0, 
                    row.len() as isize - 1, 
                    |p| column.cmp( & row[p as usize].0 ) 
                ).map(|x| row[ x as usize ].1.clone() )        
    }  

    // row lookup
    fn row(                     &self,  index: &Self::RowIndex    )   -> Self::Row   { 
        self.vec_of_vec[*index].iter()
    }
    // fn row_result(                 &   self, index: & Self::RowIndex   )   -> Result< Self::Row, Self::RowIndex >  {
    //     if *index >= self.vec_of_vec.len() { return None }
    //     else { return Some( self.row( index ) ) }
    // }
    fn row_reverse(             &self,  index: &Self::RowIndex    )       -> Self::RowReverse  { 
        self.vec_of_vec[*index].iter().rev()
    }
    // fn row_reverse_result(         &   self, index: & Self::RowIndex   )   -> Result< Self::RowReverse, Self::RowIndex >  {
    //     if *index >= self.vec_of_vec.len() { return None }
    //     else { return Some( self.row_reverse( index ) ) }        
    // }  
    
    // column lookup
    fn column(                  &self,  index: &Self::ColumnIndex )       -> Self::Column {
        VecOfVecMatrixColumn{
            vec_of_vec:             self,
            row_index:                 0,
            column_index:                 index.clone(),
            phantom_column_index:         PhantomData,            
        }
    }
    // fn column_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex > {
    //     Some( self.column(index) )
    // }   
    fn column_reverse(          &self,  index: &Self::ColumnIndex )       -> Self::ColumnReverse  { 
        VecOfVecMatrixColumnReverse{
            vec_of_vec:             self,
            row_index:                 self.vec_of_vec.len(),
            column_index:                 index.clone(),
            phantom_column_index:         PhantomData,            
        }
    }            
    // fn column_reverse_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex >{
    //     Some( self.column_reverse(index) )
    // }
    
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool {
        *index < self.vec_of_vec.len()
    }
    
    /// this data structure does not specify a maximum column index, so every column index is allowed.
    fn has_column_for_index(  &   self, _index: & Self::ColumnIndex)   -> bool {
        true
    }
}   



/// Represents a column of a `VecOfVec`, with entries appearing in descending order of index.
/// 
/// This has a column for every column index of type `ColumnIndex` (all but finitely many of these columns are zero).
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::MatrixOracle;
///         
/// // define a row-major sparse matrix representing the array
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] ).ok().unwrap();
/// itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).column_reverse( &0 ) );
/// itertools::assert_equal( vec![ (0,5) ], (& matrix).column_reverse( &1 ) );
/// itertools::assert_equal( vec![ (1,7), (0,5) ], (& matrix).column_reverse( &2 ) );  
/// ```
pub struct VecOfVecMatrixColumnReverse< 'a, ColumnIndex, Coefficient > 
    where   ColumnIndex:     Clone,    
            Coefficient:     Clone,  
{
    vec_of_vec:             &'a VecOfVec< ColumnIndex, Coefficient >,
    row_index:                 usize,
    column_index:                 ColumnIndex,
    phantom_column_index:         PhantomData< ColumnIndex >,
}

// Iterator
impl < 'a, ColumnIndex, Coefficient > 

        Iterator for 

        VecOfVecMatrixColumnReverse< 'a, ColumnIndex, Coefficient >

    where   ColumnIndex:     Clone + PartialEq,    
            Coefficient:     Clone,          
{
    type Item = (usize, Coefficient);

    fn next( &mut self ) -> Option< Self::Item > {
        // extract the underlying vector of vectors
        let vecvec_data =  self.vec_of_vec.vec_of_vec_borrowed();                
        // (assuming the matrix is row-major) scan the rows of the matrix from bottom to top
        while self.row_index > 0 {
            // drop the row number by 1
            self.row_index -= 1;            
            // get the row
            let row = & vecvec_data[ self.row_index ];
            // scan the row to see if it contains an entry of form ( my_column_index, snzval )
            for ( column_index, snzval ) in row {
                // if it does, then return ( row_number, snzval )
                if column_index != & self.column_index { continue }
                else { return Some( ( self.row_index, snzval.clone() ) ) }
            }
        }
        None
    }
}        

/// Represents a column of a `VecOfVec`, with entries appearing in ascending order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::MatrixOracle;
///         
/// // define a row-major sparse matrix representing the array
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] ).ok().unwrap();
/// itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).column( &0 ) );
/// itertools::assert_equal( vec![ (0,5) ], (& matrix).column( &1 ) );
/// itertools::assert_equal( vec![ (0,5), (1,7) ], (& matrix).column( &2 ) );  
/// ```
pub struct VecOfVecMatrixColumn< 'a, ColumnIndex, Coefficient > 
    where   ColumnIndex:     Clone,    
            Coefficient:     Clone,  
{
    vec_of_vec:             &'a VecOfVec< ColumnIndex, Coefficient >,
    row_index:                 usize,
    column_index:                 ColumnIndex,
    phantom_column_index:         PhantomData< ColumnIndex >,
}

// Iterator
impl < 'a, ColumnIndex, Coefficient > 

        Iterator for 

        VecOfVecMatrixColumn< 'a, ColumnIndex, Coefficient >

    where   ColumnIndex:     Clone + PartialEq,    
            Coefficient:     Clone,          
{
    type Item = (usize, Coefficient);

    fn next( &mut self ) -> Option< Self::Item > {
        // extract the underlying vector of vectors
        let vecvec_data =  self.vec_of_vec.vec_of_vec_borrowed();                
        // (assuming the matrix is row-major) scan the rows of the matrix from bottom to top
        while self.row_index < vecvec_data.len() {
            // get the row
            let row = & vecvec_data[ self.row_index ];
            // scan the row to see if it contains an entry of form ( my_column_index, snzval )
            for ( column_index, snzval ) in row {
                // if it does, then return ( row_number, snzval )
                if column_index != & self.column_index { continue }
                else { 
                    // grow the row number by 1 (this is one of two branches where we do this)
                    self.row_index += 1;
                    
                    // return the entry
                    return Some( ( self.row_index - 1, snzval.clone() ) ) 
                }
            }
            // grow the row number by 1 (this is one of two branches where we do this)
            self.row_index += 1;              
        }      

        // in this case the iterator is exhausted
        None
    }
}            








#[cfg(test)]
mod tests {    

    

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_vec_of_vec_from_iterable_of_iterables() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;
        
        let iter = (0..2).map( |x| vec![(x,x)] );
        let vec_of_vec = VecOfVec::from_iterable_of_iterables( iter ).ok().unwrap();
        itertools::assert_equal( (& vec_of_vec).row(&0), vec![(0,0)]);
        itertools::assert_equal( (& vec_of_vec).row(&1), vec![(1,1)])   
    }

    #[test]    
    fn test_vec_of_vec_simple_descending_column() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;
        
        // define a row-major sparse matrix representing the array
        //  0   5   5
        //  0   0   7
        let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] ).ok().unwrap();
        itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).column_reverse( & 0 ) );
        itertools::assert_equal( vec![ (0,5) ], (& matrix).column_reverse( & 1 ) );
        itertools::assert_equal( vec![ (1,7), (0,5) ], (& matrix).column_reverse( & 2 ) );        
    }

    #[test]     
    fn test_vec_of_vec_simple_ascending_column() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;

        // define a row-major sparse matrix representing the array
        //  0   5   5
        //  0   0   7
        let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] ).ok().unwrap();
        println!("waytpoint 1");
        itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).column( & 0 ) );
        println!("waytpoint 2");        
        itertools::assert_equal( vec![ (0,5) ], (& matrix).column( & 1 ) );
        println!("waytpoint 3");        
        itertools::assert_equal( vec![ (0,5), (1,7) ], (& matrix).column( & 2 ) ); 
    }

    #[test]    
    fn test_vec_of_vec_simple_antitranspose_deep() {
        
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;
        
        use itertools::Itertools;
        
        let matrix     
            =   & VecOfVec::new(
                    vec![
                        vec![ (0,"0,0"), (1,"0,1"), (2,"0,2"), ],
                        vec![ (0,"1,0"), (1,"1,1"), (2,"1,2"), ],
                    ]
                ).ok().unwrap();

        let antitranspose   
            =   & VecOfVec::new(
                    vec![
                        vec![ (0,"1,2"), (1,"0,2"), ],
                        vec![ (0,"1,1"), (1,"0,1"), ],
                        vec![ (0,"1,0"), (1,"0,0"), ],                        
                    ]
                ).ok().unwrap();

        let antitranspose_deep  =  matrix.antitranspose_deep(3).unwrap();
        let antitranspose_deep  =  & antitranspose_deep;

        // check that calculation is correct
        for p in 0..3 {
            let a   =   antitranspose.row(&p).collect_vec();
            let b   =   antitranspose_deep.row(& p ).collect_vec();
            assert_eq!(a, b)
        }

        // check that matrix returns None when too few rows are specified
        assert_eq!( matrix.antitranspose_deep(2), None );
    }



    #[test]
    fn test_vec_of_vec_simple_bimajor_data() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;
        
        use itertools::Itertools;
        
        let matrix   
            =   VecOfVec::new(
                    vec![
                        vec![ (0,"a"), (1,"b"), (2,"c"), ],
                        vec![ (0,"d"), (1,"e"), (2,"f"), ],
                    ]
                ).ok().unwrap();
            
        let bimajor =   matrix.clone().bimajor(3).unwrap();
            
        for p in 0..2 {
            assert_eq!(     ( & matrix  ).row( & p ).collect_vec(),
                            ( & bimajor ).row( & p ).collect_vec(),            );
            assert_eq!(     ( & matrix  ).row_reverse( & p ).collect_vec(),
                            ( & bimajor ).row_reverse( & p ).collect_vec(),    );
        }
        for p in 0..3 {
            assert_eq!(     ( & matrix  ).column( & p ).collect_vec(),
                            ( & bimajor ).column( & p ).collect_vec(),     );
            assert_eq!(     ( & matrix  ).column_reverse( & p ).collect_vec(),
                            ( & bimajor ).column_reverse( & p ).collect_vec(),    );
        }        
    }

    #[test]
    fn test_misc() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        
        // This line should not panic.
        let matrix = VecOfVec::new( vec![ vec![ (0,5), (1,6) ] ] ); // no panic
        
        // This line should panic.
        // let matrix = VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ] );
        // And here is the test that confirms it.
        let result = std::panic::catch_unwind(|| {VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ] )} );
        assert!(result.is_err()); 
    }

}
