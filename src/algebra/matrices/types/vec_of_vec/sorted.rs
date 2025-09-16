//! Variant of vec-of-vec where internal vectors are stored in sorted order
//! 
//! See the [parent module](super) for details on what vector-of-vectors format means.

use crate::algebra::matrices::operations::umatch::row_major::Umatch;
use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::matrices::types::packet::MatrixAlgebraPacket;
use crate::algebra::rings::traits::{DivisionRingOperations, SemiringOperations};
use crate::algebra::vectors::entries::{KeyValGet};
use crate::algebra::matrices::query::MatrixOracle;

use crate::utilities::binary_search::{find_sorted_binary_oracle};
use crate::utilities::functions::evaluate::EvaluateFunction;
use crate::utilities::order::{is_sorted_strictly, JudgePartialOrder, OrderOperatorAuto, OrderOperatorByKey };


use derive_new::new;
use ordered_float::OrderedFloat;
use rand::Rng;                                          // we use this module to generate random elements
use rand::distributions::{Bernoulli, Distribution};
use serde::{Deserialize, Serialize};     // we use this module to generate random elements

use std::collections::HashMap;
use std::iter::{Rev, Cloned};
use std::fmt::Debug;
use std::marker::PhantomData;
use std::mem;
use std::slice::Iter;
use itertools::Itertools;

use super::super::bimajor::MatrixBimajorData;





/// Vector-of-vectors sparse matrix format, with internal vectors stored in sorted order
/// 
/// See the [parent module](super) for details on what vector-of-vectors format means.
/// 
/// - Entries in each internal vector are stored in *strictly* ascending order of index (no repeat indices).
/// - Order of incides is determined by Rust's `PartialOrd` trait.
/// 
/// # Details
/// 
/// This data structure will report that it has a column for every column index of type `ColumnIndex`.
/// It will report that it has a row for every nonzero index `i < N` where `N` is the length of the outer
/// vector in the vector-of-vectors.
/// 
/// # See also
/// 
/// If you need a customized order for entries, see the [variant for custom orders](super::sorted_custom).
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// 
/// // Create a new matrix with three rows and at most one entry per row
/// let matrix  =   VecOfVec::new(
///                     vec![ 
///                         vec![(1,2.0)],  // row 0 has an entry with value 2.0 in column 1
///                         vec![],         // row 1 has no nonzero entries
///                         vec![(3,5.0)]   // row 2 has an entry with value 5.0 in column 3
///                     ],
///                 ).ok().unwrap();
/// ```
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
// #[serde(bound = "ColumnIndex: Serialize, Coefficient: Serialize")]
pub struct VecOfVec
            < ColumnIndex, Coefficient >

{
    vec_of_vec:         Vec< Vec< ( ColumnIndex, Coefficient ) > >,
    order_operator:     OrderOperatorByKey,
}

impl < ColumnIndex, Coefficient >
        
        VecOfVec
            < ColumnIndex, Coefficient > 

{

    /// Number of rows of the matrix
    /// 
    /// Calculated as the number of vectors contained in the vector-of-vectors.
    pub fn number_of_rows( &self ) -> usize {
        self.vec_of_vec.len()
    }

    /// Returns the number of structural nonzero entries in the matrix.
    pub fn number_of_structural_nonzeros( &self ) -> usize {
        self.vec_of_vec.iter().map( |v| v.len() ).sum()
    }

    /// Returns the number of structural nonzero entries in row `index`.
    /// 
    /// Returns `None` if `index` is out of bounds.
    pub fn number_of_structural_nonzeros_in_row( &self, index: usize ) -> Option< usize > {
        if index + 1 > self.number_of_rows() {
            None
        } else {
            Some( self.vec_of_vec[ index ].len() )
        }
    }    

    /// Make a new [VecOfVec]  
    /// 
    /// Returns `Ok(VecOfVec)` if the input is valid, meaning that each vector is sorted in increasing order, by index.
    /// Otherwise returns `Err(input_data)`.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// // This produces an Ok(VecOfVec)
    /// let data = vec![ vec![ (0,5), (1,6) ] ];
    /// let matrix = VecOfVec::new( data ); 
    /// 
    /// // This line returns an Err(input_data)
    /// let data = vec![ vec![ (1,6), (0,5) ] ];
    /// let err = VecOfVec::new( data.clone() );
    /// assert_eq!( err, Err( data.clone() ) ); 
    /// ```
    pub fn new( 
                vecvec: Vec < Vec < ( ColumnIndex, Coefficient ) > > 
        ) -> 
        Result< 
            Self, 
            Vec < Vec < ( ColumnIndex, Coefficient ) > > 
        >

        where   ColumnIndex: Ord + Debug,
                Coefficient: Debug,
                ( ColumnIndex, Coefficient ):     KeyValGet < Key = ColumnIndex, Val = Coefficient >, // if we comment this out then  we get an error sayint that `OrderOperatorByKey` doesn't implement the `JudgePartialOrder` trait; probably this has something to do with the fact that `OrderOperatorByKey` autoimplements `OrderOperatorByKey` for structs that implement `KeyValGet`
    {
        for (row_number, vec) in vecvec.iter().enumerate() {
            if ! is_sorted_strictly( vec, &mut OrderOperatorByKey::new() ) {
                println!("Attempt to construct `VecOfVec` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.  ");
                print!("The first unsorted vector occurs in row {:?}; it is:\n{:?}\n", row_number, &vec);
                return Err(vecvec)
            }
        }        

        Ok( VecOfVec{   
                    vec_of_vec:         vecvec,   
                    order_operator:   OrderOperatorByKey::new(),                 
                    // phantom_lifetime:   PhantomData,                  
                } )
    }

    /// Returns a clone of the internally stored order comparator.
    pub fn order_operator( &self ) -> OrderOperatorByKey 
        where OrderOperatorByKey :   Clone
    { self.order_operator.clone() }

    /// Returns an immutable reference to the `Vec< Vec< (ColumnIndex, Coefficient) > >` that stores the entries of of the matrix, internally.
    pub fn inner_vec_of_vec_ref( &self ) -> & Vec< Vec< (ColumnIndex, Coefficient) > > { & self.vec_of_vec }

    /// Returns the `Vec< Vec< (usize, Coefficient) > >` that stores the entries of of the matrix, internally.
    /// 
    /// This method consumes `self`.
    pub fn inner_vec_of_vec( self ) -> Vec< Vec< (ColumnIndex, Coefficient) > > { self.vec_of_vec } 

    /// Get a reference to the `i`th vector in the internally stored vector of vectors
    /// 
    /// This operation is different from `self.row(index)`, which would return an interator instead of a reference to a vector.
    pub fn row_ref( &self, index: usize )  -> & Vec< (ColumnIndex, Coefficient) > {
        & self.vec_of_vec[ index ]
    }


    /// Returns an iterator that runs over all `(row, column, coefficient)` triplets of the nonzero entries of the matrix.
    /// 
    /// Entries are sorted by row index, then by column index.
    /// 
    /// # Example
    /// 
    ///  ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    ///
    /// let matrix = VecOfVec::new(
    ///     vec![
    ///         vec![ (0, 1.0), (1, 2.0) ], // row 0 has entries in columns 0 and 2
    ///         vec![           (1, 3.0) ],           // row 1 has an entry in column 1
    ///     ]
    /// ).ok().unwrap();
    /// 
    /// let triplets =  matrix.triplets();
    /// assert!( triplets.eq(
    ///     vec![ (0, 0, 1.0 ),
    ///           (0, 1, 2.0 ),
    ///           (1, 1, 3.0 )
    ///     ])
    /// );
    /// ```
    pub fn triplets( &self ) -> VecOfVecTriplets< '_, ColumnIndex, Coefficient >
        where ColumnIndex: Clone,
              Coefficient: Clone,
    {
        VecOfVecTriplets::new( self )
    }


    /// Returns `M[p]`, where `p` is any sequence of indices.
    /// 
    /// Concretely, the result is a [VecOfVec] matrix `N` such that `N[i] = M[p[i]]` for all `i` in `0 .. p.len()`.
    /// 
    /// This operation is out of place -- it does not change `self`.
    pub fn permute_rows_out_of_place< Indices >( &self, indices: Indices ) -> Self 
        where 
            Indices:                IntoIterator< Item = usize >,
            ColumnIndex:            Clone + Debug + Ord,
            Coefficient:            Clone + Debug,
    {
        let row_iterator    =
        indices
            .into_iter()
            .map( |i| self.row_ref(i).clone() );

        VecOfVec::from_iterable_of_iterables( row_iterator ).ok().unwrap()
    }

    /// Assigns a new column index to each nonzero entry
    /// 
    /// We replace each entry `(i,a)` with an entry `(j,a)`, where `j` is the index provided by `column_index_map`.
    /// Concretely, we have `j = column_index_map.evaluate_function(i)`.
    /// 
    /// After re-assigning each column index, each row is sorted.
    /// 
    /// # Error handling
    /// 
    /// It is possible that `column_index_map` assigns the same value `j` to two different column indices, `p` and `q`.
    /// In that case, the entries in one of the new rows might be sorted, but not **strictly** sorted.
    /// Strict sorting is important for a variety of matrix operations to function correctly, and it's a requirement of
    /// the [VecOfVec] data structure. Therefore, if the strict sorting requirement is violated, this function will
    /// return an `Err(Vec<Vec< (NewColumnIndex, Coefficient) >>)`, which contains the vector of reindexed vectors.
    /// 
    /// # Examples
    /// 
    /// Here is an example of a successful re-indexing:
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    ///  
    /// let matrix              =    VecOfVec::new(
    ///                                  vec![
    ///                                      vec![ (0,"a"), (1,"b"), ],
    ///                                      vec![ (0,"c"), (1,"d"), ],
    ///                                  ]
    ///                              ).ok().unwrap();
    /// 
    /// // Vectors implement the EvaluateFunction trait. This vector will map 0 to 3 and 1 to 2.
    /// let column_index_map    =   vec![ 3, 2 ]; 
    /// let matrix_reindexed    =   matrix.reassign_column_indices_out_of_place( & column_index_map );
    /// 
    /// // This is what the permuted matrix *should* be
    /// let matrix_reindexed_ground_truth              
    ///                         =    VecOfVec::new(
    ///                                  vec![
    ///                                      vec![ (2,"b"), (3,"a"), ],
    ///                                      vec![ (2,"d"), (3,"c"), ],
    ///                                  ]
    ///                              ).ok().unwrap();
    /// 
    /// // Check that the permuted matrix is correct. We wrap `matrix_reindexed_ground_truth`
    /// // in an `Ok(..)` because the function `reassign_column_indices` does this. If you
    /// // ever need to get at the value inside the `Ok(..)`, you can use `unwrap` or any
    /// // of the other methods for `Result`.
    /// assert_eq!( matrix_reindexed, Ok(matrix_reindexed_ground_truth) );
    /// ```
    /// 
    /// Here is an example of a re-indexing which fails because it assigns the same index to two different columns:
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    ///  
    /// let matrix              =    VecOfVec::new(
    ///                                  vec![
    ///                                      vec![ (0,"a"), (1,"b"), ],
    ///                                      vec![ (0,"c"), (1,"d"), ],
    ///                                  ]
    ///                              ).ok().unwrap();
    /// 
    /// let reindexed_entries   =    vec![
    ///                                 vec![ (0,"a"), (0,"b"), ],
    ///                                 vec![ (0,"c"), (0,"d"), ],
    ///                              ];
    /// 
    /// // Vectors implement the EvaluateFunction trait. This vector will map 0 to 0 and 1 to 0.
    /// let column_index_map    =   vec![ 0, 0 ]; 
    /// let matrix_reindexed    =   matrix.reassign_column_indices_out_of_place( & column_index_map );
    ///
    /// assert_eq!( matrix_reindexed, Err(reindexed_entries) );
    /// ```
    pub fn reassign_column_indices_out_of_place< 'a, ColumnIndexMap, NewColumnIndex >( &'a self, column_index_map: ColumnIndexMap ) 
            -> 
            Result< 
                VecOfVec< NewColumnIndex, Coefficient >,
                Vec<Vec< (NewColumnIndex, Coefficient) >>,
            >
            where 
                ColumnIndexMap:         EvaluateFunction< & 'a ColumnIndex, NewColumnIndex >,
                ColumnIndex:            'a,
                Coefficient:            Clone + Debug,
                NewColumnIndex:         Clone + Debug + Ord,
    {
        let mut vec_of_vec  =   Vec::with_capacity( self.number_of_rows() ); // pre-allocate a list of lists
        for original_row in self.inner_vec_of_vec_ref().iter() { // for each list in the orgiinal list of lists
            let mut new_row     =   Vec::with_capacity( original_row.len() ); // pre-allocate a new inner list
            for (i,a) in original_row { // for each entry in the original list
                let j   =   column_index_map.evaluate_function( i ); // calculate the new column index
                new_row.push( (j, a.clone()) ) // push the new entry to the new vector
            }

            // now we have re-indexed all the entries in the old row. the entries in the new row may not be sorted, we we have to sort them.
            // we can sort the row using the usual Rust sort operation for vectors, because VecOfVec uses the defaul order on indices
            // that is provided by Rust.
            new_row.sort_unstable_by( |x,y| x.0.cmp( & y.0 ));

            // push the new list to the list-of-lists
            vec_of_vec.push( new_row );
        }

        // it is possible that the re-assignment function has assigned the same index to two different columns.
        // in that case the following constructor will return an error; otherwise it will return a valid VecOfVec
        VecOfVec::new( vec_of_vec )
    }


    /// Append a vector to the internally stored `Vec<Vec<EntryType>>` struct.  The new vector must be sorted in strictly ascending order; the function panics, otherwise.
    pub fn push_row( &mut self, new_sorted_vec: Vec < (ColumnIndex, Coefficient) > ) 
        where 
            OrderOperatorByKey: JudgePartialOrder< (ColumnIndex, Coefficient)>
    {
        if ! is_sorted_strictly( & new_sorted_vec, &mut self.order_operator ) {
            panic!("Attempt to append a non-strictly-sorted vector to `sorted::VecOfVec`.  To proceed, the new vector must be sorted in *strictly* ascending order, by index.") 
        }        
        self.vec_of_vec.push( new_sorted_vec );
    }  

    /// Replace a row of the matrix, returning the previous row as a vector.
    /// 
    /// If the new row is not strictly sorted by index, or if the user attempts to insert a row at index `i`
    /// in matrix that does not have a row at index `i` because it has `i` or fewer rows, the
    /// function will leave the matrix unchanged and return `Err( vector_that_we_tried_to_insert ) `.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// 
    /// // Define the matrix
    /// let mut matrix          =    VecOfVec::new(
    ///     vec![
    ///         vec![ (0,"a"), (1,"b"), ],
    ///         vec![ (0,"c"), (1,"d"), ],
    ///     ]
    /// ).ok().unwrap();
    ///
    /// // Successfully replace a row
    /// // --------------------------
    ///      
    /// let old_row             =   matrix.replace_row_and_return_old( 
    ///             1,                              // the row index
    ///             vec![ (0,"e"), (1,"f") ],       // the new row
    ///         );
    ///      
    ///      
    /// // Check that the new matrix is correct
    /// let desired_matrix      =    VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"a"), (1,"b"), ],
    ///                 vec![ (0,"e"), (1,"f"), ],
    ///             ]
    ///         ).ok().unwrap();
    /// 
    /// assert_eq!( matrix, desired_matrix );
    ///      
    /// // The update function returned the old row:
    /// assert_eq!( old_row, Ok( vec![ (0,"c"), (1,"d"), ] ) );
    ///      
    ///      
    /// // Attempt to insert a row at an index which is too high
    /// // -----------------------------------------------------
    ///      
    /// let old_row             =   matrix.replace_row_and_return_old( 
    ///             2,                              // the row index
    ///             vec![ (0,"g"), (1,"h") ],       // the new row
    ///         );
    ///      
    /// // The update returns the new row without inserting it into the matrix:
    /// assert_eq!( old_row, Err( vec![ (0,"g"), (1,"h"), ] ) );
    ///      
    /// // Attempt to insert a row that is not sorted
    /// // -----------------------------------------------------
    ///      
    /// let old_row             =   matrix.replace_row_and_return_old( 
    ///             1,                              // the row index
    ///             vec![ (1,"g"), (0,"h") ],       // the new row
    ///         );
    ///      
    /// // The update returned the new row without inserting it into the matrix:
    /// assert_eq!( old_row, Err( vec![ (1,"g"), (0,"h"), ] ) ); 
    /// ```
    pub fn replace_row_and_return_old( &mut self, row_index: usize, new_sorted_vec: Vec < (ColumnIndex, Coefficient) > ) 
            -> 
            Result<
                Vec < (ColumnIndex, Coefficient) >, 
                Vec < (ColumnIndex, Coefficient) >,
            >
        where 
            OrderOperatorByKey: JudgePartialOrder< (ColumnIndex, Coefficient)>    
    {
        // make sure the index is not too high
        if row_index >=  self.number_of_rows() {
            println!("Attempted to replace a row at index {:?}, but the matrix only has {:?} rows.", row_index, self.number_of_rows() );
            return Err( new_sorted_vec )
        }        

        // make sure the new row is sorted
        if ! is_sorted_strictly( & new_sorted_vec, &mut self.order_operator ) {
            println!("Attempted to replace a row of a `sorted::VecOfVec` with a non-strictly-sorted vector.  To proceed, the new vector must be sorted in *strictly* ascending order, by index.");
            return Err( new_sorted_vec )
        }        


        // swap the new and old rows
        let old_row     =   mem::replace( 
                                &mut self.vec_of_vec[ row_index ], 
                                new_sorted_vec
                            );

        // return the old row
        return Ok( old_row )
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

    /// Creates a [`VecOfVec`] from an iterable of iterables that run over tuples.
    /// 
    /// This function first turns the iterable of iterables into a variable `vec_of_vec` which has type `Vec<Vec<(index,coefficient)>>`. Then it
    /// calls [VecOfVec::new] on the vector of vectors to create a [VecOfVec]. 
    /// 
    /// **NB** The function [VecOfVec::new] returns `Ok(VecOfVec)` if
    /// the vectors inside `vec_of_vec` are sorted in strictly increasing order according to index, and returns `Err(vec_of_vec)` otherwise.
    /// 
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    ///     
    /// let iter = (0..2).map( |x| vec![(x,x)] );
    /// let vec_of_vec = VecOfVec::from_iterable_of_iterables( iter ).ok().unwrap();
    /// let ground_truth =  VecOfVec::new(
    ///                         vec![ 
    ///                             vec![ (0,0)        ],
    ///                             vec![        (1,1) ],
    ///                         ]
    ///                     ).ok().unwrap();
    /// assert_eq!( vec_of_vec, ground_truth  );
    /// ```
    pub fn from_iterable_of_iterables< I >( iter: I ) 
            -> 
            Result<
                VecOfVec< ColumnIndex, Coefficient >,
                Vec<Vec<(ColumnIndex, Coefficient)>>,
            >
        where
            I:          IntoIterator,
            I::Item:    IntoIterator< Item = (ColumnIndex, Coefficient) >,
            ColumnIndex:     Clone + Debug + Ord,
            Coefficient:     Clone + Debug,
    {
        let vecvec =    iter.into_iter().map( |x| x.into_iter().collect_vec() ).collect_vec();
        VecOfVec::new( vecvec )
    }

    /// Wraps `self` in a [MatrixAlgebraPacket] with a user-specific ring operator
    /// 
    /// Wrapping a vec-of-vec in a [MatrixAlgebraPacket] allows you to perform algebra operations,
    /// such as multiplying with a vector or another matrix.
    /// 
    /// Ordinarily one has to specify [order operators](crate::utilities::order) for row and column entries and indices.
    /// However, the [vec_of_vec::sorted::VecOfVec](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec) data
    /// structure always uses the default order on indices provided by Rust, so we don't have to specify order operators
    /// explicitly in this case.
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::matrices::operations::MatrixOracleOperations;
    /// use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
    /// 
    /// // define the matrix
    /// let matrix              =    VecOfVec::new(
    ///                                  vec![
    ///                                      vec![ (0,true), (1,true),  ],
    ///                                      vec![ (0,true), (1,false), ],
    ///                                  ]
    ///                              ).ok().unwrap();
    /// 
    /// // wrap it in a packet that uses the two element field {true, false}
    /// let ring_operator       =   BooleanField; // the two element field
    /// let algebraic_matrix    =   matrix.matrix_algebra_packet( ring_operator );
    /// 
    /// // multiply with a row vector
    /// let v                   =   algebraic_matrix
    ///                                 .multiply_with_row_vector( vec![(0,true), (1,true)]  )        // this line computes the product of the matrix with a vector
    ///                                 .collect::<Vec<_>>();                                   // this line collects the elements of the iterator into a Rust Vec
    /// 
    /// // check that the product is correct
    /// assert_eq!( v, vec![ (1,true) ] )
    /// ```    
    pub fn matrix_algebra_packet< RingOperator >( &self, ring_operator: RingOperator ) 
            ->
            MatrixAlgebraPacket< 
                &Self, 
                RingOperator, 
                OrderOperatorByKey, 
                OrderOperatorAuto, 
                OrderOperatorByKey, 
                OrderOperatorAuto 
            >
    {
        MatrixAlgebraPacket::with_default_order( &self, ring_operator)
    }




    /// Save a copy of `self` to file
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// // Create an instance of MyStruct
    /// let matrix = VecOfVec::new( vec![ vec![ (0,1)], vec![ (1,1) ] ] ).ok().unwrap();
    /// 
    /// // Create a temporary file in the system's temp directory
    /// let mut file_path: std::path::PathBuf = std::env::temp_dir();
    /// file_path.push("test_data.json");
    /// let file_path_str = file_path.to_str().unwrap();
    /// 
    /// // Save the struct to the file
    /// matrix.save_to_json(file_path_str).expect("Failed to save file");
    /// 
    /// // Load the struct from the file
    /// let loaded_matrix = VecOfVec::load_from_json(file_path_str).expect("Failed to load file");
    /// 
    /// // Assert that the loaded data is the same as the original
    /// assert_eq!(matrix, loaded_matrix);
    /// 
    /// // Clean up: delete the file after the test
    /// std::fs::remove_file(file_path_str).expect("Failed to delete file");
    /// ```
    pub fn save_to_json(
        & self,
        file_path: &str
    ) -> std::io::Result<()> 
    where
        ColumnIndex:        Serialize,
        Coefficient:        Serialize,

    {
        use std::io::Write;

        let file = std::fs::File::create(file_path)?;
        let json = serde_json::to_string( &self ).expect("Failed to serialize");
        writeln!(&file, "{}", json)?;
        Ok(())
    }



    /// Load a `VecOfVec` from a JSON file
    /// 
    /// This method only works for files created by [VecOfVec::save_to_json]
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// // Create an instance of VecOfVec
    /// let matrix = VecOfVec::new( vec![ vec![ (0,1)], vec![ (1,1) ] ] ).ok().unwrap();
    /// 
    /// // Create a temporary file in the system's temp directory
    /// let mut file_path: std::path::PathBuf = std::env::temp_dir();
    /// file_path.push("test_data.json");
    /// let file_path_str = file_path.to_str().unwrap();
    /// 
    /// // Save the struct to the file
    /// matrix.save_to_json(file_path_str).expect("Failed to save file");
    /// 
    /// // Load the struct from the file
    /// let loaded_matrix = VecOfVec::load_from_json(file_path_str).expect("Failed to load file");
    /// 
    /// // Assert that the loaded data is the same as the original
    /// assert_eq!(matrix, loaded_matrix);
    /// 
    /// // Clean up: delete the file after the test
    /// std::fs::remove_file(file_path_str).expect("Failed to delete file");
    /// ```
    pub fn load_from_json(file_path: &str) -> std::io::Result<Self> 
    where
        ColumnIndex:    for<'a> Deserialize<'a>,
        Coefficient:    for<'a> Deserialize<'a>,    
    {
        use std::io::Read;

        let mut file = std::fs::File::open(file_path)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        let my_struct: Self = serde_json::from_str(&contents).expect("Failed to deserialize");
        Ok(my_struct)
    }
        



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


    /// Row and column where the maximum column index appears in an explicitly stored entry (or `None``, if no entries are stored)
    /// 
    /// Concretely, returns `(i,j)` such that row `i`, column `j` is a structural nonzero element of the matrix and `j` is as large as possible.
    /// 
    /// If there are multiple rows where the maximum column index appear, then this function returns the first.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// // example with a nonempty matrix
    /// let a: Vec<Vec<(usize,f64)>>    =   vec![
    ///                                         vec![          ],
    ///                                         vec![          ],
    ///                                         vec![ (0,1.)   ], 
    ///                                     ];   // matrix with 1 entry 
    /// let matrix = VecOfVec::new(a).ok().unwrap();        
    /// assert_eq!( matrix.row_and_column_with_max_column_index(), Some((2,0)) );
    ///
    /// // example with an empty matrix        
    /// let b: Vec<Vec<(usize,f64)>>    =   Vec::new();   
    /// let matrix = VecOfVec::new(b).ok().unwrap();                      
    /// assert_eq!( matrix.row_and_column_with_max_column_index(), None      );  
    /// ```
    pub fn row_and_column_with_max_column_index( &self ) -> Option< (usize,usize) > {
        let mut max_col = None;
        let mut max_row = None;
        let mut potential_new_max_col;
        for (row_index, row) in self.vec_of_vec.iter().enumerate() {
            potential_new_max_col = row.iter().map(|x| x.0).max();
            if potential_new_max_col > max_col {
                max_col                 =   potential_new_max_col.clone();
                max_row                 =   Some(row_index);
            }
        }     
        match (max_row, max_col) {
            (Some(a),Some(b)) => Some((a,b)),
            _ => None
        }
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




    /// Returns an `n x n` diagonal matrix
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    ///     
    ///     
    /// // Construct the matrix
    /// let diagonal_element    =   1usize;
    /// let size                =   2usize;
    /// let matrix              =   VecOfVec::diagonal_matrix( diagonal_element, size );
    /// 
    /// // Verify the construction
    /// let ground_truth: VecOfVec<usize,usize>         =   VecOfVec::new( 
    ///                                                         vec![ 
    ///                                                             vec![(0,1)        ], 
    ///                                                             vec![       (1,1) ],
    ///                                                         ] 
    ///                                                     ).ok().unwrap(); 
    /// assert_eq!( matrix, ground_truth )
    /// ```
    pub fn diagonal_matrix( diagonal_element: Coefficient, size: usize )  -> Self 
        where 
            Coefficient:              Clone + Debug,
    {
        let mut outer       =   Vec::with_capacity(size);

        for p in 0 .. size {
            let mut inner       =   Vec::with_capacity(1);
            inner.push( (p, diagonal_element.clone()) );
            outer.push( inner );
        }
        return VecOfVec::new(outer).ok().unwrap()
    }







    /// Reverses the order of columns.
    /// 
    /// Concretely, we replace each row `[(i1, a1), .., (in, an)]` with `[ (max_column_index - in, an), .., (max_column_index - i1, a1) ]`
    /// where `max_column_index` is a parameter passed by the user.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use std::slice::Iter;
    ///  
    /// let a = VecOfVec::new(
    ///                             vec![
    ///                                 vec![   (0, "a"),  (1, "b"),  (2, "c")   ],
    ///                                 vec![              (1, "d"),  (2, "e"),  ],
    ///                             ]
    ///                         ).ok().unwrap();
    /// let mut b = VecOfVec::new(
    ///                             vec![
    ///                                 vec![   (0, "c"),  (1, "b"),  (2, "a")   ],
    ///                                 vec![   (0, "e"),  (1, "d"),             ],
    ///                             ]
    ///                         ).ok().unwrap();
    /// b.reverse_the_sequence_of_columns_in_place( 2 );
    /// assert_eq!( a, b );
    /// ```
    pub fn reverse_the_sequence_of_columns_in_place( & mut self, max_column_index: usize ) -> Result<(),()> {
        for row in self.vec_of_vec.iter_mut() {
            for (i,_a) in row.iter_mut() {
                if *i > max_column_index {
                    return Err(())
                } else {
                    *i = max_column_index - i.clone();
                }
            }
            row.reverse();
        }
        return Ok(())
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
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
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
    ///         ).ok().unwrap();
    /// 
    /// let transpose   
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"0,0"), (1,"1,0"), ],
    ///                 vec![ (0,"0,1"), (1,"1,1"), ],
    ///             ]
    ///         ).ok().unwrap();
    /// 
    /// let transpose_lazy  =   matrix.transpose();
    /// let transpose_deep  =   & matrix.transpose_deep(2).unwrap();
    ///
    ///  
    /// for p in 0..2 {
    ///     let a   =   transpose_deep.row(&p).collect_vec();
    ///     let b   =   transpose_lazy.row(&p).collect_vec();
    ///     let c   =   transpose.row(&p).collect_vec();
    ///     assert_eq!(&a, &b);
    ///     assert_eq!(&a, &c);
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


    /// Generates a new [VecOfVec] representing the anti-transpose, with copied (not borrowed) data.
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
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
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
    ///         ).ok().unwrap();
    ///                 
    /// let antitranspose   
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"1,2"), (1,"0,2"), ],
    ///                 vec![ (0,"1,1"), (1,"0,1"), ],
    ///                 vec![ (0,"1,0"), (1,"0,0"), ],                        
    ///             ]
    ///         ).ok().unwrap();
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
        let mut transpose_vecvec    =   transpose.inner_vec_of_vec();
        transpose_vecvec.reverse();
        let num_rows_in_transpose = transpose_vecvec.len();
        for p in 0 .. num_rows_in_transpose {
            (* transpose_vecvec[ p ] ).reverse();
        }        
        for row in transpose_vecvec.iter_mut() {  
            row.iter_mut()
                .for_each(
                    |(j,_v)|
                    *j = num_rows_original - *j -1
                );      
            println!("row after subtracting: {:?}", row);                 
        }
        return VecOfVec::new( transpose_vecvec ).ok();
    }    

    /// Returns a [MatrixBiajorData](crate::algebra::matrices::types::bimajor::MatrixBimajorData); this doubles memory use but gives fast access to columns
    /// 
    /// Concretely, the [MatrixBiajorData](crate::algebra::matrices::types::bimajor::MatrixBimajorData) stores a copy of the matrix and its
    /// transpose; then the user needs to look up a column, it looks up the corresponding row of the transpose, which is much faster than
    /// looking up a column of a row-major matrix.
    /// 
    /// This doubles memory usage; however it allows efficient access to **both rows and columns**.
    /// 
    /// # Arguments
    /// 
    /// If `num_rows_of_transpose` is less than 1 + the maximum column index of a
    /// any explicitly stored entry in `self`, then `None` is returned.  
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






    /// Compute the product `self * other_vec_of_vec_matrix` and write the result to another [VecOfVec]
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
    /// 
    /// // define the matrix
    /// let matrix              =    VecOfVec::new(
    ///                                  vec![
    ///                                      vec![ (0,true), (1,true),  ],
    ///                                      vec![ (0,true),            ],
    ///                                  ]
    ///                              ).ok().unwrap();
    ///                          
    /// // define an ring operator for the two element field {true, false}
    /// let ring_operator       =   BooleanField::new(); // the two element field
    ///                          
    /// // multiply the matrix with itself
    /// let product             =   matrix.multiply_on_the_left_and_write_the_product_to_a_vec_of_vec( 
    ///                                 & matrix, 
    ///                                 ring_operator 
    ///                             );
    ///                          
    /// // check the calculation
    /// let ground_truth        =    VecOfVec::new(
    ///                                  vec![
    ///                                      vec![            (1,true),  ],
    ///                                      vec![ (0,true),  (1,true),  ],
    ///                                  ]
    ///                              ).ok().unwrap();
    /// assert_eq!( product, Ok( ground_truth ) );   
    /// ```  
    pub fn multiply_on_the_left_and_write_the_product_to_a_vec_of_vec< RighthandMatrixColumnIndex, RingOperator >( 
            &self, 
            right_matrix:       & VecOfVec< RighthandMatrixColumnIndex, Coefficient >,
            ring_operator:      RingOperator,
        )
        -> 
        Result<
            VecOfVec< RighthandMatrixColumnIndex, Coefficient >,
            HashMap< &str, usize >,
        >

        where
            RighthandMatrixColumnIndex:        Clone + Debug + Ord + Eq,                
            Coefficient:                        Clone + Debug + PartialEq,
            RingOperator:                       Clone + SemiringOperations< Element = Coefficient >,            
    {
        // check to make sure that none of the column indices of the lefthand matrix exceeds the number of rows of the righthand matrix
        let right_matrix_number_of_rows                =   right_matrix.number_of_rows();
        let max_row_and_column_index            =   self.row_and_column_with_max_column_index();

        if let Some( (i,j) ) = max_row_and_column_index {
            if j >= right_matrix_number_of_rows {
                let mut err                     =   HashMap::new();
                err.insert("Row of the lefthand matrix: ", i );
                err.insert("Max column index of the specified row of the lefthand matrix: ", j );
                err.insert("Number of rows of righthand matrix: ", right_matrix_number_of_rows );
                return Err( err )
            }
        }

        // compute the product
        let algebra_packet      =   right_matrix.matrix_algebra_packet(ring_operator).into_peekable_matrix();
        let mut product_matrix         =   Vec::with_capacity( self.number_of_rows() );
        for row in self.inner_vec_of_vec_ref() {
            let product_row     =   algebra_packet.multiply_with_row_vector(row).collect::<Vec<_>>();
            product_matrix.push( product_row );
        }

        Ok(VecOfVec::new( product_matrix ).ok().unwrap())
    }





    /// Returns a generalized inverse formatted as a `VecOfVec`
    /// 
    /// The generalized inverse is computed as `S M^- T^{-1}`, where `(T,M,D,S)` is a proper U-match factorization
    /// and `M^-` is the generalized inverse of `M` computed by transposing `M` and inverting its nonzero entries.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
    /// 
    /// // define the matrix
    /// let matrix              =    VecOfVec::new(
    ///                                  vec![
    ///                                      vec![ (0,true), (1,true ),  ],
    ///                                      vec![ (0,true), (1,false),  ],
    ///                                  ]
    ///                              ).ok().unwrap();
    ///                          
    /// // the inverse o fthe matrix is the (unique) generalized inverse
    /// let inverse             =    VecOfVec::new(
    ///                                  vec![
    ///                                      vec![            (1,true ), ],
    ///                                      vec![ (0,true ), (1,true ), ],
    ///                                  ]
    ///                              ).ok().unwrap();
    ///                          
    /// // define an ring operator for the two element field {true, false}
    /// let ring_operator       =   BooleanField::new(); // the two element field
    ///                          
    /// // multiply the matrix with itself
    /// let number_of_columns   =   2;
    /// let generalized_inverse =   matrix.generalized_inverse( ring_operator, number_of_columns );
    ///                          
    /// // check the calculation
    /// assert_eq!( generalized_inverse, inverse );
    /// ```
    pub fn generalized_inverse< RingOperator >( & self, ring_operator: RingOperator, number_of_columns: usize ) -> Self 
        where
            RingOperator:       Clone + DivisionRingOperations< Element = Coefficient >,   
            Coefficient:        Clone + Debug + PartialEq,    
    
    {
        let packet                  =   MatrixAlgebraPacket::with_default_order(
                                                                & self, 
                                                                ring_operator.clone(),
                                                            );
        let umatch                                      =   Umatch::new( 
                                                                        packet,  
                                                                        ( 0 .. self.number_of_rows() ).rev()
                                                                    );

        let t_inv                                       =   umatch.target_comb_inverse();
        let s                                           =   umatch.source_comb();
        let m_inv                                       =   umatch.generalized_matching_matrix_ref().generalized_inverse( ring_operator.clone() );
        let m_inv_packet                                =   MatrixAlgebraPacket::with_default_order(m_inv, ring_operator.clone());

        let generalized_inverse                         =   s.multiply_on_the_left_of( m_inv_packet )
                                                                                                    .multiply_on_the_left_of( t_inv );
        let row_iterator                                =   ( 0 .. number_of_columns )
                                                                            .map( |i|  generalized_inverse.row( & i ) );
        
        VecOfVec::from_iterable_of_iterables( row_iterator ).ok().unwrap()
    }





    /// Returns a new [VecOfVec] with all entries that are strictly greater than `threshold` removed.
    /// 
    /// This is an out-of-place operation, i.e. it does not modify `self`.
    pub fn strip_values_above_threshold(
            &self,
            threshold: Coefficient
        ) -> Self
        where
            Coefficient: PartialOrd + Clone + Debug,
    {
        let mut vec_of_vec = Vec::with_capacity(self.vec_of_vec.len());
        for row in self.vec_of_vec.iter() {
            let new_row = row.iter()
                .filter(|(_, value)| value <= &threshold)
                .cloned()
                .collect_vec();
            vec_of_vec.push(new_row);
        }
        VecOfVec::new(vec_of_vec).ok().unwrap()
    }


    /// Constructs a [VecOfVec] from a ragged matrix
    /// 
    /// # Arguments
    /// 
    /// The input to this function is a vector of vectors of coefficients. Not all vectors need to have the same length.
    /// 
    /// ```text
    /// [
    ///     [a0, ],
    ///     [b0, b1, b2, ],
    ///     [c0, c1, ]
    /// ]
    /// ```
    /// 
    /// # Returns
    /// 
    ///  A [VecOfVec] where each row is a vector of tuples `(column_index, coefficient)`.  For the example above,
    ///   the output would be:
    /// ```text
    /// [
    ///     [(0, a0), ],
    ///     [(0, b0), (1, b1), (2, b2)],
    ///     [(0, c0), (1, c1)],
    /// ]
    /// ```
    /// 
    /// # Example
    /// 
    /// ```rust
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// 
    /// // Generate a matrix from ragged data
    /// let ragged_data: Vec< Vec< isize > > = vec![
    ///       vec![ 1 ],
    ///       vec![ 2, 3 ],
    ///       vec![ 4, 5, 6 ],
    /// ];
    /// let ragged_construction = VecOfVec::from_ragged(&ragged_data);
    ///  
    /// // Generate the same matrix from standard data
    /// let standard_data: Vec<Vec<(usize, isize)>> = vec![
    ///       vec![(0, 1)],
    ///       vec![(0, 2), (1, 3)],
    ///       vec![(0, 4), (1, 5), (2, 6)],
    /// ];
    /// let standard_construction = VecOfVec::new(standard_data).ok().unwrap();
    /// 
    /// // Check that the two constructions yield the same result    
    ///  assert_eq!(ragged_construction, standard_construction);
    /// ```
    pub fn from_ragged(
            ragged_matrix: & Vec<Vec< Coefficient >>
        ) ->
        Self

        where 
            Coefficient: Clone + Debug,
    {
        let mut vec_of_vec = Vec::with_capacity(ragged_matrix.len());
        for row in ragged_matrix.iter() {
            let mut new_row = Vec::with_capacity(row.len());
            for entry in row.iter().cloned().enumerate() {
                new_row.push(entry);
            }
            vec_of_vec.push(new_row);
        }
        VecOfVec::new(vec_of_vec).ok().unwrap()
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

        VecOfVec::new(vec_of_vec).ok().unwrap() // formally wrap the matrix in a VecOfVec struct
    }

    /// Generate a random invertible upper triangular matrix with integer coefficients.
    /// 
    /// - Entries on the diagonal are equal to 1.
    /// - Each entry strictly above the diagonal is either (i) strucutrally zero, 
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

        VecOfVec::new( vec_of_vec ).ok().unwrap() // formally wrap the matrix in a VecOfVec struct
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

        VecOfVec::new(vec_of_vec).ok().unwrap()// formally wrap the matrix in a VecOfVec struct
    }

    /// Generate a random `VecOfVec< usize, usize >`
    /// 
    /// For each entry we flip a coin. If the value is true we draw a (possibly zero) entry in an iid fasion from [0, .., `modulus`).
    /// Otherwise the entry is structurally zero.
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

        VecOfVec::new( vecvec ).ok().unwrap()

    }    
}









impl VecOfVec< usize, OrderedFloat<f64> > {



    /// Generates a random symmetric `number_of_points x number_of_points` matrix
    /// 
    /// - the process begins by generating a random symmetric matrix with entries draw in an iid fashion from the uniform distribution on `[0,1]`.
    /// - diagonal entries are set to zero
    /// - values above the enclosing radius threshold removed
    pub fn random_symmetric_zero_diagonal_with_enclosing_radius_threshold(
            number_of_points: usize,
        ) -> Self

    {
        use crate::utilities::iterators::general::minmax;   

        // initialize a random dissimilarity matrix
        let dissimilarity_matrix_data = crate::utilities::random::random_symmetric_matrix_zero_diag(number_of_points);
        
        // get enclosing radius
        let dissimilarity_value_max = 
            minmax( 
                    (0..number_of_points).map(
                            |x| 
                            dissimilarity_matrix_data[x].iter()
                        ) 
                ).unwrap(); 

        // strip large values and stroe in VecOfVec                
        let dissimilarity_matrix_vecvec = VecOfVec::from_ragged(
            &dissimilarity_matrix_data // this data isn't ragged but the construction works just the same
        ).strip_values_above_threshold(*dissimilarity_value_max);              

        // check that the matrix is symmetric
        for i in 0 .. number_of_points {
            for j in i .. number_of_points {
                assert_eq!( 
                    (&dissimilarity_matrix_vecvec).structural_nonzero_entry(&i,&j), 
                    (&dissimilarity_matrix_vecvec).structural_nonzero_entry(&j,&i) 
                );
            }
        }

        return dissimilarity_matrix_vecvec
    }




}





































//  ORACLE IMPLEMENTATIONS FOR &'a VecOfVec
//  ---------------------------------------------


impl < 'a, ColumnIndex, Coefficient > 

    MatrixOracle for 
    &'a VecOfVec < ColumnIndex, Coefficient >

    where   ColumnIndex:        Clone + Debug + Ord + Eq,    
            Coefficient:        Clone + Debug + PartialEq,    

{ 
    type Coefficient        =   Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   usize;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   ColumnIndex;       // The type of column indices
    
    type RowEntry           =   ( ColumnIndex, Coefficient );          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   ( usize, Coefficient );       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                =   Cloned< Iter< 'a, (ColumnIndex, Coefficient) > >;  // What you get when you ask for a row.
    type RowReverse         =   Cloned<Rev<std::slice::Iter<'a, (ColumnIndex, Coefficient)>>>;  // What you get when you ask for a row with the order of entries reversed
    
    type Column             =   VecOfVecMatrixColumn< 'a, ColumnIndex, Coefficient >; // What you get when you ask for a column
    type ColumnReverse      =   VecOfVecMatrixColumnReverse< 'a, ColumnIndex, Coefficient >; // What you get when you ask for a column with the order of entries reversed 

    // index validation
    fn has_column_for_index(  &   self, _index: & Self::ColumnIndex)   -> bool { true } // always returns true
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool {
        index < & self.vec_of_vec.len() // return true iff index is strictly less than the number of rows
    }

    // entry lookup
    fn structural_nonzero_entry( & self, row: &Self::RowIndex, column: &Self::ColumnIndex )   ->  Option< Self::Coefficient > {
        let row = & self.vec_of_vec[ *row ];
        find_sorted_binary_oracle( 
                    0, 
                    row.len() as isize - 1, 
                    |p| column.cmp( & row[p as usize].0 ) 
                ).map(|x| row[ x as usize ].1.clone() )        
    }  

    // row lookup
    fn row(                     &self,  index: &Self::RowIndex    )   -> Self::Row   { 
        self.vec_of_vec[*index].iter().cloned() 
    }
    fn row_reverse(             &self,  index: &Self::RowIndex    )       -> Self::RowReverse  { 
        self.vec_of_vec[*index].iter().rev().cloned()
    }
    
    fn column(                  &self,  index: &Self::ColumnIndex )       -> Self::Column {
        VecOfVecMatrixColumn{
            vec_of_vec:                     self,
            row_index:                      0,
            column_index:                   index.clone(),
            phantom_column_index:           PhantomData,            
        }
    }
    fn column_result(              &   self, index: & Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex > {
        Ok( self.column(index) )
    }
   
    fn column_reverse(          &self,  index: &Self::ColumnIndex )       -> Self::ColumnReverse  { 
        VecOfVecMatrixColumnReverse{
            vec_of_vec:                     self,
            row_index:                      self.vec_of_vec.len(),
            column_index:                   index.clone(),
            phantom_column_index:           PhantomData,            
        }
    }            
    fn column_reverse_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex > {
        Ok( self.column_reverse(index) )
    }
}   










//  MATRIX ORACLE OPERATIONS
//  ---------------------------------------------------------------------


impl < 'a, ColumnIndex, Coefficient > 

    MatrixOracleOperations for 

    &'a VecOfVec < ColumnIndex, Coefficient >
{}    













/// Represents a column of a `VecOfVec`, with entries appearing in descending order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::{VecOfVec, VecOfVecMatrixColumnReverse};
/// use oat_rust::algebra::matrices::query::MatrixOracle;
///         
/// // define a sparse matrix representing the array
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] ).ok().unwrap();
/// itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).column_reverse( &0 ) );
/// itertools::assert_equal( vec![ (0,5) ], (& matrix).column_reverse( &1 ) );
/// itertools::assert_equal( vec![ (1,7), (0,5) ], (& matrix).column_reverse( &2 ) );  
/// ```
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, )]
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
        let vecvec_data =  self.vec_of_vec.inner_vec_of_vec_ref();                
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
/// // define a sparse matrix representing the array
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] ).ok().unwrap();
/// itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).column( &0 ) );
/// itertools::assert_equal( vec![ (0,5) ], (& matrix).column( &1 ) );
/// itertools::assert_equal( vec![ (0,5), (1,7) ], (& matrix).column( &2 ) );  
/// ```
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, )]
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
        let vecvec_data =  self.vec_of_vec.inner_vec_of_vec_ref();                
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










/// Runs over the `(row, column, coefficient)` triplets of the structural nonzero entries of a `VecOfVec`.
pub struct VecOfVecTriplets< 'a, ColumnIndex, Coefficient >
    where   ColumnIndex:     Clone,
            Coefficient:     Clone,
{
    vec_of_vec_struct_ref:      &'a VecOfVec< ColumnIndex, Coefficient >,
    row_index:                  usize,
    column_pointers:            std::ops::Range< usize >,
}







impl < 'a, ColumnIndex, Coefficient >

    VecOfVecTriplets< 'a, ColumnIndex, Coefficient >

    where   ColumnIndex:     Clone,
            Coefficient:     Clone,
{
    pub fn new( vec_of_vec: &'a VecOfVec< ColumnIndex, Coefficient > ) -> Self {

        let max_column_index = vec_of_vec.number_of_structural_nonzeros_in_row(
            0
        ).unwrap_or(0);
        let column_pointers = 0 .. max_column_index;


        VecOfVecTriplets {
            vec_of_vec_struct_ref: vec_of_vec,
            row_index: 0,
            column_pointers: column_pointers,
        }
    }
}








impl < 'a, ColumnIndex, Coefficient >

        Iterator for

        VecOfVecTriplets< 'a, ColumnIndex, Coefficient >

    where   ColumnIndex:     Clone,
            Coefficient:     Clone,        
{
    type Item = (usize, ColumnIndex, Coefficient);

    fn next( &mut self ) -> Option< Self::Item > {

        loop {

            let inner_vec_of_vecref = self.vec_of_vec_struct_ref.inner_vec_of_vec_ref();

            if self.row_index >= inner_vec_of_vecref.len() {
                // if we are at the end of the matrix, then return None
                return None;
            }

            match self.column_pointers.next() {
                Some( column_pointer ) => {
                    // return the triplet for this row and column
                    let (column_index, coefficient) = inner_vec_of_vecref[self.row_index][column_pointer].clone();
                    return Some( (self.row_index, column_index, coefficient) )
                } None => {
                    match self.row_index + 1 >= self.vec_of_vec_struct_ref.number_of_rows()  {
                        true => {
                            // if we are at the end of the matrix, then return None
                            return None;
                        },
                        false => {
                            // otherwise, increment the row index and reset the column pointer
                            self.row_index += 1;
                            self.column_pointers = 0 .. inner_vec_of_vecref[self.row_index].len();
                        }
                    }
                }
            }
        }

    }
}














#[cfg(test)]
mod tests {    

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
            let a   =   antitranspose.row( & p ).collect_vec();
            let b   =   antitranspose_deep.row( & p ).collect_vec();
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
            assert_eq!(     ( & matrix  ).row( & p ).collect_vec(),
                            ( & bimajor ).row( & p ).collect_vec(),     );
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

}
