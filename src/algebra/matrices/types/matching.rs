//! (Generalized) matching matrices
//! 
//! A generalized matching matrix is a matrix with at most one non-zero entry in each row and column.
//! This module provides several different data types for representing generalized matching matrices,



use crate::algebra::rings::traits::DivisionRingOperations;
use crate::algebra::vectors::entries::{KeyValGet, KeyValNew, ExtractKey, ExtractVal};
use crate::algebra::matrices::query::MatrixOracle;
use crate::algebra::matrices::operations::MatrixOracleOperations;

use crate::utilities::functions::compose::ComposeFunctions;
use crate::utilities::functions::evaluate::{LogicalNot, Map};
use crate::utilities::iterators::general::{TwoTypeIterator, Filter, FilterOutMembers, FilterOutNonmembers,};
use crate::utilities::sequences_and_ordinals::BijectiveSequence;
use crate::utilities::sets::MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper;
use std::collections::HashMap;
use std::hash::Hash;
use std::iter::{Empty, FromIterator, Zip};
use std::marker::PhantomData;
use std::fmt::Debug;


use derive_getters::Dissolve;

















//  ================================================================================
//  UNORDERED MATCHING MATRIX
//  ================================================================================





/// A generalized matching matrix
#[derive(Debug, Clone)]
pub struct GeneralizedMatchingMatrix
                < ColumnIndex, RowIndex, EntryMaj, EntryMin, Coefficient >
    where   ColumnIndex:     Clone + Hash + std::cmp::Eq,  // cloning seems unavoidable for push/insert operations, where we first need a row/column index as a value, then, later, as a key 
            RowIndex:     Clone + Hash + std::cmp::Eq,  // cloning seems unavoidable for push/insert operations, where we first need a row/column index as a value, then, later, as a key 
            EntryMin:   Clone + KeyValNew < Key = ColumnIndex, Val = Coefficient >, // clone is ncessary in order for HashMap< RowIndex, EntryMin > to implement `EvaluateFunction< RowIndex, EntryMin >`
            EntryMaj:   Clone + KeyValNew < Key = RowIndex, Val = Coefficient >, // clone is ncessary in order for HashMap< ColumnIndex, EntryMaj > to implement `EvaluateFunction< ColumnIndex, EntryMaj >`
            Coefficient:     Clone,                       
{  
    column_index_to_entry:    HashMap< ColumnIndex, EntryMaj >,
    row_index_to_entry:    HashMap< RowIndex, EntryMin >,
    phantom_coefficient:     PhantomData< Coefficient >
}


impl < 'a, ColumnIndex, RowIndex, EntryMin, EntryMaj, Coefficient > 

    GeneralizedMatchingMatrix
        < ColumnIndex, RowIndex, EntryMaj, EntryMin, Coefficient >

    where   ColumnIndex:     Clone + Hash + std::cmp::Eq, // cloning seems unavoidable for push/insert operations, where we first need a row/column index as a value, then, later, as a key 
            RowIndex:     Clone + Hash + std::cmp::Eq, // cloning seems unavoidable for push/insert operations, where we first need a row/column index as a value, then, later, as a key 
            EntryMin:   Clone + KeyValNew < Key = ColumnIndex, Val = Coefficient >,
            EntryMaj:   Clone + KeyValGet < Key = RowIndex, Val = Coefficient > + KeyValNew < Key = RowIndex, Val = Coefficient >,   
            Coefficient:     Clone, 
{
    /// Create a new, empty matching array.
    /// 
    /// See [`GeneralizedMatchingMatrixWithSequentialOrder`] for usage.
    pub fn new() -> GeneralizedMatchingMatrix< ColumnIndex, RowIndex, EntryMaj, EntryMin, Coefficient > 
    {
        GeneralizedMatchingMatrix{    
            column_index_to_entry:    HashMap::new(),
            row_index_to_entry:    HashMap::new(),
            phantom_coefficient:     PhantomData,
        }        
    }

    /// Push a new structural nonzero entry to the array.    
    pub fn push( 
                    &mut self,
                    column_index:     ColumnIndex, 
                    row_index:     RowIndex, 
                    coefficient:     Coefficient 
                )
    {
        let a = self.column_index_to_entry.insert( column_index.clone(), EntryMaj::new( row_index.clone(), coefficient.clone() ) ); 
        let b = self.row_index_to_entry.insert( row_index, EntryMin::new( column_index, coefficient ) );               
        if a.is_some()  { panic!("attempted to over-write the entry matched to a column index in a generalized matching matrix")  }
        if b.is_some()  { panic!("attempted to over-write the entry matched to a row index in a generalized matching matrix")  }
    }

        
    /// Returns an object `X` such that `X.evaluate_function( column_index ) = row_index`, where `column_index` is a matched column index
    /// and `row_index` is the associated row index; panics when the user calls `X.evaluate_function( column_index )` on
    /// an unmatched index `column_index`.    
    /// 
    /// # Design note
    /// 
    /// The function `f = function_row_index_to_column_index( &self )` works by cloning the entry associated with `row_index` then
    /// extracting the associated index.  This is not always the most efficient way to get the column index matched
    /// to a row index (some data structures allow more efficient work-arounds), but (i) it is as efficient as we have
    /// thus far been able to make it, without making limiting assumptions about the data structures involved, and (ii) 
    /// we don't currently foresee that this function will be used heavily.
    pub fn function_column_index_to_row_index( &self ) 
            -> 
            ComposeFunctions< 
                    ColumnIndex, &EntryMaj, RowIndex, 
                    & HashMap< ColumnIndex, EntryMaj >,                     
                    ExtractKey
                > 
    {
        ComposeFunctions::new( & self.column_index_to_entry, ExtractKey::new() )
    }

    /// Returns an object `X` such that `X.evaluate_function( row_index ) = column_index`, where `row_index` is a matched row index
    /// and `column_index` is the associated column index; panics when the user calls `X.evaluate_function( row_index )` on
    /// an unmatched index `row_index`.    
    /// 
    /// # Design note
    /// 
    /// The function `f = function_row_index_to_column_index( &self )` works by cloning the entry associated with `row_index` then
    /// extracting the associated index.  This is not always the most efficient way to get the column index matched
    /// to a row index (some data structures allow more efficient work-arounds), but (i) it is as efficient as we have
    /// thus far been able to make it, without making limiting assumptions about the data structures involved, and (ii) 
    /// we don't currently foresee that this function will be used heavily.
    pub fn function_row_index_to_column_index( &self ) 
            -> 
            ComposeFunctions< 
                    RowIndex, & EntryMin, ColumnIndex, 
                    & HashMap< RowIndex, EntryMin >,                     
                    ExtractKey
                > 
    {
        ComposeFunctions::new( & self.row_index_to_entry, ExtractKey::new() )
    }    

    /// Returns an object `X` such that `X.evaluate_function( column_index ) = coefficient`, where `column_index` is a matched column index
    /// and `coefficient` is the scalar value of the associated nonzero entry; panics when the user calls `X.evaluate_function( column_index )` on
    /// an unmatched index `column_index`.    
    /// 
    /// # Design note
    /// 
    /// The function `f = function_row_index_to_column_index( &self )` works by cloning the entry associated with `row_index` then
    /// extracting the associated scalar.  This is not always the most efficient way to get the column index matched
    /// to a row index (some data structures allow more efficient work-arounds), but (i) it is as efficient as we have
    /// thus far been able to make it, without making limiting assumptions about the data structures involved, and (ii) 
    /// we don't currently foresee that this function will be used heavily.
    pub fn function_column_index_to_coefficient( &self ) 
            -> 
            ComposeFunctions< 
                    ColumnIndex, & EntryMaj, Coefficient, 
                    & HashMap< ColumnIndex, EntryMaj >,                     
                    ExtractVal
                > 
    {
        ComposeFunctions::new( & self.column_index_to_entry, ExtractVal::new() )
    }

    /// Returns an object `X` such that `X.evaluate_function( row_index ) = coefficient`, where `row_index` is a matched row index
    /// and `coefficient` is the scalar value of the associated nonzero entry; panics when the user calls `X.evaluate_function( row_index )` on
    /// an unmatched index `row_index`.    
    /// 
    /// # Design note
    /// 
    /// The function `f = function_row_index_to_column_index( &self )` works by cloning the entry associated with `row_index` then
    /// extracting the associated scalar.  This is not always the most efficient way to get the scalar associated
    /// to a row index (some data structures allow more efficient work-arounds), but (i) it is as efficient as we have
    /// thus far been able to make it, without making limiting assumptions about the data structures involved, and (ii) 
    /// we don't currently foresee that this function will be used heavily.
    pub fn function_row_index_to_coefficient( &self ) 
            -> 
            ComposeFunctions< 
                    &RowIndex, & EntryMin, Coefficient, 
                    &HashMap< RowIndex, EntryMin >,                     
                    ExtractVal
                > 
    {
        ComposeFunctions::new( & self.row_index_to_entry, ExtractVal::new() )
    }    




    
    /// Returns an immutable reference `& HashMap< ColumnIndex, EntryMaj >` representing the partial bijection from 
    /// column indices to pairs of form `(associated_row_key, associated_nonzero_scalar)`.
    pub fn hashmap_ref_column_index_to_structural_nonzero_entry( &self ) -> & HashMap< ColumnIndex, EntryMaj > { & self.column_index_to_entry }

    /// Returns an immutable reference `& HashMap< RowIndex, EntryMaj >` representing the partial bijection from 
    /// row indices to pairs of form `(associated_minor_key, associated_nonzero_scalar)`.
    pub fn hashmap_ref_row_index_to_structural_nonzero_entry( &self ) -> & HashMap< RowIndex, EntryMin > { & self.row_index_to_entry }




    /// Returns the nonzero entry that corresponds to the row index, or returns `None` if the row index is unmatched.
    pub fn row_index_to_structural_nonzero_entry_refopt( &self, row_index: &RowIndex ) -> Option< & EntryMin > { 
        self.row_index_to_entry.get( row_index )
    }

    /// Returns the nonzero entry that corresponds to the column index, or returns `None` if the column index is unmatched.
    pub fn column_index_to_structural_nonzero_entry_refopt( &self, column_index: &ColumnIndex ) -> Option< & EntryMaj > { 
        self.column_index_to_entry.get( column_index )
    }    
 
    /// The nonzero coefficient associated with this `row_index` index.
    pub fn coefficient_opt_for_row_index( &self, row_index: & RowIndex ) -> Option< Coefficient >
        // where Coefficient: Clone    
    {
        self.row_index_to_entry.get( row_index ).map( |x| x.val() )
    }

    /// The nonzero coefficient associated with this `column_index` index.
    pub fn coefficient_opt_for_column_index( &self, column_index: & ColumnIndex ) -> Option< Coefficient >
        // where Coefficient: Clone    
    {
        self.column_index_to_entry.get( column_index ).map( |x| x.val() )
    }    

    /// Returns true iff the given column index is matched.
    pub fn has_a_match_for_column_index(  &self, column_index: &ColumnIndex ) -> bool {
        self.column_index_to_entry.contains_key( column_index )
    }

    /// Returns true iff the given row index is matched.
    pub fn has_a_match_for_row_index(  &self, row_index: &RowIndex ) -> bool {
        self.row_index_to_entry.contains_key( row_index )
    }    



}

























//  ================================================================================
//  GENERALIZED MATCHING ARRAY STRUCT WITH ORDINALS
//  ================================================================================

/// A generalized matching matrix (with an order imposed on entries)
/// 
/// 
/// The [GeneralizedMatchingMatrixWithSequentialOrder] struct encodes a generalized matching matrix, which is
/// a matrix that has at most one nonzero entry per row and column.
/// It also provides encodes an sequential ordering
/// `(r0, c0, x0), .. (rm, cm, xm)` on the nonzero entries of the matrix, where each entry is a
/// `(row_index, column_index, coefficient)` triplet. This ordering **does not respect** any ordering
/// on row or column indices. It is not generally true, for example, that `r0 < .. < r_m` or that `c0 < .. < cm`.
/// 
/// The special feature of this data structure is that it makes it easy to map back and forth
/// between each row index, its corresponding column index, and its **ordinal**, meaning
/// its order in the sequence. For example, it's easy to take a row index `ri` and find the 
/// corresponding column index `ci`, or the corresponding ordinal `i`.  
/// 
/// An unordered matching array can be encoded as a bijection {row indices} <-> {column indices}
/// together with a function {row indices} -> {structural nonzero values}.  By contrast, a 
/// `GeneralizedMatchingMatrixWithSequentialOrder` encodes bijections
/// {row indices} <-> {0, .., N} <-> {column indices} PLUS a function {0, ..., N} -> {structural nonzero values}.
///
/// Concretely, the struct stores these data in three different private fields
/// - a [BijectiveSequence](`crate::utilities::sequences_and_ordinals::BijectiveSequence`) 
///  object representing a bijection from row indices to {0, .., N}.
/// - a [BijectiveSequence](`crate::utilities::sequences_and_ordinals::BijectiveSequence`) 
///  object representing the corresponding bijection from column indices to {0, .., N}.
/// - a vector representing the sequence of structural nonzero values
/// 
/// # Technical details
/// 
/// This data structure will report that it has a row and column for every possible index.
/// That is, calling `self.has_row_for_index( &index )` and `self.has_column_for_index( &index )`
/// will return `true` for any index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder;
/// 
/// // initialize new matching array
/// let mut matching_array  =   GeneralizedMatchingMatrixWithSequentialOrder::new();
/// 
/// // push some nonzero entries
/// matching_array.push( 0, 1, 3 );
/// matching_array.push( 1, 2, 4 );
/// matching_array.push( 2, 0, 5 );    
/// 
/// assert_eq!( matching_array.structural_nonzero_values_in_sequence(), &vec![ 3, 4, 5 ] );
/// 
/// assert_eq!( matching_array.column_index_for_row_index( &0 ), Some(1) );  
/// assert_eq!( matching_array.row_index_for_column_index( &0 ), Some(2) );   
/// 
/// assert_eq!( matching_array.ordinal_for_row_index( &1 ), Some(1) ); 
/// assert_eq!( matching_array.ordinal_for_column_index( &1 ), Some(0) );
/// 
/// assert_eq!( matching_array.row_index_for_ordinal( 0 ), 0 );     
/// assert_eq!( matching_array.column_index_for_ordinal( 0 ), 1 );
/// 
/// assert_eq!( matching_array.coefficient_for_row_index( &0 ), 3 );       
/// assert_eq!( matching_array.coefficient_for_column_index( &0 ), 5 );
/// 
/// assert_eq!( matching_array.coefficient_for_ordinal( 0 ), 3 );
/// ```
#[derive(Clone,Debug,Dissolve,Eq,PartialEq)] // we don't include Ord or PartialOrd because the struct contains hashmaps, which don't implement those traits
pub struct GeneralizedMatchingMatrixWithSequentialOrder< ColumnIndex, RowIndex, Coefficient >
    where   ColumnIndex:        Clone + Hash + std::cmp::Eq,       // without cloning, it is difficult to write a method that returns the column_index at a given index (one gets "cannot move" errors)
            RowIndex:           Clone + Hash + std::cmp::Eq,       // without cloning, it is difficult to write a method that returns the row_index at a given index (one gets "cannot move" errors)            
{  
    bimap_col:  BijectiveSequence< ColumnIndex >,
    bimap_row:  BijectiveSequence< RowIndex >,
    vec_of_coefficients: Vec< Coefficient >,
}

impl < 'a, ColumnIndex, RowIndex, Coefficient > 

    GeneralizedMatchingMatrixWithSequentialOrder
    
    < ColumnIndex, RowIndex, Coefficient >
    where   ColumnIndex:    Clone + Hash + std::cmp::Eq,            // without cloning, it is difficult to write a method that returns the column_index at a given index (one gets "cannot move" errors)
            RowIndex:       Clone + Hash + std::cmp::Eq,            // without cloning, it is difficult to write a method that returns the row_index at a given index (one gets "cannot move" errors)            
    
{
    /// Create a new, empty matching array.
    /// 
    /// See [`GeneralizedMatchingMatrixWithSequentialOrder`] for usage.
    pub fn new() -> GeneralizedMatchingMatrixWithSequentialOrder< ColumnIndex, RowIndex, Coefficient > 
    {
        GeneralizedMatchingMatrixWithSequentialOrder{    
            bimap_row:  BijectiveSequence::new(),
            bimap_col:  BijectiveSequence::new(),
            vec_of_coefficients: Vec::new(),
        }        
    }

    /// Push a new structural nonzero entry to the array.    
    pub fn push( 
                    &mut self,
                    row_index:     RowIndex,                     
                    column_index:     ColumnIndex, 
                    coefficient:     Coefficient 
                )
            {
                self.bimap_col.push( column_index ); 
                self.bimap_row.push( row_index );               
                self.vec_of_coefficients.push( coefficient );
            }

    /// Reverse the order of the ordinals
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder;
    /// 
    /// let mut matching = GeneralizedMatchingMatrixWithSequentialOrder::new();
    /// matching.push(1, 1, 1.0);
    /// matching.push(2, 2, 2.0);
    /// assert_eq!( matching.bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order(), &vec![1, 2] ); // the vector sending ordinals to column indices
    /// assert_eq!( matching.bijection_row_indices_to_ordinals_and_inverse().vec_elements_in_order(), &vec![1, 2] ); // the vector sending ordinals to row indices
    /// assert_eq!( matching.structural_nonzero_values_in_sequence(), &vec![1.0, 2.0] ); // the vector sending ordinals to linear coefficients
    /// matching.reverse_order_of_matches();
    /// assert_eq!( matching.bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order(), &vec![2, 1] ); // the vector sending ordinals to column indices
    /// assert_eq!( matching.bijection_row_indices_to_ordinals_and_inverse().vec_elements_in_order(), &vec![2, 1] ); // the vector sending ordinals to row indices
    /// assert_eq!( matching.structural_nonzero_values_in_sequence(), &vec![2.0, 1.0] ); // the vector sending ordinals to linear coefficients
    /// ```
    pub fn reverse_order_of_matches( &mut self ) {
        self.bimap_col.reverse_ordinals();    
        self.bimap_row.reverse_ordinals();
        self.vec_of_coefficients.reverse();             
    }     



    /// Returns the transpose of `self`, leaving `self` unchanged
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder;
    /// 
    /// let mut matching = GeneralizedMatchingMatrixWithSequentialOrder::new();
    /// matching.push(0, 1, 25);
    /// matching.push(2, 3, 26);
    /// 
    /// 
    /// let mut transpose = GeneralizedMatchingMatrixWithSequentialOrder::new();
    /// transpose.push(1, 0, 25);
    /// transpose.push(3, 2, 26);
    /// 
    /// assert_eq!( transpose, matching.transpose_out_of_place()  );
    /// ```
    pub fn transpose_out_of_place( & self ) -> GeneralizedMatchingMatrixWithSequentialOrder< RowIndex, ColumnIndex, Coefficient > 
        where
            Coefficient:    Clone
    {
        GeneralizedMatchingMatrixWithSequentialOrder{
            bimap_row:              self.bimap_col.clone(),
            bimap_col:              self.bimap_row.clone(),
            vec_of_coefficients:    self.vec_of_coefficients.clone(),
        }
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
    pub fn generalized_inverse< RingOperator >( & self , ring_operator: RingOperator) 
        -> GeneralizedMatchingMatrixWithSequentialOrder< RowIndex, ColumnIndex, Coefficient > 
        where
            Coefficient:        Clone,
            RingOperator:       DivisionRingOperations< Element = Coefficient >        
        
        {

            let vec_of_coefficients    =   self.vec_of_coefficients
                                            .iter()
                                            .cloned()
                                            .map(|x| ring_operator.invert(x) )
                                            .collect::<Vec<_>>();

        GeneralizedMatchingMatrixWithSequentialOrder{
            bimap_row:          self.bimap_col.clone(),
            bimap_col:          self.bimap_row.clone(),
            vec_of_coefficients
        }
    }    

    
    /// Number of structural nonzero entries in the matching matrix
    /// 
    /// This can also be thought of as the number of matched pairs `(i,j)`, if we think of
    /// a generalized matchiing matrix as a weighted matching between the set of row indices
    /// and the set of column indices.
    pub fn number_of_structural_nonzeros(&self) -> usize { self.bimap_row.len() }

    /// Returns an iterator that runs over all `(&RowIndex, &ColumnIndex)` pairs; pairs occur in the usual `ord_to_val` order.
    pub fn support( &self ) -> std::iter::Zip<std::slice::Iter<'_, RowIndex>, std::slice::Iter<'_, ColumnIndex>> {
        let b = self.bimap_col.vec_elements_in_order().iter();                
        let a = self.bimap_row.vec_elements_in_order().iter();
        a.zip(b)
    }

    /// Takes an iterator `I` as input and returns an iterator `J` that skips over matched column indices
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder;
    /// 
    /// let mut matching = GeneralizedMatchingMatrixWithSequentialOrder::new();
    /// matching.push(0, 1, 5.0);
    /// matching.push(2, 3, 7.0);
    /// 
    /// let filtered = matching.filter_out_matched_column_indices( vec![0,1,2,3] );
    /// assert!(  filtered.eq( vec![0,2] )  );
    /// ```
    pub fn filter_out_matched_column_indices< I: IntoIterator< Item = ColumnIndex > >( &self, iter: I ) 
        -> 
        FilterOutMembers< I::IntoIter, & HashMap< ColumnIndex, usize > >
    {
        Filter::new( 
            iter.into_iter(), 
            LogicalNot::new(
                    MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper::new(
                        self.bimap_col.hashmap_element_to_ordinal()
                    )
                ) 
            )
    }    

    /// Takes an iterator `I` as input and returns an iterator `J` that skips over matched row indices
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder;
    /// 
    /// let mut matching = GeneralizedMatchingMatrixWithSequentialOrder::new();
    /// matching.push(0, 1, 5.0);
    /// matching.push(2, 3, 7.0);
    /// 
    /// let filtered = matching.filter_out_matched_row_indices( vec![0,1,2,3] );
    /// assert!(  filtered.eq( vec![1,3] )  );
    /// ```
    pub fn filter_out_matched_row_indices< I: IntoIterator< Item = RowIndex > >( &self, iter: I ) 
        -> 
        FilterOutMembers< I::IntoIter, & HashMap< RowIndex, usize > >
    {
        Filter::new( 
            iter.into_iter(), 
            LogicalNot::new(
                        MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper::new(
                                self.bimap_row.hashmap_element_to_ordinal()
                            )
                    ) 
                )
    }   


    /// Takes an iterator `I` as input and returns an iterator `J` that skips over **unmatched** column indices
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder;
    /// 
    /// let mut matching = GeneralizedMatchingMatrixWithSequentialOrder::new();
    /// matching.push(0, 1, 5.0);
    /// matching.push(2, 3, 7.0);
    /// 
    /// let filtered = matching.filter_out_unmatched_column_indices( vec![0,1,2,3] );
    /// assert!(  filtered.eq( vec![1,3] )  );
    /// ```
    pub fn filter_out_unmatched_column_indices< I: IntoIterator< Item = ColumnIndex > >( &self, iter: I ) 
        ->
        FilterOutNonmembers< I::IntoIter, & HashMap< ColumnIndex, usize > >
    {
        Filter::new( 
            iter.into_iter(), 
            MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper::new(
                        self.bimap_col.hashmap_element_to_ordinal()
                    )
                )
    }    

    /// Takes an iterator `I` as input and returns an iterator `J` that skips over **unmatched** row indices
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder;
    /// 
    /// let mut matching = GeneralizedMatchingMatrixWithSequentialOrder::new();
    /// matching.push(0, 1, 5.0);
    /// matching.push(2, 3, 7.0);
    /// 
    /// let filtered = matching.filter_only_matched_row_indices( vec![0,1,2,3] );
    /// assert!(  filtered.eq( vec![0,2] )  );
    /// ```
    pub fn filter_only_matched_row_indices< I: IntoIterator< Item = RowIndex > >( &self, iter: I ) 
        -> 
        FilterOutNonmembers< I::IntoIter, & HashMap< RowIndex, usize > >
    {
        Filter::new( 
            iter.into_iter(), 
                MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper::new(
                    self.bimap_row.hashmap_element_to_ordinal()
                ) 
            )
    }           
        
    /// Returns an immutable reference to a `BijectiveSequence< ColumnIndex >` object that stores the 
    /// bijection between column indices and [ordinals](crate::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder).
    /// This object stores both the map from column index to ordinal and its inverse.
    pub fn bijection_column_indices_to_ordinals_and_inverse( &self ) -> & BijectiveSequence< ColumnIndex > { & self.bimap_col }

    /// Returns an immutable reference to a `BijectiveSequence< RowIndex >` object that stores the 
    /// bijection between row indices and [ordinals](crate::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder).
    /// This object stores both the map from row index to ordinal and its inverse.
    pub fn bijection_row_indices_to_ordinals_and_inverse( &self ) -> & BijectiveSequence< RowIndex > { & self.bimap_row }    

    /// A struct that implements [EvaluateFunction], mapping matched column indices bijectively onto matched row indices.
    /// 
    /// Concretely, if the structural nonzero entries of the matrix are `(r0,c0,v0) .. (rN,cN,vN)`, then this object maps `ci` to `ri`.
    pub fn bijection_column_indices_to_row_indices( &self ) -> ComposeFunctions< ColumnIndex, usize, RowIndex, &HashMap<ColumnIndex, usize>, &Vec<RowIndex> > {
        ComposeFunctions::new( 
                self.bijection_column_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal(),
                self.bijection_row_indices_to_ordinals_and_inverse().vec_elements_in_order()
            )
    }

    /// A struct that implements [EvaluateFunction] matrix_to_factor matched row indices bijectively onto matched column indices.
    /// 
    /// Concretely, if the structural nonzero entries of the matrix are `(r0,c0,v0) .. (rN,cN,vN)`, then this object maps `ri` to `ci`.
    pub fn bijection_row_indices_to_column_indices( &self ) -> ComposeFunctions< RowIndex, usize, ColumnIndex, &HashMap<RowIndex, usize>, &Vec<ColumnIndex> > {
        ComposeFunctions::new( 
                self.bijection_row_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal(),
                self.bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order()
            )
    } 

    /// A struct that implements [EvaluateFunction] matrix_to_factor matched column indices bijectively onto matched row indices.
    /// 
    /// Concretely, if the structural nonzero entries of the matrix are `(r0,c0,v0) .. (rN,cN,vN)`, then this object maps `ci` to `Some(ri)`;
    /// every unmatched column index maps to `None`.
    pub fn bijection_column_indices_to_row_indices_opt( &self ) 
        -> ComposeFunctions< 
                ColumnIndex, 
                Option< usize >, 
                Option< RowIndex >, 
                &HashMap<ColumnIndex, usize>, 
                Map< &Vec<RowIndex> >,
            > 
    {
        ComposeFunctions::new( 
                self.bijection_column_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal(),
                Map{ mapping_rule: self.bijection_row_indices_to_ordinals_and_inverse().vec_elements_in_order() }
            )
    }

    /// A struct that implements [EvaluateFunction] matrix_to_factor matched row indices bijectively onto matched column indices.
    /// 
    /// Concretely, if the structural nonzero entries of the matrix are `(r0,c0,v0) .. (rN,cN,vN)`, then this object maps `ri` to `Some(ci)`.
    /// All unmatched row indices map to `None`.
    pub fn bijection_row_indices_to_column_indices_opt( &self ) 
        -> ComposeFunctions< 
                RowIndex, 
                Option< usize >, 
                Option< ColumnIndex >, 
                &HashMap<RowIndex, usize>, 
                Map< &Vec<ColumnIndex> >,
            > {
        ComposeFunctions::new( 
                self.bijection_row_indices_to_ordinals_and_inverse().hashmap_element_to_ordinal(),
                Map{ mapping_rule: self.bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order() }
            )
    }       

    /// A vector `v` containing the matched row indices
    /// 
    /// Concretely, if the ordered sequence of structural nonzero entries of the matrix is
    /// `(r0,c0,v0) .. (rN,cN,vN)`,
    /// then `v` is the sequence of row indices `r0, .., rN`.
    /// 
    /// # Caution
    /// 
    /// The sequence `r0, .., rN` is not necessarily sorted in
    /// increasing order. It only represents the order in which `(row,column,value)`
    /// triplets were inserted in the matrix.
    pub fn matched_row_indices_in_sequence( &self ) -> & Vec< RowIndex > { 
        self.bijection_row_indices_to_ordinals_and_inverse().vec_elements_in_order()
    } 

    /// A vector `v` containing the matched column indices
    /// 
    /// Concretely, if the ordered sequence of structural nonzero entries of the matrix is
    /// `(r0,c0,v0) .. (rN,cN,vN)`,
    /// then `v` is the sequence of column indices `c0, .., cN`.
    /// 
    /// # Caution
    /// 
    /// The sequence `c0, .., cN` is not necessarily sorted in
    /// increasing order. It only represents the order in which `(row,column,value)`
    /// triplets were inserted in the matrix.
    pub fn matched_column_indices_in_sequence( &self ) -> & Vec< ColumnIndex > { 
        self.bijection_column_indices_to_ordinals_and_inverse().vec_elements_in_order()
    }

    /// A vector `v` containing the structural nonzero values
    /// 
    /// Concretely, if the ordered sequence of structural nonzero entries of the matrix is
    /// `(r0,c0,v0) .. (rN,cN,vN)`,
    /// then `v` is the sequence of column indices `v0, .., vN`.
    pub fn structural_nonzero_values_in_sequence( &self ) -> & Vec< Coefficient > { 
        & self.vec_of_coefficients 
    }      


    // DEPRECATED AS UNSAFE    
    // /// Returns an immutable reference to the sequence of structural nonzero values associated with the array.
    // pub fn vec_structural_nonzeros_refmut( &mut self ) -> &mut Vec< Coefficient > { &mut self.vec_of_coefficients }      

    /// Returns the [ordinal](crate::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder) of the row index, or returns `None` if the row index is unmatched.
    pub fn ordinal_for_row_index( &self, row_index: &RowIndex ) -> Option< usize > { 
        self.bimap_row.ordinal_for_element( row_index )
    }

    /// Returns the [ordinal](crate::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder) of the column index, or returns `None` if the column index is unmatched.
    pub fn ordinal_for_column_index( &self, column_index: &ColumnIndex ) -> Option< usize > { 
        self.bimap_col.ordinal_for_element( column_index )
    }       

    /// Returns the row index corresponding to the [ordinal](crate::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder).
    pub fn row_index_for_ordinal( &self, ord: usize ) -> RowIndex { 
        self.bimap_row.element_for_ordinal( ord )
    }

    /// Returns the column index corresponding to the [ordinal](crate::algebra::matrices::types::matching::GeneralizedMatchingMatrixWithSequentialOrder).
    pub fn column_index_for_ordinal( &self, ord: usize ) -> ColumnIndex { 
        self.bimap_col.element_for_ordinal( ord )
    }           

    /// Returns the column index of the row index, or returns `None` if the row index is unmatched.
    pub fn column_index_for_row_index( &self, row_index: &RowIndex ) -> Option< ColumnIndex > { 
        self.ordinal_for_row_index( row_index ).map(|ord| self.column_index_for_ordinal( ord ))
    }

    /// Returns the row index of the column index, or returns `None` if the row index is unmatched.
    pub fn row_index_for_column_index( &self, column_index: &ColumnIndex ) -> Option< RowIndex > { 
        self.ordinal_for_column_index( column_index ).map(|ord| self.row_index_for_ordinal( ord ))
    }    
 
    /// The nonzero coefficient associated with this `row_index` index.
    pub fn coefficient_for_row_index( &self, row_index: & RowIndex ) -> Coefficient 
        where Coefficient: Clone    
    {
        self.vec_of_coefficients[ 
                self.ordinal_for_row_index( row_index ).unwrap()                
            ]            
            .clone()
    }

    /// The structural nonzero coefficient associated with this `row_index` (or `None`, if this row contains no structural nonzero entry)
    pub fn coefficient_opt_for_row_index( &self, row_index: & RowIndex ) -> Option< Coefficient >
        where Coefficient: Clone    
    {
        self.ordinal_for_row_index( row_index )
            .map( 
                | ord | 
                self.vec_of_coefficients[ ord ].clone()
            )
    }      

    /// The nonzero coefficient associated with this `column_index` index.
    pub fn coefficient_for_column_index( &self, column_index: & ColumnIndex ) -> Coefficient 
        where Coefficient: Clone    
    {
        self.vec_of_coefficients[ 
                self.ordinal_for_column_index( column_index ).unwrap() 
            ]
            .clone()
    }    


    /// The structural nonzero coefficient associated with this `column_index` (or `None`, if this column contains no structural nonzero entry)
    pub fn coefficient_opt_for_column_index( &self, column_index: & ColumnIndex ) -> Option< Coefficient >
        where Coefficient: Clone    
    {
        self.ordinal_for_column_index( column_index )
            .map( 
                | ord | 
                self.vec_of_coefficients[ ord ].clone()
            )
    }          


    /// The nonzero coefficient associated with the `ord`th row index.
    pub fn coefficient_for_ordinal( &self, ord: usize ) -> Coefficient 
        where Coefficient: Clone
    {
        self.vec_of_coefficients[ ord ]
            .clone()
    }    

    /// Returns true iff the given column index is matched.
    pub fn has_a_match_for_column_index(  &self, column_index: &ColumnIndex ) -> bool {
        self.bimap_col.has_ordinal_for_element( column_index )
    }

    /// Returns true iff the given row index is matched.
    pub fn has_a_match_for_row_index(  &self, row_index: &RowIndex ) -> bool {
        self.bimap_row.has_ordinal_for_element( row_index )
    }   
    
    /// Returns true iff the given column index is unmatched.
    pub fn lacks_a_match_for_column_index(  &self, column_index: &ColumnIndex ) -> bool {
        ! self.bimap_col.has_ordinal_for_element( column_index )
    }
 
    /// Returns true iff the given row index is unmatched.
    pub fn lacks_a_match_for_row_index(  &self, row_index: &RowIndex ) -> bool {
        ! self.bimap_row.has_ordinal_for_element( row_index )
    }      


    /// Iterates over matched `(row_index, column_index)` pairs, in order.
    /// 
    /// Concretely, if the ordered sequence of structural nonzero entries of the matrix is
    /// `(r0,c0,v0) .. (rN,cN,vN)`, then this iterator will return the pairs
    /// `(r0,c0), (r1,c1), .., (rN,cN)`, in that order.
    /// 
    /// # Caution
    /// 
    /// The sequence `r0, .., rN` is not necessarily sorted in
    /// increasing order. It only represents the order in which `(row,column,value)`
    /// triplets were inserted in the matrix. The same holds for `c0, .., cN`.
    pub fn iter_index_pairs( & self ) 
            ->  
            Zip< std::slice::Iter< RowIndex >, std::slice::Iter< ColumnIndex >  > 
    {
        let row_iterator    =   self.bimap_row.vec_elements_in_order().iter();
        let col_iterator    =   self.bimap_col.vec_elements_in_order().iter();
        return row_iterator.zip( col_iterator )
    }

    /// Iterates over `((row_index, column_index), coefficient)` tuples.
    /// 
    /// The order in which entries appear is consistent with the ordering of entries used everywhere else with
    /// this data structure.
    pub fn iter_entries( &self ) -> 
        Zip<
            Zip< std::slice::Iter< RowIndex >, std::slice::Iter< ColumnIndex >  > ,
            std::slice::Iter< Coefficient >,
        >
    {
        let row_iterator    =   self.bimap_row.vec_elements_in_order().iter();
        let col_iterator    =   self.bimap_col.vec_elements_in_order().iter();
        let coefficient_iterator   =   self.vec_of_coefficients.iter();
        return row_iterator.zip( col_iterator ).zip( coefficient_iterator )
    }

    // THE FOLLOWING METHODS ARE DEPRECATED; HOWEVER WE MAY REVIVE THEM IF WE EVER NEED AN ORDINAL MATRCHING ARRAY WITH ADDITIONAL STRUCTURE        

    // /// Returns the ordinal of the column index associated with the given row index.
    // /// Panicks if the column index is not associated with a row index.    
    // pub fn row_index_to_ordinalmin( &self, row_index: &RowIndex ) -> usize { 
    //     self.ord_to_ordmin[
    //             * self.ordinal_for_row_index.get( row_index ).unwrap()
    //         ]
    // }

    // /// Returns the ordinal of the row index associated with the given column index.
    // /// Panicks if the row index is not associated with a column index.
    // pub fn ordinal_for_column_index( &self, column_index: &ColumnIndex ) -> usize { 
    //     self.ordmin_to_ord[
    //             * self.column_index_to_ordinalmin.get( column_index ).unwrap()
    //         ]
    // }    

    // /// Returns either `Some(k)`, where `k` is the ordinal of the row index associated with the given column index,
    // /// or `None`, if there exists no associated row index.
    // pub fn column_index_to_ordinal_opt( &self, column_index: &ColumnIndex ) -> Option< usize > { 
    //     match self.column_index_to_ordinalmin.get( column_index ) {
    //         Some( & ordmin ) => { Some( self.ordmin_to_ord[ ordmin ].clone() ) },
    //         None => None
    //     }
    // } 
    
    // /// Returns either `Some(k)`, where `k` is the ordinal of the row index associated with the given column index,
    // /// or `None`, if there exists no associated row index.
    // pub fn row_index_to_ordinalmin_opt( &self, row_index: &RowIndex ) -> Option< usize > { 
    //     match self.ordinal_for_row_index.get( row_index ) {
    //         Some( & ord ) => { Some( self.ord_to_ordmin[ ord ].clone() ) },
    //         None => None
    //     }
    // }    

    // pub fn column_index_for_row_index( &self, row_index: & RowIndex ) -> & ColumnIndex { 
    //     & self.ordmin_to_column_index[ 
    //             self.row_index_to_ordinalmin( row_index )
    //         ]        
    // }

    // pub fn row_index_for_column_index( &self, column_index: & ColumnIndex ) -> & RowIndex { 
    //     & self.row_index_for_ordinal[ 
    //             self.ordinal_for_column_index( column_index )
    //         ]
    // }    

    // /// The nonzero coefficient associated with this `row_index` index.
    // pub fn coefficient_for_row_index( &self, row_index: & RowIndex ) -> & Coefficient {
    //     & self.coefficient_for_ordinal[ * self.ordinal_for_row_index.get( row_index ).unwrap() ]
    // }

    // /// The nonzero coefficient associated with this `column_index` index.
    // pub fn coefficient_for_column_index( &self, column_index: & ColumnIndex ) -> & Coefficient {
    //     & self.coefficient_for_ordinal[ self.ordinal_for_column_index( column_index ) ]
    // }  


}


                          


impl < 'a, KeyBoth, Coefficient > 

    GeneralizedMatchingMatrixWithSequentialOrder
        < KeyBoth, KeyBoth, Coefficient > where   
    
    KeyBoth:     Clone + Hash + std::cmp::Eq, // without cloning, it is difficult to write a method that returns the column_index at a given index (one gets "cannot move" errors)
{
    /// Returns `true` iff `key` is a matched row or column index.
    /// 
    /// This method is only available if row and column indices have the same type.
    pub fn has_match_for_index( &self, key: &KeyBoth ) -> bool { self.has_a_match_for_column_index(key) || self.has_a_match_for_row_index(key) }

    /// Returns `true` iff `key` is an unmatched row index **and** an unmatched column index.
    /// 
    /// This method is only available if row and column indices have the same type.
    pub fn lacks_match_for_index( &self, key: &KeyBoth ) -> bool { self.lacks_a_match_for_column_index(key) && self.lacks_a_match_for_row_index(key) } 
    
    /// A binned count of unmatched items.
    /// 
    /// Concretely, this means vector `v` such that `v[i]` is the number of umatched items 
    /// in the iterator with rank `i` (where rank is determined by `rank_fn`).
    /// 
    /// This method is only available if row and column indices have the same type.
    pub fn histogram_of_unmatched_indices< I, F >( &self, iter_keyboth: I, mut rank_fn: F ) 
            -> 
            Vec< usize > 
        where
            I: Iterator<Item=KeyBoth>, 
            F: FnMut(KeyBoth)->usize        
            
        {
        let mut histo = Vec::new();
        let mut rank;
        for key in iter_keyboth { 
            if self.lacks_match_for_index( &key ) { 
                rank = rank_fn(key);
                while histo.len() < rank + 1 { histo.push(0) }
                histo[rank] +=1
            } 
        }
        histo
    }

}


//  EXTEND TRAIT IMPLEMENTATION
//  ------------------------------------------------------------------------------

impl    < 'a, ColumnIndex, RowIndex, Coefficient > 
        Extend< ( RowIndex, ColumnIndex, Coefficient ) > 
        for 
        GeneralizedMatchingMatrixWithSequentialOrder < ColumnIndex, RowIndex, Coefficient > 
    where   ColumnIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{
    /// Add a nonzero entry to a matching array.
    fn extend< I: IntoIterator<Item= ( RowIndex, ColumnIndex, Coefficient ) >>( &mut self, iter: I) {
        for x in iter.into_iter() { self.push( x.0, x.1, x.2 ) }
    }
}


//  FROM_ITER TRAIT IMPLEMENTATION
//  ------------------------------------------------------------------------------

impl    < 'a, ColumnIndex, RowIndex, Coefficient > 
        
        FromIterator
        < ( RowIndex, ColumnIndex, Coefficient ) > for 
        
        GeneralizedMatchingMatrixWithSequentialOrder 
        < ColumnIndex, RowIndex, Coefficient > 

    where   ColumnIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{
    fn from_iter< I: IntoIterator<Item= ( RowIndex, ColumnIndex, Coefficient ) >>(iter: I) -> Self {

        let mut matching_array  =   GeneralizedMatchingMatrixWithSequentialOrder::new();

        for x in iter.into_iter() { matching_array.push( x.0, x.1, x.2 ) }
        
        matching_array
    }
}


//  ORACLE TRAIT IMPLEMENTATIONS
//  ------------------------------------------------------------------------------






impl < 'a, ColumnIndex, RowIndex, Coefficient >

    MatrixOracle for 

    GeneralizedMatchingMatrixWithSequentialOrder
        < ColumnIndex, RowIndex, Coefficient >

    where   ColumnIndex:        Clone + Debug + Eq + Hash,
            RowIndex:           Clone + Debug + Eq + Hash,
            Coefficient:        Clone + Debug + PartialEq,        
{
    type Coefficient            =   Coefficient;

    type RowIndex               =   RowIndex;

    type ColumnIndex            =   ColumnIndex;

    type RowEntry               =   ( ColumnIndex, Coefficient );

    type ColumnEntry            =   ( RowIndex, Coefficient );

    type Row                    =   TwoTypeIterator<
                                        Empty< ( ColumnIndex, Coefficient ) >,
                                        std::iter::Once< ( ColumnIndex, Coefficient ) >,
                                    > ;

    type RowReverse             =   Self::Row;

    type Column                 =   TwoTypeIterator<
                                        Empty< ( RowIndex, Coefficient ) >,
                                        std::iter::Once< ( RowIndex, Coefficient ) >,
                                    > ;

    type ColumnReverse          =   Self::Column;

    fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
        match self.ordinal_for_row_index( &index ) {
            None => {
                TwoTypeIterator::Version1( std::iter::empty() )
            }
            Some( ord ) => {
                TwoTypeIterator::Version2( 
                        std::iter::once(        (   
                                                    self.column_index_for_ordinal( ord ),
                                                    self.coefficient_for_ordinal( ord ),  
                                                )          
                            )
                    )
            }
        }
    }

    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
        self.row( index )
    }

    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
        match self.ordinal_for_column_index( &index ) {
            None => {
                TwoTypeIterator::Version1( std::iter::empty() )
            }
            Some( ord ) => {
                TwoTypeIterator::Version2( 
                        std::iter::once(        (   
                                                    self.row_index_for_ordinal( ord ),
                                                    self.coefficient_for_ordinal( ord ),  
                                                )          
                            )
                    )
            }
        }
    }

    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
        self.column( index )
    }

    fn has_row_for_index(     &   self, _index: & Self::RowIndex   )   -> bool {
        true
    }

    fn has_column_for_index(  &   self, _index: & Self::ColumnIndex)   -> bool {
        true
    }

    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        match self.ordinal_for_row_index( row ) {
            None => { None }
            Some( ord ) => {
                let matched_column_index    =   self.column_index_for_ordinal( ord );

                match matched_column_index == * column {
                    true    =>  Some( self.coefficient_for_ordinal( ord ) ),
                    false   =>  None
                }
            }
        }        
    }
}









//  MATRIX ORACLE OPERATIONS
//  ---------------------------------------------------------------------



impl < 'a, ColumnIndex, RowIndex, Coefficient >

    MatrixOracleOperations for 

    GeneralizedMatchingMatrixWithSequentialOrder
        < ColumnIndex, RowIndex, Coefficient >

    where   ColumnIndex:        Clone + Hash + std::cmp::Eq,       // required by type definition
            RowIndex:           Clone + Hash + std::cmp::Eq,       // required by type definition
{}   








//  =========================================================================================================
//  TESTS
//  =========================================================================================================


//  ---------------------------------------------------------------------
//  Docstring tests (draft)
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod docstring_tests_draft {
    
    #[test]
    fn test_matching_array() {
        use super::GeneralizedMatchingMatrixWithSequentialOrder;

        let mut matching_array  =   GeneralizedMatchingMatrixWithSequentialOrder::new();
        matching_array.push( 1, 0, 3 );
        matching_array.push( 2, 1, 4 );
        matching_array.push( 0, 2, 5 );    
        
        assert_eq!( matching_array.structural_nonzero_values_in_sequence(), &vec![ 3, 4, 5 ] );

        assert_eq!( matching_array.row_index_for_column_index( &0 ), Some(1) );
        assert_eq!( matching_array.column_index_for_row_index( &0 ), Some(2) );     

        assert_eq!( matching_array.ordinal_for_column_index( &1 ), Some(1) );
        assert_eq!( matching_array.ordinal_for_row_index( &1 ), Some(0) ); 

        assert_eq!( matching_array.column_index_for_ordinal( 0 ), 0 );
        assert_eq!( matching_array.row_index_for_ordinal( 0 ), 1 );     

        assert_eq!( matching_array.coefficient_for_column_index( &0 ), 3 );
        assert_eq!( matching_array.coefficient_for_row_index( &0 ), 5 );       
        
        assert_eq!( matching_array.coefficient_for_ordinal( 0 ), 3 );
    }    

}


























