//! This is a deprecated version of the gimbled U-match factorization.
//! 
//! It wraps its COMB's in `MatrixAlgebraPacket` to remove the double-reversal from order-operators.
//! 
//! We've opted out of this design choice for now, but we're keeping a copy of the code here for reference, in case it becomes important to reverse this choice.


use ndarray::Order;

use crate::{algebra::{
    matrices::{
        operations::umatch::row_major::{comb::{SourceComb, SourceCombInverse, TargetComb, TargetCombInverse}, Umatch}, 
        query::{MatrixAlgebra, MatrixOracle}, 
        types::{matching::GeneralizedMatchingMatrixWithSequentialOrder, packet::MatrixAlgebraPacket, transpose::OrderAntiTranspose, two_type::TwoTypeMatrix}
    }, 
    rings::traits::DivisionRingOperations, 
    vectors::entries::KeyValPair
}, utilities::order::{JudgeOrder, OrderOperatorByKeyCustom, ReverseOrder}};



use std::hash::Hash;













#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GimbledUmatch< MatrixToFactor >
    where   
        MatrixToFactor:                     MatrixAlgebra,
        MatrixToFactor::ColumnIndex:        Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
{
    Standard( 
        Umatch< MatrixToFactor > 
    ),
    AntiTranspose( 
        (
            Umatch< 
                OrderAntiTranspose<
                    MatrixToFactor 
                >
            >,
            GeneralizedMatchingMatrixWithSequentialOrder< // a copy of the un-antitransposed generalized matching matrix
                MatrixToFactor::ColumnIndex,
                MatrixToFactor::RowIndex,
                MatrixToFactor::Coefficient 
            >, 
        ),
    ),
}






impl < MatrixToFactor >

    GimbledUmatch
        < MatrixToFactor >

    // these are required for the underlying `Umatch` struct to implement the COMB lookup methods
    where   
        MatrixToFactor:                             MatrixAlgebra,    
        MatrixToFactor::RingOperator:               DivisionRingOperations,        
        MatrixToFactor::ColumnIndex:                Hash, 
        MatrixToFactor::RowIndex:                   Hash, 
        MatrixToFactor::RowEntry:                   KeyValPair,
        MatrixToFactor::ColumnEntry:                KeyValPair,        
{

    /// Returns `true` if the U-match is in column-major form, `false` otherwise.
    /// 
    /// "Column-major" means that the inner data structure is a [Umatch] decomposition of the antitranspose of the matrix to factor.
    /// "Row-major" means that the inner data structure is a [Umatch] decomposition of the matrix to factor (no antitranspose).
    /// 
    /// We typically expect the column-major to be faster for looking up columns of the differential COMB,
    /// and row-major to be faster for looking up rows of the inverse differential COMB.
    pub fn is_column_major( &self ) -> bool {
        match self {
            GimbledUmatch::Standard( _ ) => false,
            GimbledUmatch::AntiTranspose(_) => true,
        }
    }

    /// Returns a reference to the matrix we wish to factor.
    pub fn matrix_to_factor_ref( &self ) -> & MatrixToFactor {
        match self {
            GimbledUmatch::Standard( umatch ) => {
                umatch.matrix_to_factor_ref()
            },
            GimbledUmatch::AntiTranspose( (umatch, _) ) => {
                umatch
                    .matrix_to_factor_ref()
                    .matrix_to_antitranspose()
            },
        }
    }

    /// Returns a reference to the generalized matching matrix of the U-match.
    pub fn generalized_matching_matrix_ref( &self ) -> 
        & GeneralizedMatchingMatrixWithSequentialOrder<
            MatrixToFactor::ColumnIndex,
            MatrixToFactor::RowIndex,
            MatrixToFactor::Coefficient
        > 
    {
        match self {
            GimbledUmatch::Standard( umatch ) => {
                umatch.generalized_matching_matrix_ref()
            },
            GimbledUmatch::AntiTranspose( (_, gmm) ) => {
                & gmm
            },
        }
    }


    /// Returns a reference to the matching array of the internally stored  U-match factorization,
    /// wrapped in a convenient convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket)
    pub fn generalized_matching_matrix_ref_packet( &self ) -> MatrixAlgebraPacket< 
            & GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient >,
            MatrixToFactor::RingOperator, 
            OrderOperatorByKeyCustom < MatrixToFactor::OrderOperatorForColumnIndices >, // order operator for row entries
            MatrixToFactor::OrderOperatorForRowIndices, // order operator for column indices            
            OrderOperatorByKeyCustom< MatrixToFactor::OrderOperatorForRowIndices >, // order operator for column entries            
            MatrixToFactor::OrderOperatorForColumnIndices, // order operator for column indices
        >
    {
        MatrixAlgebraPacket{ 
            matrix: self.generalized_matching_matrix_ref(), 
            ring_operator: self.ring_operator(), 
            order_operator_for_row_entries:     OrderOperatorByKeyCustom::< MatrixToFactor::OrderOperatorForColumnIndices >::new(  // note: we have to use this instead of `matrix_to_factor_ref().order_operator_for_row_entries()` because the order operator for row entries is specific to the type of row entries in the matrix
                                                    self.matrix_to_factor_ref().order_operator_for_column_indices() 
                                                ),
            order_operator_for_row_indices:     self.matrix_to_factor_ref().order_operator_for_row_indices(),
            order_operator_for_column_entries:  OrderOperatorByKeyCustom::< MatrixToFactor::OrderOperatorForRowIndices >::new(  // note: we have to use this instead of `matrix_to_factor_ref().order_operator_for_column_entries()` because the order operator for column entries is specific to the type of row entries in the matrix
                                                    self.matrix_to_factor_ref().order_operator_for_row_indices() 
                                                ),            
            order_operator_for_column_indices:  self.matrix_to_factor_ref().order_operator_for_column_indices(),
        }
    }     
    //     -> MatrixAlgebraPacket< 
    //         & GeneralizedMatchingMatrixWithSequentialOrder< MatrixToFactor::ColumnIndex, MatrixToFactor::RowIndex, MatrixToFactor::Coefficient >,
    //         MatrixToFactor::RingOperator, 
    //         MatrixToFactor::OrderOperatorForRowEntries, // order operator for row entries
    //         MatrixToFactor::OrderOperatorForRowIndices, // order operator for column indices            
    //         MatrixToFactor::OrderOperatorForColumnEntries, // order operator for column entries            
    //         MatrixToFactor::OrderOperatorForColumnIndices, // order operator for column indices
    //     >
    // {
    //     MatrixAlgebraPacket{ 
    //         matrix: self.generalized_matching_matrix_ref(), 
    //         ring_operator: self.ring_operator(), 
    //         order_operator_for_row_entries:     self.matrix_to_factor_ref().order_operator_for_row_entries() ,
    //         order_operator_for_row_indices:     self.matrix_to_factor_ref().order_operator_for_row_indices(),
    //         order_operator_for_column_entries:  self.matrix_to_factor_ref().order_operator_for_column_entries(),      
    //         order_operator_for_column_indices:  self.matrix_to_factor_ref().order_operator_for_column_indices(),
    //     }
    // }     


    /// Returns the source comb of the U-match.
    pub fn source_comb< 'a >( &'a self ) 
        -> TwoTypeMatrix<
            SourceComb<'a, MatrixToFactor >,
            MatrixAlgebraPacket<
                OrderAntiTranspose<
                    TargetCombInverse<
                        'a, 
                        OrderAntiTranspose< MatrixToFactor>
                    >
                >,
                MatrixToFactor::RingOperator,
                MatrixToFactor::OrderOperatorForRowEntries,
                MatrixToFactor::OrderOperatorForColumnIndices,
                MatrixToFactor::OrderOperatorForRowEntries,
                MatrixToFactor::OrderOperatorForColumnIndices
            >
        > {
        match self {
            GimbledUmatch::Standard( umatch ) => {
                TwoTypeMatrix::Version1(
                    umatch.source_comb(),
                )
            },
            GimbledUmatch::AntiTranspose( (umatch, _) ) => {
                TwoTypeMatrix::Version2(
                    MatrixAlgebraPacket{                    // WE USE THIS PACKET TO ENSURE THAT ORDER OPERATORS HAVE THE CORRECT TYPE
                        matrix:                             OrderAntiTranspose::new( umatch.target_comb_inverse() ),    
                        ring_operator:                      self.ring_operator(),
                        order_operator_for_row_entries:     self.order_operator_for_row_entries(),
                        order_operator_for_row_indices:     self.order_operator_for_column_indices(),
                        order_operator_for_column_entries:  self.order_operator_for_row_entries(),
                        order_operator_for_column_indices:  self.order_operator_for_column_indices(),
                    }                    
                )
            },
        }
    }


    /// Returns the inverse source comb of the U-match.
    pub fn source_comb_inverse< 'a >( &'a self ) 
        -> TwoTypeMatrix<
            SourceCombInverse<'a, MatrixToFactor >,
            MatrixAlgebraPacket<
                OrderAntiTranspose<
                    TargetComb<
                        'a, 
                        OrderAntiTranspose< MatrixToFactor>
                    >
                >,
                MatrixToFactor::RingOperator,
                MatrixToFactor::OrderOperatorForRowEntries,
                MatrixToFactor::OrderOperatorForColumnIndices,
                MatrixToFactor::OrderOperatorForRowEntries,
                MatrixToFactor::OrderOperatorForColumnIndices,
            > 
        >
    {
        match self {
            GimbledUmatch::Standard( umatch ) => {
                TwoTypeMatrix::Version1(
                    umatch.source_comb_inverse(),
                )
            },
            GimbledUmatch::AntiTranspose( (umatch, _) ) => { 
                TwoTypeMatrix::Version2(
                    MatrixAlgebraPacket{                    // WE USE THIS PACKET TO ENSURE THAT ORDER OPERATORS HAVE THE CORRECT TYPE
                        matrix:                             OrderAntiTranspose::new( umatch.target_comb() ),    
                        ring_operator:                      self.ring_operator(),
                        order_operator_for_row_entries:     self.order_operator_for_row_entries(),
                        order_operator_for_row_indices:     self.order_operator_for_column_indices(),
                        order_operator_for_column_entries:  self.order_operator_for_row_entries(),
                        order_operator_for_column_indices:  self.order_operator_for_column_indices(),
                    }                    
                )
            },
        }
    }    


    /// Returns the target comb of the U-match.
    pub fn target_comb< 'a >( &'a self ) -> TwoTypeMatrix<
        TargetComb<'a, MatrixToFactor>,
        MatrixAlgebraPacket<
            OrderAntiTranspose<
                SourceCombInverse<
                    'a, 
                    OrderAntiTranspose< MatrixToFactor>
                >
            >,
            MatrixToFactor::RingOperator,
            MatrixToFactor::OrderOperatorForColumnEntries,
            MatrixToFactor::OrderOperatorForRowIndices,
            MatrixToFactor::OrderOperatorForColumnEntries,
            MatrixToFactor::OrderOperatorForRowIndices,
        > 
    > {
        match self {
            GimbledUmatch::Standard( umatch ) => {
                TwoTypeMatrix::Version1(
                    umatch.target_comb(),
                )
            },
            GimbledUmatch::AntiTranspose( (umatch, _) ) => { 
                TwoTypeMatrix::Version2(
                    MatrixAlgebraPacket{                    // WE USE THIS PACKET TO ENSURE THAT ORDER OPERATORS HAVE THE CORRECT TYPE
                        matrix:                             OrderAntiTranspose::new( umatch.source_comb_inverse() ),    
                        ring_operator:                      self.ring_operator(),
                        order_operator_for_row_entries:     self.order_operator_for_column_entries(),
                        order_operator_for_row_indices:     self.order_operator_for_row_indices(),
                        order_operator_for_column_entries:  self.order_operator_for_column_entries(),
                        order_operator_for_column_indices:  self.order_operator_for_row_indices(),
                    }                    
                )
            },
        }
    }



    /// Returns the target comb of the U-match.
    pub fn target_comb_inverse< 'a >( &'a self ) 
        -> TwoTypeMatrix<
            TargetCombInverse<'a, MatrixToFactor>,
            MatrixAlgebraPacket<
                OrderAntiTranspose<
                    SourceComb<
                        'a, 
                        OrderAntiTranspose< MatrixToFactor>
                    >
                >,
                MatrixToFactor::RingOperator,
                MatrixToFactor::OrderOperatorForColumnEntries,
                MatrixToFactor::OrderOperatorForRowIndices,
                MatrixToFactor::OrderOperatorForColumnEntries,
                MatrixToFactor::OrderOperatorForRowIndices,
            > 
        > 
    {
        match self {
            GimbledUmatch::Standard( umatch ) => {
                TwoTypeMatrix::Version1(
                    umatch.target_comb_inverse(),
                )
            },
            GimbledUmatch::AntiTranspose( (umatch, _) ) => { 
                TwoTypeMatrix::Version2(
                    MatrixAlgebraPacket{                    // WE USE THIS PACKET TO ENSURE THAT ORDER OPERATORS HAVE THE CORRECT TYPE
                        matrix:                             OrderAntiTranspose::new( umatch.source_comb() ),    
                        ring_operator:                      self.ring_operator(),
                        order_operator_for_row_entries:     self.order_operator_for_column_entries(),
                        order_operator_for_row_indices:     self.order_operator_for_row_indices(),
                        order_operator_for_column_entries:  self.order_operator_for_column_entries(),
                        order_operator_for_column_indices:  self.order_operator_for_row_indices(),
                    }                    
                )
            },
        }
    }    

    /// Returns the ring operator for the matrix to factor.
    pub fn ring_operator( &self ) -> MatrixToFactor::RingOperator {
        self.matrix_to_factor_ref().ring_operator()
    }


    /// Rank of the factored matrix
    /// 
    /// Equivalently, 
    /// - the dimension of the image of the linear map represented by the matrix
    /// - the number of nonzero entries in the generalized matching matrix of the U-match factorization
    pub fn rank( &self ) -> usize        
    {
        self.generalized_matching_matrix_ref().number_of_structural_nonzeros()
    }  


    /// Returns a copy of the order comparator for (column-index, coefficient) pairs
    pub fn order_operator_for_row_entries( &self ) -> MatrixToFactor::OrderOperatorForRowEntries 
    { self.matrix_to_factor_ref().order_operator_for_row_entries() } 
    
    /// Returns a copy of the **inverted** order comparator for (column-index, coefficient) pairs
    pub fn order_operator_for_row_entries_reverse( &self ) -> ReverseOrder< MatrixToFactor::OrderOperatorForRowEntries >
    { ReverseOrder::new(self.matrix_to_factor_ref().order_operator_for_row_entries()) }    

    /// Returns a copy of the order comparator for row indices
    pub fn order_operator_for_row_indices( &self ) -> MatrixToFactor::OrderOperatorForRowIndices
    { self.matrix_to_factor_ref().order_operator_for_row_indices() } 

    /// Returns a copy of the **inverted** order comparator for row indices
    pub fn order_operator_for_row_indices_reverse( &self ) -> ReverseOrder< MatrixToFactor::OrderOperatorForRowIndices >
    { ReverseOrder::new(self.matrix_to_factor_ref().order_operator_for_row_indices()) }                

    /// Returns a copy of the order comparator for (row-index, coefficient) pairs
    pub fn order_operator_for_column_entries( &self ) -> MatrixToFactor::OrderOperatorForColumnEntries 
    { self.matrix_to_factor_ref().order_operator_for_column_entries() }  

    /// Returns a copy of the **inverted** order comparator for (row-index, coefficient) pairs
    pub fn order_operator_for_column_entries_reverse( &self ) -> ReverseOrder< MatrixToFactor::OrderOperatorForColumnEntries >
    { ReverseOrder::new( self.matrix_to_factor_ref().order_operator_for_column_entries() ) } 

    /// Returns a copy of the order comparator for (row-index, coefficient) pairs
    pub fn order_operator_for_column_indices( &self ) -> MatrixToFactor::OrderOperatorForColumnIndices 
    { self.matrix_to_factor_ref().order_operator_for_column_indices() }  

    /// Returns a copy of the **inverted** order comparator for (row-index, coefficient) pairs
    pub fn order_operator_for_column_indices_reverse( &self ) -> ReverseOrder< MatrixToFactor::OrderOperatorForColumnIndices >
    { ReverseOrder::new( self.matrix_to_factor_ref().order_operator_for_column_indices() ) }      

    
    /// The sequence of matched row indices in *ascending order*
    /// 
    /// Concretely, this is the sequence of matched row indices `r_0 < .. < r_k`, where
    /// order is deteremined by the order operator for row indices associated with the factored matrix.
    pub fn matched_row_indices_in_ascending_order( &self ) -> &Vec< MatrixToFactor::RowIndex > {
        self.generalized_matching_matrix_ref().matched_row_indices_in_sequence()
    }

    /// The sequence of matched column indices, ordered according to the associated row indices
    /// 
    /// Concretely, this is the sequence of matched column indices `c_0, .., c_k`, obtained by 
    /// ordering the sequence of matched row-column index pairs `(r0,c0), .., (rk,ck)` 
    /// such that `r_0 < .. < r_k`.
    /// 
    /// **In particular, there is no guarantee that `c_0 < .. < c_k`**.
    pub fn matched_column_indices_in_matched_row_order( &self ) -> &Vec< MatrixToFactor::ColumnIndex > {
        self.generalized_matching_matrix_ref().matched_column_indices_in_sequence()
    }  

    /// The sequence of matched column indices in *ascending order*
    /// 
    /// Concretely, this is the sequence of matched column indices `c_0 < .. < c_k`, where
    /// order is deteremined by the order operator for column indices associated with the factored matrix.
    /// 
    /// # Performance
    /// 
    /// The U-match data structure stores matched column indices in a different order. Thus to obtain this
    /// sequence, we must copy the stored data, and sort the column indices according to the order operator.
    /// If all you need is the sequence of matched column indices, use [Umatch::matched_column_indices_in_matched_row_order] instead.
    pub fn matched_column_indices_in_ascending_order( &self ) -> Vec< MatrixToFactor::ColumnIndex > {
        let mut indices = self.generalized_matching_matrix_ref().matched_column_indices_in_sequence().clone();
        let order_operator = self.matrix_to_factor_ref().order_operator_for_column_indices();
        indices.sort_by( |a,b| order_operator.judge_cmp( a, b ) );
        indices
    }      
    
}






























#[cfg(test)]
mod test {
    use itertools::Itertools;

    use super::*;
    use crate::algebra::{matrices::{debug::{matrix_oracle_is_internally_consistent, matrix_order_operators_are_internally_consistent, product_is_identity_matrix}, types::{packet::MatrixAlgebraPacket, product::ProductMatrix, vec_of_vec::sorted::VecOfVec}}, rings::types::field_prime_order::PrimeOrderField};












    /// Checks that Umatch decomposition is correct (using a random example matrix, D) in the following sense:
    /// T^{-1} * T = I
    /// S^{-1} * S = I
    /// T^{-1} * D * S = M   
    /// And the rows of T, T^{-1}, S, and S^{-1} appear in strictly ascending order 
    fn test_umatchrowmajor_comprehensive_overall() {   

        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
        use crate::algebra::matrices::types::product::ProductMatrix;
        use crate::algebra::matrices::query::MatrixOracle;               

        let num_indices_row             =   10;
        let num_indices_col             =   20;
        let approximate_density           =   0.2;
        let modulus                     =   17;
        let allow_nonstructural_zero     =   true;

        let ring_operator           =   PrimeOrderField::new( modulus );
        let matrix_to_factor_data       =   VecOfVec::random_mod_p_with_density( num_indices_row, num_indices_col, approximate_density, modulus, allow_nonstructural_zero );
        let matrix_to_factor = MatrixAlgebraPacket::with_default_order( & matrix_to_factor_data, ring_operator );
     

        let umatch 
            =   Umatch::new( 
                    matrix_to_factor, 
                    (0..num_indices_row).rev(), 
                );




        let matching = umatch.generalized_matching_matrix_ref();
  

        let comb_target = umatch.target_comb();
        let comb_target_inv = umatch.target_comb_inverse();        
        let comb_source = umatch.source_comb();        
        let comb_source_inv = umatch.source_comb_inverse();  

        let comb_target_ref         =   & comb_target;
        let comb_target_inv_ref         =   & comb_target_inv;
        let comb_source_ref         =   & comb_source;
        let comb_source_inv_ref         =   & comb_source_inv;            
        

        let product_source = ProductMatrix::new( comb_source_ref, comb_source_inv_ref );
        let product_target = ProductMatrix::new( comb_target_ref, comb_target_inv_ref );        
        let product_target_comb_inv_times_matrix_to_factor = ProductMatrix::new( comb_target_inv_ref, matrix_to_factor );      
        let product_target_comb_inv_times_matrix_to_factor_times_source_comb = ProductMatrix::new( product_target_comb_inv_times_matrix_to_factor, comb_source_ref );                        


      

        // ----------------------------------------------------------------------------------------------------------------



        // println!("matrix_to_factor:");
        // print_indexed_rows( & matrix_to_factor, 0 .. num_indices_row );
        // println!("matching:");
        // print_indexed_rows( & matching, 0 .. num_indices_row );        
        // println!("comb_source:");
        // print_indexed_rows( & comb_source, 0 .. num_indices_col );        
        // println!("comb_source_inv:");
        // print_indexed_rows( & comb_source_inv, 0 .. num_indices_col );     
        // println!("comb_target:");
        // print_indexed_rows( & comb_target, 0 .. num_indices_row );        
        // println!("comb_target_inv:");
        // print_indexed_rows( & comb_target_inv, 0 .. num_indices_row );                        
        // println!("comb_target_inv * matrix_to_factor * comb_source:");
        // print_indexed_rows( & product_target_comb_inv_times_matrix_to_factor_times_source_comb, 0 .. num_indices_row );                                        
        for column_index in 0 .. num_indices_col {
            println!("row: {:?}", column_index );
            println!("row {:?}: {:?}", column_index, product_source.row( & column_index ).collect_vec() );
            itertools::assert_equal( product_source.row( & column_index ), std::iter::once( (column_index, 1) )   );
        }


        // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
        for column_index in 0 .. num_indices_col { 
            assert_eq!(
                product_source.row( & column_index ).collect_vec(),
                vec![ (column_index, 1) ]
            ) 
        }
      

        // check that the product of the target COMB with its inverse is identity T * T^{-1} = I
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target.row( & row_index ).collect_vec(),
                vec![ (row_index, 1) ]
            ) 
        }   
       
        
        // check the factorization T^{-1} * D * S = M
        for row_index in 0 .. num_indices_row { 
            assert_eq!(
                product_target_comb_inv_times_matrix_to_factor_times_source_comb.row( &row_index ).collect_vec(),
                matching.row( & row_index ).collect_vec()
            ) 
        }   


        // ----------------------------------------------------------------------------------------------------------------
        // check that all four (inverse) COMB's are internally valid
        // see documentation for `matrix_oracle_is_internally_consistent`, for details

        assert!(
            matrix_oracle_is_internally_consistent(
                comb_source_ref, 
                0..num_indices_col, 
                0..num_indices_col
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_source_inv_ref, 
                0..num_indices_col, 
                0..num_indices_col
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_target_ref, 
                0..num_indices_row, 
                0..num_indices_row
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_target_inv_ref, 
                0..num_indices_row, 
                0..num_indices_row
            )                               
        );   

        // ----------------------------------------------------------------------------------------------------------------
        // check that all four (inverse) COMB's return entries in the proper order
        // see documentation for `matrix_order_operators_are_internally_consistent`, for details

        assert!(
            matrix_order_operators_are_internally_consistent(
                comb_source_ref, 
                0..num_indices_col, 
                0..num_indices_col
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_source_inv_ref, 
                0..num_indices_col, 
                0..num_indices_col
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_target_ref, 
                0..num_indices_row, 
                0..num_indices_row
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_target_inv_ref, 
                0..num_indices_row, 
                0..num_indices_row
            ).is_ok()                                 
        );        
              
    }










    fn test< MatrixToFactor >(
        umatch:                 GimbledUmatch< MatrixToFactor >,
        sorted_row_indices:     Vec<MatrixToFactor::RowIndex>,
        sorted_column_indices:  Vec<MatrixToFactor::ColumnIndex>,
    ) 
        where   
            MatrixToFactor:         MatrixAlgebra<
                                        ColumnIndex:            Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
                                        RowIndex:               Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct      
                                        RingOperator:           DivisionRingOperations,
                                        RowEntry:               KeyValPair,
                                        ColumnEntry:            KeyValPair,        
                                    >,             
    {





        let matching = umatch.generalized_matching_matrix_ref();
        let matrix_to_factor = umatch.matrix_to_factor_ref();

        let comb_target = umatch.target_comb();
        let comb_target_inv = umatch.target_comb_inverse();        
        let comb_source = umatch.source_comb();        
        let comb_source_inv = umatch.source_comb_inverse();  

        let comb_target_ref         =   & comb_target;
        let comb_target_inv_ref         =   & comb_target_inv;
        let comb_source_ref         =   & comb_source;
        let comb_source_inv_ref         =   & comb_source_inv;            
        

        let product_source = ProductMatrix::new( comb_source_ref, comb_source_inv_ref );
        let product_target = ProductMatrix::new( comb_target_ref, comb_target_inv_ref );        
        let product_target_comb_inv_times_matrix_to_factor = ProductMatrix::new( comb_target_inv_ref, matrix_to_factor );      
        let product_target_comb_inv_times_matrix_to_factor_times_source_comb = ProductMatrix::new( product_target_comb_inv_times_matrix_to_factor, comb_source_ref );                        



        let a = product_target_comb_inv_times_matrix_to_factor.order_operator_for_row_entries();

        // check that the product of the source COMB with its inverse is identity: S * S^{-1} = I
        assert!(
            product_is_identity_matrix(
                comb_source_ref,
                comb_source_inv_ref,
                sorted_column_indices.iter().cloned()
            )
        );

        assert!(
            product_is_identity_matrix(
                comb_target_ref,
                comb_target_inv_ref,
                sorted_row_indices.iter().cloned()
            )
        );        
       
        
        // check the factorization T^{-1} * D * S = M
        for row_index in sorted_row_indices.iter().cloned() { 
            assert_eq!(
                product_target_comb_inv_times_matrix_to_factor_times_source_comb.row( &row_index ).collect_vec(),
                matching.row( & row_index ).collect_vec()
            ) 
        }   


        // ----------------------------------------------------------------------------------------------------------------
        // check that all four (inverse) COMB's are internally valid
        // see documentation for `matrix_oracle_is_internally_consistent`, for details

        assert!(
            matrix_oracle_is_internally_consistent(
                comb_source_ref, 
                sorted_column_indices.iter().cloned(), 
                sorted_column_indices.iter().cloned()
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_source_inv_ref, 
                sorted_column_indices.iter().cloned(), 
                sorted_column_indices.iter().cloned()
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_target_ref, 
                sorted_row_indices.iter().cloned(), 
                sorted_row_indices.iter().cloned()
            )
            &&
            matrix_oracle_is_internally_consistent(
                comb_target_inv_ref, 
                sorted_row_indices.iter().cloned(), 
                sorted_row_indices.iter().cloned()
            )    
            // &&   
            // matrix_oracle_is_internally_consistent(
            //     product_target_comb_inv_times_matrix_to_factor, 
            //     sorted_row_indices.iter().cloned(), 
            //     sorted_column_indices.iter().cloned()
            // )                                     
        );   

        // ----------------------------------------------------------------------------------------------------------------
        // check that all four (inverse) COMB's return entries in the proper order
        // see documentation for `matrix_order_operators_are_internally_consistent`, for details

        assert!(
            matrix_order_operators_are_internally_consistent(
                comb_source_ref, 
                sorted_column_indices.iter().cloned(), 
                sorted_column_indices.iter().cloned()
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_source_inv_ref, 
                sorted_column_indices.iter().cloned(), 
                sorted_column_indices.iter().cloned()
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_target_ref, 
                sorted_row_indices.iter().cloned(), 
                sorted_row_indices.iter().cloned()
            ).is_ok()
            &&
            matrix_order_operators_are_internally_consistent(
                comb_target_inv_ref, 
                sorted_row_indices.iter().cloned(), 
                sorted_row_indices.iter().cloned()
            ).is_ok()                                 
        );        
    }








}