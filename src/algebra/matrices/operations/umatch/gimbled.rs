use ndarray::Order;

use crate::{algebra::{
    matrices::{
        operations::{
            combine_rows_and_columns::{LinearCombinationOfColumns, LinearCombinationOfRows},
            umatch::row_major::{comb::{SourceComb, SourceCombInverse, TargetComb, TargetCombInverse}, Umatch},
            MatrixOracleOperations,
        },
        query::{MatrixAlgebra, MatrixOracle}, 
        types::{matching::GeneralizedMatchingMatrixWithSequentialOrder, packet::MatrixAlgebraPacket, transpose::OrderAntiTranspose, two_type::TwoTypeMatrix}
    }, 
    rings::traits::DivisionRingOperations, 
    vectors::{entries::{KeyValGet, KeyValPair}, operations::VectorOperations, },
}, utilities::order::{JudgeOrder, OrderOperatorByKeyCustom, ReverseOrder}};



use std::hash::Hash;













#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GimbledUmatch< MatrixToFactor >
    where   
        MatrixToFactor:                     MatrixAlgebra,
        MatrixToFactor::ColumnIndex:        Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
        MatrixToFactor::RowIndex:           Hash, // required by the `GeneralizedMatchingMatrixWithSequentialOrder` struct
{
    RowMajor( 
        Umatch< MatrixToFactor > 
    ),
    ColumnMajor( 
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

    /// Returns a column-major version of the U-match, or `None`, if `self` is already in column-major form.
    /// 
    /// "Column-major" means that the inner data structure is a [Umatch] decomposition of the antitranspose of the matrix to factor.
    pub fn column_major( & self ) -> Option< Self >

        where
           MatrixToFactor:  Clone
    {
        match self {
            GimbledUmatch::RowMajor( umatch ) => {

                let column_indices = umatch.matched_column_indices_in_ascending_order();
                let matrix_to_factor_antitranspose = OrderAntiTranspose::new(
                    ( * umatch.matrix_to_factor_ref() ).clone()
                );

                let umatch_antitranspose = Umatch::new(
                    matrix_to_factor_antitranspose,
                    column_indices.into_iter()
                );

                return Some( GimbledUmatch::ColumnMajor(
                    (
                        umatch_antitranspose,
                        umatch.generalized_matching_matrix_ref().clone(),
                    )
                ));
            },
            GimbledUmatch::ColumnMajor( _ ) => {
                return None
            },
        }

    }


    /// Returns `true` if the U-match is in column-major form, `false` otherwise.
    /// 
    /// "Column-major" means that the inner data structure is a [Umatch] decomposition of the antitranspose of the matrix to factor.
    /// "Row-major" means that the inner data structure is a [Umatch] decomposition of the matrix to factor (no antitranspose).
    /// 
    /// We typically expect the column-major to be faster for looking up columns of the differential COMB,
    /// and row-major to be faster for looking up rows of the inverse differential COMB.
    pub fn is_column_major( &self ) -> bool {
        match self {
            GimbledUmatch::RowMajor( _ ) => false,
            GimbledUmatch::ColumnMajor(_) => true,
        }
    }

    /// Returns a reference to the matrix we wish to factor.
    pub fn matrix_to_factor_ref( &self ) -> & MatrixToFactor {
        match self {
            GimbledUmatch::RowMajor( umatch ) => {
                umatch.matrix_to_factor_ref()
            },
            GimbledUmatch::ColumnMajor( (umatch, _) ) => {
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
            GimbledUmatch::RowMajor( umatch ) => {
                umatch.generalized_matching_matrix_ref()
            },
            GimbledUmatch::ColumnMajor( (_, gmm) ) => {
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


    /// Returns the source comb of the U-match.
    pub fn source_comb< 'a >( &'a self ) -> 
        TwoTypeMatrix<
            SourceComb<'a, MatrixToFactor >,
            OrderAntiTranspose<
                TargetCombInverse<
                    'a, 
                    OrderAntiTranspose< MatrixToFactor>
                >
            >
        > 
    {
        match self {
            GimbledUmatch::RowMajor( umatch ) => {
                TwoTypeMatrix::Version1(
                    umatch.source_comb(),
                )
            },
            GimbledUmatch::ColumnMajor( (umatch, _) ) => {
                TwoTypeMatrix::Version2(
                    OrderAntiTranspose::new( umatch.target_comb_inverse() )
                )
            },
        }
    }


    /// Returns the inverse source comb of the U-match.
    pub fn source_comb_inverse< 'a >( &'a self ) -> 
        TwoTypeMatrix<
            SourceCombInverse<'a, MatrixToFactor >,
            OrderAntiTranspose<
                TargetComb<
                    'a, 
                    OrderAntiTranspose< MatrixToFactor>
                >
            >
        > 
    {
        match self {
            GimbledUmatch::RowMajor( umatch ) => {
                TwoTypeMatrix::Version1(
                    umatch.source_comb_inverse(),
                )
            },
            GimbledUmatch::ColumnMajor( (umatch, _) ) => {
                TwoTypeMatrix::Version2(
                    OrderAntiTranspose::new( umatch.target_comb() )
                )
            },
        }
    }    


    /// Returns the target comb of the U-match.
    pub fn target_comb< 'a >( &'a self ) -> 
        TwoTypeMatrix<
            TargetComb<'a, MatrixToFactor>,
            OrderAntiTranspose<
                SourceCombInverse<'a, OrderAntiTranspose< MatrixToFactor > >
            >
        >
    {
        match self {
            GimbledUmatch::RowMajor( umatch ) => {
                TwoTypeMatrix::Version1(
                    umatch.target_comb(),
                )
            },
            GimbledUmatch::ColumnMajor( (umatch, _) ) => {
                TwoTypeMatrix::Version2(
                    OrderAntiTranspose::new( umatch.source_comb_inverse() )
                )
            },
        }
    }



    /// Returns the target comb of the U-match.
    pub fn target_comb_inverse< 'a >( &'a self ) -> 
        TwoTypeMatrix<
            TargetCombInverse<'a, MatrixToFactor>,
            OrderAntiTranspose<
                SourceComb<'a, OrderAntiTranspose< MatrixToFactor > >
            >
        > 
    {
        match self {
            GimbledUmatch::RowMajor( umatch ) => {
                TwoTypeMatrix::Version1(
                    umatch.target_comb_inverse(),
                )
            },
            GimbledUmatch::ColumnMajor( (umatch, _) ) => {
                TwoTypeMatrix::Version2(
                    OrderAntiTranspose::new( umatch.source_comb() )
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
    







    /// Solve `Dx = b`, where `D` is the factored matrix.
    /// 
    /// Returns `None` if there is no solution.
    /// 
    /// # Arguments
    /// 
    /// The iterable `b` can iterate over entries in any order. 
    /// If `b` contains multiple entries with the same index, these entries will be summed.
    /// 
    /// # Calculation
    /// 
    /// The U-match equation implies `TM = DS` implies that `TMS^{-1} = D`. It can therefore be shown that
    /// `SM^{^-1}T^{-1}` is a generalized inverse of `D`, where `M^{-1}` is the generalized
    /// inverse of `M` obtained by transposing and inverting nonzero entries. Therefore the solution `x`
    /// can be computed as `SM^{^-1}T^{-1}b`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::operations::umatch::row_major::Umatch;
    /// use oat_rust::algebra::matrices::operations::umatch::gimbled::GimbledUmatch;
    /// use oat_rust::algebra::matrices::operations::MatrixOracleOperations;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;   
    /// use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket; 
    /// use itertools::Itertools;    
    /// 
    /// // DEFINE THE MATRIX
    /// // ===============================
    /// let matrix          =   VecOfVec::new( 
    ///                                         vec![   
    ///                                                     vec![(0,true), (1,true), (2,true)],
    ///                                                     vec![                            ], 
    ///                                                     vec![                    (2,true)], 
    ///                                         ] 
    ///                                     ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( &matrix );
    ///                                 
    /// // COMPUTE U-MATCH
    /// // ===============================
    ///                                 
    /// let umatch
    ///     =   Umatch::new(
    ///             & matrix,  // the matrix we wish to factor
    ///             (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
    ///         );   
    /// 
    /// // WRAP IN A GIMBLE
    /// // ===============================
    ///                                 
    /// let gimbled_umatch = GimbledUmatch::RowMajor( umatch );       
    /// 
    /// // SOLVE Dx = b FOR x
    /// // ===============================
    /// 
    /// let b   =   [ (0,true), (2,true) ]; 
    /// let x   =   gimbled_umatch.solve_dx_equals_b( b.clone() ).unwrap();
    /// let dx  =   matrix.multiply_with_column_vector(x);
    /// assert!( dx.eq( b ) );  
    /// 
    /// // SOLVE Dx = b FOR x (WHEN NO SOLUTION EXISTS)
    /// // ===============================
    /// 
    /// let b   =   [ (1,true) ]; 
    /// assert!( gimbled_umatch.solve_dx_equals_b( b ).is_none() ); // no solution exists
    /// 
    /// // REPLACE THE INNER U-MATCH WITH A COLUMN-MAJOR VERSION
    /// // ===============================
    /// 
    /// let gimbled_umatch_column_major = gimbled_umatch.column_major().unwrap();
    /// 
    /// // SOLVE Dx = b FOR x
    /// // ===============================
    /// 
    /// let b   =   [ (0,true), (2,true) ]; 
    /// let x   =   gimbled_umatch_column_major.solve_dx_equals_b( b.clone() ).unwrap();
    /// let dx  =   matrix.multiply_with_column_vector(x);
    /// assert!( dx.eq( b ) );  
    /// 
    /// // SOLVE Dx = b FOR x (WHEN NO SOLUTION EXISTS)
    /// // ===============================
    /// 
    /// let b   =   [ (1,true) ]; 
    /// assert!( gimbled_umatch_column_major.solve_dx_equals_b( b ).is_none() ); // no solution exists
    /// ```
    pub fn solve_dx_equals_b< Vector >( &self, b: Vector ) 
        -> 
        Option< Vec< MatrixToFactor::RowEntry > >
    where
        Vector:     IntoIterator<Item=MatrixToFactor::ColumnEntry>,                            
    {

        match self {
            GimbledUmatch::RowMajor( umatch ) => {
                umatch.solve_dx_equals_b( b )
            },
            GimbledUmatch::ColumnMajor( (umatch, _) ) => {

                // the method umatch.solve_xd_equals_b requires b to be SORTED in the order of the underlying matrix (which in this case is an antitranspose)
                // therefore we must sort
                let mut b: Vec<_> = b.into_iter().collect();
                let order_operator = umatch.matrix_to_factor_ref().order_operator_for_row_entries();
                b.sort_by( |a,b| order_operator.judge_cmp( a, b ) ); // sort the entries in reverse order, so that we can use the anti-transposed version of the problem

                // the vector b must also be simplified
                let ring_operator = umatch.ring_operator();
                let b = b.into_iter().peekable().simplify( ring_operator );

                // here we solve the anti-transposed version of the problem
                let mut solution = umatch.solve_xd_equals_b( b )?;
                solution.reverse(); // the solution is in reverse order, so we need to reverse it to get the correct order
                return Some( solution );
            },
        }
    }



    /// Solve `xD = b`, where `D` is the factored matrix.
    /// 
    /// Returns `None` if there is no solution.
    /// 
    /// # Arguments
    /// 
    /// The iterable `b` can iterate over entries in any order. 
    /// If `b` contains multiple entries with the same index, these entries will be summed.
    /// 
    /// # Calculation
    /// 
    /// The U-match equation implies `TM = DS` implies that `TMS^{-1} = D`. It can therefore be shown that
    /// `SM^{^-1}T^{-1}` is a generalized inverse of `D`, where `M^{-1}` is the generalized
    /// inverse of `M` obtained by transposing and inverting nonzero entries. Therefore the solution `x`
    /// can be computed as `SM^{^-1}T^{-1}b`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::operations::umatch::row_major::Umatch;
    /// use oat_rust::algebra::matrices::operations::umatch::gimbled::GimbledUmatch;
    /// use oat_rust::algebra::matrices::operations::MatrixOracleOperations;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;   
    /// use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket; 
    /// use itertools::Itertools;    
    /// 
    /// // DEFINE THE MATRIX
    /// // ===============================
    /// let matrix          =   VecOfVec::new( 
    ///                                         vec![   
    ///                                                     vec![(0,true), (1,true), (2,true)],
    ///                                                     vec![                            ], 
    ///                                                     vec![                    (2,true)], 
    ///                                         ] 
    ///                                     ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( &matrix );
    ///                                 
    /// // COMPUTE U-MATCH
    /// // ===============================
    ///                                 
    /// let umatch
    ///     =   Umatch::new(
    ///             & matrix,  // the matrix we wish to factor
    ///             (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
    ///         );   
    /// 
    /// // WRAP IN A GIMBLE
    /// // ===============================
    ///                                 
    /// let gimbled_umatch = GimbledUmatch::RowMajor( umatch );       
    /// 
    /// // SOLVE Dx = b FOR x
    /// // ===============================
    /// 
    /// let b   =   [ (0,true), (1,true) ]; 
    /// let x   =   gimbled_umatch.solve_xd_equals_b( b.clone() ).unwrap();
    /// let xd  =   matrix.multiply_with_row_vector(x);
    /// assert!( xd.eq( b ) );  
    /// 
    /// // SOLVE Dx = b FOR x (WHEN NO SOLUTION EXISTS)
    /// // ===============================
    /// 
    /// let b   =   [ (1,true) ]; 
    /// assert!( gimbled_umatch.solve_xd_equals_b( b ).is_none() ); // no solution exists
    /// 
    /// // REPLACE THE INNER U-MATCH WITH A COLUMN-MAJOR VERSION
    /// // ===============================
    /// 
    /// let gimbled_umatch_column_major = gimbled_umatch.column_major().unwrap();
    /// 
    /// // SOLVE Dx = b FOR x
    /// // ===============================
    /// 
    /// let b   =   [ (0,true), (1,true) ]; 
    /// let x   =   gimbled_umatch_column_major.solve_xd_equals_b( b.clone() ).unwrap();
    /// let xd  =   matrix.multiply_with_row_vector(x);
    /// assert!( xd.eq( b ) );  
    /// 
    /// // SOLVE Dx = b FOR x (WHEN NO SOLUTION EXISTS)
    /// // ===============================
    /// 
    /// let b   =   [ (1,true) ]; 
    /// assert!( gimbled_umatch_column_major.solve_xd_equals_b( b ).is_none() ); // no solution exists
    /// ```
    pub fn solve_xd_equals_b< Vector >( &self, b: Vector ) 
        -> 
        Option< Vec< MatrixToFactor::ColumnEntry > >
    where
        Vector:     IntoIterator<Item=MatrixToFactor::RowEntry>,                            
    {

        match self {
            GimbledUmatch::RowMajor( umatch ) => {

                // the method umatch.solve_xd_equals_b requires b to be SORTED in the order of the underlying matrix (which in this case is an antitranspose)
                // therefore we must sort
                let mut b: Vec<_> = b.into_iter().collect();
                let order_operator = umatch.matrix_to_factor_ref().order_operator_for_row_entries();
                b.sort_by( |a,b| order_operator.judge_cmp( a, b ) ); // sort the entries in reverse order, so that we can use the anti-transposed version of the problem

                // the vector b must also be simplified
                let ring_operator = umatch.ring_operator();
                let b = b.into_iter().peekable().simplify( ring_operator );

                umatch.solve_xd_equals_b( b )
            },
            GimbledUmatch::ColumnMajor( (umatch, _) ) => {

                // here we solve the anti-transposed version of the problem
                let mut solution = umatch.solve_dx_equals_b( b )?;
                solution.reverse(); // the solution is in reverse order, so we need to reverse it to get the correct order
                return Some( solution );
            },
        }
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
    #[test]
    fn comprehensive_test() {   

        use crate::algebra::matrices::operations::umatch::row_major::Umatch;            

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
        let umatch = GimbledUmatch::RowMajor(umatch);

        validate_umatch( 
            & umatch, 
            & (0..num_indices_row).collect_vec(), 
            & (0..num_indices_col).collect_vec() 
        );

        
        validate_umatch( 
            & umatch.column_major().unwrap(), 
            & (0..num_indices_row).collect_vec(), 
            & (0..num_indices_col).collect_vec() 
        );

  
              
    }









    #[allow(dead_code)]
    fn validate_umatch< MatrixToFactor >(
        umatch:                 & GimbledUmatch< MatrixToFactor >,
        sorted_row_indices:     & Vec<MatrixToFactor::RowIndex>,
        sorted_column_indices:  & Vec<MatrixToFactor::ColumnIndex>,
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
        let matching_packet = umatch.generalized_matching_matrix_ref_packet();
        let matrix_to_factor = umatch.matrix_to_factor_ref();

        let comb_target = umatch.target_comb();
        let comb_target_inv = umatch.target_comb_inverse();        
        let comb_source = umatch.source_comb();        
        let comb_source_inv = umatch.source_comb_inverse();  

        let comb_target_ref         =   & comb_target;
        let comb_target_inv_ref         =   & comb_target_inv;
        let comb_source_ref         =   & comb_source;
        let comb_source_inv_ref         =   & comb_source_inv;            
         
        let tinvd = ProductMatrix::new( comb_target_inv_ref, matrix_to_factor );
        let msinv = ProductMatrix::new( matching_packet, comb_source_inv_ref );


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
       
        
        // check the factorization T^{-1} * D = M * S^{-1}
        for row_index in sorted_row_indices.iter().cloned() { 
            assert!(
                tinvd.row( & row_index )
                    .eq(
                        msinv.row( & row_index )
                    )
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
        );
        assert!(
            matrix_oracle_is_internally_consistent(
                comb_source_inv_ref, 
                sorted_column_indices.iter().cloned(), 
                sorted_column_indices.iter().cloned()
            )
        );
        assert!(            
            matrix_oracle_is_internally_consistent(
                comb_target_ref, 
                sorted_row_indices.iter().cloned(), 
                sorted_row_indices.iter().cloned()
            )
        );
        assert!(
            matrix_oracle_is_internally_consistent(
                comb_target_inv_ref, 
                sorted_row_indices.iter().cloned(), 
                sorted_row_indices.iter().cloned()
            )                                       
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