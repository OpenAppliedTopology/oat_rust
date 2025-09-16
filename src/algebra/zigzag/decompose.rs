//! Computes an interval decomposition of a zigzag module
//! 
//! 


use derive_getters::{Getters, Dissolve};
use derive_new::new;

use crate::{algebra::{matrices::{operations::{invert::InverseUpperTriangularMatrix, umatch::row_major::Umatch}, query::MatrixOracle, types::{packet::MatrixAlgebraPacket, vec_of_vec::sorted::VecOfVec}}, vectors::operations::VectorOperations}, utilities::order::{OrderOperatorAuto, OrderOperatorByKey}};

use crate::algebra::rings::traits::DivisionRingOperations;

use std::{collections::HashMap, fmt::Debug};
use std::time::Instant;


#[derive(new,Clone,Debug,Getters,Dissolve,Eq,PartialEq)]
pub struct Diagonalization< RingElement >{
    list_of_single_bar_basis_vector_index_ledgers:      Vec< SingleBarBasisVectorIndexLedger >,
    bases_encoded_as_rows_of_invertible_matrices:       Vec< VecOfVec< usize, RingElement > >
}

impl < RingElement >

    Diagonalization
        < RingElement >
    {

    /// Number of bars in the diagonalization
    pub fn number_of_bars( &self ) -> usize {
        self.list_of_single_bar_basis_vector_index_ledgers.len()
    }


    /// Returns the index of the basis vector over vertex `v` that intersects bar number `b``
    /// 
    /// Returns `None` if `b ≥ {number of bars}` or `v` lies outside the interval for bar `b`
    pub fn basis_vector_index_for_bar_b_over_vertex_v(
            &self,
            b: usize, 
            v: usize 
        ) 
        -> Option< usize > {
        let bar                                         =   self.bar( b )?;  
        return bar.vertex_to_basis_vector_index_opt( v )
    }


    /// Returns the the basis vector over vertex `v` that intersects bar number `b``
    /// 
    /// Returns `None` if `b ≥ {number of bars}` or `v` lies outside the interval for bar `b`
    pub fn basis_vector_for_bar_b_over_vertex_v(
            & self,
            b: usize, 
            v: usize 
        ) 
        -> Option< & Vec< (usize, RingElement ) > > {        
        let basis_vector_index_in_v                     =   self.basis_vector_index_for_bar_b_over_vertex_v( b, v )?;
        let basis_vector                                =   self
                                                                .bases_encoded_as_rows_of_invertible_matrices[ v ]
                                                                .row_ref( basis_vector_index_in_v );
        return Some( basis_vector )        
    }


    /// The sequence of bases that diagonalize the representation
    /// 
    /// The sequence is formatted as a `& Vec< VecOfVec< usize, RingElement > >`
    pub fn bases( &self ) -> & Vec< VecOfVec< usize, RingElement > > {
        &self.bases_encoded_as_rows_of_invertible_matrices
    }


    /// Retuns an object which records every basis vector intersected by bar number `i`
    /// 
    /// Returns `None` if the barcode contains `i` or fewer bars.
    pub fn bar< 'a >( &'a self,  i: usize ) -> Option< & SingleBarBasisVectorIndexLedger > {
        if i    >=    self.list_of_single_bar_basis_vector_index_ledgers.len() {
            None
        } else {
            Some(
                & self.list_of_single_bar_basis_vector_index_ledgers[ i ]
            )
        }
    }


    /// Returns a reference to the list of bars
    pub fn bars( &self ) -> & Vec< SingleBarBasisVectorIndexLedger > {
        & self.list_of_single_bar_basis_vector_index_ledgers
    }
}




/// A list of the indices of the basis vectors through which a bar passes
#[derive(new,Clone,Debug,Dissolve,Eq,PartialEq,Ord,PartialOrd)]
pub struct SingleBarBasisVectorIndexLedger {
    leftmost_vertex:        usize,
    basis_vector_indices:   Vec< usize >,
}

impl SingleBarBasisVectorIndexLedger {


    /// Returns the pair `(a,b)` such that the interval passes through spaces `a, a+1, .., b-1`
    /// 
    /// Concretely, the interval includes `b-1` but not `b`.
    pub fn supporting_interval( &self ) -> (usize,usize) {
        return ( self.leftmost_vertex, self.leftmost_vertex + self.basis_vector_indices.len() )
    }




    /// A `Range<usize>` that iterats over vertices in the associated interval
    pub fn supporting_vertices( &self ) -> std::ops::Range<usize> {
        let (a,b) = self.supporting_interval();
        a .. b
    }


    /// Iterates over all `(vertex, basis_vector_index)` pairs, in ascending order of vertex
    pub fn iter( &self ) -> std::iter::Zip<std::ops::Range<usize>, std::iter::Cloned<std::slice::Iter<'_, usize>>> {
        self.supporting_vertices().zip( self.basis_vector_indices.iter().cloned() )
    }



    /// A `Range<usize>` that iterates over the indices of the directed edges in the interval
    /// 
    /// Concretely, we call the arrow `p <--> p+1` the `p`th edge.
    /// 
    /// We say that arrow `p` is in the interval if both endpoints are in the interval.
    pub fn supporting_edge_indices( &self ) -> std::ops::Range<usize> {
        let (a,b) = self.supporting_interval();
        a .. ( b - 1 )
    }


    /// Length of the bar
    /// 
    /// Calculated as the length of `self.basis_vector_indices`
    pub fn bar_length( &self ) -> usize {
        self.basis_vector_indices.len()
    }


    /// Returns `Some(i)` if this bar passes through the `i` the basis vector in space `p`
    /// 
    /// Retuns `None` if the bar does not pass through vector space `p` at all
    pub fn vertex_to_basis_vector_index_opt( &self, n: usize ) -> Option< usize > {

        if n < self.leftmost_vertex { return None }
        if n >= self.leftmost_vertex + self.basis_vector_indices.len() { return None }

        return Some( self.basis_vector_indices[ n - self.leftmost_vertex ] )
    }


    /// Returns the left endpoint of the bar
    pub fn leftmost_vertex( &self ) -> usize { self.leftmost_vertex.clone() }


    /// The pair `(a,b)`, where `[a,b)` is the half-open interval of the bar
    pub fn inclusive_left_and_exclusive_right_endpoints( &self ) -> (usize, usize) {
        (
            self.leftmost_vertex(),
            self.leftmost_vertex() + self.bar_length()
        )
    }
}






/// A quiver of type `An` together with a representation
/// 
/// Concretely, the quiver is a directed graph on vertex set `0, .., n` have one directed edge of form `(i, i+1)` or `(i+1,i)` for all `i < n`.
/// The representation consists of
/// - A copy of `F^{i_k}` to each vertex `k`, where `F` is a vector space, `i_k` is a nonnegative integer
/// - A matrix to each directed edge, where matrix dimensions match the dimensions of the vector spaces assigned to each of the incident vertices.
///   This means that the matrix `M` for a directed edge `(p,q)` should have size `i_p x i_q`.
///   **We think of this matrix as a linear map on row vectors sending `r` to `r * M`, not a map of column vectors sending `c` to `M * c`.**
///   This convention is a bit unusual, but we adopt it because allows us to use some powerful computational tools.
///   In particular, it allows us to use a U-match factorization package which works most efficiently with row-major matrices.
/// 
/// 
/// This struct holds three vectors:
/// 
/// - `vector_space_dimensions` is the dimension of the vector space sitting over vertex `i`
/// - `arrow_directions[i]` equals `true` iff arrow `i <--> i+1` points forward, i.e. from `i` to `i+1`
/// - `matrix[i]` is the matrix representation of the map between vector spaces `i` and `i+1`
/// 
/// It also stores a
/// 
/// - `ring_operator`, which is an object that can perform the basic algebraic operations of the coefficient field (addition, multiplication, division, etc.)
#[derive(new,Clone,Debug,Dissolve,Eq,PartialEq)]
pub struct QuiverReprsentation
                < RingOperator > 
    where
        RingOperator:               DivisionRingOperations,
{
    arrow_directions:               Vec< bool >,
    matrices:                       Vec< VecOfVec<usize, RingOperator::Element> >,
    vector_space_dimensions:        Vec< usize >,
    ring_operator:                  RingOperator,
}

impl < RingOperator > 

    QuiverReprsentation
        < RingOperator > 

    where 
        RingOperator:               Clone + DivisionRingOperations,
{

    /// Number of arrows in the quiver
    pub fn number_of_arrows( &self ) -> usize {
        self.arrow_directions.len()
    }

    // Number of vertices in the quiver
    pub fn number_of_vertices( &self ) -> usize {
        self.vector_space_dimensions.len()
    }

    /// The list of dimensions of the vector spaces
    pub fn vector_space_dimensions( &self ) -> & Vec< usize > {
        & self.vector_space_dimensions
    }

    /// Reference to the internally stored vector of arrow directions
    /// 
    /// If `v` is this vector then `v[p] == true` implies that arrow `p` points right; otherwise arrow `p` points left.
    pub fn arrow_directions( &self ) -> & Vec< bool > {
        & self.arrow_directions
    }


    /// Returns a reference to the sequence of matrices corresponding to the arrows in the directed graph
    pub fn arrow_matrices( &self ) -> & Vec< VecOfVec<usize, RingOperator::Element> > {
        & self.matrices
    }


    /// Returns the dimension of the `n`th space in the sequence
    pub fn dimension_of_space_over_vertex( & self, n: usize ) -> Option< usize > {
        if n <= self.number_of_arrows() {
            Some( self.vector_space_dimensions[ n ].clone() )
        } else {
            None
        }
    }

    /// Returns the ring operator for the coefficient field.
    /// 
    /// A "ring operator" is an object that performs the basic algebraic operations of a
    /// ring on the elements of the ring. For example, you could use a ring operator to 
    /// multiply or add two elements. See the OAT documentation for [rings](crate::algebra::rings) for details.
    pub fn ring_operator( & self ) -> RingOperator {
        self.ring_operator.clone()
    }


    /// Returns a matrix algebra packet for the given vertex, or `None` if the vertex is out of bounds
    pub fn matrix_packet_for_vertex( &self, vertex: usize ) 
        -> Option< MatrixAlgebraPacket< 
            & VecOfVec< usize, RingOperator::Element >,
            RingOperator,
            OrderOperatorByKey,
            OrderOperatorAuto,
            OrderOperatorByKey,
            OrderOperatorAuto,                        
        > > 
    {
        if vertex >= self.number_of_vertices() {
            return None
        } else {
            Some( MatrixAlgebraPacket::with_default_order(
                & self.matrices[ vertex ],
                self.ring_operator.clone()
            ))
        }
    }


    /// Checks that the user input data is a valid quiver representation
    /// 
    /// Specifically, it checks that 
    /// - `n_arrows = n_matrices = n_vertices - 1` if `n_vertices > 0`
    /// - for each arrow `p <--> p + 1`, the size of the associated matrix agrees with the reported dimensions for `p` and `p+1`
    /// 
    /// If these tests pass then the function returns `Ok(())`. Otherwise it returns `Err(hash)`, where `hash` is a dictionary containing
    /// information about the error.
    pub fn validate_representation( & self ) -> Result< (), HashMap< &str, usize > >    
    {
        let start = Instant::now(); // Start the timer




        let n_arrows                        =   self.arrow_directions.len();
        let n_matrices                      =   self.matrices.len();
        let n_vertices                        =   self.vector_space_dimensions.len();
        
        let directions                          =   & self.arrow_directions;
        let dimensions                                  =   self.vector_space_dimensions();

        // ensure the number of arrows matches the number of matrices and agrees with the number of space dimensions
        if  ( n_arrows      !=  n_matrices  )
            ||
            ( n_arrows + 1  !=  n_vertices    ) 
        {  
            let mut err                     =   HashMap::new();
            err.insert("Error: there's a disagreement between at least two of the following three numbers: number of arrows, number of vertices, and number of matrices.", 0);
            err.insert("number of arrows", n_arrows);
            err.insert("number of matrices", n_matrices);
            err.insert("number of vertices", n_vertices);                        
            return Err( err )
        }

        // in this case there are no arrows, so there is no way for matrix versus versus vector space dimensions to disagree
        if n_vertices <= 1 {
            return Ok(())
        }

        let mut n_rows;
        let mut min_num_columns;
        let mut expected_n_rows;
        let mut expected_n_columns;

        // for each vertex p < n_arrows, check that the reported dimension of space p agrees with the size of the matrix associated with the arrow p <--> p+1        
        for p in 0 .. n_arrows {

            let source                      =   if directions[p] { p   } else { p+1 };
            let target                      =   if directions[p] { p+1 } else { p   };   

            n_rows                          =   self.matrices[p].number_of_rows(); // number of rows of associated matrix
            min_num_columns                 =   self.matrices[p].max_column_index().map(|x| x + 1 ).unwrap_or(0); // maximum index of any structural nonzero entry

            if directions[p] {
                expected_n_rows             =   dimensions[p];
                expected_n_columns          =   dimensions[p+1];
            } else {
                expected_n_rows             =   dimensions[p+1];
                expected_n_columns          =   dimensions[p];              
            }
            
            // Check that vector space dimensions agree with the number of rows of each matrix
            if n_rows != expected_n_rows {
                let mut err                 =   HashMap::new();
                err.insert("Edge source: ", source );
                err.insert("Edge target: ", target );                    
                err.insert("Vector space dimension: ", expected_n_rows );
                err.insert("Number of matrix rows: ", n_rows );                    
                return Err(err)
            }

            // Check that vector space dimensions agree with the maximum column indices of each matrix
            if expected_n_columns < min_num_columns {
                let mut err                 =   HashMap::new();             
                err.insert("Edge source: ", source );
                err.insert("Edge target: ", target );                    
                err.insert("Vector space dimension: ",    expected_n_columns );
                err.insert("Minimum number of columns required: ",  min_num_columns );                    
                return Err(err)
            }     
        }       

        println!("Time to validate representation: {:?}", start.elapsed()); // Output the duration        

        return Ok(())  

    }



    /// Decompose the quiver representation as a direct sum of interval modules
    pub fn diagonalize( &self ) -> 
        Result< 
            Diagonalization< RingOperator::Element >, 
            HashMap< &str, usize > ,
        > {            

        //  ENSURE THE REPRESENATION IS VALID
        // ------------------------------------------------------------------------        
        self.validate_representation()?; 


        let start = Instant::now(); // Start the timer        
        

        // INITIALIZE SOME VARIABLES
        // ------------------------------------------------------------------------

        let n_arrows                        =   self.number_of_arrows();
        let n_vertices                        =   self.number_of_vertices();
        let dimensions                      =   self.vector_space_dimensions();
        let directions                      =   self.arrow_directions();

        // initialize list of bases
        let mut B                           =   Vec::with_capacity(n_vertices);   // a sequence of bases (one for each space)        
        let unity                           =   RingOperator::one();
        B.push( VecOfVec::diagonal_matrix( unity.clone(), dimensions[0] ) ); // push an identity matrix representing the standard basis for vector space 0


        // this is a variable we will use to hold the inverse of certain matrices; we will update it multiple times with a for loop, 
        // so we initialize it here to ensure the memory persists
        let mut Binv                        =   HashMap::new();
        Binv.insert(
            0,
            VecOfVec::diagonal_matrix( unity, dimensions[0] ),
        );


        // initialize a list of bar ledgers for vector space p = 0

        let mut bars_anchored_at                        =   HashMap::new();
        let mut bars_anchored_at_0                =   Vec::with_capacity(dimensions[0]);   // a sequence of ledgers (one for each basis vector)
        for basis_vector_index_in_0 in 0 .. dimensions[0] {
            let new_ledger                  =   SingleBarBasisVectorIndexLedger{
                                                    leftmost_vertex:            0,
                                                    basis_vector_indices:       vec![ basis_vector_index_in_0 ],
                                                };
            bars_anchored_at_0.push( new_ledger );
        }
        bars_anchored_at.insert( 0, bars_anchored_at_0 );


        // compute a weight vector w such that we can add an interval module with support [a,d] to an interval module with support [b,d] iff w[a] ≤ w[b]
        let mut min                         =   0isize;
        let mut max                         =   0isize;
        let mut left_endpoint_weights       =   Vec::with_capacity( n_vertices );
        
        left_endpoint_weights.push( 0 ); // the weight of space 0 is 0
        for p in 0 .. n_arrows {
            if directions[ p ] {
                max += 1;
                left_endpoint_weights.push( max );  // in this case we think about adding a quiver of form v --> v to a quiver of form 0 --> u, which we can do by replacing v --> v with 0 --> v
            } else {
                min -= 1;
                left_endpoint_weights.push( min );  // in this case we think about adding a quiver of form 0 <-- v to a quiver of form u <-- u      
            }
        }
        
        // initialize a list of bars which have been fully constructed
        let mut completed_bars              =   Vec::new();


        // DIAGONALIZE THE RERPRESENTATION
        // ------------------------------------------------------------------------

        for p in 0 .. n_arrows {


            // get a handle on bars whose right endpoint is p
            let mut bars_anchored_at_p                  =   bars_anchored_at.remove( &p ).unwrap();            

            // get the direction of the edge that joins vertices p and p + 1 (true = forward, false = backward)
            let arrow_points_right          =   directions[p];

            // right now we are on iteration p; it is possible that on iteration p+1 we will need to have a copy of
            // matrix B[ p + 1 ]^{-1}. specifically, we might need to multiply A[ p + 1 ] on the left with B[ p + 1 ]^{-1}, 
            // if the (p + 1)th arrow points left.
            // in that case it will help to do some extra operations during the current 
            // iteration as preparation. to that end, we'll define a variable that records whether or not to perform
            // those preparatory steps.
            let we_will_need_a_copy_of_Bp1_inverse 
                                            =   ( ( p + 1 ) < directions.len() ) && ( ! directions[ p + 1 ] );

            

            // compute X := [A,B]_p
            let X                           =   if arrow_points_right {
                                                    // B[p] * A[p]
                                                    B[p].multiply_on_the_left_and_write_the_product_to_a_vec_of_vec( 
                                                        & self.matrices[p],
                                                        self.ring_operator.clone(),
                                                    )
                                                    .unwrap()
                                                } else {                                                    
                                                    // A[p] * B[p]^{-1}
                                                    // |  Bp_inverse standards for `B[p] inverse`
                                                    // |  We only use the label Bp1_inverse because we compute Bp1_inverse on iteration p, so the name makes sense at the time it is constructed.
                                                    // |  The reason we compute Bp1_inverse on iteration p is because it's a little easier there (there is some stuff involving permutations that can make it more complicated, later)
                                                    self
                                                        .matrices[p]
                                                        .multiply_on_the_left_and_write_the_product_to_a_vec_of_vec( 
                                                            & Binv.remove( & p ).unwrap(),
                                                            self.ring_operator.clone()
                                                        )
                                                        .unwrap()
                                                };   

            let X = X.matrix_algebra_packet(self.ring_operator());

            // compute U-match factorization  TM = DS
            // equivalently, M Sinv = Tinv D
            // this determines a partial matching of vectors
            // row_i( Tinv ) |--->  M[i,j] * row_j(Sinv)
            // for every (i,j) such that M[i,j] is nonzero
            let umatch                      =   Umatch::new(
                                                            & X, // matrix to decompose
                                                            ( 0 .. X.matrix_ref().number_of_rows() ).rev(), // row indices in reverse order
                                                        );

            // compute the matrix U obtained from Sinv by replacing row_j(Sinv) with M[i,j] * row_j(Sinv) for all (i,j) such that M[i,j] is nonzero.
            // store this matrix in a VecOfVec
            let get_row                                 =   | row_index: usize | -> Vec<(usize,RingOperator::Element)> {
                                                                let mut row    =   umatch
                                                                    .source_comb_inverse()
                                                                    .row( & row_index )
                                                                    .collect::<Vec<_>>();
                                                                if let Some( scalar )   =   umatch.generalized_matching_matrix_ref().coefficient_opt_for_column_index( & row_index ) {
                                                                    for ( _column_index, coefficient ) in row.iter_mut() {
                                                                        *coefficient      =   self.ring_operator.multiply( scalar.clone(),  coefficient.clone() );
                                                                    }
                                                                }
                                                                row
                                                            };
            let row_index_iterator                      =   if arrow_points_right{  0 .. dimensions[ p + 1 ]  }   else  {  0 .. dimensions[ p ]  };
            let Sinv_scaled                             =   row_index_iterator  
                                                                .map(|i| get_row(i) )
                                                                .collect::<Vec<_>>();
            let Sinv_scaled                             =   VecOfVec::new ( Sinv_scaled ).ok().unwrap();



            // store the matrix Tinv in a VecOfVec
            let get_row                                 =   | row_index: usize | -> Vec<(usize,RingOperator::Element)> {
                                                                umatch
                                                                    .target_comb_inverse()
                                                                    .row( & row_index )
                                                                    .collect::<Vec<_>>()
                                                            };
            let row_index_iterator                      =   if arrow_points_right{  0 .. dimensions[ p ]  }   else  {  0 .. dimensions[ p + 1 ]  };
            let Tinv                                    =   row_index_iterator  
                                                                .map(|i| get_row(i) )
                                                                .collect::<Vec<_>>();
            let Tinv                                    =   VecOfVec::new ( Tinv ).ok().unwrap();            
            


            // compute the corresponding basis update for spaces p and p + 1
            // the current basis for space p is given by the rows of B[p]. the basis update will replace B[p] with `basis_update_for_space_p * B[p]`.
            // (at this point in the procedure, the basis for space p + 1 is the standard basis of unit vectors, so the "basis update" becomes the new basis)            
            let ( basis_update_for_space_p, basis_update_for_space_p1 ) =   if arrow_points_right {
                                                                                ( Tinv, Sinv_scaled )
                                                                            } else {
                                                                                ( Sinv_scaled, Tinv )
                                                                            };





            // // compute the corresponding basis update for space p
            // // the current basis for space p is given by the rows of B[p]. the basis update will replace B[p] with `basis_update_for_space_p * B[p]`.
            // let get_row                     =   | i | {
            //                                         match arrow_points_right {
            //                                             // if the arrow points right, then the dimension of space p matches the number of rows of the matrix; this matches the size of the target COMB
            //                                             true    =>  umatch
            //                                                             .target_comb_inverse()
            //                                                             .row( i )
            //                                                             .collect::<Vec<_>>(),
            //                                             // if the arrow points left, then the dimension of space p matches the number of columns of the matrix; this matches the size of the source COMB                                                                        
            //                                             false   =>  umatch
            //                                                             .comb_domain_inv()
            //                                                             .row( i )
            //                                                             .collect::<Vec<_>>()
            //                                         }
            //                                     };
            // let basis_update_for_space_p    =   ( 0 .. dimensions[p] )
            //                                         .map(|i| get_row(i) )
            //                                         .collect::<Vec<_>>();
            // let basis_update_for_space_p    =   VecOfVec::new( basis_update_for_space_p );

            // // compute the corresponding basis update for space p + 1
            // // (at this point in the procedure, the basis for space p + 1 is the standard basis of unit vectors, so the "basis update" becomes the new basis)
            // let get_row                     =   | i | {
            //                                         match arrow_points_right {
            //                                             // if the arrow points right, then the dimension of space p matches the number of columns of the matrix; this matches the size of the source COMB
            //                                             true    =>  umatch
            //                                                             .comb_domain_inv()
            //                                                             .row( i )
            //                                                             .collect::<Vec<_>>(),
            //                                             // if the arrow points right, then the dimension of space p matches the number of rows of the matrix; this matches the size of the target COMB                                                                        
            //                                             false   =>  umatch
            //                                                             .target_comb_inverse()
            //                                                             .row( i )
            //                                                             .collect::<Vec<_>>()
            //                                         }
            //                                     };
            // let basis_update_for_space_p1   =   ( 0 .. dimensions[ p + 1 ] )
            //                                         .map(|i| get_row(i) )
            //                                         .collect::<Vec<_>>();
                   

            // UPDATE BASES
            // ---------------------------------------------------

            // iterate over every bar that passes through both p and p-1
            // for each bar, we will update any basis vectors that correspond to this bar in B[0] .. B[p-1]
            // (this includes B[p-1] but excludes B[p])
            for ( right_basis_index, bar_ledger ) in bars_anchored_at_p.iter().enumerate() {

                // the left endpoint of the bar
                let leftmost_vertex                       =   bar_ledger.leftmost_vertex();

                // skip ahead if the bar doesn't extend strictly to the left of p
                if leftmost_vertex == p {
                    continue
                }                

                // the coefficients used to combine a set of old basis vectors to make a new basis vector in vector space p
                let linear_combination_right            =   ( & basis_update_for_space_p )
                                                            .row( & right_basis_index )
                                                            .collect::<Vec<_>>();
                
                // for each vertex strictly to the left of p that is touched by the bar
                for precursor_vertex in leftmost_vertex .. p {

                    // translate the family of coefficients for vertex p into a family of coefficients for precursor_vertex
                    let linear_combination_left         =   linear_combination_right
                                                                .iter()
                                                                .filter_map( 
                                                                    |(n,a)|                                                         // for each (n,a)
                                                                    {                                                               
                                                                        bars_anchored_at_p[ *n ]                                     // get the bar, I, that passes through basis vector n in space p
                                                                            .vertex_to_basis_vector_index_opt( precursor_vertex )   // get the (index, i, of the) basis vector for I over `precursor_vertex`; if I does not have a basis vector over precursor_vertex, then omit this term
                                                                            .map( |i|  ( i, a.clone() )  )                          // place i into a tupe with the corresponding coefficient
                                                                    }
                                                                );

                    // the vector produced by that linear combination in the precursor space
                    let new_row                         =   linear_combination_left
                                                                .multiply_self_as_a_row_vector_with_matrix(
                                                                    B[ precursor_vertex ]
                                                                        .matrix_algebra_packet(self.ring_operator())
                                                                )
                                                                .collect::<Vec<_>>();    

                    // insert this new vector in the precursor basis
                    let left_basis_index                =   bar_ledger
                                                                .vertex_to_basis_vector_index_opt( precursor_vertex )
                                                                .unwrap();
                    B[ precursor_vertex ]
                        .replace_row_and_return_old(
                            left_basis_index,
                            new_row,
                        ).unwrap();                                                                                                                                         

                }

            }


            let Bp                                      =   basis_update_for_space_p
                                                                .multiply_on_the_left_and_write_the_product_to_a_vec_of_vec( 
                                                                    & B[p],
                                                                    self.ring_operator(),
                                                                ).unwrap();
            B[ p     ]                                  =   Bp;            
            B.push(
                // set  B[ p + 1 ] =   basis_update_for_space_p1;                
                basis_update_for_space_p1
            );

            // store a copy of B[p + 1]^{-1}    (only if we will need it later)
            if  we_will_need_a_copy_of_Bp1_inverse {
                
                // get a lazy copy of the inverse; this is just an oracle, not a vec-of-vec
                let Bp1_inv_lazy                        =   InverseUpperTriangularMatrix::new(
                                                                B[ p + 1 ].matrix_algebra_packet(self.ring_operator()),
                                                            );
                // get an iterator that runs over the rows in order
                let Bp1_inv_rows                        =   ( 0 .. dimensions[ p + 1] ).map( 
                                                                | i |
                                                                Bp1_inv_lazy.row( &i )
                                                            );
                // store the rows in a vec-of-vec; now we don't need the oracle any more
                let Bp1_inv                             =   VecOfVec::from_iterable_of_iterables( 
                                                                Bp1_inv_rows
                                                            ).ok().unwrap();
                Binv.insert(
                    p + 1,
                    Bp1_inv
                );
            }

            // GET FUNCTIONS TO EVALUATE PARTIAL MATCHINGS BETWEEN INDEX SETS
            // ---------------------------------------------------        


            let matching_matrix                         =   umatch.generalized_matching_matrix_ref();  

            let basis_vector_index_in_p1_to_basis_vector_index_in_p                 
                                                        =   | i: usize | -> Option< usize > {
                                                                if directions[ p ] {
                                                                    matching_matrix
                                                                        .row_index_for_column_index(& i )
                                                                } else {
                                                                    matching_matrix
                                                                        .column_index_for_row_index( & i )                                                                    
                                                                }
                                                            };  

            let basis_vector_in_p_is_matched_in_p1      =   | i: usize | ->  bool  {
                                                                if directions[ p ] {
                                                                    matching_matrix.has_a_match_for_row_index( & i )
                                                                } else {
                                                                    matching_matrix.has_a_match_for_column_index( & i )                                                                
                                                                }
                                                            };                                                                           


            // COLLECT A LIST OF BAR LEDGERS FOR VERTEX p + 1
            // ---------------------------------------------------
                                    
            let mut bars_anchored_at_p1                 =   Vec::with_capacity( dimensions[ p + 1] );
            
            for basis_vector_index_in_p1 in 0 .. dimensions[p+1]  {
                if let Some( basis_vector_index_in_p )  =   basis_vector_index_in_p1_to_basis_vector_index_in_p(    
                                                                basis_vector_index_in_p1      
                                                            ) 
                {
                    let mut ledger                      =   bars_anchored_at_p[ basis_vector_index_in_p ].clone();
                    ledger
                        .basis_vector_indices
                        .push( basis_vector_index_in_p1 );
                    bars_anchored_at_p1.push(ledger);
                } else {
                    let new_ledger                      =   SingleBarBasisVectorIndexLedger{ 
                                                                leftmost_vertex: p + 1,  
                                                                basis_vector_indices: vec![ basis_vector_index_in_p1 ] 
                                                            };
                    bars_anchored_at_p1.push( new_ledger );
                }
            }

            //  DRAIN THE BAR LEDGER FOR SPACE p
            // ---------------------------------------------------            
            for ( basis_vector_index_in_p, ledger ) in bars_anchored_at_p.drain(..).enumerate() {

                if  !   basis_vector_in_p_is_matched_in_p1( basis_vector_index_in_p ) {
                    completed_bars.push( ledger )
                }
                // if the bar does not end at p then we discard it, since we already made a copy and pushed it to the list of bars achored at p+1

            }

            // (1)  PERMUTE THE BAR LEDGERS FOR VERTEX p + 1, AND THE BASIS VECTORS FOR VERTEX p + 1
            // (2)  IF NECESSARY ALSO PERMUTE THE COLUMNS OF Binv[p+1] accordingly, so tthat B[p+1] * Binv[p+1] = identity
            // (3)  PLACE THE BAR LEDGERS FOR VERTEX p + 1 INTO `bars_anchored_at_p`
            // ---------------------------------------------------            

            // get a list of form 0 .. dimensions[p+1]
            let mut permutation             = (0 .. dimensions[p+1]).collect::<Vec<_>>();
            
            // permute the list according the the weight of the left endpoint of each index
            permutation.sort_by_key(
                |&i| 
                {
                    let leftmost_vertex                 =   bars_anchored_at_p1[i].leftmost_vertex.clone();
                    let weight                          =   left_endpoint_weights[ leftmost_vertex ];
                    -  weight   // !! NOTE: we sort in reverse order of weight because we want vectors in the smallest filtration space to appear at the end of the list
                }                                
            );

            // permute the list of ledgers
            let bars_anchored_at_p1_in_old_order                =   bars_anchored_at_p1.clone();
            bars_anchored_at_p1.clear();
            for basis_vector_index_in_p1 in permutation.iter() {
                bars_anchored_at_p1.push(  // here we are re-using the (empty) vector bars_anchored_at_p because it is pre-allocated
                    bars_anchored_at_p1_in_old_order[ *basis_vector_index_in_p1 ].clone()
                );
            }
            
            // update the ledgers to reflect the fact that we have permuted the basis
            for ( basis_vector_index_in_p1, ledger ) in bars_anchored_at_p1.iter_mut().enumerate() {
                let bar_length                          =   ledger.bar_length();
                ledger
                    .basis_vector_indices[ bar_length - 1 ] // this entry, ledger.basis_vector_indices[ bar_length - 1 ], is the basis vector index that our ledger stores for vector space p+1
                                                        =   basis_vector_index_in_p1;
            }


            // permute the basis vectors for space p+1 (represented by rows of Bp1)
            B[ p + 1 ]                                   =   B[ p + 1 ]
                                                                .permute_rows_out_of_place( 
                                                                    permutation.iter().cloned() 
                                                                );
            
            // if we need the inverse of B[ p + 1 ] for next iteration, then permute its columns accordingly
            // note that if we move row `i` of Bp1 to position `p`, then we should change the index of
            // column `i` in `Bp1_inverse` to `p`
            if we_will_need_a_copy_of_Bp1_inverse {

                let mut inverse_permutation             =   vec![ 0; permutation.len() ];
                for ( counter, n ) in permutation.iter().enumerate() {
                    inverse_permutation[ *n ]           =   counter;
                }

                Binv.entry(p+1).and_modify(| m | {
                    *m                                  =   (*m).reassign_column_indices_out_of_place( & inverse_permutation ).unwrap();
                });
            }

            // put bars anchored at p + 1 into our book-keeping device
            bars_anchored_at.insert( p+1, bars_anchored_at_p1 );

            // delete Binv[ p ], if it exists
            Binv.remove( & p ); 

        }

        // add in all the bars that terminate at the last vertex
        let bars_anchored_at_rightmost_vertex           =   bars_anchored_at.remove( & (n_vertices - 1) ).unwrap();
        completed_bars
            .extend( bars_anchored_at_rightmost_vertex );


        println!("Time to decompose representation: {:?}", start.elapsed()); // Output the duration                


        return  Ok(
                    Diagonalization::new( completed_bars, B  )
                )

    }











    /// Verify a diagonalization
    /// 
    /// This will check that
    /// - every basis vector in the diagonalization maps either to zero or to another basis vector
    /// - no two basis vectors map to the same basis vector
    /// - the information recorded for each bar in the barcode accurately tracks the chain of associated bases vectors in the representation
    pub fn validate_diagonalization(
            & self,
            diagonalization:    & Diagonalization< RingOperator::Element >,
        )
        ->
        Result<
            (),
            HashMap< &str, usize >
        >
        where 
            RingOperator::Element:    PartialEq
    {
        let directions                                  =   self.arrow_directions();
        let bases                                                   =   diagonalization.bases();
        let arrow_matrices                              =   self.arrow_matrices();

        // Check that the representation homomorphisms carry basis vectors to basis vectors in the prescribed manner
        for (bar_number, bar) in diagonalization.bars().iter().enumerate()  {
            for p in bar.supporting_edge_indices() {

                let vec_p                               =   (&bases[p]).row( & bar.vertex_to_basis_vector_index_opt(p).unwrap() );
                let vec_p1                              =   (&bases[p+1]).row( & bar.vertex_to_basis_vector_index_opt(p+1).unwrap() );
                

                if directions[ p ] {
                    let pushforward                     =   vec_p
                                                                .multiply_self_as_a_row_vector_with_matrix(
                                                                    arrow_matrices[p]
                                                                        .matrix_algebra_packet(self.ring_operator())
                                                                )
                                                                .collect::<Vec<_>>();
                    
                    if    !      pushforward.iter().cloned().eq( vec_p1 )  {
                        let mut err                     =   HashMap::new();
                        err.insert( "Diagonalization is not valid: basis vector a over vertex p does not map to vertex b over vertex p+1. This contradicts the information stored for bar number n, where n = ", bar_number );
                        err.insert("basis vector index over p",  bar.vertex_to_basis_vector_index_opt(p).unwrap()  );
                        err.insert("basis vector index over p+1",  bar.vertex_to_basis_vector_index_opt(p+1).unwrap()  );
                        err.insert("vertex p ",  p  );
                        err.insert("bar number", bar_number );
                        return Err( err )
                    }  
                } else {
                    let pushforward                     =   vec_p1
                                                                .multiply_self_as_a_row_vector_with_matrix( 
                                                                    arrow_matrices[p]
                                                                        .matrix_algebra_packet(self.ring_operator()),
                                                                );
                    if    !      pushforward.eq( vec_p )  {
                        let mut err                     =   HashMap::new();
                        err.insert( "Diagonalization is not valid: basis vector a over vertex p+1 does not map to vertex b over vertex p. This contradicts the information stored for bar number n, where n = ", bar_number );
                        err.insert("basis vector index a",  bar.vertex_to_basis_vector_index_opt(p+1).unwrap()  );
                        err.insert("basis vector index b",  bar.vertex_to_basis_vector_index_opt(p).unwrap()  );
                        err.insert("vertex p",  p  );  
                        err.insert("bar number", bar_number );                                              
                        return Err( err )
                    }                     
                }
            }
        }

        // Check that basis vectors at either end of an interval map to zero (if their arrows point outward)
        for ( bar_number, bar ) in diagonalization.bars().iter().enumerate() {
            let (a,b) = bar.inclusive_left_and_exclusive_right_endpoints();
            
            // in this case the bar starts at vertex a and there is an arrow ( a - 1 ) <-- a
            if ( a > 0 ) && (   ! directions[ a - 1 ]   ) { 
                let i                                   =   bar.vertex_to_basis_vector_index_opt( a ).unwrap();
                let vec_a                               =   (&bases[a]).row( &i );
                let pushforward                         =   vec_a
                                                                .multiply_self_as_a_row_vector_with_matrix(
                                                                    arrow_matrices[ a - 1 ]
                                                                        .matrix_algebra_packet(self.ring_operator()),
                                                                );
                if    !      pushforward.count() == 0  {
                    let mut err                     =   HashMap::new();
                    err.insert( "Diagonalization is not valid: basis vector i over vertex a does not map to zero. . This contradicts the information stored for bar number n, where n = ", bar_number );
                    err.insert("vertex a", a );
                    err.insert("basis vector index i",  i  );
                    err.insert("bar number",  bar_number  );                        
                    return Err( err )
                }      
            }


            // in this case the bar starts at vertex a and there is an arrow ( b - 1 ) --> b
            if ( b < self.number_of_vertices() ) &&  ( b > 0 )  &&  (  directions[ b - 1 ]   ) { 
                let i                                   =   bar.vertex_to_basis_vector_index_opt( b - 1 ).unwrap();
                let vec_a                               =   (&bases[b-1]).row( &i );
                let pushforward                         =   vec_a
                                                                .multiply_self_as_a_row_vector_with_matrix( 
                                                                    arrow_matrices[ b - 1 ]
                                                                        .matrix_algebra_packet(self.ring_operator())
                                                                );
                if    !      pushforward.count() == 0  {
                    let mut err                     =   HashMap::new();
                    err.insert( "Diagonalization is not valid: basis vector i over vertex b does not map to zero. . This contradicts the information stored for bar number n, where n = ", bar_number );
                    err.insert("vertex b", b );
                    err.insert("basis vector index i",  i  );
                    err.insert("bar number",  bar_number  );                        
                    return Err( err )
                }      
            }            


        }

        // Check that the list of recorded bars "covers" all the basis vectors
        let mut covered_vector_indices                      =   Vec::new();
        for vertex in 0 .. self.number_of_vertices() {
            
            covered_vector_indices.clear(); // re-set our list of covered indices
            covered_vector_indices.extend(
                diagonalization
                    .bars()
                    .iter()
                    .filter_map(
                        | bar |
                        bar.vertex_to_basis_vector_index_opt( vertex )
                    )
            );
            covered_vector_indices.sort();

            let ground_truth                            =   ( 0 .. self.dimension_of_space_over_vertex(vertex).unwrap()  );
            let ok                                      =   covered_vector_indices.iter().cloned().eq( ground_truth );
            


            if    !   ok {
                let mut err         =       HashMap::new();
                err.insert(  "Error: There should be a bijection between basis vectors over each vertex v and bars in the barcode that pass over vertex v: however this is not true for vertex  v = ", vertex  );
                return Err( err )
            }
        }

        return Ok(())
    }

















}










































#[cfg(test)]
mod tests {    

    use num::Integer;

    use crate::algebra::rings::types::field_prime_order::PrimeOrderField;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    /// Identity map zigzag
    /// 
    /// Tests the diagonalization method in the special case where all vector spaces have equal dimension and all maps are identity.
    #[test]
    fn test_zigzag_identity() {

        for n_arrows in vec![ 0, 1, 2, 3, 10, 20, 100 ] { // 0, 1, 2, 3, 10, 20, 100
            for uniform_dimension in vec![ 0, 1, 2, 3, 10, 20, 100 ] { // 0, 1, 2, 3, 10, 20, 100
                
                let n_vertices                                  =   n_arrows + 1;

                let unity                                       =   1;
                let matrices                                    =   ( 0 .. n_arrows )
                                                                        .map(|_i| VecOfVec::diagonal_matrix(unity, uniform_dimension) )
                                                                        .collect::<Vec<_>>();

                let arrow_directions                                =   ( 0 .. n_arrows )
                                                                        .map(|i|  i.is_even() )
                                                                        .collect::<Vec<_>>();
                let vector_space_dimensions                     =   vec![ uniform_dimension; n_vertices ];
        
                let ring_operator                                   =   PrimeOrderField::new(7);
        
                let quiver_representation                       =   QuiverReprsentation{
        
                                                                        arrow_directions,
                                                                        matrices,
                                                                        vector_space_dimensions,
                                                                        ring_operator,
                                                                    };                
        
                let diagonalization                             =   quiver_representation.diagonalize().unwrap(); 

                // println!("diagonalization: {:#?}", & diagonalization );
        
                let result                                      =   quiver_representation.validate_diagonalization( & diagonalization );
                // println!("{:?}", result );
                assert!( result.is_ok() )
            }
        }
    }













    /// Single fixed matrix diagonstic
    /// 
    /// Here we just investigate a single fixed matrix.
    #[test]
    fn test_specific_matrix_2_x_2() {


        let modulus                                     =   2;
        let ring_operator                                   =   PrimeOrderField::new( modulus );




        for n_arrows in vec![ 1, ] { // 0, 1, 2, 3, 10, 20, 100
            for uniform_dimension in vec![ 2,  ] { // 0, 1, 2, 3, 10, 20, 100
                
                let n_vertices                                  =   n_arrows + 1;

                let matrices                                    =   ( 0 .. n_arrows )
                                                                        .map(
                                                                            |_i| 
                                                                            VecOfVec::random_mod_p_with_density(
                                                                                uniform_dimension, // num_indices_major: usize,
                                                                                uniform_dimension, // num_indices_minor: usize,
                                                                                0.5, // approximate_density: f64,
                                                                                2, // modulus: usize,
                                                                                false, // allow_nonstructural_zero: bool
                                                                            ) 
                                                                        )
                                                                        .collect::<Vec<_>>();

                let structure_map                       =   VecOfVec::new(
                                                                vec![
                                                                    vec![ (0,1),  ( 1,1) ],
                                                                    vec![ (0,1),  ( 1,1) ],
                                                                ]
                                                            ).ok().unwrap();
                let matrices                            =   vec![ structure_map ];

                let arrow_directions                                =   ( 0 .. n_arrows )
                                                                        .map(|i|  i.is_even() )
                                                                        .collect::<Vec<_>>();
                let vector_space_dimensions                     =   vec![ uniform_dimension; n_vertices ];
        
                let quiver_representation                       =   QuiverReprsentation{
        
                                                                        arrow_directions,
                                                                        matrices,
                                                                        vector_space_dimensions,
                                                                        ring_operator: ring_operator.clone(),
                                                                    };                
        
                let diagonalization                             =   quiver_representation.diagonalize().unwrap(); 

                let basis_0                             =   & diagonalization.bases()[0];
                let basis_1                             =   & diagonalization.bases()[1];
                let basis_1_inverse                     =   & basis_1.generalized_inverse( ring_operator.clone(), uniform_dimension );
                let structure_map                       =   & quiver_representation.arrow_matrices()[0];

                
                let matching                            =   basis_0
                                                                .multiply_on_the_left_and_write_the_product_to_a_vec_of_vec(
                                                                    structure_map,
                                                                ring_operator.clone()
                                                            )
                                                            .unwrap()
                                                            .multiply_on_the_left_and_write_the_product_to_a_vec_of_vec(
                                                                basis_1_inverse,
                                                                ring_operator.clone()
                                                            )
                                                            .unwrap();
                println!("matching: {:#?}", & matching );

                // -------------------------------------
                // check generalized inverse


                let x                                   =   basis_1.multiply_on_the_left_and_write_the_product_to_a_vec_of_vec(basis_1_inverse, ring_operator.clone() );
                println!("basis_1 times basis_1_inverse: {:#?}", & x );                


                // -------------------------------------                
                // check bases

                println!("bases: {:#?}: ", diagonalization.bases() );         


                // -------------------------------------                
                // print bars

                println!("BARS: ==================  {:#?}", diagonalization.bars() );                       



                // -------------------------------------                
                // check umatch decomposition



                // -------------------------------------   

                // println!("diagonalization: {:#?}", & diagonalization );
        
                let result                                      =   quiver_representation.validate_diagonalization( & diagonalization );
                
                println!("{:#?}", result );
                // assert!( result.is_ok() )
            }
        }
    }




















    /// Single fixed matrix diagonstic
    /// 
    /// Here we just investigate a single fixed matrix.
    #[test]
    fn test_specific_matrix_3_x_3() {


        let modulus                                     =   2;
        let ring_operator                                   =   PrimeOrderField::new( modulus );
                
        let n_arrows                                            =   1;
        let n_vertices                                  =   n_arrows + 1;
        let uniform_dimension                                   =   3;

        let matrices                                    =   ( 0 .. n_arrows )
                                                                .map(
                                                                    |_i| 
                                                                    VecOfVec::random_mod_p_with_density(
                                                                        uniform_dimension, // num_indices_major: usize,
                                                                        uniform_dimension, // num_indices_minor: usize,
                                                                        0.5, // approximate_density: f64,
                                                                        2, // modulus: usize,
                                                                        false, // allow_nonstructural_zero: bool
                                                                    ) 
                                                                )
                                                                .collect::<Vec<_>>();
                                                        

        let structure_map                       =   VecOfVec::new(
                                                        vec![
                                                            vec![ (0,1),          (2,1) ],
                                                            vec![         (1,1),        ],
                                                            vec![         (1,1),        ],                                                                    
                                                        ]                                                                
                                                    ).ok().unwrap();
        let matrices                            =   vec![ structure_map ];

        let arrow_directions                                =   ( 0 .. n_arrows )
                                                                .map(|i|  i.is_even() )
                                                                .collect::<Vec<_>>();
        let vector_space_dimensions                     =   vec![ uniform_dimension; n_vertices ];

        let quiver_representation                       =   QuiverReprsentation{

                                                                arrow_directions,
                                                                matrices,
                                                                vector_space_dimensions,
                                                                ring_operator: ring_operator.clone(),
                                                            };                

        let diagonalization                             =   quiver_representation.diagonalize().unwrap(); 

        let basis_0                             =   & diagonalization.bases()[0];
        let basis_1                             =   & diagonalization.bases()[1];
        let basis_1_inverse                     =   & basis_1.generalized_inverse( ring_operator.clone(), uniform_dimension );
        let structure_map                       =   & quiver_representation.arrow_matrices()[0];

        
        let matching                            =   basis_0
                                                        .multiply_on_the_left_and_write_the_product_to_a_vec_of_vec(
                                                            structure_map,
                                                        ring_operator.clone()
                                                    )
                                                    .unwrap()
                                                    .multiply_on_the_left_and_write_the_product_to_a_vec_of_vec(
                                                        basis_1_inverse,
                                                        ring_operator.clone()
                                                    )
                                                    .unwrap();
        println!("matching: {:#?}", & matching );

        // -------------------------------------
        // check generalized inverse


        let x                                   =   basis_1.multiply_on_the_left_and_write_the_product_to_a_vec_of_vec(basis_1_inverse, ring_operator.clone() );
        println!("basis_1 times basis_1_inverse: {:#?}", & x );                


        // -------------------------------------                
        // check bases

        println!("bases: {:#?}: ", diagonalization.bases() );         


        // -------------------------------------                
        // print bars

        println!("BARS: ==================  {:#?}", diagonalization.bars() );                       



        // -------------------------------------                
        // check umatch decomposition



        // -------------------------------------   

        // println!("diagonalization: {:#?}", & diagonalization );

        let result                                      =   quiver_representation.validate_diagonalization( & diagonalization );
        
        println!("{:#?}", result );
        assert!( result.is_ok() )
    }






















    /// Random matrices uniform dimension
    /// 
    /// Tests the diagonalization method in the special case where all vector spaces have equal dimension and all maps are identity.
    #[test]
    fn test_random_mod_p_matrices() {


        let modulus                                     =   89;
        let ring_operator                                   =   PrimeOrderField::new( modulus );




        for n_arrows in vec![ 0, 1, 2, 3, 10, 20, ] { // 0, 1, 2, 3, 10, 20, 30, 40,
            for uniform_dimension in vec![  0, 1, 2, 3, 10, 20,  ] { // 0, 1, 2, 3, 10, 20, 30, 40,
                
                let n_vertices                                  =   n_arrows + 1;

                let matrices                                    =   ( 0 .. n_arrows )
                                                                        .map(
                                                                            |_i| 
                                                                            VecOfVec::random_mod_p_with_density(
                                                                                uniform_dimension, // num_indices_major: usize,
                                                                                uniform_dimension, // num_indices_minor: usize,
                                                                                0.3, // approximate_density: f64,
                                                                                modulus, // modulus: usize,
                                                                                false, // allow_nonstructural_zero: bool
                                                                            ) 
                                                                        )
                                                                        .collect::<Vec<_>>();

                let arrow_directions                                =   ( 0 .. n_arrows )
                                                                        .map(|i|  
                                                                            i.is_even() 
                                                                            // false
                                                                        )
                                                                        .collect::<Vec<_>>();
                let vector_space_dimensions                     =   vec![ uniform_dimension; n_vertices ];
        
                let quiver_representation                       =   QuiverReprsentation{
        
                                                                        arrow_directions,
                                                                        matrices,
                                                                        vector_space_dimensions,
                                                                        ring_operator: ring_operator.clone(),
                                                                    };                
        
                let diagonalization                             =   quiver_representation.diagonalize().unwrap(); 

        
                let result                                      =   quiver_representation.validate_diagonalization( & diagonalization );

                
                // println!("{:#?}", result );
                assert!( result.is_ok() )
            }
        }
    }    

















    /// Random matrices increasing dimension
    /// 
    /// Tests the diagonalization method in the special case where vector space p has dimension p
    #[test]
    fn test_random_mod_p_matrices_increasing_dimension() {


        let modulus                                     =   89;
        let ring_operator                                   =   PrimeOrderField::new( modulus );




        for n_arrows in vec![ 0, 1, 2, 3, 10, 20, ] { // 0, 1, 2, 3, 10, 20, 30, 40,
                
            let n_vertices                                  =   n_arrows + 1;

            let arrow_directions                                =   ( 0 .. n_arrows )
                                                                    .map(|i|  
                                                                        i.is_even() 
                                                                        // false
                                                                    )
                                                                    .collect::<Vec<_>>();
            let vector_space_dimensions                     =   (0 .. n_vertices).collect::<Vec<_>>();


            let matrices                                    =   ( 0 .. n_arrows )
                                                                    .map(
                                                                        | matrix_number|
                                                                        {
                                                                            let ( n_rows, n_cols )      =   if arrow_directions[ matrix_number ] { 
                                                                                                                ( vector_space_dimensions[ matrix_number ], vector_space_dimensions[ matrix_number + 1 ] )
                                                                                                            } else {
                                                                                                                ( vector_space_dimensions[ matrix_number + 1 ], vector_space_dimensions[ matrix_number ] )
                                                                                                            };
                                                                            VecOfVec::random_mod_p_with_density(
                                                                                n_rows, // num_indices_major: usize,
                                                                                n_cols, // num_indices_minor: usize,
                                                                                0.3, // approximate_density: f64,
                                                                                modulus, // modulus: usize,
                                                                                false, // allow_nonstructural_zero: bool
                                                                            ) 
                                                                        } 
                                                                    )
                                                                    .collect::<Vec<_>>();

    
            let quiver_representation                       =   QuiverReprsentation{
    
                                                                    arrow_directions,
                                                                    matrices,
                                                                    vector_space_dimensions,
                                                                    ring_operator: ring_operator.clone(),
                                                                };                
    
            let diagonalization                             =   quiver_representation.diagonalize().unwrap(); 

    
            let result                                      =   quiver_representation.validate_diagonalization( & diagonalization );

            
            // println!("{:#?}", result );
            assert!( result.is_ok() )
        }
    }       




    /// Random matrices decreasing dimension
    /// 
    /// Tests the diagonalization method in the special case where vector space p has dimension p
    #[test]
    fn test_random_mod_p_matrices_decreasing_dimension() {


        let modulus                                     =   89;
        let ring_operator                                   =   PrimeOrderField::new( modulus );




        for n_arrows in vec![ 0, 1, 2, 3, 10, 20, ] { // 0, 1, 2, 3, 10, 20, 30, 40,
                
            let n_vertices                                  =   n_arrows + 1;

            let arrow_directions                                =   ( 0 .. n_arrows )
                                                                    .map(|i|  
                                                                        i.is_even() 
                                                                    )
                                                                    .collect::<Vec<_>>();
            let mut vector_space_dimensions                     =   (0 .. n_vertices).collect::<Vec<_>>();
            vector_space_dimensions.reverse();


            let matrices                                    =   ( 0 .. n_arrows )
                                                                    .map(
                                                                        | matrix_number|
                                                                        {
                                                                            let ( n_rows, n_cols )      =   if arrow_directions[ matrix_number ] { 
                                                                                                                ( vector_space_dimensions[ matrix_number ], vector_space_dimensions[ matrix_number + 1 ] )
                                                                                                            } else {
                                                                                                                ( vector_space_dimensions[ matrix_number + 1 ], vector_space_dimensions[ matrix_number ] )
                                                                                                            };
                                                                            VecOfVec::random_mod_p_with_density(
                                                                                n_rows, // num_indices_major: usize,
                                                                                n_cols, // num_indices_minor: usize,
                                                                                0.3, // approximate_density: f64,
                                                                                modulus, // modulus: usize,
                                                                                false, // allow_nonstructural_zero: bool
                                                                            ) 
                                                                        } 
                                                                    )
                                                                    .collect::<Vec<_>>();

    
            let quiver_representation                       =   QuiverReprsentation{
    
                                                                    arrow_directions,
                                                                    matrices,
                                                                    vector_space_dimensions,
                                                                    ring_operator: ring_operator.clone(),
                                                                };                
    
            let diagonalization                             =   quiver_representation.diagonalize().unwrap(); 

    
            let result                                      =   quiver_representation.validate_diagonalization( & diagonalization );

            
            // println!("{:#?}", result );
            assert!( result.is_ok() )
        }
    }            






}
