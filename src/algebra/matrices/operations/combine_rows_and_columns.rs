//! Type aliases for linear combinations of rows and columns
//! 
//! This module defines some [type aliases](https://doc.rust-lang.org/beta/reference/items/type-aliases.html) to help code from getting too verbose.
//! 


use crate::{algebra::{matrices::query::MatrixAlgebra, vectors::operations::{Scale, Simplify}}, utilities::{iterators::merge::hit::IteratorsMergedInSortedOrder, order::ReverseOrder}};


/// Represents a linear combination of rows of a matrix
/// 
/// This is a [type alias](https://doc.rust-lang.org/beta/reference/items/type-aliases.html) for a
/// data structure that represents a linear combination of rows of a matrix of type `Matrix`
/// 
/// # Example
/// 
/// Here we construct a linear combination of rows, `u`, which has type [LinearCombinationOfRows].
/// 
/// ```
/// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
/// use oat_rust::algebra::rings::types::native::FieldFloat64;
/// use oat_rust::algebra::vectors::operations::VectorOperations;
/// 
/// // define inputs
/// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
/// let data            =   vec![   
///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
///                             vec![          (1, 1.), (2, 1.)     ],
///                             vec![                   (2, 1.)     ],  
///                         ];
/// let data            =   VecOfVec::new( data ).ok().unwrap();
/// let matrix          =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &data );
/// 
/// // multiply
/// let u               =   v.multiply_self_as_a_row_vector_with_matrix( matrix );
///
/// // verify
/// assert!( u.eq(  vec![ (0, 1.), (1, 2.), (2, 3.) ]  ) );
/// ```
pub type LinearCombinationOfRows< 
        Matrix: MatrixAlgebra 
    > 
    =   Simplify<
            IteratorsMergedInSortedOrder<
                Scale< Matrix::Row, Matrix::RingOperator, >,
                Matrix::OrderOperatorForRowEntries,
            >,
            Matrix::RingOperator,
        >;

/// Represents a linear combination of rows of a matrix, with entries appearing in **REVERSE** order
///         
/// This is a [type alias](https://doc.rust-lang.org/beta/reference/items/type-aliases.html) for a
/// data structure that represents a linear combination of rows of a matrix of type `Matrix`;
/// the entries of the 
/// 
/// # Example
/// 
/// Here we construct a linear combination of (reversed) rows, `u`, which has type [LinearCombinationOfRowsReverse].
/// 
/// ```
/// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
/// use oat_rust::algebra::rings::types::native::FieldFloat64;
/// use oat_rust::algebra::vectors::operations::VectorOperations;
/// use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};
/// 
/// // define inputs
/// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
/// let data            =   vec![   
///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
///                             vec![          (1, 1.), (2, 1.)     ],
///                             vec![                   (2, 1.)     ],  
///                         ];
/// let data            =   VecOfVec::new( data ).ok().unwrap();
/// let matrix          =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &data );
/// 
/// // multiply
/// let u               =   v.multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order( matrix );
///
/// // verify
/// assert!( u.eq(  vec![ (2, 3.), (1, 2.), (0, 1.) ]  ) );
/// ```
pub type LinearCombinationOfRowsReverse< 
        Matrix: MatrixAlgebra 
    > 
    =   Simplify<
            IteratorsMergedInSortedOrder<
                Scale< Matrix::RowReverse, Matrix::RingOperator >,
                ReverseOrder< Matrix::OrderOperatorForRowEntries >,
            >,
            Matrix::RingOperator,
        >;        


/// Represents a linear combination of columns of a matrix
/// 
/// This is a [type alias](https://doc.rust-lang.org/beta/reference/items/type-aliases.html) for a
/// data structure that represents a linear combination of columns of a matrix of type `Matrix`
/// 
/// # Example
/// 
/// Here we construct a linear combination of columns, `u`, which has type [LinearCombinationOfColumns].
/// 
/// ```
/// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
/// use oat_rust::algebra::rings::types::native::FieldFloat64;
/// use oat_rust::algebra::vectors::operations::VectorOperations;
/// use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};
/// 
/// // define inputs
/// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
/// let data            =   vec![   
///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
///                             vec![          (1, 1.), (2, 1.)     ],
///                             vec![                   (2, 1.)     ],  
///                         ];
/// let data            =   VecOfVec::new( data ).ok().unwrap();
/// let matrix          =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &data );
/// 
/// // multiply
/// let u               =   v.multiply_self_as_a_column_vector_with_matrix( matrix );
///
/// // verify
/// assert!( u.eq(  vec![ (0, 3.), (1, 2.), (2, 1.) ]  ) );
/// ```
pub type LinearCombinationOfColumns< 
        Matrix: MatrixAlgebra 
    > 
    =   Simplify<
            IteratorsMergedInSortedOrder<
                Scale< Matrix::Column, Matrix::RingOperator >,
                Matrix::OrderOperatorForColumnEntries,
            >,
            Matrix::RingOperator,
        >;

/// Represents a linear combination of columns of a matrix, with entries appearing in **REVERSE** order
///         
/// This is a [type alias](https://doc.rust-lang.org/beta/reference/items/type-aliases.html) for a
/// data structure that represents a linear combination of columns of a matrix of type `Matrix`;
/// the entries of the resulting vector in **REVERSE** order.
/// 
/// # Example
/// 
/// Here we construct a linear combination of (reversed) columns, `u`, which has type [LinearCombinationOfRowsReverse].
/// 
/// ```
/// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
/// use oat_rust::algebra::rings::types::native::FieldFloat64;
/// use oat_rust::algebra::vectors::operations::VectorOperations;
/// use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};
/// 
/// // define inputs
/// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
/// let data            =   vec![   
///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
///                             vec![          (1, 1.), (2, 1.)     ],
///                             vec![                   (2, 1.)     ],  
///                         ];
/// let data            =   VecOfVec::new( data ).ok().unwrap();
/// let matrix          =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &data );
/// 
/// // multiply
/// let u               =   v.multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order( matrix );
///
/// // verify
/// assert!( u.eq(  vec![ (2, 1.), (1, 2.), (0, 3.) ]  ) );
/// ```
pub type LinearCombinationOfColumnsReverse< 
        Matrix: MatrixAlgebra 
    > 
    =   Simplify<
            IteratorsMergedInSortedOrder<
                Scale< Matrix::ColumnReverse, Matrix::RingOperator >,
                ReverseOrder< Matrix::OrderOperatorForColumnEntries >,
            >,
            Matrix::RingOperator,
        >;  