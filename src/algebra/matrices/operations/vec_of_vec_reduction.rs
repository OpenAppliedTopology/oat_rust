//! Only valid for `vec_of_vec` matrices (not of general iterest)

use crate::algebra::rings::traits::DivisionRingOperations;
use crate::algebra::vectors::operations::VectorOperations;
use crate::algebra::vectors::entries::KeyValGet;
use std::collections::HashMap;
use std::fmt::Debug;


//  CLEAR A VECTOR
//  --------------

/// Clear an entry of the `clearee` sparse vector by adding a scalar multiple of the `clearor`.
/// 
/// # Important notes
/// 
/// - Does nothing if the entry to be cleared in clearee is either a structural zero or 
///   structurally zero but equal to 0.
/// - If a nonzero multiple of `clearor` is required, then zero entries will be dropped. 
/// - Assumes that `pivot_entry` is and entry of the clearor vector.  
/// - Panicks if `pivot_entry` is zero but the corresponding entry of `clearee` is nonzero.
///
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
/// use oat_rust::algebra::matrices::operations::vec_of_vec_reduction::clear_if_in;;
///
/// let     clearor         =   vec![ (0, 1.), (1, 1.)          ];
/// let mut clearee         =   vec![          (1, 1.), (2, 1.) ];
/// let mut buffer          =   Vec::new();
/// let     pivot_entry     =   (1, 1.);
/// let     ring_operator   =   RingOperatorForNativeRustNumberType::<f64>::new();
///
/// clear_if_in(
///     &       clearor,
///     &mut    clearee,
///     &mut    buffer,
///     &       pivot_entry,
///             ring_operator.clone()
/// );
///
/// assert_eq!( &clearee, &vec![(0, -1.), (2, 1.)]);
///
/// let     clearor     =   vec![ (0, 1.), (1, 1.), (2, 0.) ];
/// let mut clearee     =   vec![ (0, 0.), (1, 0.), (2, 1.) ];
///
/// clear_if_in(
///     &       clearor,
///     &mut    clearee,
///     &mut    buffer,
///     &       pivot_entry,
///             ring_operator.clone()
/// );
///
/// assert_eq!( &clearee, &vec![ (0, 0.), (1, 0.), (2, 1.) ]);
///
/// let     clearor     =   vec![ (0, 1.), (1, 1.)          ];
/// let mut clearee     =   vec![                   (2, 1.) ];
///
/// clear_if_in(
///     &       clearor,
///     &mut    clearee,
///     &mut    buffer,
///     &       pivot_entry,
///     ring_operator.clone()
/// );
///
/// assert_eq!( &clearee, &vec![ (2, 1.) ]);
/// ```
/// 
pub fn  clear_if_in< Key, Val, RingOperator > (
    clearor:        &    Vec< (Key, Val) >,
    clearee:        &mut Vec< (Key, Val) >,
    buffer:         &mut Vec< (Key, Val) >,
    pivot_entry:    &         (Key, Val),
    ring_operator:                RingOperator 
)
where   RingOperator: DivisionRingOperations< Element = Val > + Clone,
        Key: Clone + Debug + PartialEq + PartialOrd,
        Val: Clone + Debug +PartialOrd

{
    let entry_to_clear_opt  =   clearee
                                .iter()
                                .find( |&x| x.key() == pivot_entry.key() );

    if let Some(entry_to_clear) = entry_to_clear_opt 
    {
        if ring_operator.is_0( entry_to_clear.val()) { return }              // short circuit if the entry to be cleared is zero

        let scalar          =   ring_operator.divide( 
                                    ring_operator.negate( entry_to_clear.val() ),
                                    pivot_entry.val()
                                );

        let merged          =   itertools::merge(                   // merge iterators, preserving
                                    clearee.iter().cloned(),
                                    clearor
                                        .iter()
                                        .cloned()
                                        .scale_by( scalar, ring_operator.clone() )
                                )
                                .peekable()                         // make peekable (necessary to gather coefficients)
                                .gather( ring_operator.clone() )             // gather coefficients
                                .drop_zeros( ring_operator );                // drop zeros
        buffer.clear();
        buffer.extend( merged );

        clearee.clear();
        clearee.append( buffer);    // note that buffer is already declared to be a mutable reference
    }
}



/// Reduce the specified columns of the `clearee_matrix` using the `clearor` column.
/// 
/// This is achieved by applying the function [`clear_if_in`] to each of the columns
/// specified.  See the documentation for that function for important notes about 
/// how the clearing is performed.
pub fn clear_cols< RingOperator, Key, Val, IndexIter: IntoIterator< Item = usize > >(
    clearor:        &    Vec< (Key, Val) >,
    clearee_matrix: &mut Vec< Vec< (Key, Val) > >,
    col_ind_2clear:      IndexIter,      
    pivot_entry:    &         (Key, Val),
    ring_operator:                RingOperator,     
    )
    where   RingOperator: DivisionRingOperations< Element = Val > + Clone,
            Key: Clone + Debug + PartialEq + PartialOrd,
            Val: Clone + Debug +PartialOrd    
{
    let mut buffer  =   Vec::new();
    for col_ind in col_ind_2clear {
        clear_if_in(
            clearor, 
            &mut clearee_matrix[ col_ind ],
            &mut buffer,
            pivot_entry, 
            ring_operator.clone()
        )
    }
}




//  RIGHT REDUCE
//  ------------

type Key = usize;

/// Compute the right-reduced matrix of input `matrix`
/// 
/// Important assumptions:
///     * all zero entries are also structurally nonzero.
///     * the entries in each column are SORTED
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
/// use oat_rust::algebra::matrices::operations::vec_of_vec_reduction::right_reduce;
/// use std::iter::FromIterator;
///
/// /// Input matrix
/// let mut matrix      =   vec![
///                             vec![                   (2, 1.), (3, -1.)   ],
///                             vec![                   (2, 1.), (3, 1.)    ],                                    
///                             vec![          (1, 1.), (2, 1.)             ],
///                             vec![ (0, 1.), (1, 1.)                      ],
///                             vec![ (0, 1.),                              ],
///                         ];
///
/// /// Correctly reduced matrix
/// let reduced_correct =   vec![
///                             vec![                   (2, 1.), (3, -1.)   ],
///                             vec![                   (2, 2.),            ],                                    
///                             vec![          (1, 1.)                      ],
///                             vec![ (0, 1.),                              ],
///                             vec![                                       ],
///                         ];                                                                        
///
/// /// Compute the actual matrix and (sorted sequence of) pivot pairs
/// let hash = right_reduce( 
///                 &mut matrix, 
///                 RingOperatorForNativeRustNumberType::<f64>::new() 
///             );            
/// let mut pivot_pairs = Vec::from_iter( hash );        
/// pivot_pairs.sort();
///
/// // Check
/// assert_eq!( pivot_pairs, vec![ (0,3), (1,2), (2,1), (3,0)] );        
/// assert_eq!( reduced_correct, matrix );   
/// ```

pub fn right_reduce 
    < Val, RingOperator > 
    
    ( 
    matrix:     &mut Vec< Vec< (Key, Val) > >,
    ring_operator:       RingOperator
    )
    ->
    HashMap::<Key, Key>

    where   RingOperator: DivisionRingOperations< Element = Val > + Clone,
            Key: Clone + Debug + PartialEq + PartialOrd + Eq + std::hash::Hash,
            Val: Clone + Debug +PartialOrd

{
    let mut pivot_hash        =   HashMap::< Key, Key >::new();
    let mut buffer          =   Vec::new();

    for clearee_count in 0..matrix.len() {

        let mut clearee     =   matrix[ clearee_count ].clone();
        
        //  REDUCE THE CLEAREE
        while let Some( clearee_entry ) = clearee.last(){
            if let Some( clearor_index ) = pivot_hash.get( &clearee_entry.key() ) {

                let  clearor        =   matrix[ *clearor_index ].clone();
                let  clearor_entry  =   clearor.last().unwrap();
                let  scalar         =   ring_operator.divide( 
                                            ring_operator.negate(clearee_entry.val()),
                                            clearor_entry.val()
                                        );                                              

                let merged          =   itertools::merge(                   // merge iterators, preserving
                                            clearee.iter().cloned(),
                                            clearor
                                                .iter()
                                                .cloned()
                                                .scale_by( scalar, ring_operator.clone() )
                                        )
                                        .peekable()                         // make peekable (necessary to gather coefficients)
                                        .gather( ring_operator.clone() )             // gather coefficients
                                        .drop_zeros( ring_operator.clone() );        // drop zeros

                buffer.clear();
                buffer.extend( merged );
        
                clearee.clear();
                clearee.append( &mut buffer);
            } else {
                break;
            }
        }

        //  UPDATE MATRIX + HASHMAP

        matrix[ clearee_count ].clear();                             // clear this column's slot in the matrix
        if let Some( pivot_entry ) = clearee.last() {
            pivot_hash.insert( pivot_entry.key(), clearee_count );      // update hashmap
            matrix[ clearee_count ].append( &mut clearee );          // write in the nonzero reduced column
        } 
    }

    pivot_hash
}







//  ---------------------------------------------------------------------------
//  TESTS
//  ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use crate::algebra::matrices::operations::vec_of_vec_reduction::*;
    use crate::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
    use std::iter::FromIterator;

    #[test]
    fn test()
    {

        // Input matrix
        let mut matrix      =   vec![
                                    vec![                   (2, 1.), (3,-1.)   ],
                                    vec![                   (2,-1.), (3, 2.)    ],                                    
                                    vec![          (1, 1.), (2, 1.)             ],
                                    vec![ (0, 1.), (1, 1.)                      ],
                                    vec![ (0, 1.),                              ],
     

  
                                ];

        // Correctly reduced matrix
        let reduced_correct =   vec![
                                    vec![                   (2, 1.), (3, -1.)   ],
                                    vec![                   (2, 1.),            ],                                    
                                    vec![          (1, 1.)                      ],
                                    vec![ (0, 1.),                              ],
                                    vec![                                       ],
                                ];                                                                        

        // Compute the actual matrix and (sorted sequence of) pivot pairs
        let hash = right_reduce( 
                        &mut matrix, 
                        RingOperatorForNativeRustNumberType::<f64>::new() 
                    );            
        let mut pivot_pairs = Vec::from_iter( hash );        
        pivot_pairs.sort();

        // Check
        assert_eq!( pivot_pairs, vec![ (0,3), (1,2), (2,1), (3,0)] );        
        assert_eq!( reduced_correct, matrix );                
    }     

    #[test]
    fn test_clear_cols()
    {
        let matrix          =   vec![
                                    vec![   (0, 1.),    (1, 1.),    (2, 1.)     ],
                                    vec![               (1, 1.),    (2, 1.)     ],
                                    vec![   (0, 1.),    (1, 1.),                ],
                                    vec![   (0, 1.),                (2, 1.)     ],                            
                                    vec![   (0, 1.),    (1, 0.),    (2, 1.)     ],
                                    vec![   (0, 1.),    (1, 1.),    (2, 1.)     ],                                       
                                    vec![   (0, 1.),    (1, 2.),    (2, 1.)     ],                                                                                            
                                    vec![                                       ],
                                    vec![   (0, 1.),    (1, 1.),    (2, 1.)     ],                            
                                ];

        let mut clearee_matrix  =   matrix;                                
        
        let clearor         =       vec![   (0, 1.),    (1, 1.),                ];

        let pivot_entry     =   *clearor.last().unwrap();

        let col_ind_2clear  =   1..7;
        
        clear_cols(
            &       clearor,
            &mut    clearee_matrix,
                    col_ind_2clear,
            &       pivot_entry,
                    RingOperatorForNativeRustNumberType::<f64>::new() ,
        );  
        
        let target_matrix   =   vec![
                                    vec![   (0, 1.),    (1, 1.),    (2, 1.)     ],
                                    vec![   (0, -1.),               (2, 1.)     ],
                                    vec![                                       ],
                                    vec![   (0, 1.),                (2, 1.)     ],                        
                                    vec![   (0, 1.),    (1, 0.),    (2, 1.)     ],                                                        
                                    vec![                           (2, 1.)     ],    
                                    vec![   (0, -1.),               (2, 1.)     ],                                                                                                                
                                    vec![                                       ],
                                    vec![   (0, 1.),    (1, 1.),    (2, 1.)     ],                            
                                ];
        assert_eq!( clearee_matrix, target_matrix );
    
    }

}
