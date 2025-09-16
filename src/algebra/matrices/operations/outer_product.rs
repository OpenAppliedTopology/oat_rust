use std::collections::HashMap;

use crate::algebra::{matrices::query::{MatrixAlgebra, MatrixOracle}, rings::traits::SemiringOperations, vectors::entries::KeyValGet};






pub trait OuterProduct {

    fn outer_product_row( & self, index: & Self::RowIndex ) 
        -> HashMap< Self::RowIndex, Self::Coefficient > 

        where
            Self:                   MatrixOracle + MatrixAlgebra + Sized,
            Self::RowIndex:         std::hash::Hash,
    
    {

        let mut result = HashMap::new();
        let ring_operator =  self.ring_operator();

        let row = self.row( index );
        for row_entry in row {
            let column_index = row_entry.key();
            let coefficient = row_entry.val();
            for column_entry in self.column( & column_index ) {
                let row_index = column_entry.key();
                let product_coefficient = ring_operator.multiply(coefficient.clone(), column_entry.val());
                if let Some(existing_coefficient) = result.get_mut( &row_index ) {
                    *existing_coefficient = ring_operator.add(*existing_coefficient, product_coefficient);
                } else {
                    result.insert( row_index, product_coefficient );
                }
            }
        }
        result

    }
}









// pub struct OuterProduct< Matrix >{
//     pub matrix: Matrix
// }




// impl < Matrix > MatrixOracle for
//     OuterProduct< Matrix >

//     where
//         Matrix: MatrixAlgebra
// {
//     type Coefficient            Matrix::Coefficient;

//     type RowIndex               Matrix::RowIndex;

//     type ColumnIndex            Matrix::RowIndex;

//     type RowEntry               Matrix::ColumnEntry;

//     type ColumnEntry            Matrix::ColumnEntry;

//     type Row                    ;

//     type RowReverse             ;

//     type Column                 ;

//     type ColumnReverse          ;

//     fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
//         todo!()
//     }

//     fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
//         todo!()
//     }

//     fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
//         todo!()
//     }

//     fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
//         todo!()
//     }

//     fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
//         todo!()
//     }

//     fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
//         todo!()
//     }

//     fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
//         todo!()
//     }
// }