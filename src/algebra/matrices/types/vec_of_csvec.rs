
use std::iter;


pub struct VecCsv< ColumnIndex, Coefficient >
{
    major_dimension: MajorDimension, 
    min_ind: Vec< Vec< ColumnIndex > > ,
    snz_val: Vec< Vec< Coefficient > >
}


impl    < ColumnIndex, Coefficient >
        VecOfVec 
        < ColumnIndex, Coefficient > 
{
    // Make new (empty) VecOfVec. 
    pub fn new( major_dimension: MajorDimension ) -> Self  
    {
        atrixOracle { scalar: scalar,
                             major_dimension: major_dimension,
                             phantom: PhantomData 
                            }
    }
}



