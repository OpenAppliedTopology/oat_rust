//! Tools to decompose a zigzag module, including zigzag homology and quivers of type A_n.

pub mod cospans;
pub mod decompose;
pub mod hypergraph_pipeline;
pub mod cospan_pipeline;
pub mod spans;
pub mod span_pipeline;




// data schema
//
// INPUT
// file for all hyperedges
// file for all hypergraphs
//
// INTERMEDIATE
// separate json file for each arrow:
//   ( direction, matrix, )
// separate json file for each space:
//   grading on basis vectors
//
// OUTPUT
// basis for each space written in sequence of vectors of simplices (becomes a dataframe)
// bar objects + grading 
// --- contained in a dataframe with one row per bar, one column for bar objects, and another column for grading]



// 1: user creates a list of data objects in python
// 2: convert this to a list of rust data objects
// 3: in embarassingly parallel fashion, 