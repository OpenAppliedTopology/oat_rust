//! Notes on how to do this
//! 
//! 
//! - CONSIDER SETTING UP LIBRARY FOR SETS AS SORTED LISTS, BUT KEEP IN MIND THAT THERE ARE DETAILS OF COMPUTER IMPLEMENTATION THAT COULD MAKE THE STANDARD RUST TOOLS FASTER THAN YOU'D THINK ... SO MAYBE GO WITH STANDARD, HEAVILY UNIT TESTED, CONCISE CODE FIRST, AND TRY FOR THIS OPTIMIZATION LATER (ALSO AFTER PRIORITIZING OTHER THINGS LIKE CORE USER API'S)
//! - start with a matrix M
//! - define a matrix N such that N[i,:] = sortperm(M[i,:])
//! - to enumerate all k-simplices with birth time B
//!   - for each row R of N define an iterator J_R as a kmerge on 
//!     [I_1, I_2, ..], where I_p is an iterator that runs over all simplices 
//!     with (i) at least p vertices that are born at time p in row R, and (ii)
//!     all other vertices born strictly before p.
//!     - NB: will probbably just have to concatenate the two sets and sort them
//!       as a vec [**this has the advantage of conceptual simplicity, even if
//!       it's not optimal**]
//!     - CONJECTURE: if you enumerate elements from the sets (i) and (ii) in 
//!       lexicographic order then perhaps you can find a sneaky way to combine
//!       them into a lazy iterator that runs over things in order (rather than
//!       placing in a vector and sorting)
//!   - do a kmerge on the set of all J_R, deleting duplicates, in order to get all
//!     simplices of birth time B in sorted lexicographic order
//! - to enumerate all k-simplices in (birth x vertex-lexicographic) order:
//!   - use the approach above to enumerate all k-simplices with birth time B,
//!     for B = 0, 1, 2, ..., and concatenate the iterators
//! - to enumerate cofacets
//!   - define a hash or vec representing a map: (new vertex) -> birth time
//!   - for each row R, check if the row contains the simplex, then enumerate 
//!     all other vertices in that row and their birth times; either add them to
//!     the hash/vec and update birth times if necessary
//!   - return an iterator that implicitly runs over these "addable" vertices, sorted in 
//!     order of birth, adding each vertex to the current simplex and updating birth time accordingly