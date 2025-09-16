
//! # Open Applied Topology
//! 
//! Open Applied Topology (OAT) is software for fast, user-friendly algebra and topology.
//! 
//! - [Welcome!](#welcome)
//! - [Community](#community)
//! - [Values](#values)
//! - [Mission](#features)
//! - [Get started](#get-started)
//! 
//! # Welcome!
//! 
//! Welcome!  This is the the documentation homepage for OAT-Rust, which is the Rust component
//! of the Open Applied Topology (OAT) package. For Python tools, check out [oat_python](https://pypi.org/project/oat-python/),
//! and for a project overview, check out the [OAT homepage](https://openappliedtopology.github.io).
//! OAT-Rust is registered as [oat_rust](https://crates.io/crates/oat_rust) on Crates.io.  It provides powerful tools for applied topology, including
//! 
// //! - Persistent homology
// //! - Simplicial complexes
// //! - Homological algebra
// //! - Cycle optimization
// //! - Interactive 2d and 3d visualization
//! - [Algebra](crate::algebra)
//!   - [Rings](crate::algebra::rings)
//!   - [Vectors](crate::algebra::vectors)
//!   - [Matrices](crate::algebra::matrices)
//! - [Topology](crate::topology)
//!   - [Simplicial complexes](crate::topology::simplicial)
//!   - [Point clouds](crate::topology::point_cloud)
//!   - [Chain complexes](crate::algebra::chain_complexes)
//!   - [Barcodes and persistence diagrams](crate::algebra::chain_complexes::barcode)
//!   - [(Persistent) (relative) (co)homology](crate::tutorials::persistent_homology)
//! - [Utilities](crate::utilities)
//!   - [Iterators](crate::utilities::iterators)
//!   - [Optimization](crate::utilities::optimization)
//!   - [Order](crate::utilities::order)
//!   - [Sets](crate::utilities::sets)
//!   - [Indexing](crate::utilities::indexing_and_bijection)
//!   - and more!
//! - Python Interactive 2d and 3d visualization
//!   - See our sister library, [oat_python](https://pypi.org/project/oat-python/)
//! 
//! 
//! 
//! # Community
//! 
//! OAT is by and for the open source community.  <span style="color: SteelBlue;">Reach out to the developers</span> if you
//! - Need help getting started
//! - Wish for a missing feature
//! - Want to try coding
//! 
//! A collaboration of 20 research centers at colleges, universities, private, and public organizations support OAT's
//! development. The founding developers are Princton University, Macalester College, and the University of Delaware
//! The National Science Foundation granted seed funding for OAT in
//! 2019 under the [ExHACT]((https://www.nsf.gov/awardsearch/showAward?AWD_ID=1854748&HistoricalAwards=false))
//! project, and [Pacific Northwest National Laboratory](https://www.pnnl.gov/) (PNNL) provides continuing financial
//! support.  PNNL now coordinates development.  See <span style="color: SteelBlue;">[here](./ATTRIBUTIONS.md)</span>
//! for further details.
//! 
//! # Values
//! 
//! Our [shared values](https://github.com/OpenAppliedTopology/oat_python/blob/main/CODE_OF_CONDUCT.md) are
//! 
//! - Inclusion
//! - Respect, and 
//! - Expanding human knowledge through algebraic topology
//! 
//! 
//! # Mission
// //! 
// // //! OAT offers powerful features for beginners through advanced users.
// //!  
// //! <details>
// //! <summary>Performance</summary>
// //! <br>
// //! OAT is a first-class engine for cutting-edge applications.  It is specially suited to large, sparse // matrices.
// //!     The core library is written in Rust, a low-level systems programming language designed for safety and // performance.
// //!     High-level wrappers are available in Python. 
// //! </details>
// //! 
// //! <details>
// //! <summary>Reliability</summary>
// //! <br>
// //!   - Verification and test coverage: the OAT library is extensively tested.  In addition, the modular design // of the library makes it easy for users to generate their own certificates of correctness.
// //!   - Safety: OAT inherits strong safety guarantees from the features of the Rust compiler, especially in the // open source development process
// //! </details>
// //! 
// //! <details>
// //! <summary>Transparency</summary>
// //! <br>
// //!   - Documentation:  emphasizes clarity and accessibility for users with all backgrounds. OAT docs provide // explicit descriptions of both code *and* the underlying mathematical concepts. 
// //!   - [Tutorials](crate::tutorials) offer examples and helpful tips for beginners through advanced Rust users.
// //!   - Indexing: is one of the most pervasive challenges to writing transparent, interpretable code in // computational topology.  OAT's matrices and vectors can be indexed by simplices, cubes, and other user-defined // data structures, in addition to integers.
// //! </details>
// //! 
// //! <details>
// //! <summary>Modularity</summary>
// //! <br>
// //! 
// // //! Creative recombination of building blocks is among the most important ways we innovate.  
// //! 
// //!   OAT breaks problems into the same basic building blocks that topologists use when writing on a chalk // board. Users can mix and match those blocks with confidence, with a simple, streamlined interface.  They can // even create new components that work seemlessly with the rest of the library, including coefficient rings, // sparse matrix data structures, and customized filtrations on simplicial complexes.
// //! </details>
//! 
//! 
//! 
//! **Performance**
//!     
//! OAT is a first-class solver for cutting-edge applications.  It is ideally suited to large, sparse data sets.
//!     The core library is written in Rust, a low-level systems programming language designed for safety and performance.
//!     High-level wrappers are available in Python. 
//! 
//! 
//! **Reliability**
//! 
//! OAT has more unit tests than type definitions and function definitions, combined.
//!   Its modular design enables end users to write their own checks for correctness, with ease.
//!   The library inherits strong safety guarantees from the the Rust compiler.
//! 
//! 
//! **Transparency**
//! 
//! OAT documentation emphasizes clarity and accessibility for users with all backgrounds.  It includes more than 180 working examples, and describes both code and underlying mathematical concepts in detail.
//! [Tutorials](crate::tutorials) illustrate how to combine multiple tools into larger applications.
//! The platform's modular design breaks large solvers into basic components, which it exposes to the user for inspection.  In addition, the library provides powerful methods to inspect and analyze objects, consistent with the way humans naturally think about problems; for example, you can look up rows and columns of boundary matrices using *cubes*, *simplices*, or *cells* as keys.
//!   
//! 
//! **Modularity**
//! 
// //! Creative recombination of building blocks is among the most important ways we innovate.  
//! 
//!   OAT reduces complex problems to the same basic building blocks that topologists use when writing on a chalk board. Users can mix and match those blocks with confidence, using a simple, streamlined interface.  They can even create new components that work seemlessly with the rest of the library, including coefficient rings, sparse matrix data structures, and customized filtrations on simplicial complexes.
//! 
//! 
//! 
//! 
//! # Get Started
//! 
//! 
//! **Try the tutorials**
//! 
//! - Check out the Python tutorials in [oat_python](https://oat_python.readthedocs.io).
//! - Try the [Rust tutorials](crate::tutorials), including a
//! [quick start guide to writing your first program with OAT](crate::tutorials::oat_quick_start).
//! 
//! **Find answers in the documentation**
//! 
//! Rust documentation places lists (of objects, functions, etc.) in two places: either the bottom of a page, or in the
//! menu bar on the left. You can also use the search bar at the top to pull up a list of related terms.
//! The question mark button to the right of the bar gives other helpful tips (for example, you can search for functions based on their type signature).
//! 
//! 
//! 
//!



pub mod algebra;
pub mod topology;
pub mod utilities;
pub mod tutorials;


