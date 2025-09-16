//! Sparse vectors represented by iterators.
//! 
//! OAT has some helpful tools for sparse vectors, organized into two groups: 
//! - [entries]: operations on entries (read, write, etc.)
//! - [operations]: operations on vectors (add, multiply, etc.)
//!
//! 
//!  # What is a sparse vector in OAT?
//! 
//! OAT typically represents a sparse vector as an iterator that runs over the structural nonzero index-coefficient pairs `(j1, u1), (j2, u2), ..`. 
//! - See below for a detailed description of what that means, mathematically!
//! - In practice, an index-coefficient pair doesn't have to be a tuple `(j, u)`. It can be any data structure that represents a key-value pair.
//!   OAT lets you treat an object as a key-value pair as long as it implements the [KeyValGet] trait. See the documentation page for
//!   [vector entries](crate::algebra::vectors::entries), for more details!
//! 
//!  # What is a sparse vector, mathematically?
//! 
//! A *vector* in linear algebra can be thought of as a sequence of coefficients `(v1, .., vN)`. We can equivalently think
//! about a vector as
//! 
//! - a function `v: {1,..,N} -> R`, where `R` is the set of possible coefficients
//! - a collection of index-coefficient pairs `(1,v1), .., (N, vN)`
//! 
//! We call `{1,..,N}` the *index set*. These different perspectives help us to think about vectors in more general terms.
//! Concretely, if we have a finite set `I` then we can define a vector with index set `I = {i1, .., iN}` as either
//! 
//! - a function `v: I -> R`, or
//! - a collection of index-coefficient pairs `(i1,v1), .., (iN, vN)`
//! 
//! Suppose we want to store a vector `v` in computer memory. If `v` contains many zero coefficients, then
//! it can be efficient to represent `v` as a collection of index-coefficient pairs `(i1,v1), .., (iN, vN)`, and *delete*
//! all the pairs with a coefficient equal to zero. This leaves a collection of pairs `(j1, u1), (j2, u2), ..`
//! 
//! In OAT, we think of a  a sparse vector, mathematically, as an index set `I = {i1, .., iN}` and a collection of index-coefficient pairs `(j1, u1), (j2, u2), ..`, where
//! `j1, j2, ..` runs over a subset of `I`.   We call the coefficients
//! *structural nonzero entries*. However, we don't actually require them to be nonzero, in practice!
//! 
//! 
//! 
//! 
// //! Recall that an *entry* is an object that stores the data of a `(key, value)` pair.
// //! Not every entry is a tuple; an entry can be 
// //! any struct that implements the [KeyValGet](crate::algebra::vectors::entries::KeyValGet) trait.
// //! 
// 
// //! 
// //!
// //! **Example**
// //! 
// //! ```
// //! let v = vec![ (1, 1.5), (2, 2.5) ]; // this is a standard Rust command; it's not special to OAT
// //! 
// //! for (index, coeff) in v {
// //!     println!("(index, coeff) = {:?}", (index, coeff));
// //! }
// //! 
// //! // Should display:
// //! // (index, coeff) = (1, 1.5)
// //! // (index, coeff) = (2, 2.5)
// //! ```
//! 
//! 
//! 
//! 
// I THINK THE FOLLOWING IS NOW DEPRECATED
// //! # Working with sparse vectors
// //!
// //! -  Simplify< IteratorsMergedInSortedOrder > is the preferred method of working with a linear combination of vectors where
// //!     - each vector returns entries in sorted order
// //!     - you want the sum of the entries to be returned in sorted order
// //!     - you want to create the sum now and modify it later, e.g. by adding new vectors
// //! - The `Pair< Index, Coeff >` items may outperform tuples of form `( Index, Coeff )` (at least
// //! for now) because tuples demand a certain memory structure that may impede performance, c.f. 
// //! [rewriting in memory](https://www.reddit.com/r/rust/comments/79ry4s/tuple_performance/).




pub mod entries;
pub mod operations;












// /// Iterating over the entries of a sparse vector.
// /// 
// /// This trait doesn't provide any methods. It's just a "wrapper" that helps write simpler, more readabel code. 
// /// 
// /// Concretely, this trait auto-implements for any 
// /// object of type `T` such that `T` is an iterator and the items
// /// iterated by `T` implement the `KeyValGet` and [KeyValNew] methods. It also provides three associated types
// /// - `Entry`, which equals `Self::Item`
// /// - `Index`, which equals `Self::Item::Key`, i.e. the `Key` type from [KeyValGet]
// /// - `Coefficient`, which equals `Self::Item::Val`, i.e. the `Val` type from [KeyValGet]
// /// 
// /// This allows us to write things like:
// /// 
// /// ```
// /// use oat_rust::algebra::vectors::IteratesOverIndexCoefficientPairs;
// /// 
// /// pub struct TwoVectors< I >
// /// 
// ///     where 
// ///         I:          IteratesOverIndexCoefficientPairs
// /// { 
// ///     vector_a:       I,
// ///     vector_b:       Vec< (I::Index, I::Coefficient) >,
// /// }
// /// 
// /// ```
// /// 
// /// Instead of something more verbose:
// /// 
// /// ```
// /// use oat_rust::algebra::entries::KeyValGet;
// /// 
// /// pub struct TwoVectors< I >
// /// 
// ///     where 
// ///         I:          Iterator,
// ///         I::Item:    KeyValGet
// /// {
// ///     vector_a:       I,
// ///     vector_b:       Vec< ( < I::Item as KeyValGet >::Key, < I::Item as KeyValGet >::Val )
// /// }
// /// 
// /// ```
// /// 
// /// # How to implement
// /// 
// /// To implement this trait for a new type `T`, simply ensure that `T` implements `Itetator`, and that `T::Item` implements [KeyValGet] and [KeyValNew].
// /// This trait will then auto-implement for `T`.
// pub trait IteratesOverIndexCoefficientPairs

// where    
//     Self:                               SparseVector + 
//                                         Iterator< Item = < Self as IntoIterator >::Item >,
//     < Self as IntoIterator >::Item:     KeyValGet,

// {}

// impl < T > 

//     IteratesOverIndexCoefficientPairs for 
    
//     T

//     where    
//     T:                                  IntoIterator< IntoIter = T > +
//                                         SparseVector + 
//                                         Iterator< Item = < T as IntoIterator >::Item >,
//     < T as IntoIterator >::Item:        KeyValGet,

// {}



// /// Similar to [IteratesOverIndexCoefficientPairs], but entries must implement [KeyValSet]
// /// 
// /// This trait doesn't provide any methods. It's just a "wrapper" that provides a slighty simpler way to
// /// write code (and produce more readable code) that involves iteration over the entries of a sparse vector. In particular,
// /// this trait auto-implements for any object of type `T` such that `T` is an iterator and the items
// /// iterated by `T` implement the `KeyValGet` and [KeyValNew] methods.
// /// 
// /// **As an added convenience** this trait has two associated types: `Index` and `Coefficient`. The `Index`
// /// type always equals the `Key` type from [KeyValGet] and the `Coefficient` type always equals the `Val`
// /// type from [KeyValGet]. 
// pub trait IteratesOverChangeableIndexCoefficientPairs

// where
//     Self:           IteratesOverIndexCoefficientPairs,

// {}

// impl < T > 

//     IteratesOverChangeableIndexCoefficientPairs for 
    
//     T

//     where       T:           IteratesOverIndexCoefficientPairs,  
//                 T::Entry:     KeyValNew,    

// {}






// /// Iterating over the entries of a sparse vector.
// /// 
// /// This trait doesn't provide any methods. It's just a "wrapper" that helps write simpler, more readabel code. 
// /// 
// /// Concretely, this trait auto-implements for any 
// /// object of type `T` such that `T` is an iterator and the items
// /// iterated by `T` implement the `KeyValGet` and [KeyValNew] methods. It also provides three associated types
// /// - `Entry`, which equals `Self::Item`
// /// - `Index`, which equals `Self::Item::Key`, i.e. the `Key` type from [KeyValGet]
// /// - `Coefficient`, which equals `Self::Item::Val`, i.e. the `Val` type from [KeyValGet]
// /// 
// /// This allows us to write things like:
// /// 
// /// ```
// /// use oat_rust::algebra::vectors::SparseVector;
// /// 
// /// pub struct TwoVectors< I >
// /// 
// ///     where 
// ///         I:          IntoIterator,
// ///         I::Item:    KeyValGet,
// /// { 
// ///     vector_a:       I,
// ///     vector_b:       Vec< ( < I as SparseVector >::Index, < I as SparseVector >::Coefficient) >,
// /// }
// /// 
// /// ```
// /// 
// /// Instead of something more verbose:
// /// 
// /// ```
// /// use oat_rust::algebra::entries::KeyValGet;
// /// 
// /// pub struct TwoVectors< I >
// /// 
// ///     where 
// ///         I:          Iterator,
// ///         I::Item:    KeyValGet
// /// {
// ///     vector_a:       I,
// ///     vector_b:       Vec< ( < < I as IntoIterator>::Item as KeyValGet >::Key, < < I as IntoIterator>::Item as KeyValGet >::Val )
// /// }
// /// 
// /// ```
// /// 
// /// # How to implement
// /// 
// /// To implement this trait for a new type `T`, simply ensure that `T` implements `Itetator`, and that `T::Item` implements [KeyValGet] and [KeyValNew].
// /// This trait will then auto-implement for `T`.
// /// 
// /// # Developer comments
// /// 
// /// We originally wanted to call this trait `IteratesOverIndexCoefficientPairs`, but this turned out to have very poor readability in real examples.
// /// The term `SparseVector` is a little unfortunate; ordinarily you would want to save this type of name for something more important than a
// /// just-for-convenience trait.
// pub trait SparseVector 

//     where       Self:               IntoIterator,
//                 Self::Item:         KeyValGet,
// {
//     type Entry;
//     type Index;
//     type Coefficient;
// }

// impl < T >

//     SparseVector for

//     T

//     where       Self:                                   IntoIterator,
//                 < Self as IntoIterator >::Item:         KeyValGet,
// {
//     type Entry          =   Self::Item;
//     type Index          =   < Self::Item as KeyValGet >::Key;
//     type Coefficient    =   < Self::Item as KeyValGet >::Val;
// }


// fn eat_vector< SparseVector >( v: SparseVector ) 
//     where 
//         SparseVector:   IntoIterator<
//                             IntoIter:   Iterator< 
//                                             Item:   KeyValGet 
//                                         >,
//                         >
// {
//     let a = v.into_iter().next().unwrap().key();
// }


// fn eat_vector_2< SparseVector >( v: SparseVector ) 
//     where 
//         SparseVector:   IntoIterator<
//                             Item:   KeyValGet 
//                         >
// {
//     let a = v.into_iter().next().unwrap().key();
// }





// /// Iterating over the entries of a sparse vector.
// /// 
// /// This trait doesn't provide any methods. It's just a "wrapper" that helps write simpler, more readabel code. 
// /// 
// /// Concretely, this trait auto-implements for any 
// /// object of type `T` such that `T` is an iterator and the items
// /// iterated by `T` implement the `KeyValGet` and [KeyValNew] methods. It also provides three associated types
// /// - `Entry`, which equals `Self::Item`
// /// - `Index`, which equals `Self::Item::Key`, i.e. the `Key` type from [KeyValGet]
// /// - `Coefficient`, which equals `Self::Item::Val`, i.e. the `Val` type from [KeyValGet]
// /// 
// /// This allows us to write things like:
// /// 
// /// ```
// /// use oat_rust::algebra::vectors::EntriesAreIndexCoefficientPairs;
// /// use oat_rust::algebra::rings::traits::SemiringOperations;
// /// 
// /// pub struct SparseVector< I, R >
// /// 
// ///     where 
// ///         R:      Semiring,
// ///         I:      Iterator + EntriesAreIndexCoefficientPairs< Index = usize, Coefficient = R::Element >
// /// { 
// ///     entries:                I,
// ///     coefficient_ring:       R,
// /// }
// /// 
// /// ```
// /// 
// /// Instead of something a bit less readable:
// /// 
// /// ```
// /// use oat_rust::algebra::entries::KeyValGet;
// /// 
// /// pub struct SparseVector< I, R >
// /// 
// ///     where 
// ///         R:          Semiring,
// ///         I:          Iterator,
// ///         I::Item:    KeyValGet< Key = usize, Val = R::Element >
// /// {
// ///     entries:                I,
// ///     coefficient_ring:       R,
// /// }
// /// 
// /// ```
// /// 
// /// # How to implement
// /// 
// /// To implement this trait for a new type `T`, simply ensure that `T` implements `IntoItetator`, and that `T::Item` implements [KeyValGet].
// /// This trait will then auto-implement for `T`.
// pub trait EntriesAreIndexCoefficientPairs 

//     where       Self:               IntoIterator,
//                 Self::Item:         KeyValGet,
// {
//     type Entry;
//     type Index;
//     type Coefficient;
// }


// impl < T >

//     EntriesAreIndexCoefficientPairs for

//     T

//     where       T:               IntoIterator,
//                 T::Item:         KeyValGet,
// {
//     type Entry          =   Self::Item;
//     type Index          =   < Self::Item as KeyValGet >::Key;
//     type Coefficient    =   < Self::Item as KeyValGet >::Val;
// }


             


// /// Similar to [EntriesAreIndexCoefficientPairs], but entries must implement [KeyValNew]
// /// 
// /// This trait doesn't provide any methods. It's just a "wrapper" that provides a slighty simpler way to
// /// write code (and produce more readable code) that involves iteration over the entries of a sparse vector. In particular,
// /// this trait auto-implements for any object of type `T` such that `T` is an iterator and the items
// /// iterated by `T` implement the `KeyValGet` and [KeyValNew] methods.
// /// 
// /// **As an added convenience** this trait has two associated types: `Index` and `Coefficient`. The `Index`
// /// type always equals the `Key` type from [KeyValGet] and the `Coefficient` type always equals the `Val`
// /// type from [KeyValGet]. 
// pub trait EntriesAreChangeableIndexCoefficientPairs

// where
//     Self:           EntriesAreIndexCoefficientPairs,
//     Self::Item:     KeyValNew,

// {}

// impl < T > 

// EntriesAreChangeableIndexCoefficientPairs for 
    
//     T

//     where       T:           EntriesAreIndexCoefficientPairs,  
//                 T::Item:     KeyValNew,    

// {}