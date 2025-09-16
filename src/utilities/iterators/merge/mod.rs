//! Utilities to merge multiple iterators, while preserving sorted order.
//! 
//! Example: combine `[ 1, 2, 4 ]` and `[ 2, 3, 5 ]` into a list `[ 1, 2, 2, 4, 5 ]`.
//! 
//! # Usage
//! 
//! - to merge two iterators of different types, use
//!   - to use the default order on entries: [merge](https://docs.rs/itertools/latest/itertools/fn.merge.html)
//!   - to use customized orders using a "closure": [merge_by](https://docs.rs/itertools/latest/itertools/trait.Itertools.html#method.merge_by)
//!   - to use customized orders using a OAT order operator (useful when you need to be explicit about object types): [MergeTwoIteratorsByOrderOperator](crate::utilities::iterators::merge::two_type::MergeTwoIteratorsByOrderOperator) 
//!   
//! - to merge k ≥ 3 iterators
//!   - to use the default order on entries: [kmerge](https://docs.rs/itertools/latest/itertools/fn.kmerge.html) or [hit_merge](crate::utilities::iterators::merge::hit::hit_merge),
//!   - to use the *reverse* of the default order: [hit_merge_descend](crate::utilities::iterators::merge::hit::hit_merge_descend)
//!   - to use customized orders using a "closure": [kmerge_by](https://docs.rs/itertools/latest/itertools/fn.kmerge_by.html) or [hit_merge_by_fnmut](crate::utilities::iterators::merge::hit::hit_merge_by_fnmut)
//!   - to use customized orders using a OAT order operator (useful when you need to be explicit about object types): the [hit_merge_by_predicate](crate::utilities::iterators::merge::hit::HitMergeByPredicateTrait::hit_merge_by_predicate) method of the [HitMergeByPredicateTrait](crate::utilities::iterators::merge::hit::HitMergeByPredicateTrait); the functions , [hit_merge_by_predicate](crate::utilities::iterators::merge::hit::hit_merge_by_predicate), 
//! 
//! - to merge k ≥ 3 iterators of different types
//!   - no method currently available; might try combining types into an Rust enum, c.f. [TwoTypeIterator](crate::utilities::iterators::general).
//! 
//! 
pub mod hit;
pub mod two_type;


