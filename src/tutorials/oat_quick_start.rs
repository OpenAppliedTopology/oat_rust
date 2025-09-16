

//! # Quick start with OAT
//! 
//! Here's an example of how to write your own program using features from OAT.  
//!
//! **Step 1: Install Rust** 
//! 
//! Instructions can be found on the [Rust website](https://www.rust-lang.org/tools/install) and in the [Rust Book](https://doc.rust-lang.org/book/ch01-01-installation.html)
//! 
//! **Step 2: Make a crate** 
//! 
//! A crate is a collection of Rust files that are designed to work together.  These files all reside within
//! a folder.  For this exercise, just open an terminal/shell and cd to wherever you'd like the new folder
//! to be located (there are no restrictions on where it can be placed).  Then run
//! 
//!    ```bash
//!    cargo new oat_arithmetic
//!    ```
//!   This will make a new folder named `oat_arithmetic,` which contains several Rust files.
//!   If you're curious, more detailed instructions to build a new crate be found in the [Rust Book](https://doc.rust-lang.org/book/ch01-02-hello-world.html).
//! 
//! 
//! **Step 3: Clone the OAT repository**  
//! 
//! This will create a local copy of the repository in a folder on your computer.
//! You'll need to remember the path to this folder.
//! 
//! **Step 4: Add OAT to your crate's dependencies.** 
//! 
//!   To do this, find the `Cargo.toml` file inside the folder that houses your crate.
//!   Open the file with a text editor, and you will see something like the following:
//!
//!    ```bash
//!    [package]
//!    name = "oat_arithmetic"
//!    version = "0.1.0"
//!    edition = "2022"
//!    
//!    # See more keys and their definitions at https://doc.rust-lang.org/cargo/reference/manifest.html
//!    
//!    [dependencies]
//!    
//!    
//!    ```
//! 
//!   Modify this file by adding a line directly under `[dependencies]` as follows:
//!   
//!   ```bash
//!   [package]
//!   name = "oat_arithmetic"
//!   version = "0.1.0"
//!   edition = "2022"
//!   
//!   # See more keys and their definitions at https://doc.rust-lang.org/cargo/reference/manifest.html
//!   
//!   [dependencies]
//!   oat_rust = { path = "path/to/your/copy/of/the/oat_rust/repository" }
//!   
//!   ```
//! 
//!   **Note that the line is `oat_rust = { ..` not `oat = { ..`.
//! 
//! **Step 5: Modify the file `main.rs` to use a feature from the OAT.** 
//! 
//! Let's use the OAT feature for arithmetic over finite fields.  Open
//! `oat_arithmetic/src/main.rs`.  This should look like the following:
//! 
//!   ```
//!   fn main() {
//!       println!("Hello, world!");
//!   }
//!   
//!   ```
//! 
//!   To use a feature from OAT, you'll add an "import" statement to the
//!   top of the file.  This statment typically takes the form `use <feature you want to use>`.
//!   
//!   ```
//!   // Import the definition of a "ring operator" for prime order fields.
//!   // Ring operators are tools that OAT uses to perform arithmetic over rings.
//!   use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
//!   
//!   // Import the set of commands that ring operators use to perform basic arthimetic.
//!   use oat_rust::algebra::rings::traits::SemiringOperations;
//!   
//!   fn main() {
//!   
//!       // Define the order of the field.
//!       let modulus = 2;
//!   
//!       // Define the ring operator.
//!       let ring_operator = PrimeOrderField::new( modulus  );
//!   
//!       // Add one plus one
//!       let sum = ring_operator.add( 1, 1 );
//!   
//!       // Print the result
//!       println!(""); // create an empty line
//!       println!("Did you know that one plus one equals {:?}, mod {:?}?", sum, modulus);
//!       println!(""); // create an empty line    
//!   }
//!   
//!   ```
//! 
//! 
//! **Step 6: Run your program**
//! 
//! Make sure that the changes to `main.rs` are saved.  Open a terminal, and cd into the `oat_arithmetic` directory.  Then run
//! ```bash
//! cargo run
//! ```
//! 
//! 
//! This should print
//! ```bash
//! $ Did you know that one plus one equals 0, mod 2?
//! ```
//! 
//! Congratulations, you've just written and run a Rust program, with OAT!
//! 
//! ## Run your code *faster*
//! 
//! Rust has several different options for running a program.  The default options run programs
//! in "debug mode," which is safer and helps the user to find and correct errors.  However, debug mode is
//! *many times slower* than another option, called "release mode."  You can do this with 
//!
//! ```bash
//! cargo run --release
//! ```
//! 
//! See the official Rust website for further details.
//! 
//! ## Import packages, organize files and folders, run multiple programs, etc.
//! 
//! For more complex projects, you'll want to use additional packages, split programs into separate files,
//! compile code to executable, etc.  [This small demonstration crate](https://stash.pnnl.gov/projects/TDAC/repos/demo_rust_library/browse) can help get you started.
//! - The `README.md` has many useful tips
//! - Feel free to copy parts of the demo for your own project!
//! 
// //! ## Add dependencies
// //! 
// //! This example shows how to use OAT in a new Rust crate.  There are two essential steps in this process: (1) modify
// //! the `Cargo.toml` folder by updating the section headed by `[dependencies]`, and (2) place `use` statements in your
// //! program file to indicate what features from OAT you want to use.
// //! 
// //! You can use features from other crates as well.  Many of the most commonly used crates are available online, through
// //! a website called [crates.io](https://crates.io/).  To add one of these crates to your `Cargo.toml` file, you just
// //! need the name of the crate and a version number; see the [Rust Book](https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html)
// //! for details.
// //! 
// //! 
// //! # Organize your crate and run multiple programs
// //! 
// //! 
// //! Once you write enough lines of code, you may eventually want to split them into separate files or even separate folders.
// //! Rust offers different ways to do this, but due to the large number of options it can be hard to find
// //! simple answers to simple questions.
// //! 
// //! 
// //! [This simple demo crate](https://stash.pnnl.gov/projects/TDAC/repos/demo_rust_library/browse) can help get you started.
// //! - Make sure to read the `README.md` file in this demo, for pointers
// //! - Feel free to copy parts of this library for your own project!
// //! 
// //! 
// //! **Nota bene** there are two types of crates: library and binary.  The demo is a library crate.  The `oat_arithmetic` crate that you just created is a binary crate.
// //! 
// //! 
// //! 
// //! 
// //! General guidelines for library crates:
// //! 
// //! - files that you wish to turn into executable go in `src/bin`
// //! - other program files go in `src/`; each file shoud end in `.rs`
// //! - you can group files in `src/` by placing them into folders; each folder should contain a `mod.rs` file
// //! - the `mod.rs` file lists other files within the same folder, which you want to make "visible" to other
// //! files outside the folder
// //! - the `lib.rs` is essentially the same as a `mod.rs` file; it just gets a special name because it sits in `src/` and not
// //! a subfolder of `src/`.
//! 

//!     
// //! > Create the directory bin under src, then create any number of rs files, e.g. `src/bin/foo.rs`. Then, either use VS Code's Run button or `cargo run --bin foo`. If you need them to be multiple files, instead create `src/bin/foo/main.rs` and any other modules under `src/bin/foo`.