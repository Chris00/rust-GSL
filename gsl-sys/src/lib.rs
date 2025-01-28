//
// FFI binding for the GSL library
//

#![cfg_attr(feature = "dox", feature(doc_cfg))]
#![allow(improper_ctypes)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(clippy::redundant_static_lifetimes)]
#![allow(clippy::upper_case_acronyms)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

pub use libc::FILE;

// mod auto;

// pub use auto::*;
