mod tree;
mod hopscotch;
mod util;
mod bucket_info;
mod hash_file;

use failure::Fail;
use std::io;

#[derive(Debug, Fail)]
pub enum CacheError {
    #[fail(display = "invalid dna base: {}", base)]
    InvalidDnaBase{
        base: char
    },
    #[fail(display = "invalid dna base index: {}", index)]
    InvalidDnaBaseIndex{
        index: usize
    },
    #[fail(display = "unable to store dna sequence with {} bases", bases)]
    UnableToStoreDnaSequence{
        bases: usize
    },
    #[fail(display = "unable to find an empty slot")]
    UnableToFindEmptySlot,
    #[fail(display = "IO error: {}", error)]
    IoError { error: io::Error },
}

impl From<io::Error> for CacheError {
    fn from(error: io::Error) -> Self {
        CacheError::IoError {error}
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
