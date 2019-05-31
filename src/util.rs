use std::convert::{TryFrom, TryInto};
use crate::CacheError;
use bv::{BitSlice, BitSliceMut, Bits, BitVec};
use bv::BitsMut;
use std::str::{FromStr, Chars};
use rand::Rng;

pub fn gen_rand_seq<R: Rng>(rng: &mut R, min: usize, max: usize) -> String {
            let l = rng.gen_range(min, max);
            (0..l)
                .map(|_| {
                    let idx: usize = rng.gen_range(0, 4);
                    let c: char = NucleotideBase::try_from(idx)
                        .unwrap()
                        .into();
                    c
                }).collect()
}

pub fn gen_rand_seqs<R: Rng>(rng: &mut R, n: usize, min: usize, max: usize) -> Vec<String> {
    let n = 100;
    let k = 100;
    (0..n)
        .map(|_| {
            gen_rand_seq(rng, min, max)
        }).collect()
}

//noinspection RsExternalLinter
//noinspection RsExternalLinter
pub fn dna_seq<'a>(s: &'a str) -> DnaStrIter<'a> {
    DnaStrIter::new(s)
}

pub struct DnaStrIter<'a> {
    c: Chars<'a>
}

impl<'a> DnaStrIter<'a> {
    pub fn new(s: &'a str) -> Self {
        DnaStrIter{c: s.chars()}
    }
}

impl<'a> Iterator for DnaStrIter<'a> {
    type Item = NucleotideBase;

    fn next(&mut self) -> Option<Self::Item> {
        self.c.next()
            .map(|c| match c {
                'A' => NucleotideBase::A,
                'C' => NucleotideBase::C,
                'G' => NucleotideBase::G,
                'T' => NucleotideBase::T,
                _ => panic!("Invalid DNA base")
            })
    }
}

pub struct BitVecStorage {
    bv: BitVec<u64>
}

impl NucleotideBaseStorage for BitVecStorage {
    fn get_bit_slice(&self) -> BitSlice<u64> {
        self.bv.as_slice()
    }

    fn get_bit_slice_mut(&mut self) -> BitSliceMut<u64> {
        self.bv.as_mut_slice()
    }

    fn try_reserve(&mut self, n: usize) -> bool {
        self.bv.reserve(n as u64);
        true
    }

    fn with_base_len(bases: usize) -> Result<Self, CacheError> {
        Ok(BitVecStorage{
            bv: BitVec::<u64>::new_fill(false, 2* bases as u64)
        })
    }
}

pub trait NucleotideBaseStorage: Sized {
    fn get_bit_slice(&self) -> BitSlice<u64>;
    fn get_bit_slice_mut(&mut self) -> BitSliceMut<u64>;
    fn bases(&self) -> usize {
        (self.get_bit_slice().len()/2) as usize
    }

    fn try_reserve(&mut self, n: usize) -> bool;
    fn with_base_len(bases: usize) -> Result<Self, CacheError>;

    fn set_base(&mut self, i: usize, base: NucleotideBase) {
        let mut b = self.get_bit_slice_mut();
        let (b1, b2) = base.to_bits();

        let i = i as u64;
        b.set_bit(i*2,b1);
        b.set_bit(i*2+1, b2);
    }

    fn get_base(&self, i: usize) -> NucleotideBase {
        let b = self.get_bit_slice();
        let i = i as u64;
        NucleotideBase::from_bits(b.get_bit(i*2), b.get_bit(i*2+1))
    }
}

/// Represents a sequence of `NucleotideBase`s backed by the given bit-storage
pub struct NucleotideBaseSequence<N> {
    storage: N
}

impl<N> NucleotideBaseSequence<N>
    where N: NucleotideBaseStorage {

    pub fn base_iter(&self) -> impl Iterator<Item = NucleotideBase> + '_ {
        let n = self.bases();
        (0..n)
            .map(move |i| self.storage.get_base(i))
    }

    pub fn bases(&self) -> usize {
        self.storage.bases()
    }

    //noinspection RsExternalLinter
    pub fn complement_sequence<M>(&self) -> Result<NucleotideBaseSequence<M>, CacheError>
        where M: NucleotideBaseStorage {
        let mut storage = M::with_base_len(self.bases())?;

        for (i, b) in self.base_iter().enumerate() {
            storage.set_base(i, b.complement());
        }

        Ok(NucleotideBaseSequence{storage})
    }

    pub fn into_storage(self) -> N {
        self.storage
    }

}

impl<N> ToString for NucleotideBaseSequence<N>
    where N: NucleotideBaseStorage {
    fn to_string(&self) -> String {
        self.base_iter()
            .map(|b| {
                let c: char = b.into();
                c
            })
            .collect()
    }
}

impl<N> FromStr for NucleotideBaseSequence<N>
    where N: NucleotideBaseStorage {
    type Err = CacheError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n = s.len();
        let mut storage = N::with_base_len(n)?;

        for (i, s) in s.chars().enumerate() {
            storage.set_base(i, s.try_into()?);
        }

        Ok(NucleotideBaseSequence{
            storage
        })
    }
}

///`NucleotideBase` represents the 4 nucleotide bases
#[derive(PartialOrd, PartialEq, Ord, Eq, Copy, Clone, Debug)]
pub enum NucleotideBase {
    A,
    C,
    G,
    T,
}

impl NucleotideBase {
    /// Calculates the complement for the given Base
    /// The complements are A - T, C - G
    pub fn complement(self) -> NucleotideBase {
        use NucleotideBase::*;
        //TODO: check if rust realizes this can be written as index(base) ^ 11b

        match self {
            A => T,
            T => A,
            G => C,
            C => G
        }
    }

    /// Converts 2 Bits into a base
    pub fn from_bits(b1: bool, b2: bool) -> NucleotideBase {
        use NucleotideBase::*;

        match (b1, b2) {
            (false, false) => A,
            (false, true) => C,
            (true, false) => G,
            (true, true) => T
        }
    }

    /// Converts base into 2 Bits
    pub fn to_bits(self) -> (bool, bool) {
        use NucleotideBase::*;

        match self {
            A => (false, false),
            C => (false, true),
            G => (true, false),
            T => (true, true)
        }
    }
}

impl TryFrom<char> for NucleotideBase {
    type Error = CacheError;

    fn try_from(base: char) -> Result<Self, Self::Error> {
        Ok(match base {
            'A' => NucleotideBase::A,
            'C' => NucleotideBase::C,
            'G' => NucleotideBase::G,
            'T' => NucleotideBase::T,
            _ => return Err(CacheError::InvalidDnaBase {base})
        })
    }
}

impl Into<char> for NucleotideBase {
    fn into(self) -> char {
        use NucleotideBase::*;

        match self {
            A => 'A',
            C => 'C',
            G => 'G',
            T => 'T'
        }
    }
}

impl TryFrom<usize> for NucleotideBase {
    type Error = CacheError;

    fn try_from(index: usize) -> Result<Self, Self::Error> {
        Ok(match index {
            0 => NucleotideBase::A,
            1 => NucleotideBase::C,
            2 => NucleotideBase::G,
            3 => NucleotideBase::T,
            _ => return Err(CacheError::InvalidDnaBaseIndex {index})
        })
    }
}

impl Into<usize> for NucleotideBase {
    fn into(self) -> usize {
        use NucleotideBase::*;

        match self {
            A => 0,
            C => 1,
            G => 2,
            T => 3
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::util::{NucleotideBase, NucleotideBaseSequence, BitVecStorage};
    use std::convert::TryFrom;

    #[test]
    fn test_dna_base() {
        let bases = [NucleotideBase::A, NucleotideBase::C, NucleotideBase::G, NucleotideBase::T];

        //test vice-versa converting

        //char
        for &b in bases.iter() {
            let c: char = b.into();
            let b1 = NucleotideBase::try_from(c).
                expect("must be correct base");

            assert_eq!(b, b1)
        }

        //index
        for &b in bases.iter() {
            let i: usize = b.into();
            let b1 = NucleotideBase::try_from(i).
                expect("must be correct base");

            assert_eq!(b, b1)
        }

        //bits
        for &b in bases.iter() {
            let c = b.to_bits();
            let b1 = NucleotideBase::from_bits(c.0, c.1);
            assert_eq!(b, b1)
        }
    }

    #[test]
    fn test_base_complement() {
        let bases = [NucleotideBase::A, NucleotideBase::C, NucleotideBase::G, NucleotideBase::T];
        let bases_c = [NucleotideBase::T, NucleotideBase::G, NucleotideBase::C, NucleotideBase::A];

        for (&b, &c) in bases.iter().zip(bases_c.iter()) {
            assert_eq!(c, b.complement());
            assert_eq!(b, c.complement());
        }

    }

    #[test]
    fn test_dna_parse() {
        let seqs = ["A", "C", "G", "T", "ACGT", "TGCA"];
        for &seq in seqs.iter() {
            let dna_seq: NucleotideBaseSequence<BitVecStorage> = seq.parse()
                .expect("Must parse");

            let s = dna_seq.to_string();
            assert_eq!(seq, &s);
        }
    }
}