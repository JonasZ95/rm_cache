use bv::{BitSliceMut, BitSlice};
use nthash::NtHashForwardIterator;

use crate::util::{NucleotideBaseStorage, NucleotideBaseSequence};
use crate::hopscotch::{HopscotchHash, Empty, PartitionProvider, Partition, Bucket};
use crate::bucket_info::{BucketInfo, BucketIndexIterator};
use crate::CacheError;
use std::marker::PhantomData;
use std::io::{Seek, Read, Write, SeekFrom};
use std::cmp::Ordering;

const SEQ_LEN: usize = 4;
const MAX_SEQ: usize = 32*SEQ_LEN;
const P: usize = 12;

/// Represents a type which can represent up to 128 bases
#[derive(Default, Clone)]
struct StaticDnaSeq([u64; SEQ_LEN]);

impl NucleotideBaseStorage for StaticDnaSeq {
    fn get_bit_slice(&self) -> BitSlice<u64> {
        self.0.as_ref().into()
    }

    fn get_bit_slice_mut(&mut self) -> BitSliceMut<u64> {
        self.0.as_mut().into()
    }

    fn try_reserve(&mut self, n: usize) -> bool {
        n <= MAX_SEQ
    }

    fn with_base_len(bases: usize) -> Result<Self, CacheError> {
        if bases <= MAX_SEQ{
            Ok(StaticDnaSeq::default())
        } else {
            Err(CacheError::UnableToStoreDnaSequence { bases })
        }
    }
}

#[repr(packed)]
#[derive(Default, Clone, Copy)]
struct DnaSeqKey {
    dna_seq: [u64; SEQ_LEN],
    hash: u64,
    n: u64,
}

///Bucket to store a sequence and a key
#[repr(packed)]
#[derive(Default, Copy, Clone)]
struct DnaBucket {
    k: DnaSeqKey,
    v: u64,
    info: u64,
}

#[repr(packed)]
struct DnaPartition{
    changed: u64,
    buckets: [DnaBucket; P]
}

impl Default for DnaPartition {
    fn default() -> Self {
        DnaPartition {
            changed: 0,
            buckets: [DnaBucket::default(); P]
        }
    }
}

impl DnaBucket {
    fn bucket_info(&self) -> BucketInfo {
        BucketInfo::new(self.info)
    }
}

impl Bucket<DnaSeqKey, u64> for DnaBucket {
    type EntryIterator = BucketIndexIterator;

    fn get_key(&self) -> &DnaSeqKey {
        &self.k
    }

    fn get_key_mut(&mut self) -> &mut DnaSeqKey {
        &mut self.k
    }

    fn get_value(&self) -> &u64 {
        &self.v
    }

    fn get_value_mut(&mut self) -> &mut u64 {
        &mut self.v
    }

    fn has_entry(&self, i: usize) -> bool {
        self.bucket_info().has_entry(i)
    }

    fn set_entry(&mut self, i: usize){
        let mut bi = self.bucket_info();
        bi.set_entry(i);
        self.info =  bi.value();
    }

    fn unset_entry(&mut self, i: usize) {
        let mut bi = self.bucket_info();
        bi.unset_entry(i);
        self.info =  bi.value();
    }

    fn clear_entries(&mut self) {
        let mut bi = self.bucket_info();
        bi.clear();
        self.info =  bi.value();
    }

    fn entry_count(&self) -> usize {
        self.bucket_info().entry_count()
    }

    fn first_unset_entry(&self, h: usize) -> Option<usize> {
        self.bucket_info().first_unset_entry(h)
    }

    fn entry_iter(&self) -> Self::EntryIterator {
        self.bucket_info().index_iter()
    }

    fn is_empty(&self) -> bool {
        self.bucket_info().is_empty()
    }
}

impl Partition<DnaBucket, DnaSeqKey, u64> for DnaPartition {
    fn get_bucket(&self, i: usize) -> &DnaBucket {
        &self.buckets[i]
    }

    fn get_bucket_mut(&mut self, i: usize) -> &mut DnaBucket {
        &mut self.buckets[i]
    }

    fn mark_as_changed(&mut self) {
        self.changed = 1;
    }
}


impl Empty for DnaSeqKey {
    fn new_empty() -> Self {
        DnaSeqKey::default()
    }

    fn is_empty(&self) -> bool {
        self.n == 0
    }
}

impl HopscotchHash for DnaSeqKey {
    fn hopscotch_hash(&self) -> u64 {
        self.hash
    }
}

impl PartialEq for DnaSeqKey {
    fn eq(&self, other: &DnaSeqKey) -> bool {
        self.n == other.n && self.dna_seq == other.dna_seq
    }
}

impl DnaSeqKey {
    fn from_str(s: &str) -> Result<DnaSeqKey, CacheError> {
        let seq: NucleotideBaseSequence<StaticDnaSeq> = s.parse()?;
        let n = s.len() as u64;

        let b = s.as_bytes();
        let nt_hasher = NtHashForwardIterator::new(b, 2)
            .expect("Hash must ");

        let hash = nt_hasher.last()
            .expect("Last hash must exist");


        let dna_seq = seq.into_storage().0;
        Ok(DnaSeqKey {
            dna_seq,
            n,
            hash,
        })
    }
}

struct DnaHashFile<F>  {
    lru: lru::LruCache<usize, usize>,
    part_vec: Vec<DnaPartition>,
    h: usize,
    n: usize,
    parts: usize,
    part_size: usize,
    c: usize,
    f: F
}

impl<F> DnaHashFile<F>
    where F: Seek + Read + Write {

    pub fn new(f: F, h: usize, n: usize, parts: usize, part_size: usize, c: usize) -> Self {
        let lru = lru::LruCache::new(c);

        let part_vec = Vec::with_capacity(c);

        DnaHashFile {
            part_vec,
            lru,
            h,
            n,
            parts,
            part_size,
            c,
            f
        }
    }

    fn load_part_inner(&mut self, p: usize) -> Result<usize, CacheError> {
        if let Some(&i) = self.lru.get(&p) {
            return Ok(i);
        }

        //If vec was not filled
        let i = if self.part_vec.len() < self.c {
            self.part_vec.push(DnaPartition::default());
            self.part_vec.len() - 1
        } else {
            self.pop_part()?
        };

        self.read_part_into(i, p)?;
        self.lru.put(p, i);
        Ok(i)
    }

    fn pop_part(&mut self) -> Result<usize, CacheError> {
        let (p, i) = self.lru.pop_lru().unwrap();

        if self.part_vec[i].changed == 1 {
            self.write_part(i, p)?;
        }

        Ok(i)
    }

    fn write_part(&mut self, i: usize,  p: usize) -> Result<(), CacheError> {
        let sz = std::mem::size_of::<DnaPartition>();

        let part = &mut self.part_vec[i];
        part.changed = 0;

        unsafe  {
            let part_slice = std::slice::from_raw_parts(
                (part as *const DnaPartition) as *const u8,
                sz);

            self.f.seek(SeekFrom::Start((p *  sz) as u64))?;
            self.f.write_all(part_slice)?;
        };

        Ok(())
    }


    fn read_part_into(&mut self, i: usize, p: usize) -> Result<(), CacheError> {
        let sz = std::mem::size_of::<DnaPartition>();

        let part = &mut self.part_vec[i];
        unsafe  {
            let mut part_slice = std::slice::from_raw_parts_mut(
                part as *mut _ as *mut u8,
                sz);


            self.f.seek(SeekFrom::Start((p * sz) as u64))?;
            self.f.read_exact(&mut part_slice)?;
        };


        Ok(())
    }

}

impl<F> PartitionProvider<DnaBucket, DnaSeqKey, u64> for DnaHashFile<F>
    where F: Seek + Read + Write {
    type Partition = DnaPartition;
    type Error = CacheError;

    fn n(&self) -> usize {
        self.n
    }

    fn m(&self) -> usize {
        self.parts * self.part_size
    }

    fn h(&self) -> usize {
        self.h
    }

    fn parts(&self) -> usize {
        self.parts
    }

    fn part_size(&self) -> usize {
        self.part_size
    }

    fn load_part(&mut self, p: usize) -> Result<&Self::Partition, Self::Error> {
        self.load_part_inner(p)
            .map(move |i| &self.part_vec[i])
    }

    fn load_part_mut(&mut self, p: usize) -> Result<&mut Self::Partition, Self::Error> {
        self.load_part_inner(p)
            .map(move |i| &mut self.part_vec[i])
    }

    fn inc_items(&mut self) {
        self.n += 1;
    }
}

#[cfg(test)]
mod tests {
    use rand::prelude::SmallRng;
    use rand::Rng;
    use rand::SeedableRng;
    use crate::util::gen_rand_seqs;
    use crate::hash_file::{DnaHashFile, DnaSeqKey, DnaPartition, P};
    use crate::hopscotch::{HopscotchTable, VecPartitionProvider};
    use std::io::Cursor;

    #[test]
    fn test_vec_table() {
        let mut rng = SmallRng::from_seed([1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16]);
        let seqs = gen_rand_seqs(&mut rng, 100, 32, 100);

        let keys: Vec<_> = seqs.iter()
            .map(|s| DnaSeqKey::from_str(s).unwrap())
            .collect();


        let mut p = VecPartitionProvider::<DnaSeqKey, u64>::new(100, P, 32);
        let mut table = HopscotchTable::new(p);


        for (i, key) in keys.iter().enumerate() {
            assert!(table.insert(key.clone(), i as u64).unwrap());
        }

        for (i, key) in keys.iter().enumerate() {
            let v = table.search(key)
                .unwrap()
                .unwrap();
            assert_eq!(v, i as u64);
            println!("{}", i);
        }
    }

    #[test]
    fn test_file_table() {
        let mut rng = SmallRng::from_seed([1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16]);
        let seqs = gen_rand_seqs(&mut rng, 10, 32, 100);

        let keys: Vec<_> = seqs.iter()
            .map(|s| DnaSeqKey::from_str(s).unwrap())
            .collect();

        let sz = std::mem::size_of::<DnaPartition>();

        let mut data = vec![0u8; 100 * sz * P];
        let r = Cursor::new(&mut data);
        let file = DnaHashFile::new(r, 32, 0, 100, P, 2);

        let mut table = HopscotchTable::new(file);


        for (i, key) in keys.iter().enumerate() {
            assert!(table.insert(key.clone(), i as u64).unwrap());
        }

        for (i, key) in keys.iter().enumerate() {
            let v = table.search(key)
                .unwrap()
                .unwrap();
            assert_eq!(v, i as u64);
            println!("{}", i);
        }
    }

}

