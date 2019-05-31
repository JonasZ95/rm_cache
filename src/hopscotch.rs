use std::marker::PhantomData;
use crate::bucket_info::{BucketInfo, BucketIndexIterator};
use crate::CacheError;


pub trait Empty {
    fn new_empty() -> Self;

    fn is_empty(&self) -> bool;
}

impl Empty for u64 {
    fn new_empty() -> Self {
        u64::max_value()
    }

    fn is_empty(&self) -> bool {
        *self == Self::new_empty()
    }
}

impl HopscotchHash for u64 {
    fn hopscotch_hash(&self) -> u64 {
        *self
    }
}

pub trait HopscotchHash {
    fn hopscotch_hash(&self) -> u64;
}

pub trait Bucket<K, V> {
    type EntryIterator: Iterator<Item=usize>;

    fn get_key(&self) -> &K;
    fn get_key_mut(&mut self) -> &mut K;

    fn get_value(&self) -> &V;
    fn get_value_mut(&mut self) -> &mut V;

    fn fill(&mut self, k: K, v: V) {
        *self.get_key_mut() = k;
        *self.get_value_mut() = v;
    }

    fn has_entry(&self, i: usize) -> bool;
    fn set_entry(&mut self, i: usize);
    fn unset_entry(&mut self, i: usize);
    fn clear_entries(&mut self);
    fn entry_count(&self) -> usize;
    fn first_unset_entry(&self, h: usize) -> Option<usize>;
    fn entry_iter(&self) -> Self::EntryIterator;

    fn is_empty(&self) -> bool;
}

pub trait Partition<B, K, V>
    where B: Bucket<K, V> {
    fn get_bucket(&self, i: usize) -> &B;
    fn get_bucket_mut(&mut self, i: usize) -> &mut B;

    fn mark_as_changed(&mut self);

    fn find<F>(&self, from: usize, to: usize, f: &mut F) -> Option<usize>
        where F: FnMut(&B) -> bool {
        for i in from..=to {
            if f(self.get_bucket(i)) {
                return Some(i);
            }
        }

        None
    }
}


pub trait PartitionProvider<B, K, V>
    where B: Bucket<K, V> {
    type Partition: Partition<B, K, V> + 'static;
    type Error;

    fn n(&self) -> usize;
    fn m(&self) -> usize;
    fn h(&self) -> usize;
    fn parts(&self) -> usize;
    fn part_size(&self) -> usize;

    fn load_part(&mut self, p: usize) -> Result<&Self::Partition, Self::Error>;
    fn load_part_mut(&mut self, p: usize) -> Result<&mut Self::Partition, Self::Error>;

    fn inc_items(&mut self);

    fn load_factor(&self) -> f64 {
        let n = self.n() as f64;
        let m = self.m() as f64;
        n/m
    }

    fn next_part(&self, p: usize) -> usize {
        if (p + 1) == self.parts() {
            0
        } else {
            p + 1
        }
    }

    fn prev_part(&self, p: usize) -> usize {
        if p == 0 {
            self.parts() - 1
        } else {
            0
        }
    }


    fn part_index(&self, i: usize) -> (usize, usize) {
        let p = i / self.part_size();
        let i = i % self.part_size();

        (p, i)
    }

    fn to_index(&self, p: usize, i: usize) -> usize {
        p*self.part_size() + i
    }

    fn get_bucket(&mut self, i: usize) -> Result<&B, Self::Error> {
        let (p, i) = self.part_index(i);
        let part = self.load_part(p)?;

        Ok(part.get_bucket(i))
    }
    fn get_bucket_mut(&mut self, i: usize) -> Result<&mut B, Self::Error> {
        let (p, i) = self.part_index(i);
        let part = self.load_part_mut(p)?;

        Ok(part.get_bucket_mut(i))
    }


    fn find_empty_bucket(&mut self, from: usize, to: usize) -> Result<Option<usize>, Self::Error> {
        self.find(from, to, |b| b.is_empty())
    }

    fn place_slot(&mut self, i: usize, k: K, v: V) -> Result<(), Self::Error> {
        let (p, i) = self.part_index(i);
        let part = self.load_part_mut(p)?;

        part.mark_as_changed();
        let bucket = part.get_bucket_mut(i);
        bucket.fill(k, v);
        Ok(())
    }

    fn clear_slot(&mut self, i: usize, k: K, v: V) -> Result<(K, V), Self::Error> {
        let (p, i) = self.part_index(i);
        let part = self.load_part_mut(p)?;

        part.mark_as_changed();
        let bucket = part.get_bucket_mut(i);

        let k = std::mem::replace(bucket.get_key_mut(), k);
        let v = std::mem::replace(bucket.get_value_mut(), v);
        Ok((k, v))
    }

    fn swap_entry(&mut self, bucket_index: usize, swap_entry: usize, empty_entry: usize) -> Result<(), Self::Error> {
        let (p, i) = self.part_index(bucket_index);
        let part = self.load_part_mut(p)?;

        part.mark_as_changed();
        let bucket = part.get_bucket_mut(i);

        bucket.set_entry(empty_entry);
        bucket.unset_entry(swap_entry);
        Ok(())
    }

    fn set_entry(&mut self, bucket_index: usize, entry: usize) -> Result<(), Self::Error> {
        let (p, i) = self.part_index(bucket_index);
        let part = self.load_part_mut(p)?;

        part.mark_as_changed();
        let bucket = part.get_bucket_mut(i);

        bucket.set_entry(entry);
        Ok(())
    }

    fn find<F>(&mut self, from: usize, to: usize, mut f: F) -> Result<Option<usize>, Self::Error>
        where F: FnMut(&B) -> bool {
        let (mut p, mut from) = self.part_index(from);
        let (p_to, to) = self.part_index(to);

        loop {
            let to_r = if p == p_to && to >= from {
                to
            } else {
                self.part_size() - 1
            };

            let part = self.load_part(p)?;
            if let Some(i) = part.find(from, to_r, &mut f) {
                return Ok(Some(self.to_index(p, i)));
            }

            if p == p_to && to >= from {
                return Ok(None);
            }

            from = 0;

            p = self.next_part(p);
        }
    }
}


pub struct HopscotchTable<PP, B, K, V>
    where
        PP: PartitionProvider<B, K, V>,
        B: Bucket<K, V>,
        K: PartialEq + HopscotchHash + Empty {
    pp: PP,
    marker_b: PhantomData<B>,
    marker_k: PhantomData<K>,
    marker_v: PhantomData<V>,
}


impl<PP, B, K, V> HopscotchTable<PP, B, K, V>
    where
        PP: PartitionProvider<B, K, V>,
        B: Bucket<K, V>,
        K: PartialEq + HopscotchHash + Empty,
        V: Default + Clone {
    pub fn new(pp: PP) -> Self {
        HopscotchTable {
            pp,
            marker_k: PhantomData,
            marker_b: PhantomData,
            marker_v: PhantomData,
        }
    }

    fn find_in_h<F>(&mut self, from: usize, h: usize, mut f: F) -> Result<Option<usize>, PP::Error>
        where F: FnMut(&B) -> bool {
            debug_assert!(h < self.pp.m());
            self.pp.find(from, self.next_index(from, h), f)
    }

    /// Checks if index `j` is in the neighborhood `h` of `ì`
    fn is_in_neighborhood(&self, i: usize, j: usize) -> bool {
        let h = self.pp.h();
        let m = self.pp.m();

        let neighborhood = i + h - 1;
        if neighborhood < m {
            j >= i && j <= neighborhood
        } else {
            j >= i || j <= neighborhood - m
        }
    }
    /*
          let neighborhood = bucket_index + self.p.h() - 1;
        if neighborhood < self.p.m() {
            index >= bucket_index && index <= neighborhood
        } else {
            index >= bucket_index || index <= neighborhood - self.p.m()
        }*/


    fn next_index(&self, i: usize, k: usize) -> usize {
        (i + k) % self.pp.m()
    }


    fn prev_index(&self, i: usize, k: usize) -> usize {
        if k > i {
            self.pp.m() - (k - i)
        } else {
            i - k
        }
    }

    fn diff_index(&self, i: usize, j: usize) -> usize{
        if i > j {
            self.pp.m() - 1 + j - i
        } else {
            j-i
        }
    }


    fn hash(&self, k: &K) -> usize {
        let h = k.hopscotch_hash();
        (h % self.pp.m() as u64) as usize
    }

    pub fn search<'a>(&'a mut self, k: &K) -> Result<Option<V>, PP::Error> {
        let part_size = self.pp.part_size();
        let i = self.hash(k);

        let (p, i) = self.pp.part_index(i);
        let mut j = i;

        let mut entry_iter = {
            let part = self.pp.load_part(p)?;

            let mut entry_iter = part.get_bucket(i)
                .entry_iter();


            loop {
                j = match entry_iter.next() {
                    Some(k) => i + k,
                    None => return Ok(None)
                };

                if j >= part_size {
                    break
                }


                let b = part.get_bucket(j);
                if b.get_key().eq(k) {
                    return Ok(Some(b.get_value().clone()));
                }
            }

            entry_iter
        };

        {
            let p = self.pp.next_part(p);
            let part = self.pp.load_part(p)?;

            loop {
                let b = part.get_bucket(j%part_size);
                if b.get_key().eq(k) {
                    return Ok(Some(b.get_value().clone()));
                }

                j = match entry_iter.next() {
                    Some(k) => i + k,
                    None => return Ok(None)
                };
            }
        }


        Ok(None)
    }

    pub fn insert(&mut self, k: K, v: V) -> Result<bool, PP::Error> {
        const MAX_SEARCH: usize = 100;
        let i = self.hash(&k);
        let h = self.pp.h();

        let empty_slot = self.find_in_h(i, MAX_SEARCH.min(self.pp.m()-1), |b| b.is_empty())?;
        let mut empty_slot = match empty_slot {
            Some(s) => s,
            None => return Ok(false)
        };


        //As long as we don't
        while !self.is_in_neighborhood(i, empty_slot) {
            //Find first bucket which has an entry before the empty slot
            let mut counter = h;
            let mut swap_entry = 0;

            let start = self.prev_index(empty_slot, h-1);
            let swap_bucket = self.find_in_h(start, h - 1,
                                             |b| {
                                                 counter -= 1;
                                                 match b.entry_iter().next() {
                                                     Some(entry) if entry < counter => {
                                                         swap_entry = entry;
                                                         true
                                                     }
                                                     _ => false
                                                 }
                                             })?;

            //If no swap_bucket was found there is no swap candidate
            let swap_bucket = match swap_bucket {
                Some(s) => s,
                None => return Ok(false)
            };

            let swap_slot = self.next_index(swap_bucket, swap_entry);


            self.swap(swap_bucket,
                      swap_slot, swap_entry,
                      empty_slot, counter)?;
            empty_slot = swap_slot;
        }

        self.pp.place_slot(empty_slot, k, v)?;
        self.pp.set_entry(i, self.diff_index(i, empty_slot))?;
        self.pp.inc_items();
        Ok(true)
    }

    fn swap(&mut self, swap_bucket: usize, swap_slot: usize, swap_entry: usize, empty_slot: usize, empty_entry: usize) -> Result<(), PP::Error> {
        //Update the bucket
        {
            self.pp.swap_entry(swap_bucket, swap_entry, empty_entry)?;
        }

        //Clear the swap slot
        let k;
        let v;
        {
            let kv = self.pp.clear_slot(swap_slot, K::new_empty(), V::default())?;
            k = kv.0;
            v = kv.1;
        }

        //Place k and v into the empty slot
        {
            self.pp.place_slot(empty_slot, k, v)?;
        }

        Ok(())
    }
}


pub struct SimpleBucket<K, V> {
    k: K,
    v: V,
    bucket_info: BucketInfo,
}

pub struct VecPartition<K, V> {
    p: Vec<SimpleBucket<K, V>>
}

pub struct VecPartitionProvider<K, V> {
    parts: Vec<VecPartition<K, V>>,
    part_size: usize,
    n: usize,
    m: usize,
    h: usize,
}

impl<K, V> SimpleBucket<K, V>
    where K: Empty,
          V: Default {
    fn empty() -> Self {
        SimpleBucket {
            k: K::new_empty(),
            v: V::default(),
            bucket_info: BucketInfo::default(),
        }
    }
}

impl<K, V> Bucket<K, V> for SimpleBucket<K, V>
    where K: Empty {
    type EntryIterator = BucketIndexIterator;

    fn get_key(&self) -> &K {
        &self.k
    }

    fn get_key_mut(&mut self) -> &mut K {
        &mut self.k
    }

    fn get_value(&self) -> &V {
        &self.v
    }

    fn get_value_mut(&mut self) -> &mut V {
        &mut self.v
    }

    fn has_entry(&self, i: usize) -> bool {
        self.bucket_info.has_entry(i)
    }

    fn set_entry(&mut self, i: usize) {
        self.bucket_info.set_entry(i)
    }

    fn unset_entry(&mut self, i: usize) {
        self.bucket_info.unset_entry(i)
    }

    fn clear_entries(&mut self) {
        self.bucket_info.clear()
    }

    fn entry_count(&self) -> usize {
        self.bucket_info.entry_count()
    }

    fn first_unset_entry(&self, h: usize) -> Option<usize> {
        self.bucket_info.first_unset_entry(h)
    }

    fn entry_iter(&self) -> Self::EntryIterator {
        self.bucket_info.index_iter()
    }

    fn is_empty(&self) -> bool {
        self.k.is_empty()
    }
}

impl<K, V> Partition<SimpleBucket<K, V>, K, V> for VecPartition<K, V>
    where K: HopscotchHash + Empty + PartialEq {
    fn get_bucket(&self, i: usize) -> &SimpleBucket<K, V> {
        &self.p[i]
    }

    fn get_bucket_mut(&mut self, i: usize) -> &mut SimpleBucket<K, V> {
        &mut self.p[i]
    }

    fn mark_as_changed(&mut self) {
        //Changes are always applied directly
    }
}

impl<K, V> PartitionProvider<SimpleBucket<K, V>, K, V> for VecPartitionProvider<K, V>
    where K: HopscotchHash + Empty + PartialEq + 'static,
          V: 'static {
    type Partition = VecPartition<K, V>;
    type Error = ();

    fn n(&self) -> usize {
        self.n
    }

    fn m(&self) -> usize {
        self.m
    }

    fn h(&self) -> usize {
       self.h
    }

    fn parts(&self) -> usize {
        self.parts.len()
    }

    fn part_size(&self) -> usize {
        self.part_size
    }

    fn load_part(&mut self, p: usize) -> Result<&Self::Partition, Self::Error> {
        Ok(&self.parts[p])
    }

    fn load_part_mut(&mut self, p: usize) -> Result<&mut Self::Partition, Self::Error> {
        Ok(&mut self.parts[p])
    }

    fn inc_items(&mut self) {
        self.n += 1;
    }
}

impl<K, V> VecPartitionProvider<K, V>
    where K: Empty,
          V: Default {
    pub fn new(parts: usize, part_size: usize, h: usize) -> Self {
        let p = (0..parts)
            .map(|_| {
                let p = (0..part_size)
                    .map(|_| SimpleBucket::empty())
                    .collect();

                VecPartition{p}
            }).collect();

        VecPartitionProvider{
            m: parts*part_size,
            n: 0,
            h,
            parts: p,
            part_size
        }
    }
}


#[cfg(test)]
mod tests {
    use crate::hopscotch::{VecPartitionProvider, PartitionProvider, SimpleBucket, Bucket, HopscotchTable};
    use rand::prelude::SmallRng;
    use rand::Rng;
    use rand::SeedableRng;

    #[test]
    fn test_vec_part() {
        let mut p = VecPartitionProvider::<u64, bool>::new(2, 4, 2);

        //Properties
        assert_eq!(p.m(), 8);
        assert_eq!(p.n(), 0);
        assert_eq!(p.load_factor(), 0.0);
        assert_eq!(p.part_size(), 4);

        //Partitions
        assert_eq!(p.part_index(0), (0,0));
        assert_eq!(p.part_index(3), (0,3));
        assert_eq!(p.part_index(4), (1, 0));
        assert_eq!(p.part_index(7), (1, 3));


        //Next Part Index
        assert_eq!(p.next_part(0), 1);
        assert_eq!(p.next_part(1), 0);

        //Prev Part Index
        assert_eq!(p.next_part(1), 0);
        assert_eq!(p.next_part(0), 1);

        p.load_part(0).unwrap();
        p.load_part(1).unwrap();

        p.get_bucket(0).unwrap();
        p.get_bucket(3).unwrap();
        p.get_bucket(4).unwrap();
        p.get_bucket(7).unwrap();

        p.place_slot(1, 10, true)
            .unwrap();

        let b1 = p.get_bucket(1)
            .unwrap();
        assert_eq!(b1.k, 10);
        assert_eq!(b1.v, true);

        let (k, v) = p.clear_slot(1, 0, false)
            .unwrap();
        assert_eq!(k, 10);
        assert_eq!(v, true);

        let b1 = p.get_bucket(1)
            .unwrap();
        assert_eq!(b1.k, 0);
        assert_eq!(b1.v, false);
    }

    #[test]
    fn test_vec_part_search() {
        let mut p = VecPartitionProvider::<u64, bool>::new(2, 4, 2);
        p.place_slot(1, 10, true)
            .unwrap();

        assert_eq!(p.find(2, 1, |b| !b.is_empty()).unwrap(),
                   Some(1));

        assert_eq!(p.find(0, 7, |b| !b.is_empty()).unwrap(),
            Some(1));

        assert_eq!(p.find(1, 7, |b| !b.is_empty()).unwrap(),
                   Some(1));

        assert_eq!(p.find(2, 7, |b| !b.is_empty()).unwrap(),
                   None);

        assert_eq!(p.find(2, 2, |b| !b.is_empty()).unwrap(),
                   None);

        assert_eq!(p.find(7, 1, |b| !b.is_empty()).unwrap(),
                   Some(1));

        assert_eq!(p.find(4, 3, |b| !b.is_empty()).unwrap(),
                   Some(1));
    }

    #[test]
    fn test_simple_table() {
        let mut p = VecPartitionProvider::<u64, bool>::new(2, 2, 2);
        let mut table = HopscotchTable::new(p);

        let seq = [1,2,3,4];

        for &s in seq.iter() {
            assert!(table.insert(s, true).unwrap());
        }

        for &s in seq.iter() {
            assert_eq!(table.search(&s).unwrap(), Some(true));
        }

        for &s in seq.iter() {
            let s = s + 4;
            assert_eq!(table.search(&s).unwrap(), None);
        }

        assert!(!table.insert(5, true).unwrap());
    }

    #[test]
    fn test_simple_table1() {
        let mut p = VecPartitionProvider::<u64, bool>::new(2, 2, 2);
        let mut table = HopscotchTable::new(p);

        assert!(table.is_in_neighborhood(3, 3));
        assert!(table.is_in_neighborhood(3, 0));
        assert!(!table.is_in_neighborhood(3, 1));

        assert!(table.is_in_neighborhood(0, 0));
        assert!(table.is_in_neighborhood(0, 1));
        assert!(!table.is_in_neighborhood(0, 2));

        assert!(table.insert(4, true).unwrap());
        assert!(table.insert(8, true).unwrap());
        assert!(!table.insert(12, true).unwrap());

        assert!(table.insert(4, true).unwrap());
    }

    #[test]
    fn test_simple_table2() {
        let mut p = VecPartitionProvider::<u64, bool>::new(2, 4, 3);
        let mut table = HopscotchTable::new(p);

        let seq = [1,1, 2, 2, 3];

        for &s in seq.iter() {
            assert!(table.insert(s, true).unwrap());
        }

        for &s in seq.iter() {
            assert_eq!(table.search(&s).unwrap(), Some(true));
        }

        assert!(!table.insert(3, true).unwrap());
    }

    #[test]
    fn test_simple_table3() {
        let mut p = VecPartitionProvider::<u64, bool>::new(2, 4, 3);
        let mut table = HopscotchTable::new(p);

        let seq = [1,1, 2, 2, 3];

        for &s in seq.iter() {
            assert!(table.insert(s, true).unwrap());
        }

        for &s in seq.iter() {
            assert_eq!(table.search(&s).unwrap(), Some(true));
        }

        assert!(!table.insert(3, true).unwrap());
    }

    #[test]
    fn test_larger_table3() {
        let mut p = VecPartitionProvider::<u64, usize>::new(10, 100, 16);
        let mut table = HopscotchTable::new(p);

        let mut rng = SmallRng::from_seed([1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16]);

        let nums: Vec<u64> = (0..200)
            .map(|_| rng.gen())
                .collect();


        for (i, n) in nums.iter().enumerate() {
            assert!(table.insert(*n, i).unwrap());
        }

        for (i, n) in nums.iter().enumerate() {
            assert_eq!(table.search(&n).unwrap(), Some(i));
        }
    }

}


