#[derive(Clone, Default)]
pub struct BucketInfo(u64);
pub struct BucketIndexIterator(BucketInfo);

impl BucketInfo {
    pub fn new(n: u64) -> BucketInfo {
        BucketInfo(n)
    }

    pub fn value(&self) -> u64 {
        self.0
    }

    pub fn entry_count(&self) -> usize {
        self.0.count_ones() as usize
    }

    pub fn clear(&mut self) {
        self.0 = 0;
    }

    pub fn has_entry(&self, i: usize) -> bool {
        let i = (63 - i) as u64;
        ((1u64 << i) & self.0) != 0
    }

    pub fn set_entry(&mut self, i: usize) {
        let i = (63 - i) as u64;
        self.0 |= 1u64 << i;
    }

    pub fn unset_entry(&mut self, i: usize) {
        let i = (63 - i) as u64;
        self.0 &= !(1u64 << i)
    }

    pub fn is_empty(&self) -> bool {
        self.0 == 0
    }

    pub fn first_unset_entry(&self, h: usize) -> Option<usize> {
        let n = (!self.0).leading_zeros() as usize;
        if n < h {
            Some(n)
        } else {
            None
        }
    }

    pub fn index_iter(&self) -> BucketIndexIterator{
        BucketIndexIterator(self.clone())
    }
}

impl Iterator for BucketIndexIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.0.is_empty() {
            return None;
        }

        let next = (self.0).0.leading_zeros();
        self.0.unset_entry(next as usize);
        Some(next as usize)
    }
}

#[cfg(test)]
mod tests {
    use crate::bucket_info::{BucketInfo};

    #[test]
    fn test_bucket_info() {
        let mut bi = BucketInfo::default();
        const H: usize = 64;

        assert!(bi.is_empty());

        //Check each set
        for i in 0..H {
            assert!(!bi.has_entry(i));
            bi.set_entry(i);
            assert!(bi.has_entry(i));
            bi.unset_entry(i);
            assert!(!bi.has_entry(i));
        }

        //Check find first unset
        for i in 0..H {
            assert_eq!(bi.first_unset_entry(H), Some(i));
            bi.set_entry(i);
        }
        assert_eq!(bi.first_unset_entry(H), None);

        //Bucket iter
        let bi = BucketInfo(u64::max_value());
        for (l, r) in bi.index_iter().zip(0..H) {
            assert_eq!(r, l);
        }

        let bi = BucketInfo(0);
        assert_eq!(bi.index_iter().count(), 0);


        let mut bi = BucketInfo(0);
        bi.set_entry(0);
        bi.set_entry(63);
        let mut iter = bi.index_iter();
        assert_eq!(iter.next(), Some(0));
        assert_eq!(iter.next(), Some(63));
        assert_eq!(iter.next(), None);
    }
}