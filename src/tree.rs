use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
use id_tree::{Tree, InsertBehavior, Node, NodeId};
use std::collections::vec_deque::VecDeque;
use crate::util::NucleotideBase;
use std::convert::TryFrom;
use std::fmt::Debug;

/// `CacheTree` represents `LOUDS` tree for DNA sequences, where each node has either 0 or 4 child nodes,
///  The index of the child node relative to It's parent node represents one of the 4 nucleotide bases
pub struct CacheTree<D> {
    rs: RankSelect,
    words: usize,
    k_max: usize,
    data: Vec<Option<D>>
}

fn block_size_for_p(n: u64, p: f64) -> usize {
    //TODO: log2 for u64 instead of f64
    ((n as f64).log2() / p) as usize
}

fn round_to_k(n: usize, k: usize) -> usize {
    ((n + k - 1) / k) * k
}

impl<D> CacheTree<D> {
    /// Creates a new `CacheTree`
    /// `k0` and `k1` represent the block lenght for the 0s and 1s in the rank select data structure
    /// words represents the amount of words saved in the tree
    /// k_max represents the maximum word length which defines the depth for the tree
    pub fn new(bv: BitVec<u8>, k0: usize, k1: usize, words: usize, k_max: usize, data: Vec<Option<D>>) -> Self {
        let rs = RankSelect::new(bv, k0);//ToDo
        CacheTree { rs, words, k_max, data }
    }

    /// check if a node indexed by `n` is valid
    fn is_valid_node(&self, n: u64) -> bool {
        let b = self.rs.bits().len();
        match n {
            //0 not a valid node, just a Bit to keep the Parentheses balanced
            0 => false,
            //Maximum of b Bits
            _ if n >= b => false,
            //Previous Bit has to be a 0
            _ if n > 1 && self.rs.get(n - 1) => false,
            _ => true
        }
    }

    /// returns the root of the true, the tree always contains one node
    fn root(&self) -> usize {
        1
    }

    /// calculates the degree for the node `n`
    pub fn degree(&self, n: u64) -> usize {
        debug_assert!(self.is_valid_node(n));

        if self.is_leaf(n) {
            0
        } else {
            4
        }
    }

    pub fn is_leaf(&self, n: u64) -> bool {
        debug_assert!(self.is_valid_node(n));

        !self.rs.get(n)
    }

    fn nodes_per_level(&self, k: u32) -> u64 {
        (4_u64 * (self.words as u64)).min(4_u64.pow(k))
    }

    fn get_node_from_index(&self, i: u64) -> Option<u64> {
        if self.is_valid_node(i) {
            Some(i)
        } else {
            None
        }
    }

    /// Gets the node next to `n` in the level order
    pub fn next_node(&self, n: u64) -> Option<u64> {
        debug_assert!(self.is_valid_node(n));

        let next = self.degree(n) + 1;
        self.get_node_from_index(n + next as u64)
    }

    /// Try to get `k`-th child of node n
    /// the result may be an invalid index If n is a leaf
    pub fn unchecked_child(&self, n: u64, k: usize, level: u32) -> u64 {
        debug_assert!(self.is_valid_node(n));

        let y = self.rs.rank_1(n)
            .expect("Rank must exist for any node");
        let y = y - 1 + (k as u64);

        //let right = self.nodes_per_level(level);
        self.rs.select_0(y)
            .expect("Node must exist") + 1
    }

    /// Try to get `k`-th child of node n
    pub fn child(&self, n: u64, k: usize, level: u32) -> Option<u64> {
        debug_assert!(self.is_valid_node(n));

        if self.degree(n) <= k {
            return None;
        }

        Some(self.unchecked_child(n, k, level))
    }

    /// Traverse a path given by a sequence of bases
    /// Returns the last node on the path If one was found
    pub fn traverse_path(&self, mut bases: impl Iterator<Item=NucleotideBase>) -> Option<u64> {
        //Starting at root
        let mut n = 1;

        for (level, base) in bases.enumerate() {
            if self.is_leaf(n) {
                return None;
            }

            let ix: usize = base.into();
            n = self.unchecked_child(n, ix, level as u32);
        }

        Some(n)
    }

    pub fn get_data(&self, n: u64) -> &Option<D> {
        let index = self.rs.rank_0(n-1)
            .expect("Index must exist");

        &self.data[index as usize]
    }
}

pub struct MergeTree<D> {
    tree: Tree<(NucleotideBase, Option<D>)>
}

impl<D> MergeTree<D> {
    fn new() -> Self {
        let mut tree = Tree::new();
        tree.insert(Node::new((NucleotideBase::A, None)), InsertBehavior::AsRoot)
            .expect("Root insert should always work");

        MergeTree {
            tree
        }
    }

    fn sequence_at_node(&self, node_id: &NodeId) -> String {
        let mut result = String::new();


        let mut node_id = node_id.clone();

        loop {
            let p = self.tree.get(&node_id)
                .expect("Node must exist")
                .parent();
            if p.is_none() {
                return result;
            }

            let p = p.unwrap();

            let base = self.child_base_iter(&p)
                .find_map(|(base, child_id)| {
                    if let Some(c) = child_id {
                        if c == node_id {
                            return Some(base);
                        }
                    }

                    None
                })
                .expect("Child must exist");

            result.push(base.into());


            node_id = p.clone();
        }
    }

    fn traverse_sequence(&self, sequence: &mut impl Iterator<Item=NucleotideBase>) -> (Option<NucleotideBase>, NodeId) {
        let mut node = self.tree.root_node_id()
            .unwrap().clone();


        //while let Some(base) = sequence.next() {
        for base in sequence {
            let child = self.tree
                .children_ids(&node).unwrap()
                .find(|child_id| {
                    let n = self.tree.get(child_id)
                        .unwrap();
                    n.data().0 == base
                });

            node = match child {
                Some(c) => c.clone(),
                None => return (Some(base), node)
            };
        }


        (None, node)
    }

    fn insert(&mut self, mut sequence: impl Iterator<Item=NucleotideBase>, node_data: D) -> bool {
        let (mut base, mut node) = self.traverse_sequence(&mut sequence);


        while let Some(b) = base {
            let child_node = Node::new((b, None));
            let new_node = self.tree.insert(child_node, InsertBehavior::UnderNode(&node))
                .expect("Child insert should always work");
            self.tree.sort_children_by_key(&node, |n| n.data().0)
                .expect("Sorting Node must exist");
            node = new_node;

            base = sequence.next();
        }

        let mut data = self.tree.get_mut(&node)
            .expect("Last inserted node must exist")
            .data_mut();

        //If there is already data that doesn't work
        if data.1.is_some() {
            return false;
        }

        data.1 = Some(node_data);
        true
    }

    pub fn is_leaf(&self, node_id: &NodeId) -> bool {
        self.tree.children(node_id)
            .expect("Node must exists")
            .count() == 0
    }

    pub fn child_base_iter(&self, node_id: &NodeId) -> MergeNodeBaseIterator<D> {
        let node = self.tree.get(node_id)
            .expect("Node must exist");

        MergeNodeBaseIterator {
            merge_tree: &self.tree,
            merge_node: node,
            ix: 0,
            child_ix: 0,
        }
    }

    pub fn node_count(&self) -> usize {
        self.tree.traverse_level_order(self.tree.root_node_id().expect("Root must exist"))
            .expect("Root must exist")
            .count()
    }

    pub fn get_data(&self, node_id: &NodeId) -> &Option<D> {
        let node = self.tree.get(node_id)
            .expect("Node must exist");

        &node.data().1
    }
}


pub trait MergeWriter<D> {
    type Error;
    fn write_data(&mut self, d: Option<D>) -> Result<(), Self::Error>;
    fn write_node(&mut self, k: usize) -> Result<(), Self::Error>;
}

#[derive(Default)]
pub struct VecMergeWriter<D> {
    v: Vec<Option<D>>,
    bv: BitVec<u8>
}

impl<D> MergeWriter<D> for VecMergeWriter<D> {
    type Error = ();
    fn write_data(&mut self, d: Option<D>) -> Result<(), Self::Error> {
        self.v.push(d);
        Ok(())
    }

    fn write_node(&mut self, k: usize)  -> Result<(), Self::Error>{
        if self.bv.is_empty() {
            self.bv.push(true);
        }


        for _ in 0..k {
            self.bv.push(true);
        }
        self.bv.push(false);
        Ok(())
    }
}

pub struct MergeNodeBaseIterator<'a, D> {
    merge_tree: &'a Tree<(NucleotideBase, Option<D>)>,
    merge_node: &'a Node<(NucleotideBase, Option<D>)>,
    ix: u8,
    child_ix: u8,
}

impl<'a, D> Iterator for MergeNodeBaseIterator<'a, D> {
    type Item = (NucleotideBase, Option<NodeId>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.ix > 3 {
            return None;
        }

        let ix = self.ix as usize;
        let base = NucleotideBase::try_from(ix)
            .expect("Base must exist");
        self.ix += 1;

        let t = &self.merge_tree;
        let n = self.merge_node;

        if let Some(child_id) = n.children().get(self.child_ix as usize) {
            let child = t.get(child_id)
                .expect("Child must exist");

            if child.data().0 == base {
                self.child_ix += 1;
                let child_id = Some(child_id.clone());
                return Some((base, child_id));
            }
        }

        Some((base, None))
    }
}

struct MergeNodeContext {
    new: bool,
    index: u64,
    node_mt: Option<NodeId>,
}

pub struct MergeContext<D, W>
    where W: MergeWriter<D> {
    mt: MergeTree<D>,
    t: CacheTree<D>,
    q: VecDeque<MergeNodeContext>,
    index: u64,
    node_t: Option<u64>,
    pending_child: u64,
    w: W
}

impl<D, W> MergeContext<D, W>
    where D: Clone + Debug,
          W: MergeWriter<D>
{
    pub fn new(w: W, mt: MergeTree<D>, t: CacheTree<D>) -> Self {
        let mut q = VecDeque::new();

        //Add merge tree root
        if let Some(root_id) = mt.tree.root_node_id() {
            q.push_back(MergeNodeContext {
                new: false,
                index: 0,
                node_mt: Some(root_id.clone()),
            });
        }

        let node_t = Some(1);

        MergeContext {
            mt,
            t,
            q,
            index: 0,
            //ToDO: verify this has nodes
            node_t: Some(1),
            pending_child: 1,
            w
        }
    }


    pub fn merge(&mut self) -> Result<(), W::Error> {
        loop {
            if self.q.is_empty() && self.node_t.is_none() {
                break Ok(());
            }

            let mut t_merged = false;
            while let Some((t, mt)) = self.deque(self.index, self.node_t) {
                self.merge_node(t, mt)?;
                if t.is_some() {
                    t_merged = true;
                }
            }

            if !t_merged && self.node_t.is_some() {
                self.merge_node(self.node_t, None)?;
                t_merged = true;
            }

            if t_merged {
                self.pending_child -= 1;
            }

            if let Some(n) = self.node_t {
                self.node_t = self.t.next_node(n);
            }
            self.index += 1;
        }
    }

    pub fn output(self) -> W {
        self.w
    }

    fn deque(&mut self, index: u64, t: Option<u64>) -> Option<(Option<u64>, Option<NodeId>)> {
        if !self.q.iter()
            .any(|ctx| ctx.index == index) {
            return None;
        }

        let ctx = self.q.pop_front()
            .expect("Queue must have an element");

        if ctx.new {
            Some((None, ctx.node_mt.clone()))
        } else {
            Some((Some(t.expect("Must have t node")), ctx.node_mt.clone()))
        }
    }

    fn merge_node(&mut self, t: Option<u64>, mt: Option<NodeId>) -> Result<(), W::Error> {
        match (t, mt) {
            (None, None) => {
                self.w.write_node(0)?;
                self.w.write_data(None)
            },
            (Some(t), None) => {
                let k = self.t.degree(t);
                self.pending_child += k as u64;

                self.w.write_node(k)?;
                self.w.write_data(self.t.get_data(t).clone())
            }
            (None, Some(mt)) => {
                if self.mt.is_leaf(&mt) {
                    self.w.write_node(0)?;
                } else {
                    self.w.write_node(4)?;
                    for (_, child_id) in self.mt.child_base_iter(&mt) {
                        self.q.push_back(MergeNodeContext {
                            new: true,
                            index: self.index + self.pending_child,
                            node_mt: child_id,
                        })
                    }
                }


                self.w.write_data(self.mt.get_data(&mt).clone())
            }
            (Some(t), Some(mt)) => {
                let new = self.t.is_leaf(t);

                let nodes = if new && self.mt.is_leaf(&mt) {
                    0
                } else {
                    4
                };

                for (base, child_id) in self.mt.child_base_iter(&mt) {
                    self.q.push_back(MergeNodeContext {
                        new,
                        index: self.index + self.pending_child,
                        node_mt: child_id,
                    });


                    if !new {
                        self.pending_child += 1;
                    }
                }

                self.w.write_node(nodes)?;

                let data = match (self.t.get_data(t), self.mt.get_data(&mt)) {
                    (None, None) => None,
                    (None, Some(d)) => Some(d.clone()),
                    (Some(d), _) => Some(d.clone()),
                };
                self.w.write_data(data)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::tree::{CacheTree, MergeTree, MergeContext, VecMergeWriter};
    use bv::{bit_vec, BitVec};
    use lazy_static::lazy_static;
    use rand::prelude::SmallRng;
    use rand::Rng;
    use rand::SeedableRng;
    use std::convert::TryFrom;
    use crate::util::{NucleotideBase, DnaStrIter, dna_seq};

    fn new_cache_tree<D>(bv: BitVec<u8>, words: usize, v: Vec<Option<D>>) -> CacheTree<D> {
        CacheTree::new(bv, 1, 1, words, 100, v)
    }

    lazy_static! {
        static ref LOUDS_ROOT: CacheTree<bool>  = new_cache_tree(bit_vec![true, false], 1, vec![Some(true)]);
        static ref LOUDS_SINGLE: CacheTree<bool>  = new_cache_tree(bit_vec![true, true, false, false], 1, vec![Some(true)]);
        static ref LOUDS_DOUBLE: CacheTree<bool>  = new_cache_tree(bit_vec![true, true, true, false, false, false], 1, vec![Some(true)]);
        static ref LOUDS_2_LEVEL: CacheTree<bool>  = new_cache_tree(bit_vec![true, true, true, true, true, false, false, false, false, false], 1, vec![Some(true)]);
    }

    #[test]
    fn test_louds_child() {
        assert_eq!(LOUDS_ROOT.is_leaf(1), true);
        assert_eq!(LOUDS_ROOT.child(1, 0, 1), None);

        assert_eq!(LOUDS_SINGLE.is_leaf(1), false);
        assert_eq!(LOUDS_SINGLE.child(1, 0, 1), Some(3));

        assert_eq!(LOUDS_SINGLE.is_leaf(3), true);
        assert_eq!(LOUDS_SINGLE.child(3, 0, 1), None);

        //Two Child
        assert_eq!(LOUDS_DOUBLE.is_leaf(1), false);
        assert_eq!(LOUDS_DOUBLE.child(1, 0, 1), Some(4));
        assert_eq!(LOUDS_DOUBLE.child(1, 1, 1), Some(5));
    }

    #[test]
    fn test_louds_traverse() {
        assert!(LOUDS_2_LEVEL.traverse_path(dna_seq("A")).is_some());
        assert!(LOUDS_2_LEVEL.traverse_path(dna_seq("C")).is_some());
        assert!(LOUDS_2_LEVEL.traverse_path(dna_seq("G")).is_some());
        assert!(LOUDS_2_LEVEL.traverse_path(dna_seq("T")).is_some());
        assert!(LOUDS_2_LEVEL.traverse_path(dna_seq("AA")).is_none());
    }

    #[test]
    fn test_merge_tree() {
        let mut merge_tree = MergeTree::<usize>::new();

        assert!(merge_tree.insert(dna_seq("AA"), 1));
        assert_eq!(3, merge_tree.node_count());
        assert!(!merge_tree.insert(dna_seq("AA"), 1));

        assert!(merge_tree.insert(dna_seq("AC"), 2));
        assert_eq!(4, merge_tree.node_count());

        assert!(merge_tree.insert(dna_seq("CC"), 3));
        assert_eq!(6, merge_tree.node_count());


    }

    #[test]
    fn test_merge_tree2() {
        let mut mt = MergeTree::<usize>::new();

        let seq = [
            "AA", "CC", "GG", "TT"
        ];

        for (i, &s) in seq.iter().enumerate() {
            let mut dna_seq = dna_seq(s);
            assert!(mt.insert(&mut dna_seq, i))
        }

        for (i, &s) in seq.iter().enumerate() {
            let mut dna_seq = dna_seq(s);
            let (_, node_id) = mt.traverse_sequence(&mut dna_seq);

            let node = mt.tree.get(&node_id)
                .expect("Node must exist");

            let data = node.data().1
                .expect("Data must exist");

            assert_eq!(data, i);
        }
    }

    #[test]
    fn test_child_iter() {
        let mut mt = MergeTree::<usize>::new();

        assert!(mt.insert(dna_seq("AA"), 1));
        assert!(mt.insert(dna_seq("AC"), 2));

        let mut seq_iter = dna_seq("");
        let root = mt.traverse_sequence(&mut seq_iter).1;

        let iter = mt.child_base_iter(&root);
        for (base, child) in iter {
            if base == NucleotideBase::A {
                assert!(child.is_some());
            } else {
                assert!(child.is_none());
            }
        }
    }

    #[test]
    fn test_merge() {
        let louds: CacheTree<usize> = CacheTree::new(bit_vec![true, true, true, true, true, false, false, false, false, false],
                                   1, 1, 0, 1,
                                   vec![None, None, None, None, None]);

        let mut mt = MergeTree::<usize>::new();
        mt.insert(dna_seq("AA"), 1);
        mt.insert(dna_seq("CC"), 2);
        mt.insert(dna_seq("GG"), 3);
        mt.insert(dna_seq("TT"), 4);

        let mut merge = MergeContext::new( VecMergeWriter::default(), mt, louds);
        merge.merge().unwrap();

        let output = merge.output();
        assert_eq!(output.bv.len(), 2*(1+4+16));
    }

    #[test]
    fn test_merge2() {
        let louds: CacheTree<usize> = CacheTree::new(bit_vec![true, true, true, true, true, false, false, false, false, false],
                                                     1, 1, 0, 1,
                                                     vec![None, None, None, None, None]);

        let mut merge = MergeContext::new( VecMergeWriter::default(), MergeTree::<usize>::new(), louds);
        merge.merge().unwrap();

        let output = merge.output();
        assert_eq!(output.bv.len(), 2*(1+4));
    }


    #[test]
    fn test_large_merge() {
        let louds: CacheTree<usize> = CacheTree::new(bit_vec![true, true, true, true, true, false, false, false, false, false],
                                                     1, 1, 0, 1,
                                                     vec![None, None, None, None, None]);

        let mut rng = SmallRng::from_seed([1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16]);
        let n = 100;
        let k = 100;
        let seqs: Vec<String> = (0..n)
            .map(|_| {
                let l = rng.gen_range(32, k);
                (0..l)
                    .map(|_| {
                        let idx: usize = rng.gen_range(0, 4);
                        let c: char = NucleotideBase::try_from(idx)
                            .unwrap()
                            .into();
                        c
                    }).collect()
            }).collect();



        let mut mt = MergeTree::<usize>::new();
        for (i, s) in seqs.iter().enumerate() {
            let iter = DnaStrIter::new(&s);
            mt.insert(iter, i+1);
        }

        for (i, s) in seqs.iter().enumerate() {
            let mut dna_seq = DnaStrIter::new(&s);
            let (_, node_id) = mt.traverse_sequence(&mut dna_seq);

            let node = mt.tree.get(&node_id)
                .expect("Node must exist");

            let data = node.data().1
                .expect("Data must exist");

            assert_eq!(data, i+1);
        }

        let mut merge = MergeContext::new( VecMergeWriter::default(), mt, louds);
        merge.merge().unwrap();
        let output = merge.output();

        let louds = CacheTree::new(output.bv,
                                   1, 1, 0, 1,
                                   output.v);

        for (i, s) in seqs.iter().enumerate() {
            let mut dna_seq = DnaStrIter::new(&s);
            let n = louds.traverse_path(dna_seq)
                .unwrap();
            let v = louds.get_data(n)
                .unwrap();
            assert_eq!(v, i+1);
        }

    }
}

