#![allow(dead_code)]
use anyhow::*;
use colorsys::{Hsl, Rgb};

use newick::{Newick, NewickTree};
use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;
use std::sync::OnceLock;
use syntesuite::genebook::{FamilyID, Gene, GeneBook};

static ANCESTRAL_QUERY: OnceLock<String> = OnceLock::new();
const LEFTS_QUERY: &str = "select ancestral, direction from genomes where species=? and chr=? and start<? order by start desc limit ?";
const RIGHTS_QUERY: &str = "select ancestral, direction from genomes where species=? and chr=? and start>? order by start asc limit ?";

pub const WINDOW: usize = 15;
pub const GENE_WIDTH: f32 = 15.;
pub const GENE_SPACING: f32 = 5.;
pub const BRANCH_WIDTH: f32 = 20.;
pub const FONT_SIZE: f32 = 10.;

pub type GeneCache = HashMap<String, Gene>;
pub type ColorMap = HashMap<usize, Rgb>;

pub fn jaccard<T: std::hash::Hash + Eq>(x: &HashSet<T>, y: &HashSet<T>) -> f32 {
    x.intersection(y).count() as f32 / x.union(y).count() as f32
}

// Creates a color for a string while trying to ensure it remains readable
pub fn name2color<S: AsRef<str>>(name: S) -> Rgb {
    let bytes: [u8; 16] = md5::compute(name.as_ref().as_bytes()).into();
    let rgb = Rgb::from((bytes[0] as f32, bytes[1] as f32, bytes[2] as f32));

    let mut hsl: Hsl = rgb.into();
    hsl.set_lightness(hsl.lightness().clamp(30., 40.));

    hsl.into()
}

pub fn gene2color(id: &[u8]) -> Rgb {
    let bytes: [u8; 16] = md5::compute(id).into();
    let r = (bytes[0] as f64 / 255.).clamp(0.1, 0.9);
    let g = (bytes[1] as f64 / 255.).clamp(0.1, 0.9);
    let b = (bytes[2] as f64 / 255.).clamp(0.1, 0.9);
    Rgb::new(r, g, b, None)
}

pub fn set_reference(reference: &str) {
    ANCESTRAL_QUERY.set(format!(
        "select ancestral, species, chr, start, direction, left_tail_names, right_tail_names from genomes where {}=?", reference)).unwrap();
}

pub fn make_colormap(tree: &NewickTree, genes: &GeneCache) -> ColorMap {
    let mut colormap = ColorMap::new();
    for l in tree.leaves() {
        if let Some(g) = tree
            .name(l)
            .as_ref()
            .and_then(|name| genes.get(name.as_str()))
        {
            for tg in g.left_landscape.iter().chain(g.right_landscape.iter()) {
                colormap
                    .entry(tg.family)
                    .or_insert_with(|| gene2color(&tg.family.to_ne_bytes()));
            }
        }
    }
    colormap
}

pub fn make_genes_cache(
    t: &NewickTree,
    db_file: &str,
    id_column: &str,
) -> Result<HashMap<String, Gene>> {
    fn reorder_tails(tree: &NewickTree, node: usize, genes: &mut GeneBook) {
        fn reorder_leaves(t: &NewickTree, leave_nodes: &[usize], genes: &mut GeneBook) {
            if leave_nodes.len() < 2 {
                return;
            }
            // Chose a random leaf from the leaves featuring the longest tails as a reference
            let tails = leave_nodes
                .iter()
                .filter_map(|l| t.name(*l))
                .filter_map(|name| genes.get(name).ok())
                .map(|g| {
                    (
                        g.strand,
                        g.left_landscape
                            .iter()
                            .map(|g| g.family)
                            .collect::<Vec<_>>(),
                        g.right_landscape
                            .iter()
                            .map(|g| g.family)
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>();
            let tailsets = tails
                .iter()
                .map(|t| {
                    HashSet::<_>::from_iter(
                        t.1.iter()
                            .chain(t.2.iter())
                            .map(|tg| md5::compute(tg.to_ne_bytes())),
                    )
                })
                .collect::<Vec<_>>();
            // scores: Vec<(id: usize, scores: usize)>
            let scores = tails
                .iter()
                .enumerate()
                .map(|(i, _)| {
                    let mut s = 0;
                    for (j, _) in tails.iter().enumerate() {
                        if i != j {
                            s += (&tailsets[i] & &tailsets[j]).len();
                        }
                    }
                    s
                })
                .collect::<Vec<_>>();

            let ref_id = scores
                .iter()
                .enumerate()
                .max_by_key(|(_, s)| *s)
                .map(|(i, _)| i)
                .unwrap_or(0);
            let ref_left_tail: HashSet<_> = HashSet::from_iter(tails[ref_id].1.iter().cloned());
            let ref_right_tail: HashSet<_> = HashSet::from_iter(tails[ref_id].2.iter().cloned());

            for l_name in leave_nodes.iter().filter_map(|l| t.name(*l)) {
                if let Result::Ok(gene) = genes.get_mut(l_name) {
                    let left_tail: HashSet<FamilyID> =
                        HashSet::from_iter(gene.left_landscape.iter().map(|tg| tg.family));
                    let right_tail: HashSet<FamilyID> =
                        HashSet::from_iter(gene.right_landscape.iter().map(|tg| tg.family));

                    let direct_score =
                        jaccard(&left_tail, &ref_left_tail) + jaccard(&right_tail, &ref_right_tail);
                    let reverse_score =
                        jaccard(&left_tail, &ref_right_tail) + jaccard(&right_tail, &ref_left_tail);

                    if reverse_score > direct_score {
                        std::mem::swap(&mut gene.left_landscape, &mut gene.right_landscape);
                        gene.left_landscape.reverse();
                        gene.right_landscape.reverse();
                        gene.strand.reverse();
                    }
                }
            }
        }

        if tree.is_root(node) || tree.is_duplication(node) {
            let children = tree[node].children();
            let members = children
                .iter()
                .filter(|c| tree[**c].is_leaf())
                .cloned()
                .collect::<Vec<_>>();
            reorder_leaves(tree, &members, genes);

            for c in children.iter().filter(|c| !tree[**c].is_leaf()) {
                reorder_leaves(tree, &tree.leaves_of(*c), genes);
            }
        }

        for c in tree[node].children().iter() {
            reorder_tails(tree, *c, genes);
        }
    }

    let leaves = t.leaves().filter_map(|n| t.name(n)).collect::<Vec<_>>();
    let mut gene_book =
        GeneBook::cached(db_file, WINDOW, id_column, &leaves).map_err(|e| anyhow!(e))?;
    reorder_tails(t, t.root(), &mut gene_book);
    let r = leaves
        .into_iter()
        .map(|g| {
            gene_book.get(g).map(|mut gene| {
                // XXX: left tails are drawn right to left, so they must be reversed
                gene.left_landscape.reverse();
                (g.to_owned(), gene)
            })
        })
        .collect::<Result<HashMap<_, _>>>()?;
    Ok(r)
}

pub fn effective_losses(
    a: &HashSet<usize>,
    b: &HashSet<usize>,
    species: &NewickTree,
    actual_species: &HashSet<usize>,
) -> (usize, usize) {
    let missing_left = a.difference(b).copied().collect::<HashSet<_>>();
    let missing_right = b.difference(a).copied().collect::<HashSet<_>>();

    fn els_oneside(
        missing: &HashSet<usize>,
        species: &NewickTree,
        actual_species: &HashSet<usize>,
    ) -> (usize, usize) {
        fn id2names(xs: &[usize], s: &NewickTree) -> Vec<String> {
            xs.iter()
                .map(|x| s.name(*x).unwrap().to_owned())
                .collect::<Vec<_>>()
        }

        if missing.is_empty() {
            return (0, 0);
        }

        let log = false; // missing.contains(&"Dasypus.novemcinctus".to_string());
        if log {
            println!("\n\n\nProcessing {:?}", missing);
        }

        let mut r_large = 0;
        let mut r_all = 0;
        let mut missing = missing.iter().copied().collect::<Vec<_>>();

        while !missing.is_empty() {
            let mut mrca = missing[0];
            let mut current: Vec<usize> = vec![mrca];
            if log {
                println!("Current: {:?}", id2names(&current, species));
            }

            'goup: loop {
                let candidates = species
                    .leaves_of(mrca)
                    .into_iter()
                    .filter(|x| actual_species.contains(x))
                    .collect::<Vec<_>>();
                if log {
                    println!("Candidates: {:?}", id2names(&candidates, species));
                }

                if !candidates.is_empty() && candidates.iter().all(|x| missing.contains(x)) {
                    current = candidates.clone();
                } else {
                    break 'goup;
                }

                if let Some(new_mrca) = species.parent(mrca) {
                    mrca = new_mrca;
                } else {
                    break 'goup;
                }
            }

            if log {
                eprintln!(
                    "Loss |{}|: {:?}",
                    current.len(),
                    id2names(&current, species)
                );
            }
            r_all += 1;
            if current.len() > 1 {
                r_large += 1;
            }

            missing.retain(|x| !current.contains(x));
        }
        (r_all, r_large)
    }

    let (lr_all, lr_large) = els_oneside(&missing_left, species, actual_species);
    let (rr_all, rr_large) = els_oneside(&missing_right, species, actual_species);
    (lr_all + rr_all, lr_large + rr_large)
}

pub fn capitalize(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        None => String::new(),
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
    }
}
