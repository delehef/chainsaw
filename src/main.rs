use clap::*;
use newick::*;
use std::collections::HashSet;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::thread::current;

use anyhow::{Context, Result};

fn jaccard<T>(a: &HashSet<T>, b: &HashSet<T>) -> f32
where
    T: Eq + Hash,
{
    a.intersection(b).collect::<Vec<_>>().len() as f32 / a.union(b).collect::<Vec<_>>().len() as f32
}

fn effective_losses(
    a: &HashSet<String>,
    b: &HashSet<String>,
    species: &Tree,
    actual_species: &HashSet<usize>,
) -> (usize, usize) {
    let missing_left = a.difference(&b).collect::<HashSet<_>>();
    let missing_right = b.difference(&a).collect::<HashSet<_>>();

    fn els_oneside(
        missing: &HashSet<&String>,
        species: &Tree,
        actual_species: &HashSet<usize>,
    ) -> (usize, usize) {
        fn id2names(xs: &[usize], s: &Tree) -> Vec<String> {
            xs.iter()
                .map(|x| s[*x].name.as_ref().unwrap().to_owned())
                .collect::<Vec<_>>()
        }

        if missing.is_empty() {
            return (0, 0);
        }
        let mut r_large = 0;
        let mut r_all = 0;
        let mut missing = missing
            .iter()
            .map(|x| {
                species
                    .find_leaf(|l: &Node| l.name.as_ref().unwrap().as_str() == x.as_str())
                    .unwrap()
            })
            .collect::<Vec<_>>();
        // println!("\n\n\nProcessing {:?}", id2names(&missing, &species));

        while !missing.is_empty() {
            let mut mrca = missing[0];
            let mut current: Vec<usize> = vec![mrca];
            // println!("Current: {:?}", id2names(&current, &species));

            'goup: loop {
                let candidates = species
                    .leaves_of(mrca)
                    .into_iter()
                    .filter(|x| actual_species.contains(x))
                    .collect::<Vec<_>>();
                // println!("Candidates: {:?}", id2names(&candidates, &species));

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

            // eprintln!(
            //     "Loss |{}|: {:?}",
            //     current.len(),
            //     id2names(&current, &species)
            // );
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

fn annotate_duplications(t: &mut Tree, species_tree: &Tree, filter_species: bool) {
    let restricted_species = if filter_species {
        let present_species = t
            .leaf_names()
            .iter()
            .map(|(_, n)| n.unwrap().split('#').nth(1).unwrap().to_owned())
            .map(|name| {
                species_tree
                    .find_leaf(|n| n.name.as_ref().unwrap().as_str() == name.as_str())
                    .unwrap()
            })
            .collect::<HashSet<_>>();
        Some(present_species)
    } else {
        None
    };
    t.inners().collect::<Vec<_>>().iter().for_each(|n| {
        let children = t[*n].children.as_ref().unwrap();
        let genes: Vec<HashSet<_>> = children
            .iter()
            .map(|&c| {
                t.leaves_of(c)
                    .iter()
                    .map(|&n| {
                        t[n].name
                            .as_ref()
                            .unwrap()
                            .split('#')
                            .nth(0)
                            .unwrap()
                            .to_owned()
                    })
                    .collect()
            })
            .collect();
        let species: Vec<HashSet<_>> = children
            .iter()
            .map(|&c| {
                t.leaves_of(c)
                    .iter()
                    .map(|&n| {
                        t[n].name
                            .as_ref()
                            .unwrap()
                            .split('#')
                            .nth(1)
                            .unwrap()
                            .to_owned()
                    })
                    .collect()
            })
            .collect();
        if species.len() == 2 {
            let d = species.iter().skip(1).any(|x| !x.is_disjoint(&species[0]));
            if d {
                let dcs = jaccard(&species[0], &species[1]);
                let (elc_all, elc_large) = effective_losses(
                    &species[0],
                    &species[1],
                    species_tree,
                    restricted_species.as_ref().unwrap(),
                );
                t[*n].data.insert("D".to_string(), "Y".to_owned());
                t[*n].data.insert("DCS".to_string(), dcs.to_string());
                t[*n].data.insert("ELC".to_string(), elc_all.to_string());
                t[*n].data.insert("ELCL".to_string(), elc_large.to_string());
            } else {
                t[*n].data.insert("D".to_string(), "N".to_owned());
            }
        }
    });
}

fn main() -> Result<()> {
    let args = App::new("Chainsaw")
        .version(clap::crate_version!())
        .author(clap::crate_authors!())
        .arg(
            Arg::with_name("FILE")
                .help("Sets the input file to use")
                .required(true),
        )
        .arg(
            Arg::with_name("species-tree")
                .short("S")
                .long("species-tree")
                .help("The species tree to plot against")
                .takes_value(true),
        )
        .subcommand(SubCommand::with_name("annotate"))
        .get_matches();
    let filename = value_t!(args, "FILE", String).unwrap();
    println!("Processing {}", filename);

    match args.subcommand() {
        ("annotate", _) => {
            let mut out = String::new();
            let species_tree =
                Tree::from_filename(&value_t!(args, "species-tree", String).unwrap()).unwrap();
            for l in BufReader::new(File::open(&filename)?).lines() {
                let mut t = Tree::from_string(&l?)?;
                annotate_duplications(&mut t, &species_tree, true);
                out.push_str(&t.to_string());
                out.push_str("\n");
            }

            File::create(&filename)?
                .write_all(out.as_bytes())
                .context(format!("Cannot write to `{}`", filename))
        }
        _ => unimplemented!(),
    }
}
