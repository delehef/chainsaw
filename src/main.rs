use clap::*;
use newick::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::BufReader;

use anyhow::{anyhow, bail, Context, Result};

fn jaccard<T>(a: &HashSet<T>, b: &HashSet<T>) -> f32
where
    T: Eq + Hash,
{
    a.intersection(b).count() as f32 / a.union(b).count() as f32
}

fn effective_losses(
    a: &HashSet<String>,
    b: &HashSet<String>,
    species: &NewickTree,
    actual_species: &HashSet<usize>,
) -> (usize, usize) {
    let missing_left = a.difference(b).collect::<HashSet<_>>();
    let missing_right = b.difference(a).collect::<HashSet<_>>();

    fn els_oneside(
        missing: &HashSet<&String>,
        species: &NewickTree,
        actual_species: &HashSet<usize>,
    ) -> (usize, usize) {
        fn id2names(xs: &[usize], s: &NewickTree) -> Vec<String> {
            xs.iter()
                .map(|x| s[*x].data.name.as_ref().unwrap().to_owned())
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
        let mut missing = missing
            .iter()
            .map(|x| {
                species
                    .find_leaf(|l| l.name.as_ref().unwrap().as_str() == x.as_str())
                    .unwrap()
            })
            .collect::<Vec<_>>();

        while !missing.is_empty() {
            let mut mrca = missing[0];
            let mut current: Vec<usize> = vec![mrca];
            if log {
                println!("Current: {:?}", id2names(&current, &species));
            }

            'goup: loop {
                let candidates = species
                    .leaves_of(mrca)
                    .into_iter()
                    .filter(|x| actual_species.contains(x))
                    .collect::<Vec<_>>();
                if log {
                    println!("Candidates: {:?}", id2names(&candidates, &species));
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
                    id2names(&current, &species)
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

fn annotate_duplications(t: &mut NewickTree, species_tree: &NewickTree, filter_species: bool) {
    let restricted_species = if filter_species {
        Some(
            t.leaves()
                .filter_map(|l| t[l].data.attrs.get("S").map(|s| s.to_owned()))
                .map(|name| {
                    species_tree
                        .find_leaf(|n| n.name.as_ref().unwrap().as_str() == name.as_str())
                        .unwrap()
                })
                .collect::<HashSet<_>>(),
        )
    } else {
        None
    };
    t.inners().collect::<Vec<_>>().iter().for_each(|n| {
        let species: Vec<HashSet<_>> = t[*n]
            .children()
            .iter()
            .map(|&c| {
                t.leaves_of(c)
                    .iter()
                    .map(|&n| t[n].data.attrs.get("S").map(|s| s.to_owned()).unwrap())
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
                t[*n].data.attrs.insert("D".to_string(), "Y".to_owned());
                t[*n].data.attrs.insert("DCS".to_string(), dcs.to_string());
                t[*n]
                    .data
                    .attrs
                    .insert("ELC".to_string(), elc_all.to_string());
                t[*n]
                    .data
                    .attrs
                    .insert("ELLC".to_string(), elc_large.to_string());
            } else {
                t[*n].data.attrs.insert("D".to_string(), "N".to_owned());
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
                newick::from_filename(&value_t!(args, "species-tree", String).unwrap()).unwrap();
            let mut t = newick::from_string(&std::fs::read_to_string(&filename)?)?;
            annotate_duplications(&mut t, &species_tree, true);
            out.push_str(&t.to_newick());
            out.push('\n');

            File::create(&filename)?
                .write_all(out.as_bytes())
                .context(format!("Cannot write to `{}`", filename))
        }
        _ => unimplemented!(),
    }
}
