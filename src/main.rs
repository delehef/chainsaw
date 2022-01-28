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
    restricted_species: &HashSet<usize>,
    log: bool
) -> i32 {
    let missing_left = a.difference(&b).collect::<HashSet<_>>();
    let missing_right = b.difference(&a).collect::<HashSet<_>>();
    if log {
        println!("\n\n");
        dbg!(&missing_left);
        dbg!(&missing_right);
    }

    fn els_oneside(
        missing: &HashSet<&String>,
        species: &Tree,
        restricted_species: &HashSet<usize>,
        log: bool,
    ) -> i32 {
        if missing.is_empty() {
            return 0;
        }
        let mut r = 0;
        if log {
            println!("Processing {:?}", missing);
        }
        let mut missing = missing
            .iter()
            .map(|x| {
                species
                    .find_leaf(|l: &Node| l.name.as_ref().unwrap().as_str() == x.as_str())
                    .unwrap()
            })
            .collect::<Vec<_>>();

        while !missing.is_empty() {
            let mut mrca = Some(missing[0]);

            while mrca.is_some()
                && species
                    .leaves_of(mrca.unwrap())
                    .iter()
                    .all(|x| missing.contains(x) || !restricted_species.contains(x))
            {
                mrca = mrca.and_then(|m| species.parent(m));
            }

            if log {
                eprintln!(
                    "Loss #{}: {:?}",
                    r,
                    species
                        .leaves_of(mrca.unwrap())
                        .iter()
                        .map(|n| species[*n].name.as_ref().unwrap())
                        .collect::<Vec<_>>()
                );
            }
            if let Some(mrca) = mrca {
                    missing.retain(|x| !species.leaves_of(mrca).contains(x));
                    r += 1;
            } else {
                return r + 1;
            }
        }
        r
    }

    let lr = els_oneside(&missing_left, species, restricted_species, log);
    let ll = els_oneside(&missing_right, species, restricted_species, log);
    if log {
        dbg!(lr);
        dbg!(ll);
    }
    lr + ll
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
                let log = (dcs == 0.33846155);
                if log {
                    println!("\n\n\n\n\nDUP: {}", dcs);
                    println!("LEFT: {:#?}", &species[0]);
                    println!("RGHT: {:#?}", &species[1]);
                }
                let elc = effective_losses(
                    &species[0],
                    &species[1],
                    species_tree,
                    restricted_species.as_ref().unwrap(),
                    log,
                );
                if log {
                    dbg!(elc);
                    dbg!(elc.to_string());
                }
                t[*n].data.insert("D".to_string(), "Y".to_owned());
                t[*n].data.insert("DCS".to_string(), dcs.to_string());
                t[*n].data.insert("ELC".to_string(), elc.to_string());
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
