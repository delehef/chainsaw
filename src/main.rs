use clap::{Parser, Subcommand};
use newick::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::BufReader;

use anyhow::{anyhow, Context, Result};

mod genebook;
use genebook::GeneBook;

fn jaccard<T>(a: &HashSet<T>, b: &HashSet<T>) -> f32
where
    T: Eq + Hash,
{
    a.intersection(b).count() as f32 / a.union(b).count() as f32
}

fn effective_losses(
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

fn annotate_duplications(t: &mut NewickTree, species_tree: &NewickTree, filter_species: bool) {
    let restricted_species = if filter_species {
        Some(
            t.leaves()
                .filter_map(|l| t[l].data.attrs.get("S").map(|s| s.to_owned()))
                .map(|name| {
                    species_tree
                        .find_leaf(|n| n.name.as_ref().unwrap().as_str() == name.as_str())
                        .unwrap_or_else(|| panic!("{} not found in species tree", name))
                })
                .collect::<HashSet<_>>(),
        )
    } else {
        None
    };
    t.inners().collect::<Vec<_>>().iter().for_each(|n| {
        let species: Vec<HashSet<usize>> = t[*n]
            .children()
            .iter()
            .map(|&c| {
                t.leaves_of(c)
                    .iter()
                    .map(|&n| t[n].data.attrs.get("S").map(|s| s.to_owned()).unwrap())
                    .map(|s| {
                        species_tree
                            .find_leaf(|l| l.name.as_ref().unwrap().as_str() == s.as_str())
                            .unwrap()
                    })
                    .collect()
            })
            .collect();
        if species.len() >= 2 {
            let mrcas = species
                .iter()
                .map(|ss| species_tree.mrca(ss).unwrap())
                .collect::<Vec<_>>();
            let mut d = false;
            'find_d: for (i, &m1) in mrcas.iter().enumerate() {
                for &m2 in mrcas.iter().skip(i + 1) {
                    if species_tree.ascendance(m1).contains(&m2)
                        || species_tree.ascendance(m2).contains(&m1)
                    {
                        d = true;
                        break 'find_d;
                    }
                }
            }

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

fn annotate_mrcas(t: &mut NewickTree, species_tree: &NewickTree) -> Result<()> {
    for n in t.inners().collect::<Vec<_>>().into_iter() {
        let species: HashSet<usize> = t
            .leaves_of(n)
            .iter()
            .flat_map(|&c| {
                t.leaves_of(c)
                    .into_iter()
                    .map(|n| t[n].data.attrs.get("S").map(|s| s.to_owned()).unwrap())
                    .map(|s| {
                        species_tree
                            .find_leaf(|l| l.name.as_ref().unwrap().as_str() == s.as_str())
                            .ok_or_else(|| anyhow!(format!("{} not found in species tree", s)))
                    })
            })
            .collect::<Result<HashSet<_>>>()?;
        let mrca = species_tree.mrca(&species).unwrap();
        t[n].data.attrs.insert(
            "S".to_owned(),
            species_tree[mrca].data.name.as_ref().unwrap().to_owned(),
        );
    }
    Ok(())
}

fn geneize(t: &mut NewickTree, book: &mut GeneBook) -> Result<()> {
    t.map_leaves(&mut |n| {
        if n.data.name.is_some() {
            n.data.name = Some(
                book.get(n.data.name.as_ref().unwrap())
                    .unwrap_or_else(|_| panic!("Cannot find {:?}", n.data.name))
                    .gene,
            )
        }
    });

    Ok(())
}

fn speciesize(t: &mut NewickTree, book: &mut GeneBook) -> Result<()> {
    t.map_leaves(&mut |n| {
        if n.data.name.is_some() {
            n.data.attrs.insert(
                "S".to_owned(),
                book.get(n.data.name.as_ref().unwrap())
                    .unwrap_or_else(|_| panic!("Cannot find {:?}", n.data.name))
                    .species,
            );
        }
    });

    Ok(())
}

fn taxonize(t: &mut NewickTree, map_file: &str) -> Result<()> {
    let map = BufReader::new(File::open(map_file)?)
        .lines()
        .filter_map(|l| {
            l.ok().and_then(|l| {
                let mut s = l.split('\t');
                let src = s.next()?.to_owned().parse::<usize>().ok()?;
                let tgt = s.next()?.replace(' ', ".");
                Some((src, tgt))
            })
        })
        .collect::<HashMap<usize, String>>();

    for l in t.nodes_mut() {
        let taxon_id = l.data.attrs.get("T").and_then(|s| s.parse::<usize>().ok());
        // .map_err(|e| {
        //     println!("#children: {}", l.children().len());
        //     e
        // });
        // .context(format!("while parsing taxon {:?}", &l.data.name))?;
        if let Some(taxon_id) = taxon_id {
            if let Some(species) = map.get(&taxon_id) {
                l.data.attrs.insert("S".to_string(), species.to_owned());
            } else {
                eprintln!("`{}` has no match in the reference file", taxon_id);
            }
        } else {
            eprintln!("Node `{:?}` has no taxon specified", l.data.name);
        }
    }
    Ok(())
}

fn compress(t: &mut NewickTree) -> Result<()> {
    while t[t.root()].children().len() == 1 {
        println!("Compressing");
        t.set_root(t[t.root()].children()[0]);
    }
    Ok(())
}

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// sets the input file
    #[clap(value_parser)]
    infile: String,

    /// output file name
    #[clap(value_parser, short = 'o', long = "out")]
    outfile: Option<String>,

    #[clap(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// add duplications & ancestral nodes anntations to a tree
    Annotate {
        /// the species tree to use
        #[clap(value_parser, short = 'S')]
        species_tree: String,
    },

    /// annotate leaves in a tree with their species
    Speciesize {
        /// the database containing the id/species mapping
        #[clap(value_parser, short = 'D')]
        database: String,

        /// if set, cache the database in memory
        #[clap(value_parser, short = 'c')]
        cache_db: bool,
    },

    /// convert ids from proteins to genes
    Geneize {
        /// the database containing the proteins/genes mapping
        #[clap(value_parser, short = 'D')]
        database: String,

        /// if set, cache the database in memory
        #[clap(value_parser, short = 'c')]
        cache_db: bool,
    },

    ///
    Taxonize {
        #[clap(value_parser)]
        mapping: String,
    },

    /// compress root nodes with a single child
    Compress,
}

fn main() -> Result<()> {
    let args = Args::parse();
    println!("Parsing {}", &args.infile);
    let mut trees: Vec<NewickTree> = newick::from_filename(&args.infile)
        .with_context(|| format!("failed to parse {}", &args.infile))?;

    match args.command {
        Command::Annotate { species_tree } => {
            let mut out = String::new();
            let species_tree = newick::one_from_filename(&species_tree)
                .context(format!("while parsing {}", &species_tree))?;
            println!("Processing {} trees", trees.len());
            for t in trees.iter_mut() {
                annotate_duplications(t, &species_tree, true);
                annotate_mrcas(t, &species_tree)?;
                out.push_str(&Newick::to_newick(t));
                out.push('\n');
            }

            File::create(&args.infile)?
                .write_all(out.as_bytes())
                .context(format!("Cannot write to `{}`", &args.infile))
        }
        Command::Compress => {
            let mut out = String::new();
            for t in trees.iter_mut() {
                compress(t).map(|_| {
                    out.push_str(&Newick::to_newick(t).replace('\n', ""));
                    out.push('\n');
                })?;
            }
            let outfile = args.outfile.unwrap_or(args.infile);
            File::create(&outfile)?
                .write_all(out.as_bytes())
                .with_context(|| anyhow!("cannot write to `{}`", &outfile))
        }
        Command::Speciesize { database, cache_db } => {
            let mut out = String::new();
            let mut book = if cache_db {
                GeneBook::in_memory(&database)
            } else {
                GeneBook::inline(&database)
            }?;

            for t in trees.iter_mut() {
                speciesize(t, &mut book).map(|_| {
                    out.push_str(&Newick::to_newick(t).replace('\n', ""));
                    out.push('\n');
                })?
            }

            let outfile = args.outfile.unwrap_or(args.infile);
            File::create(&outfile)?
                .write_all(out.as_bytes())
                .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            Ok(())
        }
        Command::Taxonize { mapping } => {
            let mut out = String::new();

            for t in trees.iter_mut() {
                taxonize(t, &mapping).map(|_| {
                    out.push_str(&Newick::to_newick(t).replace('\n', ""));
                    out.push('\n');
                })?
            }

            let outfile = args.outfile.unwrap_or(args.infile);
            File::create(&outfile)?
                .write_all(out.as_bytes())
                .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            Ok(())
        }
        Command::Geneize { database, cache_db } => {
            let mut out = String::new();
            let mut book = if cache_db {
                GeneBook::in_memory(&database)
            } else {
                GeneBook::inline(&database)
            }?;

            for t in trees.iter_mut() {
                geneize(t, &mut book).map(|_| {
                    out.push_str(&Newick::to_newick(t).replace('\n', ""));
                    out.push('\n');
                })?
            }

            let outfile = args.outfile.unwrap_or(args.infile);
            File::create(&outfile)?
                .write_all(out.as_bytes())
                .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            Ok(())
        }
    }
}
