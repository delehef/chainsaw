use clap::*;
use newick::*;
use rusqlite::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::BufReader;

use anyhow::{anyhow, bail, Context, Result};

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
        let mut missing = missing.into_iter().copied().collect::<Vec<_>>();

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
                        .expect(&format!("{} not found in species tree", name))
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

fn geneize(t: &mut NewickTree, book: &mut GeneBook) -> Result<()> {
    fn get_gene(
        db: &mut Connection,
        protein: &str,
    ) -> std::result::Result<String, rusqlite::Error> {
        db.query_row(
            "SELECT gene FROM genomes WHERE protein=?",
            &[&protein],
            |r| {
                let gene: String = r.get("gene").unwrap();

                Ok(gene.to_owned())
            },
        )
    }

    t.map_leaves(&mut |n| {
        if n.data.name.is_some() {
            n.data.name = Some(book.get(n.data.name.as_ref().unwrap()).unwrap().name)
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
                let tgt = s.next()?.replace(' ', ".").to_owned();
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

fn main() -> Result<()> {
    let args = App::new("Chainsaw")
        .version(clap::crate_version!())
        .author(clap::crate_authors!())
        .arg(
            Arg::with_name("FILE")
                .help("Sets the input file to use")
                .required(true),
        )
        .subcommand(
            SubCommand::with_name("annotate").arg(
                Arg::with_name("species-tree")
                    .short("S")
                    .long("species-tree")
                    .help("The species tree to use")
                    .required(true)
                    .takes_value(true),
            ),
        )
        .subcommand(
            SubCommand::with_name("taxonize")
                .arg(Arg::with_name("mapping").short("m").takes_value(true)),
        )
        .subcommand(
            SubCommand::with_name("geneize")
                .arg(
                    Arg::with_name("database")
                        .short("D")
                        .long("database")
                        .help("The database to use")
                        .required(true)
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("cache-db")
                        .short("c")
                        .long("cache-db")
                        .help("Read the whole database at once"),
                ),
        )
        .subcommand(SubCommand::with_name("compress"))
        .get_matches();
    let filename = value_t!(args, "FILE", String).unwrap();
    println!("Processing {}", filename);
    let mut trees: Vec<NewickTree> = newick::from_filename(&filename)
        .with_context(|| format!("failed to parse {}", filename))?;

    match args.subcommand() {
        ("annotate", Some(margs)) => {
            let mut out = String::new();
            let species_tree =
                newick::one_from_filename(&value_t!(margs, "species-tree", String).unwrap())
                    .unwrap();
            println!("Processing {} trees", trees.len());
            for t in trees.iter_mut() {
                annotate_duplications(t, &species_tree, true);
                out.push_str(&Newick::to_newick(t));
                out.push('\n');
            }

            File::create(&filename)?
                .write_all(out.as_bytes())
                .context(format!("Cannot write to `{}`", filename))
        }
        ("compress", Some(margs)) => {
            let mut out = String::new();
            for t in trees.iter_mut() {
                compress(t).and_then(|_| {
                    out.push_str(&Newick::to_newick(t).replace('\n', ""));
                    out.push('\n');
                    Ok(())
                })?;
            }
            File::create("out.nhx")?
                .write_all(out.as_bytes())
                .context(format!("Cannot write to `out.nhx`"))
        }
        ("taxonize", Some(margs)) => {
            let mut out = String::new();
            let mut err = String::new();
            let map_file = value_t!(margs, "mapping", String).unwrap();

            for t in trees.iter_mut() {
                taxonize(t, &map_file)
                    .and_then(|_| {
                        out.push_str(&Newick::to_newick(t).replace('\n', ""));
                        out.push('\n');
                        Ok(())
                    })
                    .or_else::<(), _>(|_| {
                        err.push_str(&Newick::to_newick(t));
                        err.push('\n');
                        Ok(())
                    });
            }

            File::create("out.nhx")?
                .write_all(out.as_bytes())
                .context(format!("Cannot write to `out.nhx`"))?;
            if !err.is_empty() {
                File::create("err.nhx")?
                    .write_all(err.as_bytes())
                    .context(format!("Cannot write to `err.nhx`"))?
            }
            Ok(())
        }
        ("geneize", Some(margs)) => {
            let mut out = String::new();
            let mut err = String::new();
            let db_filename = margs.value_of("database").unwrap();
            let mut book = if margs.is_present("cache-db") {
                GeneBook::cached(db_filename)
            } else {
                GeneBook::inline(db_filename)
            }?;

            for t in trees.iter_mut() {
                geneize(t, &mut book)
                    .and_then(|_| {
                        out.push_str(&Newick::to_newick(t).replace('\n', ""));
                        out.push('\n');
                        Ok(())
                    })
                    .or_else::<(), _>(|_| {
                        err.push_str(&Newick::to_newick(t).replace('\n', ""));
                        err.push('\n');
                        Ok(())
                    });
            }
            File::create("out.nhx")?
                .write_all(out.as_bytes())
                .context(format!("Cannot write to `out.nhx`"))?;
            if !err.is_empty() {
                File::create("err.nhx")?
                    .write_all(err.as_bytes())
                    .context(format!("Cannot write to `err.nhx`"))?
            }
            Ok(())
        }
        _ => unimplemented!(),
    }
}
