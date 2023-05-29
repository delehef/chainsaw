use clap::{Parser, Subcommand};
use newick::*;
use std::fs::File;
use std::io::prelude::*;

use anyhow::{anyhow, bail, Context, Result};

use syntesuite::genebook::GeneBook;
mod actions;
mod utils;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// sets the input file
    #[clap(value_parser)]
    infile: String,

    /// output file name
    #[clap(value_parser, short = 'o', long = "out", global = true)]
    outfile: Option<String>,

    #[clap(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// add duplications & ancestral nodes anntations to a tree
    Annotate {
        /// the species tree to use
        #[clap(value_parser, short = 'S', long)]
        species_tree: String,
    },

    /// annotate leaves in a tree with their species
    Speciesize {
        /// the database containing the id/species mapping
        #[clap(value_parser, short = 'D', long)]
        database: String,

        /// if set, cache the database in memory
        #[clap(value_parser, long)]
        cache_db: bool,

        /// the database column corresponding to leaf IDs in the tree
        #[clap(value_parser, default_value_t = String::from("id"), long)]
        id: String,

        /// the database column containing the species
        #[clap(value_parser, default_value_t = String::from("species"))]
        species: String,
    },

    ///
    Taxonize {
        #[clap(value_parser)]
        mapping: String,
    },

    /// compress root nodes with a single child
    Compress,

    /// convert a newick-formatted tree to a phyl-formatted tree
    ToPhy,

    /// list the named leaves of the given tree
    Leaves,

    /// list the names nodes of the given tree
    Nodes,

    /// normalize a species tree according to ENSEMBL naming conventions
    Normalize,

    /// prune the specified nodes from the tree
    Prune {
        /// nodes to recursively remove
        #[clap(value_parser)]
        remove: Vec<String>,
    },

    /// ensure that the provided tree only contains binary speciations
    Binarize {},

    /// rename the leaves of a tree following the given mapping file
    Rename {
        #[clap(value_parser, short, long = "mapping")]
        mapping_file: String,

        /// if set, use as a separator in `mapping`; otherwise, split on space
        #[clap(value_parser, short, long)]
        separator: Option<String>,
    },

    /// format a newick tree in a humean-readableish way
    Format,

    /// strip the given informations from the given trees
    Strip {
        #[clap(value_enum)]
        to_strip: Vec<actions::Strippable>,
    },
    Tune {
        #[clap(value_parser, short, long = "reference")]
        reference_tree: String,
    },
}

fn main() -> Result<()> {
    let args = Args::parse();
    let mut trees: Vec<NewickTree> = newick::from_filename(&args.infile)
        .with_context(|| format!("failed to parse {}", &args.infile))?;

    match args.command {
        Command::Annotate { species_tree } => {
            let mut out = String::new();
            let species_tree = newick::one_from_filename(&species_tree)
                .context(format!("while parsing {}", &species_tree))?;
            for t in trees.iter_mut() {
                actions::annotate_mrcas(t, &species_tree)?;
                actions::annotate_duplications(t, &species_tree, true);
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
                actions::compress(t).map(|_| {
                    out.push_str(&Newick::to_newick(t));
                    out.push('\n');
                })?;
            }
            let outfile = args.outfile.unwrap_or(args.infile);
            File::create(&outfile)?
                .write_all(out.as_bytes())
                .with_context(|| anyhow!("cannot write to `{}`", &outfile))
        }
        Command::Speciesize {
            database,
            cache_db,
            id,
            species: _species,
        } => {
            let mut out = String::new();
            let mut book = if cache_db {
                GeneBook::in_memory(&database, 0, &id)
            } else {
                GeneBook::inline(&database, 0, &id)
            }?;

            for t in trees.iter_mut() {
                actions::speciesize(t, &mut book).map(|_| {
                    out.push_str(&Newick::to_newick(t));
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
                actions::taxonize(t, &mapping).map(|_| {
                    out.push_str(&Newick::to_newick(t));
                    out.push('\n');
                })?
            }

            let outfile = args.outfile.unwrap_or(args.infile);
            File::create(&outfile)?
                .write_all(out.as_bytes())
                .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            Ok(())
        }
        Command::ToPhy {} => {
            let outfile = args.outfile.unwrap_or(
                std::path::Path::new(&args.infile)
                    .with_extension("phy")
                    .to_str()
                    .with_context(|| anyhow!("invalid filename found"))?
                    .to_owned(),
            );
            let mut out = File::create(&outfile)?;

            for t in trees.iter() {
                out.write_all(actions::to_phy(t)?.as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
                out.write_all("\n".as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            }

            Ok(())
        }
        Command::Leaves {} => {
            for t in trees.iter() {
                t.leaves()
                    .filter_map(|l| t.name(l))
                    .for_each(|n| println!("{}", n));
            }
            Ok(())
        }
        Command::Nodes {} => {
            for t in trees.iter() {
                t.nodes()
                    .filter_map(|n| t.name(n))
                    .for_each(|n| println!("{}", n));
            }
            Ok(())
        }
        Command::Normalize {} => {
            let outfile = args.outfile.unwrap_or(args.infile);
            let mut out = File::create(&outfile)?;

            for t in trees.iter_mut() {
                actions::normalize(t);
                out.write_all(Newick::to_newick(t).as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
                out.write_all("\n".as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            }
            Ok(())
        }
        Command::Prune { remove } => {
            let outfile = args.outfile.unwrap_or(args.infile);
            let mut out = File::create(&outfile)?;

            for mut t in trees {
                t.delete_nodes(
                    &t.nodes()
                        .filter(|&n| t.name(n).map(|s| remove.contains(s)).unwrap_or(false))
                        .collect::<Vec<_>>(),
                );
                t.prune();
                out.write_all(Newick::to_newick(&t).as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
                out.write_all("\n".as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            }
            Ok(())
        }
        Command::Binarize {} => {
            let outfile = args.outfile.unwrap_or(args.infile);
            let mut out = File::create(&outfile)?;

            for t in trees.iter_mut() {
                actions::binarize(t);
                out.write_all(Newick::to_newick(t).as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
                out.write_all("\n".as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            }
            Ok(())
        }
        Command::Rename {
            mapping_file,
            separator,
        } => {
            let mapping = std::io::BufReader::new(
                File::open(&mapping_file)
                    .with_context(|| anyhow!("while opening `{}`", &mapping_file))?,
            )
            .lines()
            .filter_map(|l| {
                l.ok().and_then(|l| {
                    let (src, tgt) = if let Some(sep) = separator.as_ref() {
                        let mut s = l.split(sep);
                        (s.next()?.to_owned(), s.next()?.to_owned())
                    } else {
                        let mut s = l.split_whitespace();
                        (s.next()?.to_owned(), s.next()?.to_owned())
                    };
                    Some((src, tgt))
                })
            })
            .collect::<std::collections::HashMap<String, String>>();
            let outfile = args.outfile.unwrap_or(args.infile);
            let mut out = File::create(&outfile)?;

            for t in trees.iter_mut() {
                actions::rename(t, &mapping);
                out.write_all(Newick::to_newick(t).as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
                out.write_all("\n".as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            }
            Ok(())
        }
        Command::Format => {
            let outfile = args.outfile.unwrap_or(args.infile);
            let mut out = File::create(&outfile)?;

            for t in trees {
                out.write_all(Newick::to_newick(&t).as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            }
            Ok(())
        }
        Command::Strip { to_strip } => {
            println!("Stripping {:?}", to_strip);

            let outfile = args.outfile.unwrap_or(args.infile);
            let mut out = File::create(&outfile)?;

            for t in trees.iter_mut() {
                actions::strip(t, &to_strip);
                out.write_all(Newick::to_newick(t).as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
                out.write_all("\n".as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            }
            Ok(())
        }
        Command::Tune { reference_tree } => {
            if trees.len() != 1 {
                bail!("tuning is only supported on single trees");
            }
            let reference: Vec<NewickTree> = newick::from_filename(&reference_tree)
                .with_context(|| format!("failed to parse {}", &reference_tree))?;
            if reference.len() != 1 {
                bail!("tuning is only supported on single trees");
            }
            let reference = &reference[0];

            Ok(())
        }
    }
}
