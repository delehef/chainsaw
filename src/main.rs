use clap::{Parser, Subcommand};
use newick::*;
use std::fs::File;
use std::io::prelude::*;

use anyhow::{anyhow, Context, Result};

mod genebook;
use genebook::GeneBook;
mod actions;
mod utils;

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
        #[clap(value_parser, short, long)]
        species_tree: String,
    },

    /// annotate leaves in a tree with their species
    Speciesize {
        /// the database containing the id/species mapping
        #[clap(value_parser, short, long)]
        database: String,

        /// if set, cache the database in memory
        #[clap(value_parser, short)]
        cache_db: bool,

        /// the database column corresponding to leaf IDs in the tree
        #[clap(value_parser, default_value_t = String::from("protein"), long = "id")]
        id: String,

        /// the database column containing the species
        #[clap(value_parser, default_value_t = String::from("species"))]
        species: String,
    },

    /// convert ids from proteins to genes
    Convert {
        /// the database containing the proteins/genes mapping
        #[clap(value_parser, short = 'D', long = "database")]
        database: String,

        /// if set, cache the database in memory
        #[clap(value_parser, short = 'c')]
        cache_db: bool,

        /// column in the database corresponding to the current IDs
        #[clap(value_parser, long = "from")]
        from: String,

        /// column in the database corresponding to the desired IDs
        #[clap(value_parser, long = "to")]
        to: String,
    },

    ///
    Taxonize {
        #[clap(value_parser)]
        mapping: String,
    },

    /// compress root nodes with a single child
    Compress,

    /// convert a newick-formatted tree to a phyl-formatted tree
    ToPhy {},

    /// list the named leaves of the given tree
    Leaves {},

    /// list the names nodes of the given tree
    Nodes {},

    /// normalize a species tree according to ENSEMBL naming conventions
    Normalize {},
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
                actions::annotate_duplications(t, &species_tree, true);
                actions::annotate_mrcas(t, &species_tree)?;
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
                    out.push_str(&Newick::to_newick(t).replace('\n', ""));
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
            id: ids,
            species,
        } => {
            let mut out = String::new();
            let mut book = if cache_db {
                GeneBook::in_memory(&database, &ids, "gene", &species)
            } else {
                GeneBook::inline(&database, &ids, "gene", &species)
            }?;

            for t in trees.iter_mut() {
                actions::speciesize(t, &mut book).map(|_| {
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
                actions::taxonize(t, &mapping).map(|_| {
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
        Command::Convert {
            database,
            cache_db,
            from,
            to,
        } => {
            let mut out = String::new();
            let mut book = if cache_db {
                GeneBook::in_memory(&database, &from, &to, "species")
            } else {
                GeneBook::inline(&database, &from, &to, "species")
            }?;

            for t in trees.iter_mut() {
                actions::convert(t, &mut book).map(|_| {
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
                    .filter_map(|l| t[l].data.name.as_ref())
                    .for_each(|n| println!("{}", n));
            }
            Ok(())
        }
        Command::Nodes {} => {
            for t in trees.iter() {
                t.nodes()
                    .filter_map(|n| n.data.name.as_ref())
                    .for_each(|n| println!("{}", n));
            }
            Ok(())
        }
        Command::Normalize {} => {
            let outfile = args.outfile.unwrap_or(args.infile);
            let mut out = File::create(&outfile)?;

            for t in trees.iter_mut() {
                actions::normalize(t);
                out.write_all(Newick::to_newick(t).replace('\n', "").as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
                out.write_all("\n".as_bytes())
                    .with_context(|| anyhow!("cannot write to `{}`", &outfile))?;
            }
            Ok(())
        }
    }
}
