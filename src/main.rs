use clap::*;
use std::io::{BufRead, BufReader};
use newick::*;
use std::collections::HashSet;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::path::Path;

use anyhow::{Context, Result};

fn jaccard<T>(a: &HashSet<T>, b: &HashSet<T>) -> f32
where
    T: Eq + Hash,
{
    a.intersection(b).collect::<Vec<_>>().len() as f32 / a.union(b).collect::<Vec<_>>().len() as f32
}

fn annotate_duplications(t: &mut Tree) {
    t.inners().collect::<Vec<_>>().iter().for_each(|n| {
        let children = t[*n].children.as_ref().unwrap();
        let species: Vec<HashSet<_>> = children
            .iter()
            .map(|&c| {
                t.leaves_for(c)
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
        let d = species.iter().skip(1).any(|x| !x.is_disjoint(&species[0]));
        assert!(species.len() == 2);
        let dcs = jaccard(&species[0], &species[1]);
        if d {
            t[*n].data.insert("D".to_string(), "Y".to_owned());
            t[*n].data.insert("DCS".to_string(), dcs.to_string());
        } else {
            t[*n].data.insert("D".to_string(), "N".to_owned());
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
        .subcommand(SubCommand::with_name("duplications"))
        .get_matches();
    let filename = value_t!(args, "FILE", String).unwrap();
    println!("Processing {}", filename);

    match args.subcommand() {
        ("duplications", _) => {
            let mut out = String::new();

            for l in BufReader::new(File::open(&filename)?).lines() {
                let mut t = Tree::from_string(&l?)?;
                annotate_duplications(&mut t);
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
