use anyhow::*;
use clap::ValueEnum;
use newick::{Newick, NewickTree};
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
};
use syntesuite::genebook::GeneBook;

use crate::utils::{capitalize, effective_losses, jaccard};

#[derive(Debug, Clone, ValueEnum)]
pub(crate) enum Strippable {
    Ancestors,
    Attributes,
    Length,
}

pub fn annotate_duplications(t: &mut NewickTree, species_tree: &NewickTree, filter_species: bool) {
    let restricted_species = if filter_species {
        Some(
            t.leaves()
                .filter_map(|l| t.attrs(l).get("S").map(|s| s.to_owned()))
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
                    .map(|&n| t.attrs(n).get("S").map(|s| s.to_owned()).unwrap())
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
                t.attrs_mut(*n).insert("D".to_string(), "Y".to_owned());
                t.attrs_mut(*n).insert("DCS".to_string(), dcs.to_string());
                t.attrs_mut(*n)
                    .insert("ELC".to_string(), elc_all.to_string());
                t.attrs_mut(*n)
                    .insert("ELLC".to_string(), elc_large.to_string());
            } else {
                t.attrs_mut(*n).insert("D".to_string(), "N".to_owned());
            }
        }
    });
}

pub fn annotate_mrcas(t: &mut NewickTree, species_tree: &NewickTree) -> Result<()> {
    for n in t.inners().collect::<Vec<_>>().into_iter() {
        let species: HashSet<usize> = t
            .leaves_of(n)
            .iter()
            .flat_map(|&c| {
                t.leaves_of(c)
                    .into_iter()
                    .map(|n| t.attrs(n).get("S").map(|s| s.to_owned()).unwrap())
                    .map(|s| {
                        species_tree
                            .find_leaf(|l| l.name.as_ref().unwrap().as_str() == s.as_str())
                            .ok_or_else(|| anyhow!(format!("{} not found in species tree", s)))
                    })
            })
            .collect::<Result<HashSet<_>>>()?;
        let mrca = species_tree.mrca(&species).unwrap();
        t.attrs_mut(n)
            .insert("S".to_owned(), species_tree.name(mrca).unwrap().to_owned());
    }
    Ok(())
}

pub fn speciesize(t: &mut NewickTree, book: &mut GeneBook) -> Result<()> {
    t.map_leaves(&mut |n| {
        if n.data.as_ref().unwrap().name.is_some() {
            let name = n.data.as_ref().unwrap().name.as_ref().unwrap().to_owned();
            n.data.as_mut().unwrap().attrs.insert(
                "S".to_owned(),
                book.get(&name)
                    .unwrap_or_else(|_| panic!("Cannot find {:?}", name))
                    .species,
            );
        }
    });

    Ok(())
}

pub fn taxonize(t: &mut NewickTree, map_file: &str) -> Result<()> {
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
        let taxon_id = l
            .data
            .as_ref()
            .unwrap()
            .attrs
            .get("T")
            .and_then(|s| s.parse::<usize>().ok());
        if let Some(taxon_id) = taxon_id {
            if let Some(species) = map.get(&taxon_id) {
                l.data
                    .as_mut()
                    .unwrap()
                    .attrs
                    .insert("S".to_string(), species.to_owned());
            } else {
                eprintln!("`{}` has no match in the reference file", taxon_id);
            }
        } else {
            eprintln!(
                "Node `{:?}` has no taxon specified",
                l.data.as_ref().unwrap().name
            );
        }
    }
    Ok(())
}

pub fn compress(t: &mut NewickTree) -> Result<()> {
    while t[t.root()].children().len() == 1 {
        println!("Compressing");
        t.set_root(t[t.root()].children()[0]);
    }
    Ok(())
}

pub fn to_phy(t: &NewickTree) -> Result<String> {
    fn rec_to_phy(ax: &mut String, t: &NewickTree, n: usize, d: usize, id: usize) {
        ax.push_str(&format!(
            "{}{}\n",
            "\t".repeat(d),
            t.name(n).cloned().unwrap_or_else(|| String::from("UKNWN"))
        ));
        for c in t[n].children() {
            rec_to_phy(ax, t, *c, d + 1, id);
        }
    }

    let mut r = String::new();
    let root = t.root();
    rec_to_phy(&mut r, t, root, 0, 0);
    Ok(r)
}

pub fn normalize(t: &mut NewickTree) {
    let mut known_names = HashSet::new();

    for (i, n) in t.nodes_mut().enumerate() {
        if let Some(name) = n.data.as_ref().unwrap().name.clone() {
            let mut new_name = capitalize(
                &name
                    .split('_')
                    .map(|s| s.replace(|c: char| !(c.is_alphanumeric() || c == '.'), ""))
                    .collect::<Vec<String>>()
                    .join(".")
                    .to_lowercase(),
            );

            if known_names.contains(&new_name) && !n.is_leaf() {
                eprint!("/!\\ {} already known; replacing with", new_name);
                new_name = format!(" {}-{}", new_name, i);
                eprintln!("{}", new_name)
            } else {
                known_names.insert(new_name.clone());
            }
            n.data.as_mut().unwrap().name = Some(new_name);
        } else if !n.is_leaf() {
            eprintln!("/!\\ creating an ancestral name");
            n.data.as_mut().unwrap().name = Some(format!("ancestral-{}", i))
        }
    }
}

pub fn binarize(t: &mut NewickTree) {
    loop {
        let todo = t.nodes().find(|n| t.children(*n).len() > 2);

        if let Some(n) = todo {
            let n2 = t.add_node(
                None,
                Some(newick::Data {
                    name: t.name(n).map(|n| format!("{}+", n)),
                    attrs: t.attrs(n).clone(),
                }),
            );
            for c in t.children(n).to_owned().iter().skip(1) {
                t.move_node(*c, n2);
            }
            t.move_node(n2, n);
        } else {
            break;
        }
    }
}

pub fn rename(t: &mut NewickTree, mapping: &HashMap<String, String>) {
    for l in t.nodes_mut() {
        if l.data.as_ref().and_then(|d| d.name.as_ref()).is_some() {
            if let Some(tgt) = mapping.get(l.data.as_ref().unwrap().name.as_ref().unwrap()) {
                l.data.as_mut().unwrap().name = Some(tgt.clone());
            }
        }
    }
}

pub(crate) fn strip(t: &mut NewickTree, to_strip: &[Strippable]) {
    for n in t.nodes_mut() {
        for s in to_strip {
            match s {
                Strippable::Ancestors => {
                    if !n.is_leaf() {
                        if let Some(d) = n.data.as_mut() {
                            d.name = None
                        }
                    }
                }
                Strippable::Attributes => {
                    if let Some(d) = n.data.as_mut() {
                        d.attrs.clear()
                    }
                }
                Strippable::Length => n.branch_length = None,
            }
        }
    }
}
