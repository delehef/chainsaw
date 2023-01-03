use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::*;
use newick::NewickTree;

use crate::{
    genebook::GeneBook,
    utils::{effective_losses, jaccard},
};

pub fn annotate_duplications(t: &mut NewickTree, species_tree: &NewickTree, filter_species: bool) {
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

pub fn annotate_mrcas(t: &mut NewickTree, species_tree: &NewickTree) -> Result<()> {
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

pub fn convert(t: &mut NewickTree, book: &mut GeneBook) -> Result<()> {
    t.map_leaves(&mut |n| {
        if n.data.name.is_some() {
            n.data.name = Some(
                book.get(n.data.name.as_ref().unwrap())
                    .unwrap_or_else(|_| {
                        panic!(
                            "can not find {:?} in database",
                            n.data.name.as_ref().unwrap()
                        )
                    })
                    .gene,
            )
        }
    });

    Ok(())
}

pub fn speciesize(t: &mut NewickTree, book: &mut GeneBook) -> Result<()> {
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

pub fn compress(t: &mut NewickTree) -> Result<()> {
    while t[t.root()].children().len() == 1 {
        println!("Compressing");
        t.set_root(t[t.root()].children()[0]);
    }
    Ok(())
}

pub fn to_phy(t: &NewickTree) -> Result<String> {
    fn rec_to_phy(ax: &mut String, t: &NewickTree, n: usize, d: usize, id: usize) {
        let node = &t[n];
        ax.push_str(&format!(
            "{}{}\n",
            "\t".repeat(d),
            node.data
                .name
                .as_ref()
                .cloned()
                .unwrap_or_else(|| String::from("UKNWN"))
        ));
        for c in node.children() {
            rec_to_phy(ax, t, *c, d + 1, id);
        }
    }

    let mut r = String::new();
    let root = t.root();
    rec_to_phy(&mut r, t, root, 0, 0);
    Ok(r)
}
