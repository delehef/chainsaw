use std::collections::HashSet;
use std::hash::Hash;

use newick::NewickTree;

pub fn jaccard<T>(a: &HashSet<T>, b: &HashSet<T>) -> f32
where
    T: Eq + Hash,
{
    a.intersection(b).count() as f32 / a.union(b).count() as f32
}

pub fn effective_losses(
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
