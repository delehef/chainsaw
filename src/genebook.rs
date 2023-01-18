use anyhow::{anyhow, Context, Result};
use rusqlite::*;
use std::collections::HashMap;
use std::sync::Mutex;

#[allow(dead_code)]
pub enum GeneBook {
    InMemory(HashMap<String, Gene>),
    Cached(HashMap<String, Gene>),
    Inline {
        mutex: Mutex<Connection>,
        ids: String,
        species: String,
    },
}

#[derive(Clone, Debug)]
pub struct Gene {
    pub id: String,
    pub species: String,
}

impl GeneBook {
    fn get_rows<P: rusqlite::Params>(
        mut query: rusqlite::Statement,
        params: P,
    ) -> Result<HashMap<String, Gene>> {
        let genes = query
            .query_map(params, |r| {
                std::result::Result::Ok((r.get::<_, String>(0)?, r.get::<_, String>(1)?))
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(genes
            .into_iter()
            .map(|g| {
                (
                    g.0.clone(),
                    Gene {
                        id: g.0,
                        species: g.1,
                    },
                )
            })
            .collect())
    }

    pub fn in_memory(filename: &str, id: &str, species: &str) -> Result<Self> {
        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        let query = conn.prepare(&format!("SELECT {}, {} FROM genomes", id, species))?;
        let r = Self::get_rows(query, [])?;
        Ok(GeneBook::InMemory(r))
    }

    #[allow(dead_code)]
    pub fn cached<S: AsRef<str>>(filename: &str, ids: &[S]) -> Result<Self> {
        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        let query = conn.prepare(&format!(
            "SELECT id, species FROM genomes WHERE id IN ({})", // XXX
            std::iter::repeat("?")
                .take(ids.len())
                .collect::<Vec<_>>()
                .join(", ")
        ))?;
        let r = Self::get_rows(
            query,
            rusqlite::params_from_iter(ids.iter().map(|s| s.as_ref())),
        )?;
        Ok(GeneBook::Cached(r))
    }

    #[allow(dead_code)]
    pub fn inline(filename: &str, from: &str, species: &str) -> Result<Self> {
        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        Ok(GeneBook::Inline {
            mutex: Mutex::new(conn),
            ids: from.to_owned(),
            species: species.to_owned(),
        })
    }

    pub fn get(&self, g: &str) -> Result<Gene> {
        match self {
            GeneBook::InMemory(book) | GeneBook::Cached(book) => {
                book.get(g).cloned().ok_or_else(|| anyhow!("key not found"))
            }
            GeneBook::Inline {
                mutex,
                ids,
                species,
            } => {
                let conn = mutex.lock().expect("MUTEX POISONING");
                let mut query = conn.prepare(&format!(
                    "SELECT {}, {} FROM genomes WHERE {}=?",
                    ids, species, ids
                ))?;
                query
                    .query_row(&[g], |r| {
                        let gene = r.get::<_, String>(0)?;
                        let species = r.get::<_, String>(1)?;

                        rusqlite::Result::Ok(Gene { id: gene, species })
                    })
                    .with_context(|| "while accessing DB")
            }
        }
    }
}
