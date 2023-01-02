use anyhow::{anyhow, Context, Result};
use rusqlite::*;
use std::collections::HashMap;
use std::sync::Mutex;

#[allow(dead_code)]
pub enum GeneBook {
    InMemory(HashMap<String, Gene>),
    Cached(HashMap<String, Gene>),
    Inline(Mutex<Connection>),
}

#[derive(Clone, Default)]
pub struct Gene {
    pub gene: String,
    pub species: String,
}

impl GeneBook {
    fn get_rows<P: rusqlite::Params>(
        mut query: rusqlite::Statement,
        params: P,
    ) -> Result<HashMap<String, Gene>> {
        let genes = query
            .query_map(params, |r| {
                std::result::Result::Ok((
                    r.get::<_, String>(0)?,
                    r.get::<_, String>(1)?,
                    r.get::<_, String>(2)?,
                ))
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(genes
            .into_iter()
            .map(|g| {
                (
                    g.1.clone(),
                    Gene {
                        gene: g.0,
                        species: g.2,
                    },
                )
            })
            .collect())
    }

    pub fn in_memory(filename: &str) -> Result<Self> {
        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        let query = conn.prepare("SELECT gene, protein, species FROM genomes")?;
        let r = Self::get_rows(query, [])?;
        Ok(GeneBook::InMemory(r))
    }

    #[allow(dead_code)]
    pub fn cached<S: AsRef<str>>(filename: &str, ids: &[S]) -> Result<Self> {
        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        let query = conn.prepare(&format!(
            "SELECT gene, protein, species FROM genomes WHERE protein IN ({})",
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
    pub fn inline(filename: &str) -> Result<Self> {
        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        Ok(GeneBook::Inline(Mutex::new(conn)))
    }

    pub fn get(&self, g: &str) -> Result<Gene> {
        match self {
            GeneBook::InMemory(book) | GeneBook::Cached(book) => {
                book.get(g).cloned().ok_or_else(|| anyhow!("key not found"))
            }
            GeneBook::Inline(conn_mutex) => {
                let conn = conn_mutex.lock().expect("MUTEX POISONING");
                let mut query =
                    conn.prepare("SELECT gene, protein, species FROM genomes WHERE protein=?")?;
                query
                    .query_row(&[g], |r| {
                        let gene = r.get::<_, String>(0)?;
                        let species = r.get::<_, String>(5)?;

                        rusqlite::Result::Ok(Gene { gene, species })
                    })
                    .with_context(|| "while accessing DB")
            }
        }
    }
}
