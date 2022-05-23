use anyhow::{anyhow, Context, Result};
use rusqlite::*;
use std::collections::HashMap;
use std::sync::Mutex;

#[derive(Clone, Default)]
pub struct Gene {
    pub name: String,
}

pub enum GeneBook {
    Cached(HashMap<String, Gene>),
    Inline(Mutex<Connection>),
}

impl GeneBook {
    pub fn cached(filename: &str) -> Result<Self> {
        println!("Parsing the database...");

        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        let mut query = conn.prepare("SELECT gene, protein FROM genomes")?;
        let genes = query
            .query_map([], |r| {
                std::result::Result::Ok((r.get::<_, String>(0)?, r.get::<_, String>(1)?))
            })?
            .collect::<Result<Vec<_>, _>>()?;

        let r = genes
            .into_iter()
            .map(|g| (g.1.clone(), Gene { name: g.0 }))
            .collect();
        println!("Done.");
        Ok(GeneBook::Cached(r))
    }

    pub fn inline(filename: &str) -> Result<Self> {
        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        Ok(GeneBook::Inline(Mutex::new(conn)))
    }

    pub fn get(&self, g: &str) -> Result<Gene> {
        match self {
            GeneBook::Cached(book) => book
                .get(g)
                .map(|g| g.clone())
                .ok_or(anyhow::anyhow!("key not found")),
            GeneBook::Inline(conn_mutex) => {
                let conn = conn_mutex.lock().expect("MUTEX POISONING");
                let mut query = conn.prepare("SELECT gene FROM genomes WHERE protein=?")?;
                query
                    .query_row(&[g], |r| {
                        rusqlite::Result::Ok(Gene {
                            name: r.get::<_, String>(0)?,
                        })
                    })
                    .with_context(|| "while accessing DB")
            }
        }
    }
}
