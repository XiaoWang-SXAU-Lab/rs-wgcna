//! I/O for tab-delimited expression matrices and trait tables.

use faer::Mat;

use crate::error::{Result, WgcnaError};

/// Expression matrix: genes (rows) x samples (cols).
#[derive(Debug, Clone)]
pub struct ExpData {
    pub data: Mat<f64>,
    pub genes: Vec<String>,
    pub samples: Vec<String>,
}

impl ExpData {
    pub fn n_genes(&self) -> usize {
        self.data.nrows()
    }
    pub fn n_samples(&self) -> usize {
        self.data.ncols()
    }
}

/// Read a tab-delimited expression matrix.
///
/// Expected layout: first row = header (`ID<TAB>sample1<TAB>sample2...`),
/// each following row = `gene_id<TAB>v1<TAB>v2...`.
pub fn read_expression(path: &str) -> Result<ExpData> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;

    let headers = rdr.headers()?.clone();
    let samples: Vec<String> = headers.iter().skip(1).map(|s| s.to_string()).collect();
    if samples.is_empty() {
        return Err(WgcnaError::InvalidInput(
            "expression matrix has no sample columns".into(),
        ));
    }

    let mut genes = Vec::new();
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        let gene = rec.get(0).unwrap_or("").to_string();
        let mut vals = Vec::with_capacity(samples.len());
        for (i, field) in rec.iter().enumerate().skip(1) {
            if i - 1 >= samples.len() {
                break;
            }
            let v = field.parse::<f64>().map_err(|_| {
                WgcnaError::InvalidInput(format!("non-numeric value '{field}' for gene {gene}"))
            })?;
            vals.push(v);
        }
        if vals.len() != samples.len() {
            return Err(WgcnaError::InvalidInput(format!(
                "gene {gene} has {} values, expected {}",
                vals.len(),
                samples.len()
            )));
        }
        genes.push(gene);
        rows.push(vals);
    }

    if genes.is_empty() {
        return Err(WgcnaError::InvalidInput("no genes found".into()));
    }

    let g = genes.len();
    let n = samples.len();
    let mut data = Mat::<f64>::zeros(g, n);
    for (i, row) in rows.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            data[(i, j)] = v;
        }
    }

    Ok(ExpData {
        data,
        genes,
        samples,
    })
}

/// Read a tab-delimited trait table.
///
/// Expected: first column = sample id, remaining columns = traits
/// (already one-hot encoded for qualitative traits). Returns (samples, traits, data)
/// where `data` is samples x traits.
pub fn read_trait(path: &str) -> Result<(Vec<String>, Vec<String>, Mat<f64>)> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;

    let headers = rdr.headers()?.clone();
    let traits: Vec<String> = headers.iter().skip(1).map(|s| s.to_string()).collect();
    let mut samples = Vec::new();
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        samples.push(rec.get(0).unwrap_or("").to_string());
        let mut vals = Vec::with_capacity(traits.len());
        for (i, field) in rec.iter().enumerate().skip(1) {
            if i - 1 >= traits.len() {
                break;
            }
            let v = field.parse::<f64>().map_err(|_| {
                WgcnaError::InvalidInput(format!("non-numeric trait value '{field}'"))
            })?;
            vals.push(v);
        }
        rows.push(vals);
    }
    let s = samples.len();
    let t = traits.len();
    let mut data = Mat::<f64>::zeros(s, t);
    for (i, row) in rows.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            data[(i, j)] = v;
        }
    }
    Ok((samples, traits, data))
}

/// Write a genes x samples expression matrix as a tab-delimited TSV
/// (header = `ID<TAB>sample1<TAB>...`).
pub fn write_expression(
    path: &str,
    data: &Mat<f64>,
    genes: &[String],
    samples: &[String],
) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    write!(f, "ID")?;
    for s in samples {
        write!(f, "\t{s}")?;
    }
    writeln!(f)?;
    for i in 0..data.nrows() {
        write!(f, "{}", genes[i])?;
        for j in 0..data.ncols() {
            write!(f, "\t{:.6}", data[(i, j)])?;
        }
        writeln!(f)?;
    }
    Ok(())
}
