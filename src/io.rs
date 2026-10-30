//! I/O for expression matrices and trait tables: TSV/CSV (text) and
//! Excel workbooks (xlsx/xls/xlsb/ods).

use calamine::{open_workbook_auto, Data, Reader};
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

/// Auto-detect the delimiter (tab vs comma) from the first line of the file.
/// Supports TSV and CSV input.
fn detect_delimiter(path: &str) -> Result<u8> {
    use std::io::BufRead;
    let f = std::fs::File::open(path)?;
    let mut first = String::new();
    std::io::BufReader::new(f).read_line(&mut first)?;
    let tabs = first.matches('\t').count();
    let commas = first.matches(',').count();
    Ok(if commas > tabs { b',' } else { b'\t' })
}

/// Detect whether a file is an Excel workbook (xlsx/xls/xlsb/ods) by magic
/// bytes: zip `PK` header (xlsx family) or OLE compound-file header (xls).
fn is_excel(path: &str) -> bool {
    use std::io::Read;
    let mut buf = [0u8; 8];
    let n = std::fs::File::open(path)
        .ok()
        .and_then(|mut f| f.read(&mut buf).ok())
        .unwrap_or(0);
    if n >= 2 && buf[0] == b'P' && buf[1] == b'K' {
        return true;
    }
    n >= 8 && buf[..8] == [0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1]
}

fn cell_to_string(c: &Data) -> String {
    match c {
        Data::String(s) => s.clone(),
        Data::Int(i) => i.to_string(),
        Data::Float(f) => f.to_string(),
        Data::Bool(b) => b.to_string(),
        _ => String::new(),
    }
}

fn cell_to_f64(c: &Data) -> Result<f64> {
    match c {
        Data::Int(i) => Ok(*i as f64),
        Data::Float(f) => Ok(*f),
        Data::String(s) => s.parse::<f64>().map_err(|_| {
            WgcnaError::InvalidInput(format!("non-numeric value '{s}'"))
        }),
        _ => Err(WgcnaError::InvalidInput("non-numeric cell".into())),
    }
}

/// Read the first sheet of an Excel workbook as a genes × samples matrix.
fn read_expression_excel(path: &str) -> Result<ExpData> {
    let mut wb = open_workbook_auto(path)
        .map_err(|e| WgcnaError::InvalidInput(format!("cannot open Excel file: {e}")))?;
    let range = wb
        .worksheet_range_at(0)
        .ok_or_else(|| WgcnaError::InvalidInput("workbook has no sheets".into()))?
        .map_err(|e| WgcnaError::InvalidInput(format!("cannot read sheet: {e}")))?;

    let mut rows = range.rows();
    let header = rows
        .next()
        .ok_or_else(|| WgcnaError::InvalidInput("empty Excel sheet".into()))?;
    let samples: Vec<String> = header.iter().skip(1).map(cell_to_string).collect();
    if samples.is_empty() {
        return Err(WgcnaError::InvalidInput(
            "expression matrix has no sample columns".into(),
        ));
    }

    let mut genes = Vec::new();
    let mut data_rows: Vec<Vec<f64>> = Vec::new();
    for row in rows {
        if row.is_empty() {
            continue;
        }
        let gene = cell_to_string(&row[0]);
        let mut vals = Vec::with_capacity(samples.len());
        for (i, cell) in row.iter().enumerate().skip(1) {
            if i - 1 >= samples.len() {
                break;
            }
            vals.push(cell_to_f64(cell)?);
        }
        if vals.len() != samples.len() {
            return Err(WgcnaError::InvalidInput(format!(
                "gene {gene} has {} values, expected {}",
                vals.len(),
                samples.len()
            )));
        }
        genes.push(gene);
        data_rows.push(vals);
    }

    if genes.is_empty() {
        return Err(WgcnaError::InvalidInput("no genes found".into()));
    }

    let g = genes.len();
    let n = samples.len();
    let mut data = Mat::<f64>::zeros(g, n);
    for (i, row) in data_rows.iter().enumerate() {
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

/// Read an expression matrix from text (TSV/CSV) or Excel (xlsx/xls/xlsb/ods).
pub fn read_expression(path: &str) -> Result<ExpData> {
    if is_excel(path) {
        read_expression_excel(path)
    } else {
        read_expression_text(path)
    }
}

/// Read a tab/comma-delimited expression matrix.
///
/// Expected layout: first row = header (`ID<DELIM>sample1<DELIM>sample2...`),
/// each following row = `gene_id<DELIM>v1<DELIM>v2...`.
fn read_expression_text(path: &str) -> Result<ExpData> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(detect_delimiter(path)?)
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

/// Read a trait table from text (TSV/CSV) or Excel (xlsx/xls/xlsb/ods).
pub fn read_trait(path: &str) -> Result<(Vec<String>, Vec<String>, Mat<f64>)> {
    let result = if is_excel(path) {
        read_trait_excel(path)
    } else {
        read_trait_text(path)
    }?;
    validate_identifiers(&result.0, "trait sample")?;
    validate_identifiers(&result.1, "trait column")?;
    Ok(result)
}

fn validate_identifiers(ids: &[String], kind: &str) -> Result<()> {
    let mut seen = std::collections::HashSet::new();
    if ids.is_empty() {
        return Err(WgcnaError::InvalidInput(format!("no {kind} IDs")));
    }
    for id in ids {
        if id.trim().is_empty() || !seen.insert(id) {
            return Err(WgcnaError::InvalidInput(format!(
                "empty or duplicate {kind} ID: '{id}'"
            )));
        }
    }
    Ok(())
}

/// Align trait rows to expression sample IDs, rejecting ambiguous or unequal
/// sample sets. Matrix-only numerical APIs still require pre-aligned rows.
pub fn align_trait_samples(
    data: &Mat<f64>,
    trait_samples: &[String],
    expression_samples: &[String],
) -> Result<Mat<f64>> {
    validate_identifiers(trait_samples, "trait sample")?;
    validate_identifiers(expression_samples, "expression sample")?;
    if data.nrows() != trait_samples.len() || trait_samples.len() != expression_samples.len() {
        return Err(WgcnaError::InvalidInput(
            "trait/expression sample dimensions differ".into(),
        ));
    }
    let index: std::collections::HashMap<&str, usize> = trait_samples
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let order: Vec<usize> = expression_samples
        .iter()
        .map(|s| {
            index.get(s.as_str()).copied().ok_or_else(|| {
                WgcnaError::InvalidInput(format!("missing trait sample '{s}'"))
            })
        })
        .collect::<Result<_>>()?;
    Ok(Mat::from_fn(order.len(), data.ncols(), |i, j| {
        data[(order[i], j)]
    }))
}

/// Read traits and align them by sample identity rather than file order.
pub fn read_trait_aligned(
    path: &str,
    expression_samples: &[String],
) -> Result<(Vec<String>, Mat<f64>)> {
    let (samples, traits, data) = read_trait(path)?;
    let aligned = align_trait_samples(&data, &samples, expression_samples)?;
    Ok((traits, aligned))
}

/// Read the first sheet of an Excel workbook as a trait table.
fn read_trait_excel(path: &str) -> Result<(Vec<String>, Vec<String>, Mat<f64>)> {
    let mut wb = open_workbook_auto(path)
        .map_err(|e| WgcnaError::InvalidInput(format!("cannot open Excel file: {e}")))?;
    let range = wb
        .worksheet_range_at(0)
        .ok_or_else(|| WgcnaError::InvalidInput("workbook has no sheets".into()))?
        .map_err(|e| WgcnaError::InvalidInput(format!("cannot read sheet: {e}")))?;

    let mut rows = range.rows();
    let header = rows
        .next()
        .ok_or_else(|| WgcnaError::InvalidInput("empty Excel sheet".into()))?;
    let traits: Vec<String> = header.iter().skip(1).map(cell_to_string).collect();

    let mut samples = Vec::new();
    let mut data_rows: Vec<Vec<f64>> = Vec::new();
    for row in rows {
        if row.is_empty() {
            continue;
        }
        if row.len() != traits.len() + 1 {
            return Err(WgcnaError::InvalidInput(
                "trait row width differs from header".into(),
            ));
        }
        samples.push(cell_to_string(&row[0]));
        let mut vals = Vec::with_capacity(traits.len());
        for (i, cell) in row.iter().enumerate().skip(1) {
            if i - 1 >= traits.len() {
                break;
            }
            vals.push(cell_to_f64(cell)?);
        }
        data_rows.push(vals);
    }

    let s = samples.len();
    let t = traits.len();
    let mut data = Mat::<f64>::zeros(s, t);
    for (i, row) in data_rows.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            data[(i, j)] = v;
        }
    }
    Ok((samples, traits, data))
}

/// Read a tab/comma-delimited trait table.
///
/// Expected: first column = sample id, remaining columns = traits
/// (already one-hot encoded for qualitative traits). Returns (samples, traits, data)
/// where `data` is samples x traits.
fn read_trait_text(path: &str) -> Result<(Vec<String>, Vec<String>, Mat<f64>)> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(detect_delimiter(path)?)
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;

    let headers = rdr.headers()?.clone();
    let traits: Vec<String> = headers.iter().skip(1).map(|s| s.to_string()).collect();
    let mut samples = Vec::new();
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        if rec.len() != traits.len() + 1 {
            return Err(WgcnaError::InvalidInput(
                "trait row width differs from header".into(),
            ));
        }
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

/// Convert the first sheet of any supported workbook to TSV for Dashboard
/// upload validation. Values stay textual; semantic validation happens before
/// these bytes are admitted to a job.
pub fn workbook_to_tsv(input: &str, output: &str) -> Result<()> {
    let mut wb = open_workbook_auto(input)
        .map_err(|e| WgcnaError::InvalidInput(format!("cannot open Excel file: {e}")))?;
    let range = wb.worksheet_range_at(0)
        .ok_or_else(|| WgcnaError::InvalidInput("workbook has no sheets".into()))?
        .map_err(|e| WgcnaError::InvalidInput(format!("cannot read sheet: {e}")))?;
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(output)?;
    for row in range.rows() {
        writer.write_record(row.iter().map(cell_to_string))?;
    }
    writer.flush()?;
    Ok(())
}
