//! Network / TOM type enums, matching WGCNA's `.networkTypes`, `.TOMTypes`, `.TOMDenoms`.

use crate::error::{Result, WgcnaError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NetworkType {
    Unsigned,
    Signed,
    SignedHybrid,
}

impl NetworkType {
    pub fn from_str(s: &str) -> Result<Self> {
        match s {
            "unsigned" => Ok(NetworkType::Unsigned),
            "signed" => Ok(NetworkType::Signed),
            "signed hybrid" => Ok(NetworkType::SignedHybrid),
            _ => Err(WgcnaError::InvalidInput(format!(
                "unknown network type: {s}"
            ))),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TomType {
    None,
    Unsigned,
    Signed,
    SignedNowick,
    Unsigned2,
    Signed2,
    SignedNowick2,
}

impl TomType {
    pub fn from_str(s: &str) -> Result<Self> {
        match s {
            "none" => Ok(TomType::None),
            "unsigned" => Ok(TomType::Unsigned),
            "signed" => Ok(TomType::Signed),
            "signed Nowick" => Ok(TomType::SignedNowick),
            "unsigned 2" => Ok(TomType::Unsigned2),
            "signed 2" => Ok(TomType::Signed2),
            "signed Nowick 2" => Ok(TomType::SignedNowick2),
            _ => Err(WgcnaError::InvalidInput(format!("unknown TOM type: {s}"))),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TomDenom {
    Min,
    Mean,
}

impl TomDenom {
    pub fn from_str(s: &str) -> Result<Self> {
        match s {
            "min" => Ok(TomDenom::Min),
            "mean" => Ok(TomDenom::Mean),
            _ => Err(WgcnaError::InvalidInput(format!("unknown TOM denom: {s}"))),
        }
    }
}
