# Decision log

## 2026-10-06 — before modeling
- Chose erlotinib from the user's single-drug cancer-response project scope. No predictive results from other drugs were screened.
- Legacy GDSC downloads were unavailable (410/403). Used official archived versioned files on Sanger's object storage instead.
- Stopped a full expression archive download after discovering its internal structure. Retrieved just the FPKM matrix by HTTP byte ranges and validated its ZIP CRC. The partial archive is not an analysis input or shared artifact.
- Source audit found several assay concentration windows. Chose the largest exact window before modeling: 705 responses at 0.002001–2 micromolar. Four had no matching expression, leaving 701.
- Grouped all known donor/sample/parent relationships before splitting. All eligible models were in distinct known groups; no group crosses the train/test or CV boundaries.
- Fixed the split at seed 42 and feature count at 1,000. The plan was saved before fitting any model. The study is an internal benchmark, with no external validation claimed.
- Fixed gene IDs retained; ambiguous gene symbols do not cause collapsing of distinct genes.

## Verification detail
On Windows, the split checksum is computed over the actual file bytes so CRLF line endings do not cause a mismatch with a checksum computed from a normalized in-memory string. The cohort and selection records must identify the same physical split file.
