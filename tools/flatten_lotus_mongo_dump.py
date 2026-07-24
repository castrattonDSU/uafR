#!/usr/bin/env python3
"""Flatten the official LOTUS MongoDB BSON dump into CSV rows.

This tool intentionally uses only the Python standard library. It streams
`lotusUniqueNaturalProduct.bson` from the official LOTUS MongoDB ZIP download
and writes a flat table with one row per LOTUS compound/taxon/reference record.

The output is designed for:

    Rscript tools/build_lotus_index.R --input lotus_mongo_flat.csv \
      --out-file lotus_compact_index.csv --overwrite

For production plant panels, also write a manifest-backed lookup directory:

    python3 tools/flatten_lotus_mongo_dump.py --input LOTUSlatest.zip \
      --out-file lotus_mongo_flat.csv \
      --compact-index-file lotus_compact_index.csv \
      --lookup-dir lotus_lookup_index --overwrite

If the compact CSV already exists, build only the lookup directory:

    python3 tools/flatten_lotus_mongo_dump.py \
      --from-compact-index lotus_compact_index.csv \
      --lookup-dir lotus_lookup_index --overwrite

It does not fabricate taxa. Rows come from `taxonomyReferenceObjects` when that
field is present, with a conservative fallback to species-like terms in
`allTaxa`.
"""

import argparse
import csv
import io
import json
import os
import re
import shutil
import struct
import sys
import time
import zipfile
from collections import OrderedDict
from datetime import datetime, timezone


DEFAULT_MEMBER = "NPOC2021/NPOC2021/lotusUniqueNaturalProduct.bson"

TOP_LEVEL_FIELDS = {
    "allTaxa",
    "allWikidataIds",
    "chemicalTaxonomyClassyfireClass",
    "chemicalTaxonomyClassyfireDirectParent",
    "chemicalTaxonomyClassyfireKingdom",
    "chemicalTaxonomyClassyfireSuperclass",
    "chemicalTaxonomyNPclassifierClass",
    "chemicalTaxonomyNPclassifierPathway",
    "chemicalTaxonomyNPclassifierSuperclass",
    "inchi",
    "inchi2D",
    "inchikey",
    "inchikey2D",
    "iupac_name",
    "lotus_id",
    "molecular_formula",
    "smiles",
    "smiles2D",
    "synonyms",
    "taxonomyReferenceObjects",
    "traditional_name",
    "wikidata_id",
    "xrefs",
}

OUTPUT_FIELDS = [
    "lotus_id",
    "wikidata_id",
    "compound_name",
    "traditional_name",
    "iupac_name",
    "smiles",
    "inchikey",
    "molecular_formula",
    "chemicalTaxonomyNPclassifierPathway",
    "chemicalTaxonomyNPclassifierSuperclass",
    "chemicalTaxonomyNPclassifierClass",
    "chemicalTaxonomyClassyfireKingdom",
    "chemicalTaxonomyClassyfireSuperclass",
    "chemicalTaxonomyClassyfireClass",
    "allTaxa",
    "allWikidataIds",
    "reference_id",
    "doi",
    "taxonomy_provider",
    "organism_value",
    "species",
    "genus",
    "family",
    "kingdom",
    "phylum",
    "class",
    "order",
    "cleaned_organism_id",
    "taxon_wikidata_id",
    "reference_wikidata_id",
    "source_file",
]

COMPACT_FIELDS = [
    "species",
    "species_clean",
    "genus",
    "family",
    "compound_name",
    "compound_name_clean",
    "lotus_id",
    "wikidata_id",
    "cid",
    "smiles",
    "inchikey",
    "molecular_formula",
    "source_database",
    "source_record_id",
    "evidence_text",
    "evidence_url",
    "reference_id",
    "pmid",
    "doi",
    "plant_part",
    "tissue",
    "method",
    "chemical_class_pathway",
    "chemical_class_superclass",
    "chemical_class_class",
    "all_chem_classifications",
    "all_taxa",
    "retrieved_at",
    "source_file",
]

LOOKUP_FIELDS = ["index_key_type", "index_key"] + COMPACT_FIELDS


class BsonReader:
    def __init__(self, data):
        self.data = data
        self.i = 0

    def read(self, n):
        chunk = self.data[self.i : self.i + n]
        if len(chunk) != n:
            raise ValueError("short BSON read")
        self.i += n
        return chunk

    def cstring(self):
        end = self.data.index(b"\x00", self.i)
        value = self.data[self.i : end].decode("utf-8", "replace")
        self.i = end + 1
        return value

    def string(self):
        length = struct.unpack("<i", self.read(4))[0]
        raw = self.read(length)
        return raw[:-1].decode("utf-8", "replace")

    def document(self, wanted=None):
        length = struct.unpack("<i", self.read(4))[0]
        end = self.i + length - 5
        out = {}
        while self.i < end:
            element_type = self.read(1)[0]
            key = self.cstring()
            if wanted is None or key in wanted:
                out[key] = self.value(element_type)
            else:
                self.skip_value(element_type)
        self.read(1)
        return out

    def value(self, element_type):
        if element_type == 0x01:
            return struct.unpack("<d", self.read(8))[0]
        if element_type == 0x02:
            return self.string()
        if element_type == 0x03:
            return self.document()
        if element_type == 0x04:
            doc = self.document()
            return [doc[key] for key in sorted(doc.keys(), key=array_key)]
        if element_type == 0x05:
            length = struct.unpack("<i", self.read(4))[0]
            subtype = self.read(1)[0]
            self.read(length)
            return {"binary_subtype": subtype, "length": length}
        if element_type == 0x06:
            return None
        if element_type == 0x07:
            return self.read(12).hex()
        if element_type == 0x08:
            return bool(self.read(1)[0])
        if element_type == 0x09:
            return {"date_ms": struct.unpack("<q", self.read(8))[0]}
        if element_type == 0x0A:
            return None
        if element_type == 0x0B:
            return {"regex": self.cstring(), "options": self.cstring()}
        if element_type == 0x0C:
            return {"namespace": self.string(), "id": self.read(12).hex()}
        if element_type in (0x0D, 0x0E):
            return self.string()
        if element_type == 0x0F:
            length = struct.unpack("<i", self.read(4))[0]
            self.read(length - 4)
            return {"code_w_scope_length": length}
        if element_type == 0x10:
            return struct.unpack("<i", self.read(4))[0]
        if element_type == 0x11:
            inc, ts = struct.unpack("<II", self.read(8))
            return {"timestamp": ts, "increment": inc}
        if element_type == 0x12:
            return struct.unpack("<q", self.read(8))[0]
        if element_type == 0x13:
            return {"decimal128": self.read(16).hex()}
        if element_type in (0xFF, 0x7F):
            return None
        raise ValueError(f"unsupported BSON type: {element_type:#x}")

    def skip_value(self, element_type):
        if element_type == 0x01:
            self.i += 8
        elif element_type == 0x02:
            length = struct.unpack("<i", self.read(4))[0]
            self.i += length
        elif element_type in (0x03, 0x04):
            length = struct.unpack("<i", self.read(4))[0]
            self.i += length - 4
        elif element_type == 0x05:
            length = struct.unpack("<i", self.read(4))[0]
            self.i += 1 + length
        elif element_type in (0x06, 0x0A, 0xFF, 0x7F):
            return
        elif element_type == 0x07:
            self.i += 12
        elif element_type == 0x08:
            self.i += 1
        elif element_type in (0x09, 0x11, 0x12):
            self.i += 8
        elif element_type == 0x0B:
            self.cstring()
            self.cstring()
        elif element_type == 0x0C:
            length = struct.unpack("<i", self.read(4))[0]
            self.i += length + 12
        elif element_type in (0x0D, 0x0E):
            length = struct.unpack("<i", self.read(4))[0]
            self.i += length
        elif element_type == 0x0F:
            length = struct.unpack("<i", self.read(4))[0]
            self.i += length - 4
        elif element_type == 0x10:
            self.i += 4
        elif element_type == 0x13:
            self.i += 16
        else:
            raise ValueError(f"unsupported BSON type while skipping: {element_type:#x}")


def array_key(value):
    return int(value) if isinstance(value, str) and value.isdigit() else value


def read_documents(stream):
    while True:
        length_raw = stream.read(4)
        if not length_raw:
            break
        if len(length_raw) != 4:
            raise ValueError("truncated BSON document length")
        length = struct.unpack("<i", length_raw)[0]
        if length < 5:
            raise ValueError(f"invalid BSON document length: {length}")
        body = stream.read(length - 4)
        if len(body) != length - 4:
            raise ValueError("truncated BSON document")
        yield BsonReader(length_raw + body).document(wanted=TOP_LEVEL_FIELDS)


def open_bson_source(path, member):
    if zipfile.is_zipfile(path):
        archive = zipfile.ZipFile(path)
        chosen = member
        if chosen not in archive.namelist():
            candidates = [
                name
                for name in archive.namelist()
                if name.endswith("lotusUniqueNaturalProduct.bson")
            ]
            if not candidates:
                raise ValueError("No lotusUniqueNaturalProduct.bson member found in ZIP")
            chosen = candidates[0]
        return archive, archive.open(chosen), chosen
    return None, open(path, "rb"), path


def clean_text(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return "; ".join(clean_text(item) for item in value if clean_text(item))
    if isinstance(value, dict):
        return "; ".join(f"{key}={clean_text(val)}" for key, val in value.items())
    return re.sub(r"\s+", " ", str(value)).strip()


def clean_key(value):
    text = clean_text(value).lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text


def clean_name(value):
    return re.sub(r"\s+", " ", clean_text(value).lower()).strip()


def first_text(*values):
    for value in values:
        text = clean_text(value)
        if text:
            return text
    return ""


def wikidata_qid(value):
    text = first_text(value)
    match = re.search(r"\bQ\d+\b", text)
    return match.group(0) if match else text


def decode_reference_id(reference_id):
    text = first_text(reference_id)
    if not text:
        return "", ""
    decoded = text.replace("$x$x$", ".")
    doi = decoded if decoded.lower().startswith("10.") else ""
    return decoded, doi


def species_like_terms(terms):
    out = []
    for term in terms or []:
        text = first_text(term)
        if re.match(r"^[A-Z][a-z-]{2,}\s+(?:x\s+)?[a-z][a-z-]{2,}", text):
            out.append(text)
    return out


def taxonomy_rows(taxonomy_reference_objects, all_taxa):
    rows = []
    if isinstance(taxonomy_reference_objects, dict):
        for reference_id, provider_map in taxonomy_reference_objects.items():
            if not isinstance(provider_map, dict):
                continue
            for provider, records in provider_map.items():
                if isinstance(records, dict):
                    records = [records]
                if not isinstance(records, list):
                    continue
                for record in records:
                    if not isinstance(record, dict):
                        continue
                    species = first_text(record.get("species"))
                    organism = first_text(record.get("organism_value"))
                    genus = first_text(record.get("genus"))
                    family = first_text(record.get("family"))
                    if not species and organism:
                        species = first_species_like([organism])
                    if not species:
                        continue
                    decoded_ref, doi = decode_reference_id(reference_id)
                    rows.append(
                        {
                            "reference_id": decoded_ref or reference_id,
                            "doi": doi,
                            "taxonomy_provider": provider,
                            "organism_value": organism,
                            "species": species,
                            "genus": genus,
                            "family": family,
                            "kingdom": first_text(
                                record.get("kingdom"),
                                record.get("superkingdom"),
                                record.get("domain"),
                            ),
                            "phylum": first_text(record.get("phylum")),
                            "class": first_text(record.get("classx"), record.get("class")),
                            "order": first_text(record.get("order")),
                            "cleaned_organism_id": first_text(
                                record.get("cleaned_organism_id")
                            ),
                            "taxon_wikidata_id": wikidata_qid(record.get("wikidata_id")),
                            "reference_wikidata_id": wikidata_qid(
                                record.get("reference_wikidata_id")
                            ),
                        }
                    )
    if rows:
        return rows

    fallback = []
    for species in species_like_terms(all_taxa):
        parts = species.split()
        fallback.append(
            {
                "reference_id": "",
                "doi": "",
                "taxonomy_provider": "allTaxa_fallback",
                "organism_value": species,
                "species": species,
                "genus": parts[0] if parts else "",
                "family": family_from_terms(all_taxa),
                "kingdom": "",
                "phylum": "",
                "class": "",
                "order": "",
                "cleaned_organism_id": "",
                "taxon_wikidata_id": "",
                "reference_wikidata_id": "",
            }
        )
    return fallback


def first_species_like(values):
    terms = species_like_terms(values)
    return terms[0] if terms else ""


def family_from_terms(terms):
    for term in terms or []:
        text = first_text(term)
        if re.match(r"^[A-Z][A-Za-z-]*(aceae|idae|ceae)$", text):
            return text
    return ""


def base_row(doc, source_file):
    traditional = first_text(doc.get("traditional_name"))
    iupac = first_text(doc.get("iupac_name"))
    return {
        "lotus_id": first_text(doc.get("lotus_id")),
        "wikidata_id": wikidata_qid(doc.get("wikidata_id")),
        "compound_name": first_text(traditional, iupac, doc.get("lotus_id")),
        "traditional_name": traditional,
        "iupac_name": iupac,
        "smiles": first_text(doc.get("smiles"), doc.get("smiles2D")),
        "inchikey": first_text(doc.get("inchikey"), doc.get("inchikey2D")),
        "molecular_formula": first_text(doc.get("molecular_formula")),
        "chemicalTaxonomyNPclassifierPathway": first_text(
            doc.get("chemicalTaxonomyNPclassifierPathway")
        ),
        "chemicalTaxonomyNPclassifierSuperclass": first_text(
            doc.get("chemicalTaxonomyNPclassifierSuperclass")
        ),
        "chemicalTaxonomyNPclassifierClass": first_text(
            doc.get("chemicalTaxonomyNPclassifierClass")
        ),
        "chemicalTaxonomyClassyfireKingdom": first_text(
            doc.get("chemicalTaxonomyClassyfireKingdom")
        ),
        "chemicalTaxonomyClassyfireSuperclass": first_text(
            doc.get("chemicalTaxonomyClassyfireSuperclass")
        ),
        "chemicalTaxonomyClassyfireClass": first_text(
            doc.get("chemicalTaxonomyClassyfireClass")
        ),
        "allTaxa": clean_text(doc.get("allTaxa")),
        "allWikidataIds": clean_text(doc.get("allWikidataIds")),
        "source_file": source_file,
    }


def compact_row(flat_row, retrieved_at):
    compound_name = first_text(flat_row.get("compound_name"))
    lotus_id = first_text(flat_row.get("lotus_id"))
    wikidata_id = first_text(flat_row.get("wikidata_id"))
    doi = first_text(flat_row.get("doi"))
    reference_id = first_text(flat_row.get("reference_id"), lotus_id, wikidata_id)
    evidence_pieces = [
        f"LOTUS ID: {lotus_id}" if lotus_id else "",
        f"taxon: {first_text(flat_row.get('species'))}",
        f"provider: {first_text(flat_row.get('taxonomy_provider'))}",
        f"reference: {reference_id}" if reference_id else "",
        f"pathway: {first_text(flat_row.get('chemicalTaxonomyNPclassifierPathway'))}",
        f"superclass: {first_text(flat_row.get('chemicalTaxonomyNPclassifierSuperclass'))}",
        f"class: {first_text(flat_row.get('chemicalTaxonomyNPclassifierClass'))}",
        f"DOI: {doi}" if doi else "",
    ]
    return {
        "species": first_text(flat_row.get("species")),
        "species_clean": clean_name(flat_row.get("species")),
        "genus": first_text(flat_row.get("genus")),
        "family": first_text(flat_row.get("family")),
        "compound_name": compound_name,
        "compound_name_clean": clean_key(compound_name),
        "lotus_id": lotus_id,
        "wikidata_id": wikidata_id,
        "cid": "",
        "smiles": first_text(flat_row.get("smiles")),
        "inchikey": first_text(flat_row.get("inchikey")),
        "molecular_formula": first_text(flat_row.get("molecular_formula")),
        "source_database": "LOTUS",
        "source_record_id": first_text(lotus_id, wikidata_id),
        "evidence_text": "; ".join(piece for piece in evidence_pieces if piece),
        "evidence_url": "https://lotus.naturalproducts.net/download",
        "reference_id": reference_id,
        "pmid": "",
        "doi": doi,
        "plant_part": "",
        "tissue": "",
        "method": "",
        "chemical_class_pathway": first_text(
            flat_row.get("chemicalTaxonomyNPclassifierPathway")
        ),
        "chemical_class_superclass": first_text(
            flat_row.get("chemicalTaxonomyNPclassifierSuperclass")
        ),
        "chemical_class_class": first_text(
            flat_row.get("chemicalTaxonomyNPclassifierClass")
        ),
        "all_chem_classifications": "; ".join(
            piece
            for piece in [
                first_text(flat_row.get("chemicalTaxonomyNPclassifierPathway")),
                first_text(flat_row.get("chemicalTaxonomyNPclassifierSuperclass")),
                first_text(flat_row.get("chemicalTaxonomyNPclassifierClass")),
                first_text(flat_row.get("chemicalTaxonomyClassyfireKingdom")),
                first_text(flat_row.get("chemicalTaxonomyClassyfireSuperclass")),
                first_text(flat_row.get("chemicalTaxonomyClassyfireClass")),
            ]
            if piece
        ),
        "all_taxa": first_text(flat_row.get("allTaxa")),
        "retrieved_at": retrieved_at,
        "source_file": first_text(flat_row.get("source_file")),
    }


def compact_key(row):
    return (
        row["species_clean"],
        row["compound_name_clean"],
        row["source_record_id"],
        row["lotus_id"],
        row["wikidata_id"],
        row["pmid"],
        row["doi"],
    )


def lookup_prefix(index_key, prefix_length=2):
    key = clean_key(index_key)
    if not key:
        return "_"
    return key[: max(1, prefix_length)]


def compact_lookup_rows(compact):
    keys = [
        ("species", first_text(compact.get("species_clean"))),
        ("genus", clean_name(compact.get("genus"))),
        ("family", clean_name(compact.get("family"))),
    ]
    seen = set()
    rows = []
    for key_type, key in keys:
        if not key or (key_type, key) in seen:
            continue
        seen.add((key_type, key))
        row = {"index_key_type": key_type, "index_key": key}
        row.update(compact)
        rows.append(row)
    return rows


class LookupShardWriter:
    def __init__(
        self,
        lookup_dir,
        source_path,
        member,
        prefix_length=2,
        max_open=64,
        layout="exact",
    ):
        self.lookup_dir = lookup_dir
        self.source_path = source_path
        self.member = member
        self.prefix_length = prefix_length
        self.max_open = max_open
        self.layout = layout if layout in {"exact", "prefix"} else "exact"
        self.handles = OrderedDict()
        self.file_counts = {}
        self.key_counts = {}
        self.row_count = 0
        os.makedirs(self.lookup_dir, exist_ok=True)

    def _relative_path(self, row):
        prefix = lookup_prefix(row["index_key"], self.prefix_length)
        if self.layout == "prefix":
            return os.path.join("shards", f"{prefix}.csv")
        key_type = clean_key(row["index_key_type"]) or "unknown"
        key_slug = clean_key(row["index_key"]) or "_"
        return os.path.join("keys", key_type, prefix, f"{key_slug}.csv")

    def _open_writer(self, relative_path):
        if relative_path in self.handles:
            handle, writer = self.handles.pop(relative_path)
            self.handles[relative_path] = (handle, writer)
            return writer
        if len(self.handles) >= self.max_open:
            _, (old_handle, _) = self.handles.popitem(last=False)
            old_handle.close()
        path = os.path.join(self.lookup_dir, relative_path)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        is_new = not os.path.exists(path)
        handle = open(path, "a", newline="", encoding="utf-8")
        writer = csv.DictWriter(handle, fieldnames=LOOKUP_FIELDS)
        if is_new:
            writer.writeheader()
        self.handles[relative_path] = (handle, writer)
        return writer

    def writerow(self, row):
        relative_path = self._relative_path(row)
        writer = self._open_writer(relative_path)
        writer.writerow(row)
        self.row_count += 1
        self.file_counts[relative_path] = self.file_counts.get(relative_path, 0) + 1
        key_type = row["index_key_type"]
        self.key_counts[key_type] = self.key_counts.get(key_type, 0) + 1

    def close(self):
        for handle, _ in self.handles.values():
            handle.close()
        self.handles.clear()

    def write_manifest(self, doc_count, compact_count, elapsed):
        manifest = {
            "format": "uafR_lotus_lookup_index",
            "version": 1,
            "created_at": datetime.now(timezone.utc)
            .replace(microsecond=0)
            .isoformat(),
            "source_path": self.source_path,
            "bson_member": self.member,
            "lookup_layout": self.layout,
            "shard_prefix_length": self.prefix_length,
            "compact_columns": COMPACT_FIELDS,
            "lookup_columns": LOOKUP_FIELDS,
            "bson_document_count": doc_count,
            "compact_row_count": compact_count,
            "lookup_row_count": self.row_count,
            "key_type_row_counts": dict(sorted(self.key_counts.items())),
            "shards": [
                {
                    "path": relative_path,
                    "prefix": os.path.basename(relative_path).replace(".csv", ""),
                    "row_count": self.file_counts[relative_path],
                }
                for relative_path in sorted(self.file_counts)
            ],
            "elapsed_seconds": round(elapsed, 1),
            "notes": (
                "Rows are keyed by normalized species, genus, and family terms "
                "so uafR can query LOTUS without reading the complete compact CSV."
            ),
        }
        with open(
            os.path.join(self.lookup_dir, "manifest.json"), "w", encoding="utf-8"
        ) as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)


def flatten(
    path,
    out_file,
    compact_index_file,
    lookup_dir,
    lookup_layout,
    member,
    limit,
    progress_every,
):
    started = time.time()
    archive, stream, source_label = open_bson_source(path, member)
    doc_count = 0
    row_count = 0
    compact_count = 0
    retrieved_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    compact_seen = set()
    lookup_writer = None
    try:
        if lookup_dir:
            lookup_writer = LookupShardWriter(
                lookup_dir=lookup_dir,
                source_path=path,
                member=member,
                layout=lookup_layout,
            )
        with stream:
            with open(out_file, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
                writer.writeheader()
                compact_handle = None
                compact_writer = None
                try:
                    if compact_index_file:
                        compact_handle = open(
                            compact_index_file, "w", newline="", encoding="utf-8"
                        )
                        compact_writer = csv.DictWriter(
                            compact_handle, fieldnames=COMPACT_FIELDS
                        )
                        compact_writer.writeheader()
                    for doc in read_documents(stream):
                        doc_count += 1
                        base = base_row(doc, source_label)
                        if not base["compound_name"]:
                            continue
                        tax_rows = taxonomy_rows(
                            doc.get("taxonomyReferenceObjects"),
                            doc.get("allTaxa"),
                        )
                        for tax_row in tax_rows:
                            row = dict(base)
                            row.update(tax_row)
                            writer.writerow(row)
                            row_count += 1
                            if compact_writer is not None or lookup_writer is not None:
                                compact = compact_row(row, retrieved_at)
                                key = compact_key(compact)
                                if key not in compact_seen:
                                    compact_seen.add(key)
                                    if compact_writer is not None:
                                        compact_writer.writerow(compact)
                                    if lookup_writer is not None:
                                        for lookup_row in compact_lookup_rows(compact):
                                            lookup_writer.writerow(lookup_row)
                                    compact_count += 1
                        if limit and doc_count >= limit:
                            break
                        if progress_every and doc_count % progress_every == 0:
                            elapsed = time.time() - started
                            print(
                                f"processed {doc_count:,} BSON docs; "
                                f"wrote {row_count:,} flat rows; "
                                f"{compact_count:,} compact rows; "
                                f"{elapsed:.1f} sec",
                                file=sys.stderr,
                                flush=True,
                            )
                finally:
                    if compact_handle is not None:
                        compact_handle.close()
    finally:
        elapsed = time.time() - started
        if lookup_writer is not None:
            lookup_writer.close()
            lookup_writer.write_manifest(doc_count, compact_count, elapsed)
        if archive is not None:
            archive.close()
    return doc_count, row_count, compact_count, elapsed


def build_lookup_from_compact(
    compact_index_file,
    lookup_dir,
    lookup_layout,
    progress_every,
):
    started = time.time()
    row_count = 0
    lookup_writer = LookupShardWriter(
        lookup_dir=lookup_dir,
        source_path=compact_index_file,
        member="compact_csv",
        layout=lookup_layout,
    )
    try:
        with open(compact_index_file, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                compact = {field: first_text(row.get(field)) for field in COMPACT_FIELDS}
                if not compact["compound_name"] or not (
                    compact["species_clean"] or compact["genus"] or compact["family"]
                ):
                    continue
                for lookup_row in compact_lookup_rows(compact):
                    lookup_writer.writerow(lookup_row)
                row_count += 1
                if progress_every and row_count % progress_every == 0:
                    elapsed = time.time() - started
                    print(
                        f"processed {row_count:,} compact rows; "
                        f"wrote {lookup_writer.row_count:,} lookup rows; "
                        f"{elapsed:.1f} sec",
                        file=sys.stderr,
                        flush=True,
                    )
    finally:
        elapsed = time.time() - started
        lookup_writer.close()
        lookup_writer.write_manifest(None, row_count, elapsed)
    return row_count, lookup_writer.row_count, elapsed


def parse_args():
    parser = argparse.ArgumentParser(
        description="Flatten official LOTUS Mongo BSON dump to CSV."
    )
    parser.add_argument("--mongo-zip", "--input", dest="input", default="")
    parser.add_argument("--out-file", default="")
    parser.add_argument("--compact-index-file", default="")
    parser.add_argument(
        "--from-compact-index",
        default="",
        help=(
            "Build only --lookup-dir from an existing compact LOTUS CSV. This "
            "streams the compact file and avoids rewriting the flat BSON export."
        ),
    )
    parser.add_argument(
        "--lookup-dir",
        default="",
        help=(
            "Optional directory for a manifest-backed species/genus/family "
            "lookup index that uafR can query without reading the full CSV."
        ),
    )
    parser.add_argument(
        "--lookup-layout",
        choices=("exact", "prefix"),
        default="exact",
        help=(
            "Lookup index layout. 'exact' writes one file per normalized "
            "species/genus/family key and is recommended for production. "
            "'prefix' writes broader two-character shard files."
        ),
    )
    parser.add_argument("--member", default=DEFAULT_MEMBER)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=25000)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.from_compact_index:
        if not args.lookup_dir:
            raise SystemExit("--from-compact-index requires --lookup-dir")
        if not os.path.exists(args.from_compact_index):
            raise SystemExit(
                f"Compact index file does not exist: {args.from_compact_index}"
            )
        if os.path.exists(args.lookup_dir):
            if not args.overwrite:
                raise SystemExit(
                    "Lookup directory exists and --overwrite was not supplied: "
                    f"{args.lookup_dir}"
                )
            shutil.rmtree(args.lookup_dir)
        os.makedirs(os.path.abspath(args.lookup_dir), exist_ok=True)
        row_count, lookup_count, elapsed = build_lookup_from_compact(
            compact_index_file=args.from_compact_index,
            lookup_dir=args.lookup_dir,
            lookup_layout=args.lookup_layout,
            progress_every=args.progress_every,
        )
        print("LOTUS lookup index build complete.")
        print(f"Compact input: {args.from_compact_index}")
        print(f"Lookup index: {args.lookup_dir}")
        print(f"Compact rows processed: {row_count}")
        print(f"Lookup rows written: {lookup_count}")
        print(f"Elapsed seconds: {elapsed:.1f}")
        return

    if not args.input:
        raise SystemExit("--input is required unless --from-compact-index is used")
    if not args.out_file:
        raise SystemExit("--out-file is required unless --from-compact-index is used")
    if os.path.exists(args.out_file) and not args.overwrite:
        raise SystemExit(
            f"Output file exists and --overwrite was not supplied: {args.out_file}"
        )
    if (
        args.compact_index_file
        and os.path.exists(args.compact_index_file)
        and not args.overwrite
    ):
        raise SystemExit(
            "Compact index file exists and --overwrite was not supplied: "
            f"{args.compact_index_file}"
        )
    if args.lookup_dir and os.path.exists(args.lookup_dir):
        if not args.overwrite:
            raise SystemExit(
                "Lookup directory exists and --overwrite was not supplied: "
                f"{args.lookup_dir}"
            )
        shutil.rmtree(args.lookup_dir)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_file)), exist_ok=True)
    if args.compact_index_file:
        os.makedirs(
            os.path.dirname(os.path.abspath(args.compact_index_file)),
            exist_ok=True,
        )
    if args.lookup_dir:
        os.makedirs(os.path.abspath(args.lookup_dir), exist_ok=True)
    doc_count, row_count, compact_count, elapsed = flatten(
        path=args.input,
        out_file=args.out_file,
        compact_index_file=args.compact_index_file,
        lookup_dir=args.lookup_dir,
        lookup_layout=args.lookup_layout,
        member=args.member,
        limit=args.limit,
        progress_every=args.progress_every,
    )
    print("LOTUS Mongo flatten complete.")
    print(f"Input: {args.input}")
    print(f"Output: {args.out_file}")
    if args.compact_index_file:
        print(f"Compact index: {args.compact_index_file}")
    if args.lookup_dir:
        print(f"Lookup index: {args.lookup_dir}")
    print(f"BSON documents processed: {doc_count}")
    print(f"Flat rows written: {row_count}")
    if args.compact_index_file or args.lookup_dir:
        print(f"Compact rows written: {compact_count}")
    print(f"Elapsed seconds: {elapsed:.1f}")


if __name__ == "__main__":
    main()
