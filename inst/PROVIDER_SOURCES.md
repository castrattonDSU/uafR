# Public Provider Sources And Use Boundaries

Last reviewed for uafR documentation: 2026-07-16.

uafR queries or parses records from public chemical, biochemical,
natural-product, regulatory, and literature resources. The package does not
redistribute a frozen copy of those providers' databases. Small simulated test
fixtures are labeled as simulated. Users who download source exports or retain
API responses are responsible for reviewing the provider's current access,
attribution, licensing, and redistribution terms for their use case.

| Provider | Official entry point | uafR use | Interpretation boundary |
| --- | --- | --- | --- |
| PubChem | <https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest> | Compound identity, properties, structure, annotations, source links, and literature context | A PubChem record or annotation is not confirmation in a project sample. |
| KEGG | <https://www.kegg.jp/kegg/rest/keggapi.html> | Exact-normalized compound synonyms plus pathway, reaction, enzyme, module, and link context; rejected broad name candidates remain in an audit table | A broad name-search hit is not accepted as compound identity, and pathway association is not evidence of pathway activity. |
| LOTUS | <https://lotus.naturalproducts.net/download> | Manifest-backed local indexes from the frozen downloadable snapshot; bounded live discovery is retained only for small probes | Reported taxon-compound occurrence is not a complete metabolome or new sample measurement; snapshot version and checksums must remain in provenance. |
| KNApSAcK | <https://www.knapsackfamily.com/knapsack_core/top.php> | Throttled, cached organism lookup through the documented `info.php` flow with exact accepted-name or verified-alias filtering | Forward matches and partial taxon names are not accepted; raw responses are kept out of redistributable bundles pending terms review. |
| NPASS | <https://bidd.group/NPASS/downloadnpass.html> | Local manifest-backed index built from the official NPASS 3.0/NPASS-2026 general, structure, species-source, and taxonomy files | Source records require identity and taxon review; raw files remain external project resources and no live species endpoint is inferred. |
| PubMed / NCBI E-utilities | <https://www.ncbi.nlm.nih.gov/books/NBK25501/> | Literature discovery, summaries, and abstracts | Search retrieval alone does not establish chemical occurrence. |
| PubTator | <https://www.ncbi.nlm.nih.gov/research/pubtator3/api> | Batched BioC annotations for PMIDs returned by the PubMed discovery layer | Chemical/species co-mentions remain candidate evidence unless the underlying report is curated; requests stay below the documented service limit. |
| FDA Structured Product Labeling | <https://www.fda.gov/industry/fda-data-standards-advisory-board/structured-product-labeling-resources> | Source-filtered PubChem regulatory and label annotations | Label presence is not exposure, hazard probability, or risk assessment. |
| FEMA | <https://www.femaflavor.org/> | Source-filtered PubChem flavor and fragrance annotations | Source annotation is not a measured abundance or sensory response. |
| MeSH | <https://www.nlm.nih.gov/mesh/meshhome.html> | Source-filtered biomedical vocabulary annotations | Vocabulary membership is not a clinical mechanism or outcome. |

`uafRProviderContracts()` returns the machine-readable interpretation and
diagnostic contract used by the package. Provider availability, schemas, and
terms can change independently of uafR. Live workflows therefore remain
opt-in, cached, throttled, and accompanied by provider diagnostics.
