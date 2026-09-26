# Manual review of flagged tree-input proteins (27 Sep 2026)

Multi-reference check (`verify_gene_calls.py`, DIAMOND very-sensitive, hit must
cover >= 50 % of the query, best hit's NCBI name compared across 5 reference
proteomes). 4964 proteins; 6 flagged. All 6 reviewed and kept:

| Gene | Genomes | Why flagged | Evidence | Verdict |
|---|---|---|---|---|
| fliN | Hippea alviniae, Nitrosophilus labii / alvini / kaiyonis, Nitratiruptor tergarcus | best hit FliY in 3-4 refs | Kept copy is 81-92 aa (ref FliN 95 aa; FliY ~287 aa). FliY's C-terminus is a FliN-type (SpoA) domain, so short divergent FliN scores similarly against both (bitscore 45-60, ~30 % id). The genome's *other* protein named fliN (164-330 aa, separate orthogroup) matches FliY better (e.g. 130 vs 94 bits) — a FliY mis-named fliN by NCBI, correctly dropped by the orthogroup rule. | correct (FliN) |
| flhF | Nitratiruptor tergarcus | best hit Ffh in 3 refs | Direct comparison: FlhF 141 bits vs Ffh 94 bits. The FlhF hit covers 49 % of the query (FlhF's variable N-terminal region), just under the 50 % cover filter, leaving Ffh (same SRP-GTPase family) as the top hit. | correct (FlhF) |

Not flagged but uninformative in this check — covered instead by the original
OrthoFinder + eggNOG pipeline (141/141 flgJ and 102/104 flgR tree inputs are
orthology-added calls with eggNOG `confirmed` status; the other 2 flgR carry the
NCBI gene name). No separate eggNOG run on the tree inputs was done:
- **flgJ**: none of the 5 reference genomes names its FlgJ, so every hit is "unnamed".
- **flgR**: only the C. jejuni reference names FlgR; in the H. pylori reference
  FlgR is a pseudogene (no protein). Before the coverage filter this produced
  11 spurious flags from receiver-domain-only hits (CrdR, CheY, HsrA).
