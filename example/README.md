# Example dataset

`prepare_example.sh` downloads five small, complete public sequences from NCBI and
simulates PacBio HiFi-like reads (~12 kb, ~0.1 % error) from them with
`simulate_hifi.py`:

| Accession | Organism / element | Size | Depth |
|---|---|---|---|
| NC_000908.2 | *Mycoplasmoides genitalium* G37 chromosome | 580 kb | 40× |
| NC_002528.1 | *Buchnera aphidicola* APS chromosome | 641 kb | 30× |
| NC_002252.1 | *B. aphidicola* plasmid pLeu | 7.3 kb | 80× |
| NC_002253.1 | *B. aphidicola* plasmid pTrp | 7.8 kb | 80× |
| NC_001422.1 | phage phiX174 | 5.4 kb | 60× |

```bash
bash example/prepare_example.sh                      # -> example/data/reads.fq.gz
ouroassemble run -r example/data/reads.fq.gz -o example/output -t 8
```

Expected result: both chromosomes recovered as circular, junction-validated
`chromosome`s (one representative each after dereplication), plus small circular
elements classified by geNomad.
