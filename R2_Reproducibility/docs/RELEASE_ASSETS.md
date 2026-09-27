# Release assets

Upload the following files from the local `release_staging/` directory when
creating the GitHub Release for this code revision.

| Asset | Size | SHA-256 | Purpose |
|---|---:|---|---|
| `FC_R2_Results_S1-S8.zip` | 80.53 MiB | `AECD2E386F1A69F17E086F1044A244AAAFD86DBC72389AEC9617221FBD275524` | Formal revised-manuscript benchmark results |
| `FC_R2_MODDPG_Audit.zip` | 12.63 MiB | `0082C74BAC41D26EA0B12A4132ECA7EDB5D9BB02E935C09CBC245D31BED0877E` | Separated implementation-audit evidence; not a source of manuscript metrics |

The first archive contains the eight canonical result files and
`manifest.csv`. The second contains the retained pre-correction S2 record,
numerical summaries, and response-letter figures. Keeping them as separate
assets prevents accidental mixing of formal and historical data.
