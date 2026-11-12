# SPREAD-Hackathon-viz

**PlasmidLens** is an interactive R Shiny dashboard for exploring plasmid collections. It brings together the outputs of common plasmid-analysis pipelines: metadata, antimicrobial resistance (AMR) genes, mobility, mobile genetic elements, genetic context and plasmid clustering. You can browse them in one place without writing any code.

Input files are recognised by their column headers, not their file names. The app joins them by plasmid accession, so you can supply any combination of supported files. A missing file type only disables the plots that need its fields.

## Features

The app is organised into five tabs:

### 1. Overview
- Headline showing total genomes, plasmids, and replicon types
- Summary statstics for genomes, including isolate source and country
- Plasmid size and GC content
- Species and MLST distribution
- Plasmid replicon distributions

### 2. AMR gene associations
- Headline showing total plasmids, ARG types, Rep types, IS types, transposon types and virulence factor
- ARG frequency (customisable filters)
- Count plots for Rep, MOB, MPF types, IS elements, transposons, and virulence factors

### 3. Genetic context
- **Single plasmid:** select any annotated CDS from a Bakta annotation and draw its gene neighbourhood with `gggenes`. You can choose how many CDSs to show upstream and downstream of it. The tab also shows a feature table you can export (copy, CSV or Excel) and details of the selected feature.
- **Multiple plasmids:** load an interactive Clinker comparison (`plot.html`) and choose an anchor gene (homology group, e.g. *blaNDM-1*). The plot is then trimmed to N genes upstream and downstream of that gene in every cluster.

### 4. Consensus network
- Compares plasmid relatedness from different clustering tools, such as Pling and mge-cluster.
- An interactive network with edges coloured by agreement between tools.
- Adjustable thresholds: maximum containment distance, maximum DCJ-Indel distance and an option to require the same Pling community
- Nodes can be coloured by mge-cluster, Pling community or Pling subcommunity
- A heatmap and cross-tabulation of Pling subcommunities against mge-clusters
- A pairwise relationship table and input diagnostics (for example, IDs that don't match between the tools)
- Download options for edge table (CSV), node table (CSV) and network (GraphML)

### 5. Data import
- Upload one or more CSV/TSV/TXT pipeline outputs, or reload files from the data folder
- Preview the merged plasmid metadata table and the pairwise distance table
- An import summary showing recognised files, detected fields, column mappings and data notes

## Supported input files

| Tool / file | Used for | How it is recognised |
|---|---|---|
| Plasmid metadata | Isolate, species, MLST, source, country, size, GC, etc. | `accession` / `plasmid_id` plus metadata columns such as `host species`, `country`, `source type` |
| [AMRFinderPlus](https://github.com/ncbi/amr) | ARG genes, AMR drug classes, virulence genes | `Element symbol` / `Gene symbol` / `amr_genes` / `drug_classes` columns |
| [MOB-suite](https://github.com/phac-nml/mob-suite) | Replicon type, length, GC, MOB type, MPF type | `accession` plus `predicted_mobility`, `relaxase_type(s)` or `mpf_type` |
| [ISFinder](https://isfinder.biotoul.fr/) | IS elements | `sequence`, `start` and `gene` / `IS family` columns |
| [MobileElementFinder](https://pypi.org/project/MobileElementFinder/) | IS elements and transposons | `contig`, `name`, `type` and `prediction` / `mge_no` columns |
| Mash distances (optional) | Pairwise distance table | `plasmid1`, `plasmid2`, `mash_distance` |
| [Bakta](https://github.com/oschwengers/bakta) annotation table (`.tsv`/`.csv`) | Genetic context – single plasmid | `Sequence Id` header row. `.gff3` and GenBank (`.gbff`/`.gb`/`.gbk`) are also accepted when uploaded. |
| [Clinker](https://github.com/gamcil/clinker) HTML plot | Genetic context – multiple plasmids | HTML file containing the embedded `const data=` object |
| [Pling](https://github.com/iqbal-lab-org/pling) typing + distances | Consensus network | `plasmid`, `type` (typing), and `plasmid_1`, `plasmid_2`, `containment_distance`, `DCJI_distance` (distances) |
| [mge-cluster](https://gitlab.com/sirarredondo/mge-cluster) results | Consensus network | `Sample_Name`, `Standard_Cluster` |

Common alternative column names are recognised automatically (for example `organism` → `species`, `ST` → `mlst`). The applied mappings are listed in the import summary.

## Repository structure

```
.
├── app.R                  # The Shiny app
├── data/                  # Example pipeline outputs (metadata, AMRFinder, MOB-suite, ISFinder, MobileElementFinder)
├── genetic_context/       # Example Bakta table (.tsv), GenBank file and Clinker plot.html
└── plasmid_cluster/       # Example Pling typing/distances and mge-cluster results
```

## Installation

Install the required R packages once:

```r
install.packages(c("shiny", "bslib", "plotly", "DT", "dplyr", "stringr", "readr",
                   "igraph", "visNetwork", "ggplot2", "gggenes", "jsonlite"))
```

Optional, needed only for `.gff3` / GenBank annotation input (from Bioconductor):

```r
BiocManager::install(c("rtracklayer", "genbankr"))
```

## Running the app

1. Put `app.R` in the project folder and your input files in a `data/` folder next to it.
2. Start the app from R:
   ```r
   shiny::runApp("path/to/SPREAD-Hackathon-viz")
   ```
   Or open `app.R` in RStudio and click **Run App**.

The app reads every CSV/TSV/TXT file in the data folder at startup. It looks for the data folder in this order:
1. The folder set in the `PLASMID_DATA_DIR` environment variable, e.g. `Sys.setenv(PLASMID_DATA_DIR = "path/to/data")`
2. `data/` next to `app.R`
3. `data/` next to the file open in RStudio
4. `data/` in the current working directory

The **Genetic context** and **Consensus network** tabs also check the data folder (every 5 seconds) for their input files. To try them with the example files, upload the files from `genetic_context/` and `plasmid_cluster/` in those tabs, or copy them into `data/`. An uploaded file always takes priority over one found in the folder.

Uploads of up to 200 MB are allowed.
