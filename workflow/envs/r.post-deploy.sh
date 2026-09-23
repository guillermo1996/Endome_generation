#!env bash
set -o pipefail

Rscript -e 'BiocManager::install("guillermo1996/txendcutr", update = FALSE, ask = FALSE)'
Rscript -e 'BiocManager::install("guillermo1996/grpSciRmdTheme", update = FALSE, ask = FALSE)'
Rscript -e 'devtools::install_github("dzhang32/ggtranscript")'
Rscript -e 'BiocManager::install("scrapper", update = FALSE, ask = FALSE)'