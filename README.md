# Data analysis for the manuscript: "Taxonomic diversity is more sensitive than functional diversity to indicate mild mining impacts in stream fish assemblages"

## Authors
Nathália C. S. Silva, Bruno E. Soares; Fabrício B. Teresa, Érica P. Caramaschi, Miriam P. Albrecht¹

## Published in
Preprint not available. Submitted to Aquatic Ecology.

## To cite this data
Citation and doi for citing this data and coding are available by Zenodo: https://zenodo.org/record/4118631

## Description
Data and script to analyze the indirect effects of the mining activities at the National Forest of Saracá-Taquera, Trombetas River Basin, Central Amazon. These results are discussed in the manuscript "Taxonomic diversity is more sensitive than functional diversity to indicate mild mining impacts in stream fish assemblages", submitted to Aquatic Ecology.

## How to use this directory
In R, you may download all the data contained in this repository by using the download.file() function, then use the function unzip(). Example:

`
download.file(url="https://github.com/bruno-soares/MS_TaxFuncDivTrombetas/archive/master.zip", destfile = "MS_TaxFuncDivTrombetas.zip")
`

`
unzip(zipfile = "MS_TaxFuncDivTrombetas.zip")
`

R scripts are numbered in the order that the analyses are performed, from 01 to 04. All the data necessary to run analyses are availagle in /data and /results, so one can run any part of the script independently from running prior analysis. Make sure you install the necessary R packages before running the analysis.

## Data description
Data encompasses information on the fish assemblages and abiotic conditions of pristine and disturbed streams in the National Forest Saracá-Taquera (Trombetas River Basin, Brazil). Data is available in five files:

1. ab.rel_orders: average percentage of sampled taxonomic orders.

2. abiotic variables: environmental conditions at sampling sites during the expedition.

3. cpue_sRec: abundance of fishes sampled at each stream.

4. grouping variables: additional information on sampling sites that is used in the scripts.

5. traits: functional traits measured for the sampled fishes. Traits are measured by species and includes both dietary information and ecomorphological traits.

## Contact information
Please, feel free to contact me in my personal e-mail: soares.e.bruno@gmail.com.
