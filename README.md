# Short Project in Bioinformatics-Lund University
## Project scripts description

My short project entitled "Functional annotation of potential regulators of lysosomes and cell death identified in a genome-wide screen " involved scripts that carry out two major tasks.The first task primarly involves data download and analysis of genes collected from 12 databases related to autophagy,cell death,and lysosomes. However, the second task involves functional annotation of genes obtained from a genome-wide screen and basic statistical analysis of knockdown data. 
Scripts of the above tasks require the installation and import of various packages(see packages.py) and Matplotlib library(v3.1.1). <br>
Some databases were available for download(see datadownload.py), while others were not and therefore the htmls were downloaded (datadownload.py) and webscraping scripts were written to parse them(see webscraping.py). Other scripts include parse.py, mergeAutophagy.py,mergeCellDeath.py, mergeLysosome.py,overlap.py,screenCheck.py, plotting.py,MappingToUniprot.py, and MappingToLists.py.

## Dependencies

* python v3.7.4
* wget
* matplotlib v3.1.1
* pandas
* beautifulsoup4
* owlready2
* goatools

All the scripts were run using python. 

## Main Scripts
* packages.py<br>
This script imports the packages used in this short project.

* datadownload.py<br>
Involves the codes written for downloading the databases that were available for download or the htmls of databases that were not available for download.<br>
Some codes should be run a long with the correct packages from packages.py.<br>
One file was downloaded manually and therefore it is not included in this set of codes. It is from uniprot/subcellular-lysosome section(release 2020_02).

* webscraping.py<br>
Runs codes written for parsing databases that were not available for download. HTML scraping was the method used in this case. Some scraping codes were followed by another code to fix the data structure in the database.

* parse.py<br>
Involves codes written for parsing the databases that were available for download in order to extract necessary fields.

## Analysis Scripts
* mergeAutophagy.py, mergeCellDeath.py, mergeLysosome.py<br>
Run codes written for merging the databases from each category in order to collect all available information related to autophagy, cell death and lysosomes.


* overlap.py<br>
Runs codes written for finding the overlap between the end files of autophagy-cell death, autophagy-lysosome , cell death-lysosome as well as autophagy-lysosome-cell death.


* screenCheck.py<br>
Runs codes written to parse a whole genome screen of knockdown data of all known human genes. The script will split the screen file to control and sample files. These belong to genes that showed at least 15% reduction in cell count upon the knockdown in replicate1 and/or replicate2.

* plotting.py<br>
Runs codes that were used in matplotlib library to plot positive- and negative control-counts (as a part of the statistical evaluation of false positives and false negatives) as well as the codes written to get the counts that are used in generating the plots.

* MappingToUniprot.py<br>
Runs codes for the files whose entrez IDs were mapped to uniprot in order to get the accession numbers for genes in the original screen file.

* MappingToLists.py <br>
Runs codes written for for matching the screen hits with the autophagy, cell death, and lysosome lists.

## Usage
All codes were written in Jupyter notebook (v6.0.1). The name or the path of the file should be supplied and the codes can be simply run with ctrl+enter OR alt+enter.

### Step-by-step
- create new environmente with conda
- install wget: pip install wget
- install matplotlib: conda install matplotlib==3.1.1
- install pandas: conda install pandas
- install beautifulsoup4: conda install beautifulsoup4
- install owlready2: conda install -c conda-forge owlready2
- install goatools: pip install goatools
- run packages.py to import packages and then datadownload.py to get the raw data, place data in folder "rawdata"

Alternative:
- install jupyterlab or notebook, e.g. 'conda install -c conda-forge jupyterlab (if you want to look at the notebook)
- launch jupyterlab from the anaconda/miniconda prompt: jupyter lab


## Other Information
* Versions and release dates of databases were recorded if they were available, otherwise the date of data collection was written.
* All unreviewed, obsolete, and unmapped-data files that were obtained upon mapping our data to UniProt are found in "Incomplete_Mappings" folders. The final and reviewed-data files are saved to "Final_Databases" folder.
* 3 files were generated for the sample files in the screen analysis and they correspond to the 3 thresholds used in the study to assess false positives and false negatives.


### References
#### Autophagy Database

http://www.tanpaku.org/autophagy/

dx.doi.org/10.1093/nar/gkq995
