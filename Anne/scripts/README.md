# In depth comparison of somatic mutations in cancer, shows specificity of non-coding mutations near known breast cancer genes

## About the project
Somatic mutations may lead to the conversion of normal cells into cancer cells. Until now, research has focused on somatic mutations in the coding regions of the genome, and hardly at the non-coding regions of the genome. Here, we focus on creating methods for finding mutations or regions containing mutations that are non-coding, that can be linked to cancer. The focus is mainly on the following non-coding regulatory regions: promoter, transcription factor binding sites, deoxyribonuclease regions and Ultra-Conserved Non-coding Elements. We created a pipeline to compare mutation data of 2,238 donors from the International Cancer Genome Consortium, in breast versus non-breast cancer, to untangle how non-somatic mutations can contribute to cancer phenotypes. Significant regions with the nearest gene being a known breast cancer gene were discovered in this project using the Fisher-exact test and Cochran–Armitage test for trend. However, these breast cancer genes were also discovered in regions where non-breast cancer samples were noticeably more prevalent. The second goal was to build a pipeline to identify key/interesting regions and SNPs. However, further research is needed to determine whether this pipeline is efficient and whether it needs to be modified or improved to achieve better results.

### Built with
-	Python 3.8.12
-	R 4.1.2

## Getting started
Most of the scripts were run on the Gearshift cluster, some others were run on the Calculon cluster, and still others were run on my own computer. This resulted in three yaml files. It is recommended to run everything on a cluster, because some scripts can take quite a while when run on your own computer/laptop. Also, some scripts create a lot of data.
The data used for the *scripts/pipeline_variant_calling* comes from Okosun et al. (2014). It consists of six paired-paired follicular lymphoma-transformed follicular lymphoma (FL-tFL)-germline (GL) samples and two additional relapsed follicular lymphoma (FL) samples for case S2.
The rest of scripts use a dataset obtained from the ICGC (see https://dcc.icgc.org/releases/current/Projects (*simple_somatic_mutation.open.[ICGC project code].tsv.gz*)). Open-access simple mutation calls were used. These mutations calls include single and multiple base substitutions and small (<= 200 base pairs (bp)) insertions and deletions that appear in the tumor tissue but not in the normal control tissues. These mutations are present in the coding part of the genome as well as in the non-coding part.


### Prerequisites
Jupyter Notebook, Python and R should be installed and working before the scripts can be used.

### Installation 
1.	Clone the repo

``` git clone https://github.com/molgenis/non-coding-somatic-mutations-in-cancer.git ```

2.	Install all the required python packages

```conda install –-file requirements.txt ```
 
## Usage
The folder scripts contains all scripts. The scripts are divided into five folders.
- *pipeline_variant_calling*: Contains the scripts for performing variant calling over BAM files

- *database*: contains the scripts for creating and populating the database

- *add_layers*: contains the scripts that add different layers to the database

- *analysis_layers*: contains scripts that test whether certain regions/SNPs differ significantly between breast cancer donors and non breast cancer donors. Multiple testing correction is also performed.

- *chromosome*: contains scripts for making plots. These plots represent the difference in a particular region between the amount of breast cancer SNPs and non-breast cancer SNPs. This was no longer used during the research, but could be used in the future.

- Furthermore, the script contains several yaml files and python scripts that are called from multiple scripts.