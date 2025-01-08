# Summary
R code associated with Martin et al.'s "Measuring demographic stochasticity and selection due to attrition in ex-situ collections for plant conservation"

These scripts are part of a project on provenance selection and opportunity for selection in ex-situ collections for plant conservation. The project uses as a study system a living collection of Quercus arkansana (Fagaceae) maintained in the Oertli Hardy Plant Nursery of the Missouri Botanical Garden.

The "viridis" R package was used for color palettes. https://github.com/sjmgarnier/viridis 

# Data preparation
Provenance data preparation 

  * Obtain occurrence and provenance data for Q. arkansana (GBIF).
  * Map known occurrences and collection localities of the unique accessions in the common garden experiment.
  * Assign accessions to geographic regions, following Thomas et al. (2023, Biological Conservation 283: 110052).
  * Save provenance data as a text file.

Survival data preparation

  * Examine survival survey dates.
  * Compare survival data between the two surveys.
  * Examine the distribution of plantings among accessions.
  * Examine the distribution of planting deaths among accessions.
  * Save the data as a text file.

# Null model for provenance selection
Illustration of the null model, representing demographic stochasticity, used to gauge provenance selection and opportunity for selection

  * Read and examine distribution of initial plants per maternal line per provenance (i.e., geographic region), and surviving plants per maternal line per provenance.
  * Calculate a demographic stochasticity null model.

# Provenance selection and opportunity for selection

  * Calculate selection coefficients for provenances (i.e., geographic regions).
  * Examine the relationship between the number of maternal lines and the mean survival rates across provenances.
  * Calculate opportunity for selection within each provenance.
  * Compare the null and observed opportunity for selection.

# Q. arkansana provenance and occurrence map

  * Map latitudes and longitudes obtained from GBIF data and collected provenance records. Plot created with ggplot2 https://cran.r-project.org/web/packages/ggplot2/index.html 

