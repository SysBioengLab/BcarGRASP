# Transcripts Data Pre-processing
In this folder, the script processData processes the data stored in the folder rt-qpcr_data to relativize transcript quantities between conditions. To select a different reference condition, the script inputs must be modified accordingly. The scripts plotHeatmapD010 and plotHeatmapD025 create heatmaps of the relative transcript quantities once they are obtained by processData.

## additional_code
Includes functions called by processData to complete the estimation of relative transcript abundances.

## output
- ConditionsRelativeLimits.xlsx: 95% confidence intervals and mean relative transcript quantities. Each sheet is a different reference condition.
- LinearRegressionStatistics.xlsx: Statistics of the linear regression results for each of the genes' standard calibration curve. Each sheet is a gene.
- AppliedNormalizations.xlsx: Genes used as references to normalize each of the experimental conditions. Each sheet is a different reference condition.

## output/images
- linear_regressions: each image are the linear regressions of a gene standard calibration curve.
- relative_expression: each image represents the relative transcript quantities of a condition. The nomenclature of files is X_Y, where X is the reference condition and Y the relative condition.
- gene_relative_expression: each image represents the relative transcript quantities of a gene. The nomenclature is X_Y, where X is the reference condition and Y the gene relativized.
- heatmaps: heatmaps fo relative transcript quantities. The nomenclature is X_Y, where X is the dilution rate that is being represented, and Y the reference condition.