# Mammary Gland Transcriptomes - Shiny App

App for browsing mammary gland transcriptome data. 
The App is currently running on the internal server from the Ferguson-Smith group at: http://jeremy-bio.gen.private.cam.ac.uk:3838/mammary_dataviz/ 

The script `prepare_data.R` makes all the necessary files used in the App and fetches processed data from two separate folders, expected to be one level above the current project: 

- `../mammary_cellsort_hybrid_rnaseq` - these are data from hybrid mice ([code]())
- `../mammary_cellsort_zfp57_rnaseq` - data from _zfp57_ mutants ([code]())