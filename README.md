# OncoMergeAnalysis
Contains scripts for analyzing and understanding the output from OncoMerge including value added assessments.

Scripts included:

1. Assess_Value_Added.py
    - Uses R to plot OncoPlots.
2. Plot_Muts.py
4. OncoMergeWrapper.py


# Assess_Value_Added.py:
This script creates visualizations and tables to help show what OncoMerge adds to downstream analysis. 
1. Reads data from OncoMerged_file into glofmf as the dictionary of dictionaries {Cancer1:gene1:MutType1:Mutation Frequency}
    * Parses OncoMerged_file
    * Create dataframe for violinplot of individual cancers
      * Calculate delta (frequency gain with combined PAMs and CNAs)  
      * Classify as recovered mutation or reinforced mutation
  2. Create OncoPlot Table
    * Load in cancer data
      * Extract Stats for OncoPlotTable
      * Create OncoPlot Table
    * Save some data for Enrichment Analyses
 3. Define the Order of each violin plot
 4. Make Violin Plots of deltas
 5. Plot OncoPlots for top 10 genes
 6. Enrichment Analysis
    * Enrichment of Enrichment of recovered PAMs in Mutsigcv2
    * Enrichment of all OncoMerged PAMs in Mutsigcv2
    * Enrichment of recovered CNAs in sig gistic
7. Plot Heat Map of CNA/PAM cutoffs for greatest LoF/GoF enrichment
8. Plot Venn diagram of genes in MutSigCV, Gistic, and OncoMerge

# Plot_Muts.py
This script creates visualizations of OncoMerged datasets to assess OncoMerge's proper functionality. 
1. Set up
2. Plot Bar plot and raw csv of number of each CNAtype in OncoMerged file
3. Plot Histograms to plot the distribution of cancers for the number of different CNAtypes
4. Plot Violin Plots of mutaitonal frequency

# OncoMergeWrapper.py
This script is the user interface used to run OncoMerge.
1. Enter Parameters and File Paths Here
2. Run OncoMerge



    










