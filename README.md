# overwriTE
## Interrogating the influence of Transposable Element (TE) insertions on transcriptional dynamics   
  
overwriTE dependencies include: mySQL, python3, R (tidyr, ggextra)    

This figure below outlines the various components of a transcript that overwriTE will intersect with a TE insertion:  

![OverwriTE classification](images/AMG_INSERT_PLOT_v2.png)  

By intersecting overwriTE's intermediate .csv file with DEseq2 output, the potential influence of TEs can be quantified and visualized as seen below:  

Hexbin             |  Boxplot
:-------------------------:|:-------------------------:
![Hexplot demonstrating the density of TE coverage across all parts of the transcript](images/hexbin_example.png?raw=true)  |  ![Boxplots highlighting the statistically significant difference between expression and TE coverage](images/boxplot_example.png)  
