# Intraspecific Facilitation in the Context of Niche and Fitness Differences


## Files


Read me file for the article entitled 
"Intraspecific Facilitation in the Context of Niche and Fitness Differences"
by
Lisa Buche and Jürg W. Spaak
Published in 
""

Buche L. and Spaak W. J. (2025), Buchel9844/NFD_SpaakBuche2024: Initial release (V.0). Zenodo.  DOI - 10.5281/zenodo.12682869 (version V.1 will be released upon acceptance).
Contact details: Jürg W. Spaak (j.w.spaak@gmail.com) for code-related questions and Lisa Buche (buchel9844@gmail.com) for data-related questions. 

Statement of authorship:  

Abstract: 
Species interactions impact how species grow and which species persist. 
Many of these interactions, such as competition and predator-prey relationships, are well-studied.
However, there is less understanding of facilitation, particularly within the same species (i.e.,intra-specific facilitation).
The most commonly used theoretical framework for understanding species persistence, modern coexistence theory, does not apply to species with intraspecific facilitation.
We aim to reevaluate the existing methods of modern coexistence theory, which is based on the concept of niche and fitness differences, to include intra-specific facilitation.
Specifically, we reinterpret the traditional niche differences as how different inter-specific interactions are from intra-specific interactions and we reinterpret the traditional fitness differences as how the monoculutre equilibrium densities between the species compare. We apply the novel framework to an empirical community of four interacting annual plant species which can co-exist thanks to intraspecific facilitation.
The new interpretation is not only useful in the presence of intra-specific facilitation but also sheds new insights into why these can be used to interpret coexistence of natural communities.\\


Authorship of data: We want to thank Manuel Sevenello for the data collected in the Perenjori region. 

Authorship of code: R Code was written by Lisa Buche and Jürg W Spaak.

A shiny App is available. Open the webpage: https://lisa-buche.shinyapps.io/SpaakAndBuche/ or run the code the following code. First run the shiny toolbox.R to have all the functions in your environment. Then you can run the ShinyApp.R to pley. The figures match the ones from the manuscript. Top left: Diagram showing the strength and sign of species intra- and inter-specific interactions and their respective intrinsic performance. Top right: Abundance over time of species I when invading species j at its equilibrium. The three growth rates of i are shown. Bottom left: 3D visualisation of the growth rates, niche and fitness differences space. Bottom right: 2D visualisation of the niche and fitness differences.

Details of scripts: 

- **`plot_interpretation_NFD.py`**  
  Plots the interpretation of how niche and fitness differences can be separated into 3 regions and how this links to he linear interpolation of the growth rates.

- **`plot_allee_example.py`**  
  Computes niche nd fitness differences for the chosen community as shown in the main text.

- **`plot_NFD_map.py`**  
  Plots the NF-map with the corresponding lines partitioning the NF-map
