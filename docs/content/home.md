
# Pannagram

Pannagram is a package for constructing Pan-genome alignments, analyzing Structural Variants (SVs) and Mobile Element Families, as well as translating annotations between genomes and more.  
It comprises two components: command-line interface modules and an R library for downstream analysis and visualization.  
Both components are installed together through a Conda environment.  
The code is available in the [GitHub repository](https://github.com/iganna/pannagram/).  

The philosophy of Pannagram is to create a user-friendly platform that facilitates various tools for complex analysis.
The structure of Pannagram is illustrated below:

<div style="width: 70%;">
<p align="center">
  <img src="images/pannagram_scheme.png" style="width:70%; object-fit:cover;"/>
</p>
</div>

## Citation

If you use Pannagram, please cite:

- **Pannagram: unbiased pangenome alignment and Mobilome calling**  
  *Anna A. Igolkina et al.*, *bioRxiv*, 2025. [**Link**](https://doi.org/10.1101/2025.02.07.637071)

To explore Pannagram applications, we recommend:

- **A comparison of 27 *Arabidopsis thaliana* genomes and the path toward an unbiased characterization of genetic polymorphism**  
  *Anna A. Igolkina et al.*, *Nature Genetics*, 2025. [**Link**](https://doi.org/10.1038/s41588-025-02293-0)

## Acknowledgements

**Development:**
- Anna A. Igolkina - Lead Developer and Project Initiator
- Alexander D. Bezlepsky - Assistant

**Testing:**
- Anna Igolkina: Lead Tester
- Anna Glushkevich: Testing the alignment on _A. lyrata_ genomes
- Elizaveta Grigoreva: Testing the alignment on _A. thaliana_ and _A. lyrata_ genomes
- Jilong Ma: Testing the SV-graph on spider genomes
- Alexander Bezlepsky: Testing the Pannagram's functionality on Rhizobial genomes
- Gregoire Bohl-Viallefond: Testing the annotation converter on _A. thaliana_ alignment

**Utilities:**
- Parallel Processing Tool: O. Tange (2018): GNU Parallel 2018, ISBN 9781387509881, DOI [https://doi.org/10.5281/zenodo.1146014](https://doi.org/10.5281/zenodo.1146014).
