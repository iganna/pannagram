# Extract Features from Alignments

After running the `pannagram` alignment pipeline, you can extract useful features from your data.  
The features must be executed in a specific order — not all can be run independently.

The figure below shows these dependencies:
<div style="width: 40%;">
<p align="left">
  <img src="images/features_scheme.png" style="width:100%; object-fit:cover;"/>
</p>
</div>

Extractable features are organized into four main groups, each depending on the previous steps:

**`I`. Alignment Summary.** These include identifying pangenome synteny and generating consensus sequences, which serve as the foundation for all subsequent analyses.

**`II`. SNPs and Diversity.** Focus on single-nucleotide polymorphisms (SNPs) and related metrics such as nucleotide diversity (π).

**`III`. SVs and Mobilome Families.** Focus on structural variants, their properties, and mobile element families.  
This step also enables comparison of SVs with other sequence features (e.g., transposon annotations).

**`IV`. Gene annotation.** Use gene annotations from different accessions and arrange them into common annotation groups through the pangenome coordinate system.

### Run All Steps in One Command
to execute the complete feature extraction workflow:
```
features -path_project ${PATH_PROJECT} -synteny -consensus -snp -snp_pi -sv -sv_family -sv_orf
```