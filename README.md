# ProtEp
Processing of proteomic files for prolamin peptides analysis and celiac disease (CD) epitopes search for each protein fraction.

*Requisites:*
* python version 3.6 or higher
  * Library: biopython
    ```
    pip install biopython
    ```


# **Inputs:**
* Proteomic files with information of protein identification and characterization by qTOF mass spectrometry in xlsx format.
* Database with peptide sequence of: CD epitopes, innate immune response related p31-41 peptide, monoclonal antibody (moAb) recognition sites and immunoglobulin E (IgE) recognition sites. The sequence of deaminated CD epitopes can be included. The database must be in fasta format. An example is in epitopes_example.fasta. In addition, in the epitope name has to be included the category of the epitope separated by ";": e.g. g12_1;moAb. Categories: original (adaptive and innate immuno response in CD), moAb, IgE, deaminated1, deaminated2 (second deaminated variant of the original epitope sequence).

The protein fractions considered in the analysis are: alpha, omega and gamma-gliadins, High Molecular Weight (HMW) and Low Molecular Weight (LMW) glutenins, and Non Gluten Proteins (NGPs).


# **Option 1:**

Count the number of unique peptides annotated in each protein fraction per sample.

```
python unique_peptides_v2.py -p </path/to/ proteomic files directory> -o </path/to/ output directory>
```

*Output:*
* Table in tsv format.

# **Option 2:**
Count the number of epitopes, moAb recognition sites, ..., in unique peptides for each protein fraction in each sample. The output can be by epitope/moAb name or by epitope/moAb type: e.g. g12_1 moAb has g12_1 as name and g12 as type.
```
python epitopeFinder.py -p </path/to/ proteomic files directory> -o </path/to/ output directory> -e </path/to/epitope fasta file> -m <method: original | moAb | ...> -d <remove epitope duplications: yes | no> -mm <number of allowed mismatches to find epitopes: 0 | 1 | ...> -t <type of table with counts: name | type>
```
*Output:*
* Table in tsv format.
