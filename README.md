# Protein-Function-Prediction


The following codes are for extracting multiple features for protein function prediction and generating a dataset of Protein features 
to identify their function in terms of Gene ontologies assocaited with the molecular function domain.

The Dataset is prepared using the following steps.

1. Protein Data was fetched for different bacterial phyla from uniprot.org in two different formats.

.csv //comma separated file that consist of different selected columns available on Uniprot.org

.fasta //fixed format containing sequences

2. For each sequence in fasta file, consensus pattern of the motifs present in the sequence were retrieved from prosite.expasy.org/scanProsite using the BioPython package and Regular Expressions in python.

Prosite_Pattern.py yields 

sequence_pattern.csv  //file that contains sequences along with consensus patterns obtained by excluding the profiles from scan. This returns only the motifs in the proteins corresponding to which the patterns are extracted.

3. All the unique motifs were identified from sequence_pattern.csv file. 

unique_motifs.csv //file containing list of unique motifs retrieved as an output for step 2 (in sequence_pattern.csv)

4. Frequent Gene Ontology (Molecular Function) terms were extracted from uniprot file and then the samples were retrieved corresponding to those frequent Gene Ontology Terms.

modified_sample.csv //containing data downloaded from uniprot.org but only the samples which have the Gene Ontology (Molecular Function) of interest.


5. features_dataset.py is used to generate feature dataset with the following features:


A. Sequence based features

Protein length //Number of amino acids present in the protein sequence

AAC (Amino Acid Composition) // the fraction of each of the 20 standard amino acids present in the sequence

DPC (Dipeptide Composition) // the fraction of all the possible dipeptides (2 contiguous amino acids) present in the sequence

TPC (Tripeptide Composition) // the fraction of all the possible tripeptides (3 contiguous amino acids) present in the sequence

PAAC (Psedo amino acid composition) // collection of 20+𝜆 descriptor values to preserve the sequence order information of the protein sequence


B. Physicochemical based features


Molecular weight

Instability Index // numerical value showing the stability of the protein

Isoelectric point // pH at which the net charge on the protein is 0 

GRAVY (Grand average of hydropathy) // average of the hydropathy values (numerical indices corresponding to the hydrophilicity and hydrophobicity) associated with all the amino acids in a protein sequence

Extinction coeffcient // measure of the amount of light absorbed at a certain wavelength by the amino acids of the protein sequence

Secondary structure fraction // fraction of amino acids that appear to form helix, turn and sheets

GPPC (Grouped amino acid composition) // amino acids are categorically divided on the basis of their nature into 5 groups and fraction of each group is identified

Moran Autocorrealtion // numerical descriptors showing the correlation between the amino acids in terms of structural and physicochemical properties

CTD (Composition, Transition and Distribution) // the amino acid distribution pattern for different structural and physicochemical properties in a given sequence

Conjoint triad // amino acid are divided into 7 categories on the basis of properties influencing their interaction, followed by formation of all possible triads and calculating the frequency of each in the protein sequence



C. Subsequence based features


pattern count // All the unique consensus patterns obtained in step 3 are used as an attribute for the feature set being formed, with the count of every attribute in each protein sequence considered as a feature value


D. Annotation based features

Annotations regarding the subcellullar localisation, binding preference of proteins and presence of transmembrane helices are taken as a feature in the feature set

Subcellular localisation // 13 categories were prepared after analysing the unique subcellular localisation values in .tsv file, with each being treated as an attribute and its presence/absence denoted by binary value as a feature

Binding Preferences // Proteins can bind to metals, nucleotides or DNA/RNA, hence annotation regarding these (interaction denoted by a value of 1) is taken as a feature value

Transmembrane helices // annotation showing presence/absence of transmembrane helix in the form of binary values      


Final_Dataset.csv contains all these features for a protein as well as 1739 class labels of GO terms which were previously used.


# Datasets

Dataset for the 9 phyla were generated using features_dataset.py with 1,71,212 samples of protein in total.

The hypothetical Protein features file was generated the same way the final_dataset.csv is generated (removing the GO terms from the script) to make predictions. 

Complete Dataset for 9 phyla can also be found at: https://github.com/suraiyajabin/ProteinFunctionPredictionDataSet

# Training

The model is trained on full data set with 75% of 1,71,212 reviewed proteins of 9 Bacterial phyla and tested on 25% of 1,71,212 reviewed proteins 9 Bacterial phyla being considered but here due to size limitation, we show the smallest Train, Test, and Predict dataset.



# Predictions

1. Actual and Predicted Values of GO terms for the Test dataset along with uniprot ID, Entry name and Family of the Bacteria

2. Predicted values of GO terms for the Hypothetical Proteins.
