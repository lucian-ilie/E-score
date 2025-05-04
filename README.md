

<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--

-->



<br />
<div align="center">

<h1 align="center"> E-score</h1>


  <p align="center">
    Aligning Protein Sequences Using Embedding Scores
    <br />
  
   
  </p>
</div>



<!-- ABOUT THE PROJECT -->
## About The Project
The E-score project enables computation of local, global-regular, and semi-global alignments between any two protein sequences using their embedding vectors derived from state-of-the-art pre-trained models. Unlike traditional methods that rely on fixed substitution scores (e.g., BLOSUM matrices), E-score uses cosine similarity between embedding vectors of amino acids, allowing for context-dependent scoring. The alignment algorithms are adaptations of classical dynamic programming approaches tailored for protein sequence analysis.
<br />
To help you get started, we provide a simple usage example in e_score_quick_start.py. For running multiple examples from an example file within a given domain, see the full_test_example directory. Scripts for analyzing results and computing sequence distances are available in analyzing_results and distances_and_extracting_local_alignments.py, respectively.
<br />
Global alignment examples can be constructed from any multiple sequence alignment (MSA) by selecting sequence pairs within the same MSA. Creating local examples requires more care; for convenience, we provide all of our curated local alignment examples, described in the Protein Embeddings and Local Alignments paper, in the example_data directory.

## System Requirements
Recommended Python Version: 3.10
<br />
Recommended RAM: 24GB
<br />
Each of the models needs about 8-12GB of RAM and as the length of the sequences goes up, the needed RAM increases too. 

## Installation
You can install all of the needed packages using requirement.txt. Python virtual environment is recommended.

```
pip3 install -r requirement.txt
```

## Protein Language Models
The default model we found to work best is Ankh. However, we also provide instructions for using other models such as ProtT5, ProtBert, ProtAlbert, ProtXLNet, ESM1b, and ESM2—all of which are assumed to be downloaded locally from Hugging Face. This setup is ideal for server environments where downloading models at runtime may pose security risks. That said, any properly installed embedding model can be easily integrated into the pipeline. See initalizing_models_and_getting_embs.py. 



<!-- USAGE EXAMPLES -->
## Usage
The program is designed to get a fasta file as an input and compute the alignment based on chosen parameters.
<br />
By calling the alignment_file_TXT function and passing the needed parameters (which are described below) the output would be a text file containing the computed alignment, its score, and the alignment visualization.

### Parameters and Descriptions

| Parameter | Description |
| :---         |     :---:     | 
| saving_add   | the path to the directory for the output   | 
| seqs_path    | the path of the directory for the FASTA file containing two protein sequences       | 
| scoring_type  | the embedding method used to produce the embedding vectors; allowed values are: ProtT5, ESM2, ProtBert, ProtAlbert, ESM1b, ProtXLNet     | 
| alignment_type    | allowed values are Global-regular or Global-end-gap-free     | 
| gap_penalty   |  default = -0.25; Recommended Values: -1, -0.5, -0.4, -0.3, -0.2, -0.1    |
| gap_extension_penalty   |default = -0.01; Recommended Values: -0.2, -0.05, -0.04, -0.03, -0.02, -0.01      | 

The output file will be named "FastaFileName_ScoringType_AlignmentType_GapPenalty_GapExtensionPenalty_Alignment.txt".


## E-Score Web Server
We provided a web server containing the [E-score](https://e-score.csd.uwo.ca) code.


<!-- CONTACT -->
## Contact

<br />
Lucian Ilie - ilie@uwo.ca
<br />
Julia Malec - jmalec@uwo.ca
<br />
Sepehr Ashrafzadeh - sashra29@uwo.ca


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
