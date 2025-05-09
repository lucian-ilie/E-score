

<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--

-->



<br />
<div align="center">

<h1 align="center"> E-score</h1>


  <p align="center">
    Aligning Protein Sequences Using Protein Language Model Vectors
    <br />
  
   
  </p>
</div>



<!-- ABOUT THE PROJECT -->
## About The Project
The E-score project enables computation of local, global-regular, and semi-global alignments between any two protein sequences using their embedding vectors derived from state-of-the-art pre-trained models. Unlike traditional methods that rely on fixed substitution scores (e.g., BLOSUM matrices), E-score uses cosine similarity between embedding vectors of amino acids, allowing for context-dependent scoring. The alignment algorithms are adaptations of classical dynamic programming approaches tailored for protein sequence analysis.
<br />
<br />
To help you get started, we provide a simple usage example in e_score_quick_start.py. For running multiple examples from an example file, see the full_test_example directory. Scripts for analyzing results and computing sequence distances are available in analyzing_results and distances_and_extracting_local_alignments.py, respectively.
<br />
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
The model we found to work best is Ankh. However, we also provide instructions for using other models such as ProtT5, ProtBert, ProtAlbert, ProtXLNet, ESM1b, and ESM2, all of which are assumed to be downloaded locally from Hugging Face. This setup is ideal for server environments where downloading models at runtime may pose security risks. That said, any properly installed embedding model can be easily integrated into the pipeline. See initalizing_models_and_getting_embs.py. 



## Parameters and Descriptions: e_score_quick_start.py

Given a pair of protein sequences, this file can align them locally, globally, or semi-globally. Within the main function: 


| Parameter | Description |
| protein_seq1       |      | 
| protein_seq2   |  | 
| gap_open   |       | 
| gap_extension |   | 
| alignment_type    |   | 
| int_shift  |   |
<br />
| Return Value | Description |
| aligned_seq1       |      | 
| aligned_seq2   |  | 
| max_score |  | 

## Parameters and Descriptions: full_test_example

Given an example.csv file, with the same columns as those in the CSVs from the example_data directory, full_test_example demonstrates how to run all examples sequentially and save the results in formats compatible with the analysis scripts we provide.

| Parameter | Description |
| example_id       |      | 
| parms   |  | 
| model and tokenizer  |       | 
| gap_open |   | 
| gap_extension    |   | 
| shift  |   |
<br />

Each alignment test is saved to its own text file, and a single distances.csv file is generated to collect the distances from all tests.




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
