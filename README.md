

<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--

-->



<br />
<div align="center">

<h3 align="center"> E-score</h3>


  <p align="center">
    Aligning Protein Sequences Using Embedding Scores
    <br />
  
   
  </p>
</div>



<!-- ABOUT THE PROJECT -->
## About The Project
The E-score project focuses on computing Global-regular and Global-end-gap-free alignment between any two protein sequences using their embedding vectors computed by stat-of-art pre-trained models. Instead of a fixed score between two pairs of amino acids(like BLOSUM matrices), we use the cosine similarity between the embedding vectors of two amino acids and use it as the context-dependent score.

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

## Available Models
| Models | Embedding Dim | Pre-trained on
| :---         |     :---:     |  :---:     | 
| ProtT5   | 1024   | Uniref50 |
| ProtBert     | 1024       | Uniref100 |
| ProtAlbert  | 4096     | Uniref100 |
| ProtXLNet    | 1024      |  Uniref100  |
| ESM1b  | 1280     | Uniref50 |
| ESM2   | 1280      | Uniref50 |


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
| gap_penalty   |  default = -1; Recommended Values: -4, -3, -2, -1.5, -1, -0.5    |
| gap_extension_penalty   |default = -0.2; Recommended Values: -1, -0.8, -0.5, -0.3, -0.2, -0.1      | 

The output file will be named "FastaFileName_ScoringType_AlignmentType_GapPenalty_GapExtensionPenalty_Alignment.txt".
### Usage Examples

```python
saving_add =  "/content/"
seqs_path = "Test2.fasta"
scoring = "ProtT5" 
alignment_type = "Global-regular" 
gap_penalty = -1
gap_extension_penalty = -0.2

alignment_file_TXT(saving_add = saving_add , seqs_path = seqs_path, scoring = scoring, alignment_type = alignment_type,
                      gap_penalty = gap_penalty, gap_extension_penalty = gap_extension_penalty)
```

Output (Test2_ProtT5_Global-regular_-1_-0.2_Alignment.txt):

```
Seq 1 
>gi|464921|sp|P34981|TRFR_HUMAN
TILLVLIICGLGIVGNIMVVLVVMRTKHMRTPTNCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCSITAFTIERYIAICHPIKAQFLCTFSRAKKIIIFVWAFTSLYCMLWFFLLDLNISTYKDAIVISCGYKISRNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLNPIPSDPKENSKTWKNDSTHQNTNLNVNTSNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENWFLLFCRICIYLNSAINPVIYNLMS
Seq 2 
>gi|20455271|sp|Q9NSD7|R3R1_HUMAN
ISVVYWVVCALGLAGNLLVLYLMKSMQGWRKSSINLFVTNLALTDFQFVLTLPFWAVENALDFKWPFGKAMCKIVSMVTSMNMYASVFFLTAMSVTRYHSVASALKSHRTRGHGRGDCCGRSLGDSCCFSAKALCVWIWALAALASLPSAIFSTTVKVMGEELCLVRFPDKLLGRDRQFWLGLYHSQKVLLGFVLPLGIIILCYLLLVRFIADRRAAGTKGGAAVAGGRPTGASARRLSKVTKSVTIVVLSFFLCWLPNQALTTWSILIKFNAVPFSQEYFLCQVYAFPVSVCLAHSNSCLNPVLYCLVR

Alignment Type : Global-regular

Opening Gap Penalty : -1
Extension Gap Penalty : -0.2
Scoring System : ProtT5
Score : 141.1924964427947

Seq 1 : 1     TILLVLIICGLGIVGNIMVVLVVMRTK-HMRTPTNCYLVSLAVADLMVLVAAGLPNITDS    59
                      C LG  GN  V               N     LA  D               
Seq 2 : 1     ISVVYWVVCALGLAGNLLVLYLMKSMQGWRKSSINLFVTNLALTDFQFVLTLPFWAVENA    60

Seq 1 : 60    IYGSWVYGYVGCLCITYLQYLGINASSCSITAFTIERYIAICHPIKAQF-----------   108
                  W  G   C            AS    TA    RY       K              
Seq 2 : 61    LDFKWPFGKAMCKIVSMVTSMNMYASVFFLTAMSVTRYHSVASALKSHRTRGHGRGDCCG   120

Seq 1 : 109   ----LCTFSRAKKIIIFVWAFTSLYCMLWFFLLDLNISTYKDAIVISCGYKI----SRNY   160
                        AK      WA   L                          K         
Seq 2 : 121   RSLGDSCCFSAKALCVWIWALAALASLPSAIFSTTVKVMGEELCLVRFPDKLLGRDRQFW   180

Seq 1 : 161   YSPIYLMDFGVFYVVPMILATVLYGFIARILFLNPIPSDPKENSKTWKNDSTHQNTNLNV   220
                           V P       Y    R           K                   
Seq 2 : 181   LGLYHSQKVLLGFVLPLGIIILCYLLLVRFIADRR-AAGTKGG---------------AA   224

Seq 1 : 221   NTSNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVV---VNSFLSSPF------QE   271
                  R           VTK    VV  F L W P   L        F   PF        
Seq 2 : 225   VAGGRPTGASARRLSKVTKSVTIVVLSFFLCWLPNQALTTWSILIKFNAVPFSQEYFLCQ   284

Seq 1 : 272   NWFLLFCRICIYLNSAINPVIYNLMS   297
                           NS  NPV Y L  
Seq 2 : 285   VYAFPVSVCLAHSNSCLNPVLYCLVR   310
```

## Google Colab
We provided a Google colab containing the [E-score](https://colab.research.google.com/drive/1DJRt6zvfg6Qc-kaTjgRetl6Emmvc11gW?usp=sharing) notebook.


<!-- CONTACT -->
## Contact

Sepehr Ashrafzadeh - sashra29@uwo.ca
<br />
Lucian Ilie - ilie@uwo.ca






<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
