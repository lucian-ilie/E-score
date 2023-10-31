

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
| gap_penalty   |  default = -0.25; Recommended Values: -1, -0.5, -0.4, -0.3, -0.2, -0.1    |
| gap_extension_penalty   |default = -0.01; Recommended Values: -0.2, -0.05, -0.04, -0.03, -0.02, -0.01      | 

The output file will be named "FastaFileName_ScoringType_AlignmentType_GapPenalty_GapExtensionPenalty_Alignment.txt".
### Usage Examples

```python
saving_add =  "/content/"
seqs_path = "TestOne.fasta"
scoring = "ProtT5" 
alignment_type = "Global-regular" 
gap_penalty = -0.25
gap_extension_penalty = -0.01

alignment_file_TXT(saving_add = saving_add , seqs_path = seqs_path, scoring = scoring, alignment_type = alignment_type,
                      gap_penalty = gap_penalty, gap_extension_penalty = gap_extension_penalty)
```

Output (TestOne_ProtT5_Global-regular_-0.25_-0.01_Alignment.txt):

```
Seq 1 
>gi|150273654|gb|EDN00782.1|
TGVDLGTAYIVLVVLDEENNPVACEKQAAQVLRDGVVVDYTGALRIVRELKEKLEARLGTELVNCAIAMPAGTESSVRTHQYVAEGAGFEVTEILDEPSAANAIYQIENGVVVDIGGGTTGLAMLKDGVVVQTEDEPTGGTHLSLVLAGNYHISFAEAEAIKQDYARHREILPVVRPVLEKMASIVKRYVSQSDVDTIYLCGGTCCLTGIEQVFEKVTGIHTVKPANPFLVTPTGIAM
Seq 2 
>gi|168988865|pdb|3B8A|X
LAIDLGGTNLRVVLVKLSGNHTFDTTQSKYKLPHDMRTTKHQEELWSFIADSLKDFMVEQELLNTKDTLPLGFTFSYPASQNKINEGILQRWTKGFDIPNVEGHDVVPLLQNEISKRELPIEIVALINDTVGTLIASYYTDPETKMGVIFGTGVNGAFYDVVSDIEKLEGKLADDIPSNSPMAINCEYGSFDNEHLVLPRTKYDVAVDEQSPRPGQQAFEKMTSGYYLGELLRLVLLELNEKGLMLKDQDLSKLKQPYIMDTSYPARIEDDPFENLEDTDDIFQKDFGVKTTLPERKLIRRLCELIGTRAARLAVCGIAAICQKRGYKTGHIAADGSVYNKYPGFKEAAAKGLRDIYGWTGDASKDPITIVPAEDGSGAGAAV

Alignment Type : Global-regular

Opening Gap Penalty : -0.25
Extension Gap Penalty : -0.01
Scoring System : ProtT5
Score : 71.37376693785201

Seq 1 : 1     TGVDLGTAY--IVLVVLDEENNP-VACEKQAAQVLRDGVVVDYTGALRIVRELKEKLEAR    57
                 DLG      VLV L           K                               
Seq 2 : 1     LAIDLGGTNLRVVLVKLSGNHTFDTTQSKYKL-PHDMRTTKHQEELWSFIADSLKDFMVE    59

Seq 1 : 58    LGT----ELVNCAIA--MPAGTE----------------------SSVRTHQYVAEG--A    87
                                PA                           V   Q        
Seq 2 : 60    QELLNTKDTLPLGFTFSYPASQNKINEGILQRWTKGFDIPNVEGHDVVPLLQNEISKREL   119

Seq 1 : 88    GFEVTEILDEPSA---ANAIYQIENGVVVDIGGGTTGLAMLKDGVV--------------   130
                E             A      E    V  G G  G                       
Seq 2 : 120   PIEIVALINDTVGTLIASYYTDPETKMGVIFGTGVNGAFYDVVSDIEKLEGKLADDIPSN   179

Seq 1 : 131   ------------------------------------VQTEDEPTGGTHLS-------LVL   147
                                                   Q     T G  L        L L
Seq 2 : 180   SPMAINCEYGSFDNEHLVLPRTKYDVAVDEQSPRPGQQAFEKMTSGYYLGELLRLVLLEL   239

Seq 1 : 148   AGN---------------YHISFAEAEAIKQD---------------------YARHREI   171
                                Y         I  D                           I
Seq 2 : 240   NEKGLMLKDQDLSKLKQPYIMDTSYPARIEDDPFENLEDTDDIFQKDFGVKTTLPERKLI   299

Seq 1 : 172   LPVVRPVLEKMA----SIVKRYVSQSDVDT--IYLCGG-TCCLTGIEQVFEK----VTG-   219
                         A                 T  I   G       G      K      G 
Seq 2 : 300   RRLCELIGTRAARLAVCGIAAICQKRGYKTGHIAADGSVYNKYPGFKEAAAKGLRDIYGW   359

Seq 1 : 220   --------IHTVKPANPFLVTPTGIAM   238
                      I  V           G A 
Seq 2 : 360   TGDASKDPITIVPAE---DGSGAGAAV   383

```

## Google Colab
We provided a Google colab containing the [E-score](https://colab.research.google.com/drive/18F_S2sthFDZMdcp3s77uodCMMN77aoFJ?usp=sharing) notebook.


<!-- CONTACT -->
## Contact

Sepehr Ashrafzadeh - sashra29@uwo.ca
<br />
Lucian Ilie - ilie@uwo.ca






<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
