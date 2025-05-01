#1. required imports
import os
import sys
import warnings
import re
from Bio import SeqIO
from scipy import stats
import esm
import numpy as np
import pandas as pd
import torch
import shutil
from tqdm import tqdm
from transformers import T5Model, T5Tokenizer, T5EncoderModel, AutoTokenizer, AutoModelForSeq2SeqLM, T5ForConditionalGeneration
from transformers import EsmModel, EsmTokenizer
from transformers import BertModel, BertTokenizer
from transformers import AlbertModel, AlbertTokenizer
from transformers import XLNetModel, XLNetTokenizer

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu') 

#2. MODELS

#1)Ankh
def Ankh_initialize():

    ankh_model = T5ForConditionalGeneration.from_pretrained("/change/path/to/model", local_files_only=True)
    ankh_tokenizer = AutoTokenizer.from_pretrained("/change/path/to/model", local_files_only=True)
    ankh_model = ankh_model.eval()
    ankh_model = ankh_model.to(device)

    return ankh_model, ankh_tokenizer


def get_embs_Ankh(seq, ankh_model, tokenizer):
    
    inputs = tokenizer([seq], padding=True, return_tensors="pt")
    inputs = {key: value.to(device) for key, value in inputs.items()}

    with torch.no_grad(): 
        outputs = ankh_model.encoder(**inputs)
        last_hidden_states = outputs.last_hidden_state 
    sequence_embedding = last_hidden_states[0, :-1, :]
    return sequence_embedding


#2) PROT_T5
def ProtT5_initialize():

    ProtT5 = T5EncoderModel.from_pretrained("/change/path/to/model", local_files_only=True)
    ProtT5_tokenizer = T5Tokenizer.from_pretrained("/change/path/to/model", do_lower_case=False)
    ProtT5 = ProtT5.half()
    ProtT5 = ProtT5.eval()
    ProtT5 = ProtT5.to(device)

    return ProtT5, ProtT5_tokenizer

def get_embs_T5(seq, T5, ProtT5_tokenizer):

  seq = re.sub(r"[UZOB]", "X", seq)
  seq = " ".join(list(seq))
  inputs = ProtT5_tokenizer([seq], padding=True, return_tensors="pt")
  inputs = {key: value.to(device) for key, value in inputs.items()}

  with torch.no_grad(): 
    outputs = T5(**inputs) 
    last_hidden_states = outputs.last_hidden_state 
  
  sequence_embedding = last_hidden_states[0, :-1, :]
  return sequence_embedding


#3) ESM models (ESM1b, ESM2)
def Esm_initialize():
    
    esm2 = EsmModel.from_pretrained("/change/path/to/model", local_files_only=True)
    esm2_tokenizer = EsmTokenizer.from_pretrained("/change/path/to/model", do_lower_case=False)
    esm2 = esm2.half() 
    esm2 = esm2.eval()  
    esm2 = esm2 = esm2.to(device) 

    return esm2, esm2_tokenizer

def get_esm2_embedding(seq, esm2, esm2_tokenizer):

    inputs = esm2_tokenizer([seq], padding=True, return_tensors="pt")
    inputs = {key: value.to(device) for key, value in inputs.items()}

    with torch.no_grad(): 
        outputs = esm2(**inputs) 
        last_hidden_states = outputs.last_hidden_state 

    sequence_embedding = last_hidden_states[0, 1:-1, :]  

    return sequence_embedding


#4) ProtBert
def ProtBert_initialize():

    protBert = BertModel.from_pretrained("/change/path/to/model", local_files_only=True)
    protBert_tokenizer = BertTokenizer.from_pretrained("/change/path/to/model", do_lower_case=False)
    protBert = protBert.half()  
    protBert = protBert = protBert.eval() 
    protBert = protBert.to(device)

    return protBert, protBert_tokenizer


def get_protbert_embedding(seq, protBert, protBert_tokenizer):

    seq = " ".join(list(re.sub(r"[UZOB]", "X", seq)))
    inputs = protBert_tokenizer([seq], padding=True, return_tensors="pt") 
    inputs = {key: value.to(device) for key, value in inputs.items()}

    with torch.no_grad(): 
        outputs = protBert(**inputs)
        last_hidden_states = outputs.last_hidden_state 

    sequence_embedding = last_hidden_states[0, 1:-1, :] 

    return sequence_embedding


#5) ProtAlbert
def ProtAlbert_initialize():

    protAlbert = AlbertModel.from_pretrained("/change/path/to/model", local_files_only=True)
    protAlbert.half()  
    protAlbert = protAlbert.to(device)  
    protAlbert = protAlbert.eval()  
    protAlbert_tokenizer = AlbertTokenizer.from_pretrained("/change/path/to/model", do_lower_case=False)

    return protAlbert, protAlbert_tokenizer

def get_protAlbert_embedding(seq, protAlbert, protAlbert_tokenizer):

    seq = " ".join(list(re.sub(r"[UZOB]", "X", seq)))
    ids = protAlbert_tokenizer.batch_encode_plus([seq], add_special_tokens=True, padding="longest")
    input_ids = torch.tensor(ids['input_ids']).to(device)
    attention_mask = torch.tensor(ids['attention_mask']).to(device)
    with torch.no_grad():
        outputs = protAlbert(input_ids=input_ids, attention_mask=attention_mask)

    last_layer_repr = outputs.last_hidden_state  
    residue_embeddings = last_layer_repr.squeeze(0)[1:-1] 

    return residue_embeddings 


#6) XLNet
def XLNet_initialize():

    xlnet = XLNetModel.from_pretrained("/change/path/to/model", local_files_only=True)
    xlnet.half()  
    xlnet = xlnet.to(device)  
    xlnet = xlnet.eval()  
    xlnet_tokenizer = XLNetTokenizer.from_pretrained("/change/path/to/model", do_lower_case=False)

    return xlnet, xlnet_tokenizer

def get_xlnet_embedding(seq, xlnet, xlnet_tokenizer):

    seq = " ".join(list(re.sub(r"[UZOB]", "X", seq))) 

    ids = xlnet_tokenizer.batch_encode_plus([seq], add_special_tokens=True, padding="longest")
    input_ids = torch.tensor(ids['input_ids']).to(device)
    attention_mask = torch.tensor(ids['attention_mask']).to(device)
    with torch.no_grad():
        outputs = xlnet(input_ids=input_ids, attention_mask=attention_mask)

    last_layer_repr = outputs.last_hidden_state.squeeze(0) 

    num_residues = len(seq.replace(" ", ""))  
    residue_embeddings = last_layer_repr[:num_residues] 

    return residue_embeddings 




def main():

    seq1 = "AAACC"
    seq2 = "DDDFF"

    #ankh
    model, tokenizer = Ankh_initialize()
    emb1 = get_embs_Ankh(seq1, model, tokenizer).cpu().numpy()
    emb2 = get_embs_Ankh(seq2, model, tokenizer).cpu().numpy()
    print(emb1.shape)
    print(emb2.shape)
    print("\n")

    #protT5
    model, tokenizer = ProtT5_initialize()
    emb1 = get_embs_T5(seq1, model, tokenizer).cpu().numpy()
    emb2 = get_embs_T5(seq2, model, tokenizer).cpu().numpy()
    print(emb1.shape)
    print(emb2.shape)
    print("\n")

    #ESM2
    model, tokenizer = Esm_initialize()
    emb1 = get_esm2_embedding(seq1, model, tokenizer).cpu().numpy()
    emb2 = get_esm2_embedding(seq2, model, tokenizer).cpu().numpy()
    print(emb1.shape)
    print(emb2.shape)
    print("\n")

    #BERT
    model, tokenizer = ProtBert_initialize()
    emb1 = get_protbert_embedding(seq1, model, tokenizer).cpu().numpy()
    emb2 = get_protbert_embedding(seq2, model, tokenizer).cpu().numpy()
    print(emb1.shape)
    print(emb2.shape)
    print("\n")

    #ALBERT
    model, tokenizer = ProtAlbert_initialize()
    emb1 = get_protAlbert_embedding(seq1, model, tokenizer).cpu().numpy()
    emb2 = get_protAlbert_embedding(seq2, model, tokenizer).cpu().numpy()
    print(emb1.shape)
    print(emb2.shape)
    print("\n")

    #XLNET
    model, tokenizer = XLNet_initialize()
    emb1 = get_xlnet_embedding(seq1, model, tokenizer).cpu().numpy()
    emb2 = get_xlnet_embedding(seq2, model, tokenizer).cpu().numpy()
    print(emb1.shape)
    print(emb2.shape)
    print("\n")


if __name__ == "__main__":
    main()