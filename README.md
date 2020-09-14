## DeepFun
DeepFun is a deep learning based approach to predict the functional impacts of sequence alterations at single nucleotide resolution. It integrates a total of 7879 epigenomic profiles, including 1548 DNase I accessibility data sets, 1536 histone mark data sets, and 4795 transcription factor binding data sets. DeepFun implements a deep convolutional neural network (CNN) to train a prediction model based on these profiles, augmented by computational evaluation. For any query single nucleotide polymorphism (SNP), DeepFun predicts values measuring the DNA binding effects of the reference allele and the alternative allele in each individual epigenomic data set. With the pre-trained model, this web-server provides an online service to quickly assess genetic variants.

## Documentation
DeepFun is a updated model based on Basset Framework (https://github.com/davek44/Basset). 

## Tutorials
### Compute SNP Activity Difference profiles.
The most useful capability of DeepFun right now is to annotate noncoding genomic variants with their influence on functional properties like accessibility.
To run this tutorial, you'll need to either download the pre-trained model from https://www.dropbox.com/s/rguytuztemctkf8/pretrained_model.th.gz or substitute your own file here: dnacnn_best_A.zip (DNA accessibility sites, histone marks and Transcription factor CTCF binding sites) or dnacnn_best_B.zip (Transcription factors binding sites).


model_file = '../data/models/pretrained_model.th'

### Execute an in silico saturated mutagenesis

## Webserver
https://bioinfo.uth.edu/deepfun/
