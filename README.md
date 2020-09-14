## DeepFun
&#8194;&#8194;DeepFun is a deep learning based approach to predict the functional impacts of sequence alterations at single nucleotide resolution. It integrates a total of 7879 epigenomic profiles, including 1548 DNase I accessibility data sets, 1536 histone mark data sets, and 4795 transcription factor binding data sets. DeepFun implements a deep convolutional neural network (CNN) to train a prediction model based on these profiles, augmented by computational evaluation. For any query single nucleotide polymorphism (SNP), DeepFun predicts values measuring the DNA binding effects of the reference allele and the alternative allele in each individual epigenomic data set. With the pre-trained model, this web-server provides an online service to quickly assess genetic variants.

## Webserver
&#8194;&#8194;We present a user-friendly web server, DeepFun, available at https://bioinfo.uth.edu/deepfun/.

## Local analysis
&#8194;&#8194;DeepFun is an updated model based on Basset Framework (https://github.com/davek44/Basset) based on Torch7 and Python. If you want to conduct local analysis for large mounts of data, please follow Basset tutorials install dependencies packages at first. In addition, you'll need to download the pre-trained model from github here: dnacnn_best_A.zip (DNA accessibility sites, histone marks and Transcription factor CTCF binding sites) or dnacnn_best_B.zip (Transcription factors binding sites).  

## Tutorials
### Compute SNP Activity Difference profiles.
&#8194;&#8194;The most useful capability of DeepFun right now is to annotate noncoding genomic variants with their influence on functional properties like accessibility. More details can follow https://github.com/davek44/Basset/blob/master/tutorials/sad.ipynb.

### Execute an in silico saturated mutagenesis
&#8194;&#8194;Saturated mutagenesis is a powerful tool both for dissecting a specific sequence of interest and understanding what the model learned. More details can follow https://github.com/davek44/Basset/blob/master/tutorials/sat_mut.ipynb.

## Read more about the method in the manuscript here:
Pei G, Hu R, Dai Y, Manuel A M, Zhao Z, Jia P. Predicting regulatory variants using a dense epigenomic mapped CNN model elucidated the molecular basis of trait-tissue associations. Nucleic Acids Research. In submission.
