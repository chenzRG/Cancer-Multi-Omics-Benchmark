# Subtype-GAN
Subtype-GAN is a cancer subtype framework based on multi-omics profiles. The input for the framework is copy number, mRNA, miRNA, DNA methylation, and other omics data. The output is the corresponding subclass label for each sample. SubtypeGAN is mainly divided into two components: 1. GAN-based feature extraction module is used to obtain latent variables from each type of omics data. 2. Consensus clustering and the Gaussian Mixture model determine the number of subtypes and the cluster label corresponding to each sample. 
```{r}
# It is recommended to use the conda command to configure the environment:
conda env create -f sgan.yml
conda activate sgan
# the input list of BRCA omics data set is input.txt. We can use the following command to finish the subtyping process: 
python SubtypeGAN.py -m SubtypeGAN -t BRCA -i ./input/BRCA.list   
# the Clustering output file are stored in ./results/BRCA.SubtypeGAN  
```

SubtypeGAN's Consensus clustering module is used as follows:  
```{r}
python SubtypeGAN.py -m cc -t BRCA -i ./input/BRCA.list
# record the corresponding class label for each sample and the output file is ./results/BRCA.cc
```

SubtypeGAN's performance comparison module (using the Variational Auto-Encoder method as an example) is used as follows: 
```{r}
python SubtypeGAN.py -m vae -t BRCA -i ./input/BRCA.list
# record the corresponding class label for each sample and the output file is ./results/BRCA.vae
```
SubtypeGAN is based on the Python program language. The generative adversarial network's implementation was based on the open-source library scikit-learn 0.22, Keras 2.2.4, and Tensorflow 1.14.0 (GPU version). After testing, this framework has been working correctly on Ubuntu Linux release 20.04. We used the NVIDIA TESLA T4 (16G) for the model training. When the GPU's memory is not enough to support the running of the tool, we suggest simplifying the encoder's network structure.
