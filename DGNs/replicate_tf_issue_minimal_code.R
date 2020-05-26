library(reticulate)

## For setup on the renku environment:  
reticulate::use_virtualenv("/opt/conda")

## Only for xenon6 setup:
reticulate::use_virtualenv("/tungstenfs/groups/gbioinfo/sharedSoft/virtualenvs/r-reticulate-keras-2.3.0-tensorflow-2.0.0-gpu")  # TF 2.00
#reticulate::use_virtualenv("/tungstenfs/groups/gbioinfo/sharedSoft/virtualenvs/r-reticulate-keras-2.2.5-tensorflow-1.14.0-gpu") # TF 1.14
Sys.setenv("CUDA_VISIBLE_DEVICES" = "1" ) # Define visible GPU devices. Make certain that the set device is not in use by checking the output of nvidia-smi



library(keras)
use_implementation("tensorflow")
library(tensorflow)
library(Matrix)
library(SingleCellExperiment)

K <- backend() # manual add-on
tf$version$VERSION


sce <- readRDS(gzcon(url("https://github.com/fmicompbio/adv_scrnaseq_2020/blob/master/DGNs/data/SCE_MammaryGland.rds?raw=true")))
assays(sce )[["lognorm"]] <- log2(sweep( counts(sce),2,sce$library_size ,FUN="/")*1e4 +1)

combined.df.filtered <- as.matrix(assays(sce )[["lognorm"]] )  

combined.annot.study <- sce@colData$study
combined.annot.ct <- sce@colData$cell.class
combined.annot <- paste(sce@colData$study,sce@colData$cell.class,sep="_")





####### Splitting in training and validation data, converting to array
set.seed(1)
holback.fraction=0.2
holdback.samples=sample(1:ncol(sce),round(holback.fraction*ncol(sce)) ) 


##### Training Data:
study_annot.train=combined.annot.study[-holdback.samples]
ct_annot.train=combined.annot.ct[-holdback.samples]

M=combined.df.filtered[,-holdback.samples]
sc_train_x=array(M, dim= c(dim(M)[1], prod(dim(M)[-1]))) # convert to an array
sc_train_x=t(sc_train_x)                                 #Need to transpose before passing to the model
rm(M)

##### Validation Data:
study_annot.test=combined.annot.study[holdback.samples]
ct_annot.test=combined.annot.ct[holdback.samples]

M=combined.df.filtered[,holdback.samples]
sc_test_x=array( M, dim= c(dim(M)[1], prod(dim(M)[-1]))) # convert to an array
sc_test_x=t(sc_test_x)                                   # Need to transpose before passing to the model
rm(M)
###################################################################



##### Sparse variational autoencoder with one hot encodding for auxiliary input fed after the latent layer
## Ensure compatibility with both TF2 nd TF1:
if (tensorflow::tf$executing_eagerly())
  tensorflow::tf$compat$v1$disable_eager_execution()



# Parameters --------------------------------------------------------------
neck <- 64L #256L
drop_rate=0.2 #less than 0.2
gene_dim <- ncol(sc_train_x)  #Number of features (genes) in your dataset
latent_dim <- neck
intermediate_dim <- 512L #
epsilon_std <- 0.8  #Standard deviation of the prior latent distribution (def=1)
var_prior <- epsilon_std**2
log_var_prior <- log(var_prior)
kl_weight=0.1   #Weight got the kulllback leibler divergence loss (def=1 ) 

# Encoder definition --------------------------------------------------------
x <- layer_input(shape = c(gene_dim),name="gene_input")
h <- layer_dense(x, intermediate_dim, activation = "elu") #softsign +elu +linear
h <- layer_dropout(h, rate = drop_rate)
h <- layer_dense(h,256,activation="elu")
h <- layer_dropout(h, rate = drop_rate)
h <- layer_dense(h,128, activation = "elu")
h <- layer_dropout(h, rate = drop_rate)
z_mean <- layer_dense(h, latent_dim)
z_log_var <- layer_dense(h, latent_dim)

#### Sampling from the latent space:
sampling <- function(arg){
  z_mean <- arg[, 1:(latent_dim)]
  z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]
  epsilon <- K$random_normal(
    shape = c(K$shape(z_mean)[[1]]), 
    mean=0.,
    stddev=epsilon_std
  )
  z_mean + K$exp(z_log_var/2)*epsilon
}

# note that "output_shape" isn't necessary with the TensorFlow backend
z <- layer_concatenate(list(z_mean, z_log_var)) %>% 
  layer_lambda(sampling)

# we instantiate the decoder separately so as to reuse it later
decoder_h <- keras_model_sequential()
decoder_h %>%
  layer_dense(units=128,activation="elu") %>% #Start with /4 when the extra optional layer is used. Otherwise /2
  layer_dropout(rate = drop_rate) %>%
  layer_dense(units=256,activation="elu") %>% 
  layer_dropout(rate = drop_rate) %>%
  layer_dense(intermediate_dim, activation = "elu") %>%  
  layer_dropout(rate = drop_rate)
decoder_mean <- layer_dense(units = gene_dim, activation = "relu")
h_decoded <- decoder_h(z)
x_decoded_mean <- decoder_mean(h_decoded)

# end-to-end autoencoder
vae <- keras_model(x, x_decoded_mean)

# encoder, from inputs to latent space
encoder <- keras_model(x, z_mean)

# generator, from latent space to reconstructed inputs
decoder_input <- layer_input(shape = latent_dim)
h_decoded_2 <- decoder_h(decoder_input)
x_decoded_mean_2 <- decoder_mean(h_decoded_2)
generator <- keras_model(decoder_input, x_decoded_mean_2)

vae_loss <- function(x, x_decoded_mean){
  reconstruction_loss  <-  loss_mean_squared_error(x, x_decoded_mean)
  kl_loss <- -kl_weight*0.5*K$mean(1 + z_log_var-log_var_prior - K$square(z_mean)/var_prior - K$exp(z_log_var)/var_prior, axis = -1L)  # More general formula
  reconstruction_loss + kl_loss
}



#compiling the defined model with metric = accuracy and optimiser adam.
opt <-  optimizer_adam(lr =0.001,amsgrad = TRUE)# 
vae %>% compile(
  loss = vae_loss,
  optimizer = opt
  #,
  #experimental_run_tf_function=FALSE
  #run_eagerly=TRUE
)








########## Fitting the model:
batch_size <- 512 
nepochs=20 # Run 150-250 epochs with a low lr 
burn_in_lr=0.00005 # 0.000001 x 50 For GTEx-TCGA
###### Training:
k_set_value(vae$optimizer$lr, burn_in_lr )

history2 <- vae %>% fit(
  x=sc_train_x,
  y=sc_train_x, 
  shuffle = TRUE, 
  epochs = nepochs,
  initial_epoch=41,
  batch_size = batch_size, 
  validation_data=list(sc_test_x,sc_test_x)
)





