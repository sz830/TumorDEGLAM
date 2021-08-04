library(Matrix)
library(RCTD)
library(DEGLAM)
library(doParallel)
library(ggplot2)
library(xlsx)
library(dplyr)
library(fields)
library(GSA)
library(stringr)

for(file in list.files("~/DEGLAM/R")){
  source(paste0("~/DEGLAM/R/",file), echo=TRUE)
}

# Load in spatialRNA data and Reference data
pwd = getwd()
datadir <- paste0(pwd,'/data/tumor','/')
resultsdir <- paste0(pwd,'/ResultsTumor','/')

# Load RCTD which contains cropped puck, also load full puck
myRCTD = readRDS(paste0(datadir,'RCTD_merged.rds')) # 21902 genes time start 2:30pm 640 genes done at 2:53, estimated finish time 13 hours from start...
cropped_puck = myRCTD@spatialRNA
myRCTD@originalSpatialRNA = cropped_puck
full_puck = readRDS(paste0(datadir,'puck.rds'))

## Examine SpatialRNA object (optional)
print(dim(cropped_puck@counts))
hist(log(cropped_puck@nUMI,2))
print(head(cropped_puck@coords)) 
barcodes <- colnames(cropped_puck@counts) 
plot_puck_continuous(cropped_puck, barcodes, cropped_puck@nUMI, ylimit = c(0,round(quantile(cropped_puck@nUMI,0.9))),
                     title ='plot of nUMI cropped_puck')
plot_puck_continuous(full_puck, colnames(full_puck@counts) , full_puck@nUMI, ylimit = c(0,round(quantile(full_puck@nUMI,0.9))),
                     title ='plot of nUMI full_puck')

# Create explanatory.variable
cell_types = myRCTD@cell_type_info$info[[2]]
target_type = cell_types[6] # "monocyte/DC"
doublet_df = myRCTD@results$results_df
weights_doublet=myRCTD@results$weights_doublet
coords = cropped_puck@coords

# Get a list of barcodes for cells of target_type
# Filter so we have cells in the cropped puck, that are "singlets" or "certain doublets" with first or second type being the target type
target_df = dplyr::filter(doublet_df, (rownames(doublet_df) %in% barcodes) & (first_type == target_type | second_type == target_type) & (spot_class != 'reject') & (spot_class !='doublet_uncertain'))
target_barcodes = rownames(target_df)

# Names are the barcodes, value is a score computed using euclidean distance from the cells of target_type
all_barcodes = barcodes # The cropped puck subset, use rownames(doublet_df) for all barcodes
explanatory.variable = c(rep(0,length(all_barcodes)))
names(explanatory.variable) = all_barcodes

# Calculate proximity score by summing the scores across all cells of target type for each cell in cropped_puck
# Individual scores between a cell and any target cell is calculated as n_i*exp(-d_i/c) 
# n_i is the weighted nUMI of the target cell; weighted by the proportion that the pixel is the target cell type. singlets are weighted as 1.0
# d_i is the distance between the current cell and target cell
# c is a rate constant around 30-50 microns

# Create a distance table between all pairs of cells. rdist is so fast there's no need to save this.
# fields::rdist treats rows as coordinates and computes all distances between placing them in a distance matrix.
dist_matrix = fields::rdist(as.matrix(cropped_puck@coords))
rownames(dist_matrix) = rownames(cropped_puck@coords)
colnames(dist_matrix) = rownames(cropped_puck@coords)
# Precompute the exponent component of the proximity score for all pairs of cells
c = 30/.65
exponent_mat = exp(-dist_matrix/c)

# Precompute the weighted nUMI values for all target cells
weighted_nUMIs = c(rep(0,length(target_barcodes)))
names(weighted_nUMIs) = target_barcodes
for(i in 1:length(weighted_nUMIs)) {
  barcode = target_barcodes[i]
  nUMI = cropped_puck@nUMI[barcode]
  
  spot_class = doublet_df[barcode,"spot_class"]
  first_type = doublet_df[barcode,"first_type"]
  second_type = doublet_df[barcode,"second_type"]

  weight = 0.0;
  if(spot_class == "singlet") {
    weight = if (first_type == target_type) 1.0 else 0.0;
  } else {
    weight = if (first_type == target_type) weights_doublet[barcode,1] else weights_doublet[barcode,2];
  }
  weighted_nUMI = nUMI * weight
  weighted_nUMIs[i] = weighted_nUMI
}

# Use the precomputed components above to compute explanatory.variable
for(i in 1:length(all_barcodes)) {
  barcode = all_barcodes[i]
  
  exp_dists = exponent_mat[barcode,target_barcodes]
  proximity_score = weighted_nUMIs %*% exp_dists
  explanatory.variable[i]=proximity_score
}

# Normalize explanatory.variable so the values span from 0 to 1.
normalize_ev = function(ev) {
  # Threshold values over the specific 85% percentile to be 1.
  percentile = quantile(explanatory.variable,.85)
  ev = ev - min(ev)
  ev = ev / percentile
  ev[ev>1] = 1
  return(ev)
}
explanatory.variable = normalize_ev(explanatory.variable)
hist(explanatory.variable)
plot_puck_continuous(cropped_puck, colnames(cropped_puck@counts) , explanatory.variable, # ylimit = c(0,round(quantile(explanatory.variable,0.9))),
                     title ='plot of explanatory.variable cropped_puck')

# run DE. Won't run on all genes, only the highly expressed ones in both spatialRNA and reference
# DEGLAM automatically checks for what cell types to use but since ev is distance from monocytes, it doesn't make sense to measure DE of monocytes so we have a custom list
cell_types = c("CAF","hepatocyte 2","vascular smooth mc")
myRCTD <- run.de.single(myRCTD, explanatory.variable, cell_types=cell_types) # outputs plots into pwd() + "/de_plots/"

# Save/Load and plot resulting RCTDde
# saveRDS(myRCTD,file.path(resultsdir,'myRCTDde.rds'))
myRCTD = readRDS(paste0(resultsdir,'myRCTDde.rds'))
make_all_de_plots(myRCTD, resultsdir)

# Check the state of the results, Should have iterations mostly less than 50, and a proportion of converged
# genes greater than 95%, small precision values are better
# Histogram of n.iter used by DEGLAM -- myRCTD@de_results$gene_fits$n.iter
hist(myRCTD@de_results$gene_fits$n.iter)
# Check for convergence -- look at myRCTD@de_results$gene_fits$con_val and myRCTD@de_results$gene_fits$precision_val
con_val = myRCTD@de_results$gene_fits$con_val
proportion_genes_converged = sum(con_val)/length(con_val)

# Precision vals should be low
precision_val = myRCTD@de_results$gene_fits$precision_val
hist(precision_val)

# Make sure the reject cells were discarded
dim(myRCTD@internal_vars_de$my_beta)
dim(myRCTD@spatialRNA@coords)
length(myRCTD@internal_vars_de$all_barc)


########################################################################
# Dylans Gene ontoogy code, adapted to CAF for tumor analysis
de_results = myRCTD@de_results
# Get a list of the overexpressed and under expressed genes in CAF cells
over_genes <- tolower(rownames(de_results$res_gene_list$CAF[de_results$res_gene_list$CAF$log_fc > 0,]))
under_genes <- tolower(rownames(de_results$res_gene_list$CAF[de_results$res_gene_list$CAF$log_fc < 0,]))

# Load hallmark gene sets and change the genes to all lowercase so it's easier to work with.
gene_sets = GSA.read.gmt(file.path(datadir,'hallmark_genesets.gmt'))
gene_set_names = gene_sets$geneset.names
gene_set_descriptions = gene_sets$geneset.descriptions
gene_sets = gene_sets$genesets
names(gene_sets)=gene_set_names
gene_sets = lapply(gene_sets, tolower)
n_sets = length(gene_sets)

# Optional: Check for intersections between hallmark sets.
 for(i in 1:n_sets) {
  for(j in 1:50){
    if(i!=j){
      print(intersect(gene_sets[[i]], gene_sets[[j]]))
    }
  }
}
# Conclusion, there's definitely intersections between hallmark sets

# Optional: How many of the gene set genes do we have spatial counts information on
genes_list <- tolower(rownames(myRCTD@originalSpatialRNA@counts))
coverage = c()
for(i in 1:n_sets){
  coverage[i] = length(intersect(gene_sets[[i]], genes_list))/length(gene_sets[[i]])
}
hist(coverage)

# Optional: Print the percent of genes in over and under genes that are in the hallmark gene sets
for(i in 1:n_sets) {
  print(paste0(i,"-----"))
  print(length(intersect(gene_sets[[i]], over_genes))/length(gene_sets[[i]]))
  print(length(intersect(gene_sets[[i]], under_genes))/length(gene_sets[[i]]))
}

# Use binomial tests to see if the ratio of genes over/under expressed in our data is significantly different
# when compared to the ratio of over/under expressed genes in the hallmark gene sets.
# If the ratio is significantly different, we can assume that gene sets pathway is significant somehow in
# how the cells are differentially expressed.
# Given a random gene in the differentially expressed genes, it should have a p_avg chance of being over expressed.
# So for n_o + n_u random genes, we expect p_avg amount of them to be overexpressed but we actually have n_o overexpressed.
# If we find too many more overexpressed it implies our overexpressed genes seem to fit this gene set well as in this gene set is enriched.
# two.sided makes the opposite true as well, so it also tests for underexpression
p_vals = numeric(n_sets)
p_avg = length(over_genes) / (length(over_genes) + length(under_genes))
n_o_vals = numeric(n_sets)
n_u_vals = numeric(n_sets)

max_n = 0; # Used later when generating a chart to know how many rows it needs to have, rows contain the over/under expressed genes
for(i in 1:n_sets) {
  print(i)
  n_o = length(intersect(gene_sets[[i]], over_genes))
  n_u = length(intersect(gene_sets[[i]], under_genes))
  n_o_vals[i] = n_o
  n_u_vals[i] = n_u
  # arguments are heads, flips, rate. it tells you the chance that you got n_e when the fair rate was p_avg
  p_vals[i] <- binom.test(n_o, n_o + n_u, p_avg, alternative="two.sided")$p.value
  max_n = max(max_n,n_o,n_u)
}

# Thorough analysis of pvalues taking indo account false discovery rate
# https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini%E2%80%93Hochberg_procedure
fdr = 0.1 
N_sig <- max(which(p_vals[order(p_vals)] < fdr * 1:n_sets / n_sets)) # which pvalues are less than .1 scaled from 0 t .1???
sig_list <- order(p_vals)[1:N_sig]
sig_pathways = gene_set_names[sig_list]
sig_list_df <- data.frame(sig_list, p_vals[sig_list], n_o_vals[sig_list], n_u_vals[sig_list])
colnames(sig_list_df) <- c('index', 'p_val', 'n_over', 'n_under')
sig_list_df$name <- gene_set_names[sig_list]
sig_list_df$desc <- gene_set_descriptions[sig_list]
gene_list_df <- data.frame()
# 2 columns per gene set
df <- matrix(nrow = max_n, ncol = 2*length(sig_list_df$name))
colnames(df) <- c(rbind(paste("UnderExpressed_" ,sig_list_df$name,sep=""),paste("OverExpressed_" ,sig_list_df$name,sep="")))
for(i in 1:N_sig) {
  gl <- intersect(gene_sets[[sig_list[i]]], under_genes)
  if(length(gl) > 0)
    df[1:length(gl),2*i-1] <- gl
  gl <- intersect(gene_sets[[sig_list[i]]], over_genes)
  if(length(gl) > 0)
    df[1:length(gl),2*i] <- gl # i+N_sig because the over expressed gene columns are shifted over by a constant from their under expressed column counterpart
}
# TODO: create a new folder for the gene enrichment analysis results
write.csv(df,file = file.path(resultsdir, 'gene_pathways.csv'))
write.csv(sig_list_df, file = file.path(resultsdir, 'sig_list_df.csv'))

# Looking at the list
pathways = substring(tolower(gene_set_names[sig_list]),10,)
#oxidative_phosphorylation - Process of ATP formation from O2

#epithelial_mesenchymal_transition - Process in which a regular epithelial cell changes into
#a mesenchymal cell. Increasing migratory ability, invasiveness, resistance to apoptosis, and
#increase in ECM component generation. Completion of this process causes it to degrade the 
#basement membrane it came from and migrate away from its origin location. 
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2689101/

#mitotic_spindle - dynamic molecular machine that distributes duplicated genome to daughters
#during mitsos

#complement - A part of the innate immune system. 
                   
# myc_targets_v1 - Involved in cell growth, apoptosis, metabolism. Homolog of retroviral
# v-myc oncogene, found activated in various tumors. 70,000 U.S. cancer deaths 
# per year are associated with changes in the c-myc gene or its expression
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC83860/

#apical_junction - Relevant in epithelial cell polarity

#uv_response_dn - Genes down-regulated in response to ultraviolet (UV) radiation
# Of the DE genes, 5 genes in this set are over expressed and 3 are under expressed
# So the liver sample has a lack of UV exposure? this one seems odd and hard to pinpoint anything on
# https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_UV_RESPONSE_DN

# mtorc1_signaling - Genes up-regulated through activation of mTORC1 complex.
# 18 under expressed from the DE genes, no over expressed. so a downregulation of the mTORC complex?
# mtorc1 stimulates cell throwth so this implies, less cell growth. ... 
# https://pubmed.ncbi.nlm.nih.gov/28411448/

#reactive_oxygen_species_pathway - Genes upregulated by ROS to fix the cell i assume
# A lot of both over and under, but more under.

# p53_pathway -The p53 pathway responds to stresses that can disrupt the fidelity
# of DNA replication and cell division. initiates a program of cell cycle arrest, cellular senescence or apoptosis.
# So like a suicide pathway when activated
# 8 under 5 over expressed.

# Next steps, check ratio of over to under expressed for each of these.
# look at the log_fc of these genes.
# If I can find any with a heavy balance toward over or under, due to either reason
# Best to investigate further.
# 
#####################################################################

# Interesting genes
# CAF
# F13a1 


# Extra Stuff



# Look for genes with high log_fc values and low P_values
#Original cell types
# 1  "hepatocyte 1" - Main cells of the liver, produce bile, stores glycogen, filter toxins from digestion
# 2  "tumor I"            
# 3  "tumor II"            
# 4  "LSEC"-Liver sinusoidal endothelial cells, line capillaries, antigen presenting, initiates immune response              
# 5  "Kupffer"-Phagocytic, lines sunisoids, breakdown RBC, first innate immune defense           
# 6  "monocyte/DC"-Monocytes are early macrophages, large macrophages, can turn into dendritic cells, they are APC and not really relevant in fighting cancer        
# 7  "HSC"-Hematopoitic stemm cell- turns into any and all cells in blood. mostly in bone marrow                
# 8  "CAF"- cancer associated fibroblasts. fibroblasts build the connective tissues forming the ECM                 
# 9  "hepatocyte 2"        
# 10 "hepatocyte 3"       
# 11 "T" - immune adaptive response, killer or helper, kills infected cells or cancers, reads only MHC1, helper T when meeting APC will stimulate other immune cells                   
# 12 "vascular smooth mc"- vascular smooth muscle. blod vessel walls
# 13 "interferon response"-interferon is a cytokine that interferes with viral replication mechanims. covid evades IFN-1, causes inflamation 
# 14 "B" - humoral adaptive response. in blood. kills pathogens in blood by eating them or using antibodies. antibodies bind to the antigens on pathogens marking them for death.  

# Get singlet positions of cells 
cell_types = myRCTD@cell_type_info$info[[2]]
plot_cell_type(7)
# 4 and 12 are supposed to be near blood vessels, but they don't fit eachothers locations.
# 10 and 12 kind of match? 8 is the inverse of 10
# 9 and 12 match a lot
# didn't expect a lot of 7 HSC since this is liver tissue but there's a good amount
# 5 and 12 kind of match, 5 is immune cells that line 4
# The tumor 2,3 is clearly on the top left. why is the caf 8 so clustered elsewhere?
# 12 vascular smooth mc is kinda around the tumor, because of the tumor?
# 9 and 10 are kind of away from the tumor and CAFs
plot_cell_type = function(cell_type_int) {
  cell_type = cell_types[cell_type_int]
  barcodes = rownames(doublet_df)[rownames(doublet_df) %in% rownames(full_puck@coords) & doublet_df$spot_class == "singlet" & doublet_df$first_type == cell_type]
  x_vals = c(); y_vals = c();
  for(i in 1:length(barcodes)) {
    barcode = barcodes[i]
    x_vals[i] = full_puck@coords[barcode,1]
    y_vals[i] = full_puck@coords[barcode,2]
  }
  cell_type_df = data.frame(x_vals,y_vals)
  ggplot(cell_type_df,aes(x=x_vals, y=y_vals)) + geom_point(aes(stroke=.000000001))
}

log_fc_threshold = .75;
p_val_threshold = .001;

# Columns: "Z_score"   "log_fc"    "conv"      "d"         "precision" "p_val" 
hep2_data = read.csv(paste0(resultsdir,'de_summary/hepatocyte 2.csv'))
rownames(hep2_data) = hep2_data[,1]; hep2_data[,1]=NULL;
interesting_hep2_genes = rownames(hep2_data)[hep2_data$log_fc > log_fc_threshold & hep2_data$p_val < p_val_threshold]
interesting_hep2_genes
print(length(interesting_hep2_genes))
print(dim(hep2_data)[1])

caf_data = read.csv(paste0(resultsdir,'de_summary/CAF.csv'))
rownames(caf_data) = caf_data[,1]; caf_data[,1]=NULL;
interesting_caf_genes = rownames(caf_data)[caf_data$log_fc > log_fc_threshold & caf_data$p_val < p_val_threshold]
interesting_caf_genes
print(length(interesting_caf_genes))
print(dim(caf_data)[1])

smooth_mc_data = read.csv(paste0(resultsdir,'de_summary/vascular smooth mc.csv'))
rownames(smooth_mc_data) = smooth_mc_data[,1]; smooth_mc_data[,1]=NULL;
interesting_smooth_mc_genes = rownames(smooth_mc_data)[smooth_mc_data$log_fc > log_fc_threshold & smooth_mc_data$p_val < p_val_threshold]
interesting_smooth_mc_genes
print(length(interesting_smooth_mc_genes))
print(dim(smooth_mc_data)[1])

##### Understanding results
gene_fits = myRCTD@de_results$gene_fits
names(gene_fits)
# "mean_val"      "con_val"       "ll_val"        "I_val"         "I_mat"         "n.iter"        "d_vals"        "base_val"     
# "all_vals"      "precision_val" "sigma_g" 