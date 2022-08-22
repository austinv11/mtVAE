source("helper_functions.R")


# Create color palettes ----------------------------------------------------

superpws <- bb_lookup$SUPER_PATHWAY %>% unique()
brewer_pal <- brewer.pal(length(superpws), name = "Set1")
brewer_pal[[1]] <- "#FFBE0B" # Lipid superpathway
brewer_pal[[6]] <- "#EFC3E6" # Nucleotide superpathway
brewer_pal[[length(brewer_pal)]] <- "#adb5bd" # For the Unknown superpathway

superpws_color <- 
  superpws %>% 
  sapply(function(pw) brewer_pal[which(superpws == pw)],
         USE.NAMES = TRUE,
         simplify = FALSE)

anno_colors <- 
  superpws_color %>% 
  # Add opacity alpha of 0.5 (following ggplot convention).
  # This is equivalent to the hex value 80:
  # https://gist.github.com/lopspower/03fb1cc0ac9f32ef38f4
  lapply(function(val) str_c(val, "FF")) %>%
  unlist() %>%
  .[order(names(.))] %>%
  list(SUPER_PATHWAY = .)

# Alluvial plot palette
str_wrap_width <- 20
alluvial_superpws <- 
  superpws %>% 
  str_wrap(width = str_wrap_width)

alluvial_color <- 
  alluvial_superpws %>% 
  sapply(function(pw) brewer_pal[which(alluvial_superpws == pw)],
         USE.NAMES = TRUE,
         simplify = FALSE)

alluvial_color$Lipid <- "#FFBE0B"
alluvial_color$Unknown <- "#adb5bd"



# Load cosine KPCA SAGE scores -----------------------------------------------------

cosine_sage_scores <- 
  list.files("results/sage_values/cosine", full.names = TRUE) %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# metabolites
cosine_met_sage_df <- 
  cosine_sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

cosine_met_sage_wide <- 
  cosine_met_sage_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Sub-pathways
cosine_sage_subpw_df <- 
  cosine_sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


cosine_sage_subpw_wide <- 
  cosine_sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

cosine_sage_subpw <- 
  cosine_sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")


# Super-pathways
cosine_sage_superpw_df <- 
  cosine_sage_scores %>% 
  filter(sage_type == "super")

cosine_sage_superpw <- 
  cosine_sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Load sigmoid SAGE scores -----------------------------------------------------

sigmoid_sage_scores <- 
  list.files("results/sage_values/sigmoid", full.names = TRUE) %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# metabolites
sigmoid_met_sage_df <- 
  sigmoid_sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

sigmoid_met_sage_wide <- 
  sigmoid_met_sage_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Sub-pathways
sigmoid_sage_subpw_df <- 
  sigmoid_sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


sigmoid_sage_subpw_wide <- 
  sigmoid_sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

sigmoid_sage_subpw <- 
  sigmoid_sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")


# Super-pathways
sigmoid_sage_superpw_df <- 
  sigmoid_sage_scores %>% 
  filter(sage_type == "super")

sigmoid_sage_superpw <- 
  sigmoid_sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())


# Load rbf SAGE scores -----------------------------------------------------

rbf_sage_scores <- 
  list.files("results/sage_values/rbf", full.names = TRUE) %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# metabolites
rbf_met_sage_df <- 
  rbf_sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

rbf_met_sage_wide <- 
  rbf_met_sage_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Sub-pathways
rbf_sage_subpw_df <- 
  rbf_sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


rbf_sage_subpw_wide <- 
  rbf_sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

rbf_sage_subpw <- 
  rbf_sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")


# Super-pathways
rbf_sage_superpw_df <- 
  rbf_sage_scores %>% 
  filter(sage_type == "super")

rbf_sage_superpw <- 
  rbf_sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Load poly SAGE scores -----------------------------------------------------

poly_sage_scores <- 
  list.files("results/sage_values/poly", full.names = TRUE) %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# metabolites
poly_met_sage_df <- 
  poly_sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

poly_met_sage_wide <- 
  poly_met_sage_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Sub-pathways
poly_sage_subpw_df <- 
  poly_sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


poly_sage_subpw_wide <- 
  poly_sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

poly_sage_subpw <- 
  poly_sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")


# Super-pathways
poly_sage_superpw_df <- 
  poly_sage_scores %>% 
  filter(sage_type == "super")

poly_sage_superpw <- 
  poly_sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())



# Create heatmaps -----------------------------------------------------


# Reassign dimension names, instead of the dim_xx syntax
colnames(sage_subpw) <- 1:18
colnames(sage_superpw) <- 1:18
colnames(sage_met) <- 1:18

colnames(cosine_sage_subpw) <- 1:18

# cosine sub-pathway SAGE dimension-normalized
cosine_subpw_dim_normed <- 
  cosine_sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# cosine sub-pathway SAGE subpathway-normalized
cosine_subpw_pw_normed <- 
  cosine_sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))

colnames(sigmoid_sage_subpw) <- 1:18

# sigmoid sub-pathway SAGE dimension-normalized
sigmoid_subpw_dim_normed <- 
  sigmoid_sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# sigmoid sub-pathway SAGE subpathway-normalized
sigmoid_subpw_pw_normed <- 
  sigmoid_sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))

colnames(rbf_sage_subpw) <- 1:18

# rbf sub-pathway SAGE dimension-normalized
rbf_subpw_dim_normed <- 
  rbf_sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# rbf sub-pathway SAGE subpathway-normalized
rbf_subpw_pw_normed <- 
  rbf_sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))

colnames(poly_sage_subpw) <- 1:18

# poly sub-pathway SAGE dimension-normalized
poly_subpw_dim_normed <- 
  poly_sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# poly sub-pathway SAGE subpathway-normalized
poly_subpw_pw_normed <- 
  poly_sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))

# Display plots -------------------------------------------------------

# cosine
plot.new()
cosine_subpw_dim_normed
plot.new()
cosine_subpw_pw_normed

#sigmoid
plot.new()
sigmoid_subpw_dim_normed
plot.new()
sigmoid_subpw_pw_normed

#rbf
plot.new()
rbf_subpw_dim_normed
plot.new()
rbf_subpw_pw_normed

#poly
plot.new()
poly_subpw_dim_normed
plot.new()
poly_subpw_pw_normed
