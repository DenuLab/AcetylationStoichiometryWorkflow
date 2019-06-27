# Load Kegg ---------------------------------------------------------------

load_kegg <- function(fname) {
  con <- file(fname, open='r')
  kegg_map <- list()
  i <- 1
  while(length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line <- unlist(strsplit(line, split = "\t"))
    pathway <- line[1:2]
    genes <- line[-(1:2)]
    kegg_map[i] <- list(genes)
    names(kegg_map)[i] <- pathway[1]
    i <- i + 1
  }
  close(con)
  return(kegg_map)
}

# Load Fasta --------------------------------------------------------------

load_fasta <- function(fname) {
  fasta <- read.fasta(fname, 
                      seqtype="AA",set.attributes=FALSE)
  
  # Rename to just accession number.
  names(fasta) <- lapply(names(fasta), 
                         function(x) {unlist(strsplit(x,"[|]"))[2]})
  return(fasta)
}

# Trim KEGG ---------------------------------------------------------------

trim_kegg <- function(kegg, proteins) {
  # Construct background by removing proteins not detected from the kegg map
  for(i in 1:length(kegg)) {
    kegg[i] <- list(intersect(kegg[[i]], proteins))
    names(kegg)[i] <- names(kegg[i])
  }
  
  # remove kegg pathways for which no proteins were detected
  kegg_group_sizes <- sapply(kegg, length)
  kegg <- kegg[which(kegg_group_sizes != 0)]
  return(kegg)
}

# Count pathway k-sites ---------------------------------------------------

lysine_sites_per_pathway <- function(kegg, fasta) {
  ptwys <- list()
  i <- 1
  for(pathway in kegg) {
    kcnt <- 0L
    for(protein in pathway) {
      seq <- fasta[[protein]]
      kcnt <- kcnt + length(seq[which(seq =="K")])
    }
    ptwys[[i]] <- list("ptns"=pathway,"kcnt"=kcnt)
    names(ptwys)[i] <- names(kegg)[i]
    i <- i + 1L
  }
  
  sites_per_pathway <- 
    unlist(lapply(ptwys,function(x) { x[[2]] }))
  return(sites_per_pathway)
}


# z-calculation -----------------------------------------------------------

calc_z <- function(vals) {
  mean_vals <- mean(vals)
  sd_vals <- sd(vals)
  zlist <- unlist(lapply(vals, function(x) {
                  (x - mean_vals) / sd_vals
                }))
  return(zlist)
}

# Count pathway k-sites ---------------------------------------------------

detected_sites_per_pathway <- function(kegg, stoich) {
  
  sites_detected_per_protein <- stoich %>% 
    select(PG.ProteinGroups, k_site) %>%
    unique() %>% 
    group_by(PG.ProteinGroups) %>%
    summarize(num_sites = length(k_site))  
  
  detected_site_count <- unlist(lapply(names(kegg), function(x){
                  proteins <- kegg[[x]]
                  sites <- sites_detected_per_protein[
                    sites_detected_per_protein$PG.ProteinGroups %in% proteins,
                    "num_sites"]
                  sum(sites$num_sites)
                }))
  
  return(detected_site_count)
}


# Sum per pathway ---------------------------------------------------------

sum_per_pathway <- function(kegg, stoich) {
  
  stoich_sums <- unlist(lapply(names(kegg),function(x) {
    proteins <- kegg[[x]]
    stoichs <- stoich[stoich$PG.ProteinGroups %in% proteins,"stoich_mean"]
    
    # we'll start by summing up all sites in all pathways points as
    # a first pass
    sum(stoichs)
    
  }))
  
  return(stoich_sums)
}


# Mean per pathway --------------------------------------------------------

mean_per_pathway <- function(kegg, stoich) {
  
  stoich_means <- unlist(lapply(names(kegg),function(x) {
    proteins <- kegg[[x]]
    stoichs <- stoich[stoich$PG.ProteinGroups %in% proteins,"stoich_mean"]
    
    mean(stoichs)
    
  }))
  
  return(stoich_means)
}