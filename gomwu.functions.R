plot_gomwu_trees <- function(
  n,
  gomwu_plot_list,
  files_to_plot,
  goDivision,
  pcut = 0.05,
  hcut = 0.9,
  plot_tree = TRUE,
  tree_width = 5,
  gomwu_dir,
  is_WGCNA = FALSE,
  # Plot saving parameters
  width = 7,
  height = 10,
  dpi = 300,
  filetype = "png",
  name = NULL
) {
  # --- Extract relevant data from the provided lists ---
  plot_data   <- gomwu_plot_list[[n]]
  good_genes  <- plot_data$goods
  hcl         <- plot_data$hcl
  inFile      <- files_to_plot[n]
  passing_genes <- plot_data$passing_genes
  
  # --- Derive a default output name if not provided ---
  outNameBase <- if (is.null(name)) {
    paste(gsub("\\.csv$", "", basename(inFile)), goDivision, sep = "_")
  } else {
    name
  }
  
  # --- Call the function to get representative GOs (assumes it’s defined elsewhere) ---
  representative_gos <- extract_representative_GOs(
    results    = plot_data,
    goDivision = goDivision,
    input      = inFile,
    gomwu_dir  = gomwu_dir,
    pcut       = pcut,
    hcut       = hcut,
    plot_tree  = plot_tree
  )
  
  # --- Convert hierarchical cluster (HCL) object to a 'phylo' tree ---
  phy_goods <- as.phylo(hcl)
  
  # --- Prepare metadata for plotting ---
  goods_mod <- good_genes %>%
    dplyr::select(-label) %>%      # remove existing 'label' if present
    tibble::rownames_to_column(var = "label") %>%
    dplyr::mutate(
      fontface = dplyr::case_when(
        sig_cat == "p < 0.01" ~ "bold",
        sig_cat == "p < 0.05" ~ "plain",
        sig_cat == "p < 0.1"  ~ "italic"
      )
    )
  
  # --- Plot: Full (Original) GO Tree ---
  if (is_WGCNA) {
    p <- ggtree(phy_goods) %<+% goods_mod +
      geom_tiplab(aes(fontface = fontface, size = sig_cat), color = "black") +
      scale_size_manual(
        name   = "Significance",
        values = c("p < 0.01" = 4, "p < 0.05" = 3, "p < 0.1" = 3)
      ) +
      coord_cartesian(xlim = c(0, tree_width), clip = "off") +
      theme(legend.position = "right")
  } else {
    p <- ggtree(phy_goods) %<+% goods_mod +
      geom_tiplab(aes(color = direction_factor, fontface = fontface, size = sig_cat)) +
      scale_size_manual(
        name   = "Significance",
        values = c("p < 0.01" = 4, "p < 0.05" = 3, "p < 0.1" = 3)
      ) +
      scale_color_manual(
        name   = "Direction",
        values = c("down" = "blue", "up" = "red")
      ) +
      coord_cartesian(xlim = c(0, tree_width), clip = "off") +
      theme(legend.position = "right")
  }
  
  # Save the full tree if requested
  if (plot_tree) {
    ggsave(
      filename = file.path(gomwu_dir, paste0(outNameBase, "_full.", filetype)),
      plot     = p,
      width    = width,
      height   = height,
      dpi      = dpi
    )
  }
  
  # --- Plot: Reduced (Pruned) GO Tree ---
  finalLabels <- representative_gos %>%
    dplyr::left_join(good_genes, by = c("term" = "original_term")) %>%
    dplyr::pull(label)
  
  phy_goods_pruned <- drop.tip(phy_goods, setdiff(phy_goods$tip.label, finalLabels))
  goods_mod_pruned <- goods_mod %>%
    dplyr::filter(label %in% finalLabels)
  
  if (is_WGCNA) {
    p_pruned <- ggtree(phy_goods_pruned) %<+% goods_mod_pruned +
      geom_tiplab(aes(fontface = fontface, size = sig_cat), color = "black") +
      scale_size_manual(
        name   = "Significance",
        values = c("p < 0.01" = 4, "p < 0.05" = 3, "p < 0.1" = 3)
      ) +
      coord_cartesian(xlim = c(0, tree_width), clip = "off") +
      theme(legend.position = "right")
  } else {
    p_pruned <- ggtree(phy_goods_pruned) %<+% goods_mod_pruned +
      geom_tiplab(aes(color = direction_factor, fontface = fontface, size = sig_cat)) +
      scale_size_manual(
        name   = "Significance",
        values = c("p < 0.01" = 4, "p < 0.05" = 3, "p < 0.1" = 3)
      ) +
      scale_color_manual(
        name   = "Direction",
        values = c("down" = "blue", "up" = "red")
      ) +
      coord_cartesian(xlim = c(0, tree_width), clip = "off") +
      theme(legend.position = "right")
  }
  
  # Save the reduced tree if requested
  if (plot_tree) {
    ggsave(
      filename = file.path(gomwu_dir, paste0(outNameBase, "_reduced.", filetype)),
      plot     = p_pruned,
      width    = width,
      height   = height,
      dpi      = dpi
    )
  }
  
  # --- Print the plots to the current device ---
  print(p)
  print(p_pruned)
  
  # --- Return additional objects invisibly ---
  invisible(list(
    representative_gos = representative_gos,
    passing_genes      = passing_genes,
    plot_data          = plot_data,
    good_genes         = good_genes,
    hcl                = hcl,
    goods_mod_pruned   = goods_mod_pruned,
    finalLabels        = finalLabels,
    original_tree      = p,
    reduced_tree       = p_pruned
  ))
}

extract_representative_GOs <- function(results, goDivision, input, gomwu_dir, 
                                       pcut = 1e-2, hcut = 0.9, plot_tree = FALSE) {
  message("\n--- Extracting Representative GO Terms ---\n")
  
  # Set working directory to the GO-MWU directory
  if (!dir.exists(gomwu_dir)) {
    stop("Error: The specified GO-MWU directory does not exist.")
  }
  setwd(gomwu_dir)
  message("Working directory changed to: ", getwd())
  
  # Check that results contains the required 'goods' component and that it has a 'rawp' column
  if (!("goods" %in% names(results))) {
    stop("Error: The results object must contain a 'goods' element.")
  }
  if (!("rawp" %in% colnames(results$goods))) {
    stop("Error: The 'goods' element must contain a 'rawp' column.")
  }
  
  # Initialize annotation storage
  annots <- c()
  num_terms <- nrow(results$goods)
  
  if (num_terms > 1) {
    if ("hcl" %in% names(results)) {
      # Check if the hcl object has enough leaves to form a dendrogram
      if (length(results$hcl$order) < 3) {
        message("Too few terms in hcl for clustering/plotting; selecting all terms passing significance filter.")
        significant_indices <- which(results$goods$rawp <= pcut)
        if (length(significant_indices) > 0) {
          annots <- sub("\\d+/\\d+ ", "", rownames(results$goods)[significant_indices])
        }
      } else {
        # Valid hcl available: cut the tree
        ct <- tryCatch({
          cutree(results$hcl, h = hcut)
        }, error = function(e) {
          message("Error cutting dendrogram for ", input, ": ", e$message)
          return(NULL)
        })
        
        if (is.null(ct)) {
          message("Dendrogram issues for ", input, "; selecting all terms passing significance filter.")
          significant_indices <- which(results$goods$rawp <= pcut)
          if (length(significant_indices) > 0) {
            annots <- sub("\\d+/\\d+ ", "", rownames(results$goods)[significant_indices])
          }
        } else {
          # Optionally plot the dendrogram if requested
          if (plot_tree) {
            plot_title <- gsub("\\.csv$", "", input)
            if (length(results$hcl$order) >= 3) {
              tryCatch({
                plot(results$hcl, cex = 0.6, main = plot_title)
                abline(h = hcut, col = "red")
              }, error = function(e) {
                message("Error plotting dendrogram for ", input, ": ", e$message)
              })
            } else {
              message("Not enough elements to plot dendrogram for ", input, "; skipping plot.")
            }
          }
          
          # Process clusters from the cut tree
          for (ci in unique(ct)) {
            message("Processing Cluster: ", ci)
            rn <- names(ct)[ct == ci]  # Terms in this cluster
            obs <- grep("obsolete", rn)  # Exclude obsolete terms
            if (length(obs) > 0) rn <- rn[-obs]
            if (length(rn) == 0) next  # Skip empty clusters
            
            rr <- results$goods[rn, , drop = FALSE]
            bestrr <- rr[which(rr$rawp == min(rr$rawp, na.rm = TRUE)), , drop = FALSE]
            best <- 1
            
            # If there are ties, break them using a frequency-based selection
            if (nrow(bestrr) > 1) {
              nns <- sub(" .+", "", row.names(bestrr))
              fr <- numeric(0)
              for (i in seq_along(nns)) {
                fr <- c(fr, eval(parse(text = nns[i])))
              }
              best <- which(fr == max(fr))[1]
            }
            
            if (!is.na(bestrr$rawp[best]) && bestrr$rawp[best] <= pcut) {
              annots <- c(annots, sub("\\d+/\\d+ ", "", row.names(bestrr)[best]))
            }
          }
        }
      }
    } else {
      # No hcl element present: fall back to selecting overall best GO term
      message("No hcl element present for ", input, "; selecting overall best GO term.")
      best_index <- which.min(results$goods$rawp)
      if (results$goods$rawp[best_index] <= pcut) {
        annots <- sub("\\d+/\\d+ ", "", rownames(results$goods)[best_index])
      }
    }
  } else if (num_terms == 1) {
    # Single-term case: process without clustering
    message("Only one term is differentially regulated; skipping tree clustering for ", input, ".")
    single_rawp <- results$goods[1, "rawp"]
    if (!is.na(single_rawp) && length(single_rawp) > 0 && single_rawp <= pcut) {
      term_name <- rownames(results$goods)[1]
      if (is.null(term_name) || term_name == "") term_name <- "term1"
      annots <- sub("\\d+/\\d+ ", "", term_name)
    }
  } else {
    stop("Error: No terms available in the 'goods' element.")
  }
  
  # Load the MWU file and filter the results based on selected annotations
  mwus_file <- paste("MWU", goDivision, input, sep = "_")
  if (!file.exists(mwus_file)) {
    stop("Error: MWU file not found in the specified directory: ", mwus_file)
  }
  
  mwus <- read.table(mwus_file, header = TRUE)
  bestGOs <- mwus[mwus$name %in% annots, ]
  
  message("\nRepresentative GO extraction complete for ", input, ".\n")
  return(bestGOs)
}
#-----
clusteringGOs <- function(goAnnotations, goDivision, clusterCutHeight, uniqueID = "") {
  
  inname <- paste0("dissim0_", goDivision, "_", goAnnotations)
  if (uniqueID != "") {
    inname <- paste0(inname, "_", uniqueID)
  }
  outname=paste("cl_",inname,sep="")
  	if (!file.exists(outname)) {
  		diss=read.table(inname,sep="\t",header=T,check.names=F)
  		row.names(diss)=names(diss)
  		hc=hclust(as.dist(diss),method="complete")
  		cc=cutree(hc,h=clusterCutHeight)
  		ccnames=names(cc)
  		if (clusterCutHeight==0) {
  		 print("no clustering of GO terms will be done")
  		  cc=seq(1:length(cc))
  		  names(cc)=ccnames
  		}
  		write.csv(cc,file=outname,quote=F)
  #		write.csv(cc,file=paste(outname,"bkp",sep="."),quote=F)
  
  	}
}

#---------------


gomwuStats=function(input,goDatabase,goAnnotations, goDivision, Module=FALSE, Alternative="t", adjust.multcomp="BH", clusterCutHeight=0.25,largest=0.1,smallest=5,perlPath="perl", shuffle.reps=20, uniqueID = ""){
  
  extraOptions=paste("largest=",largest," smallest=",smallest," cutHeight=",clusterCutHeight, " uniqueID=", uniqueID, sep="") # Add uniqueID here
  if (Module==TRUE) { adjust.multcomp="shuffle" }
  command <- paste(perlPath,"./gomwu_a.pl",goDatabase,goAnnotations,input,goDivision,extraOptions) # uniqueID added via extraOptions
  system(command)
  
  clusteringGOs(goAnnotations,goDivision,clusterCutHeight, uniqueID) # Pass uniqueID to clusteringGOs
  
  command <- paste(perlPath,"./gomwu_b.pl",goAnnotations,input,goDivision, uniqueID) #Pass any extra options if necessary
  system(command)
  
  
  inname=paste(goDivision,"_",input,sep="")
  rsq=read.table(inname,sep="\t",header=T)
  rsq$term=as.factor(rsq$term)
  
  mwut.t=TRUE
  if (length(levels(as.factor(rsq$value)))==2) {
    cat("Binary classification detected; will perform Fisher's test\n");
    mwut.t=F
    rr=fisherTest(rsq)
  } else {
    if (Module==TRUE) {
      rsq.f=rsq
      rsq.f$value=as.numeric(rsq.f$value>0)
      rf=fisherTest(rsq.f)
      rsq.m=rsq[rsq$value>0,]
      rsq.m$term=factor(rsq.m$term,levels=unique(rsq.m$term))
      rrm=mwuTest(rsq.m,"g")
      rr0=rf[rf$term %in% rrm$term,]
      rr1=rf[!(rf$term %in% rrm$term),]
      rr0=rr0[order(rr0$term),]
      rrm=rrm[order(rrm$term),]
      rr0$pval=rr0$pval*rrm$pval
      rr=rbind(rr0,rr1)
    } else {
      cat("Continuous measure of interest: will perform MWU test\n");
      rr=mwuTest(rsq,Alternative)
    }
  }
  
  if (adjust.multcomp=="shuffle"){
    message("shuffling values to calculate FDR, ",shuffle.reps," reps")
    reps=shuffle.reps
    spv=c()
    for (i in 1:reps) {
      message("replicate ",i)
      rsqq=rsq
      rsqq$value=sample(rsq$value)
      if (Module==TRUE) {
        rsq.f=rsqq
        rsq.f$value=as.numeric(rsq.f$value>0)
        rf=fisherTest(rsq.f)
        rsq.m=rsqq[rsqq$value>0,]
        rsq.m$term=factor(rsq.m$term,levels=unique(rsq.m$term))
        rrm=mwuTest(rsq.m,"g")
        rr0=rf[rf$term %in% rrm$term,]
        rr1=rf[!(rf$term %in% rrm$term),]
        rr0=rr0[order(rr0$term),]
        rrm=rrm[order(rrm$term),]
        rr0$pval=rr0$pval*rrm$pval
        rs=rbind(rr0,rr1)
      } else {
        if (mwut.t==TRUE) { rs=mwuTest(rsqq,Alternative) } else { rs=fisherTest(rsqq) }
      }
      spv=append(spv,rs$pval)
    }
    fdr=c()
    for (p in rr$pval){
      fdr=append(fdr,(sum(spv<=p)/reps)/sum(rr$pval<=p))
    }
    fdr[fdr>1]=1
  } else {
    fdr=p.adjust(rr$pval,method=adjust.multcomp)
  }
  
  message(sum(fdr<0.1)," GO terms at 10% FDR")
  rr$p.adj=fdr
  fname=paste("MWU_",inname,sep="")
  write.table(rr,fname,row.names=F)
}

#---------------------
mwuTest=function(gotable,Alternative) {
	gos=gotable
	terms=levels(gos$term)
	gos$seq=as.character(gos$seq)
	nrg=gos[!duplicated(gos$seq),5]
	names(nrg)=gos[!duplicated(gos$seq),4]
#	nrg=nrg+rnorm(nrg,sd=0.01) # to be able to do exact wilcox test
	rnk=rank(nrg)
	names(rnk)=names(nrg)
	pvals=c();drs=c();nams=c();levs=c();nseqs=c()
	for (t in terms){
		got=gos[gos$term==t,]
		got=got[!duplicated(got$seq),]
		ngot=gos[gos$term!=t,]
		ngot=ngot[!duplicated(ngot$seq),]
		ngot=ngot[!(ngot$seq %in% got$seq),]
		sgo.yes=got$seq
		n1=length(sgo.yes)
		sgo.no=ngot$seq
		n2=length(sgo.no)
#message(t," wilcox:",n1," ",n2)
		wi=wilcox.test(nrg[sgo.yes],nrg[sgo.no],alternative=Alternative)	# removed correct=FALSE 
		r1=sum(rnk[sgo.yes])/n1
		r0=sum(rnk[sgo.no])/n2
		dr=r1-r0
		drs=append(drs,round(dr,0))
		levs=append(levs,got$lev[1])
		nams=append(nams,as.character(got$name[1]))
		pvals=append(pvals,wi$p.value)
		nseqs=append(nseqs,n1)	
	}
	res=data.frame(cbind("delta.rank"=drs,"pval"=pvals,"level"=levs,nseqs))
	res=cbind(res,"term"=as.character(terms),"name"=nams)
	res$pval=as.numeric(as.character(res$pval))
	res$delta.rank=as.numeric(as.character(res$delta.rank))
	res$level=as.numeric(as.character(res$level))
	res$nseqs=as.numeric(as.character(res$nseqs))
	return(res)
}
#------------------------
fisherTest=function(gotable) {
	gos=gotable
	terms=levels(gos$term)
	gos$seq=as.character(gos$seq)
	pft=c();nam=c();lev=c();nseqs=c()
	for (t in terms) {
		got=gos[gos$term==t,]
		got=got[!duplicated(got$seq),]
		ngot=gos[gos$term!=t,]
		ngot=ngot[!duplicated(ngot$seq),]
		ngot=ngot[!(ngot$seq %in% got$seq),]
		go.sig=sum(got$value)
		go.ns=length(got[,1])-go.sig
		ngo.sig=sum(ngot$value)
		ngo.ns=length(ngot[,1])-ngo.sig
		sig=c(go.sig,ngo.sig) # number of significant genes belonging and not belonging to the tested GO category
		ns=c(go.ns,ngo.ns) # number of not-significant genes belonging and not belonging to the tested GO category
		mm=matrix(c(sig,ns),nrow=2,dimnames=list(ns=c("go","notgo"),sig=c("go","notgo")))
		ff=fisher.test(mm,alternative="greater")
		pft=append(pft,ff$p.value)
		nam=append(nam,as.character(got$name[1]))
		lev=append(lev,got$lev[1])
		nseqs=append(nseqs,length(got[,1]))
	}
	res=data.frame(cbind("delta.rank"=rep(0),"pval"=pft,"level"=lev,nseqs,"term"=terms,"name"=nam))
	res[,1]=as.numeric(as.character(res[,1]))
	res[,2]=as.numeric(as.character(res[,2]))
	res[,3]=as.numeric(as.character(res[,3]))
	res$nseqs=as.numeric(as.character(res$nseqs))
	return(res)
}

#-------------------------
gomwuPlot=function(inFile,goAnnotations,goDivision,level1=0.1,level2=0.05,level3=0.01,absValue=-log(0.05,10),adjusted=TRUE,txtsize=1,font.family="sans",treeHeight=0.5,colors=NULL) {
	require(ape)
	
	input=inFile
	in.mwu=paste("MWU",goDivision,input,sep="_")
	in.dissim=paste("dissim",goDivision,input,goAnnotations,sep="_")
	
	cutoff=-log(level1,10)
	pv=read.table(in.mwu,header=T)
	row.names(pv)=pv$term
	in.raw=paste(goDivision,input,sep="_")
	rsq=read.table(in.raw,sep="\t",header=T)
	rsq$term=as.factor(rsq$term)
	
	if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
	heat=data.frame(cbind("pval"=pvals)) 
	row.names(heat)=pv$term
	heat$pval=-log(heat$pval+1e-15,10)
	heat$direction=0
	heat$direction[pv$delta.rank>0]=1
	if (cutoff>0) { 
		goods=subset(heat,pval>=cutoff) 
	} else {
		goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
		goods=heat[row.names(heat) %in% goods.names,]
	}
	
	if (is.null(colors) | length(colors)<4 ) {
		colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
		if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
			colors=c("black","black","grey50","grey50")
		}
	}
	goods.names=row.names(goods)
	
	# reading and subsetting dissimilarity matrix
	diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
	row.names(diss)=names(diss)
	diss.goods=diss[goods.names,goods.names]
	
	# how many genes out of what we started with we account for with our best categories?
	good.len=c();good.genes=c()
	for (g in goods.names) {
		sel=rsq[rsq$term==g,]	
		pass=abs(sel$value)>=absValue
		sel=sel[pass,]
		good.genes=append(good.genes,as.character(sel$seq))
		good.len=append(good.len,nrow(sel))
	}
	ngenes=length(unique(good.genes))
	
	#hist(rsq$value)
	totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
	row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
	row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
	row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
	
	# clustering terms better than cutoff
	GO.categories=as.dist(diss.goods)
	cl.goods=hclust(GO.categories,method="average")
	labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
	goods=goods[labs,]
	labs=sub(" activity","",labs)

	old.par <- par( no.readonly = TRUE )

	plots=layout(matrix(c(1,2,3),1,3,byrow=T),c(treeHeight,3,1),TRUE)

    par(mar = c(2,2,0.85,0))
	plot(as.phylo(cl.goods),show.tip.label=FALSE,cex=0.0000001)
	step=100
	left=1
	top=step*(2+length(labs))

    par(mar = c(0,0,0.3,0))
	plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
	ii=1
	goods$color=1
	goods$color[goods$direction==1 & goods$pval>cutoff]=colors[4]
	goods$color[goods$direction==0 & goods$pval>cutoff]=colors[3]
	goods$color[goods$direction==1 & goods$pval>(-log(level2,10))]=colors[2]
	goods$color[goods$direction==0 & goods$pval>(-log(level2,10))]=colors[1]
	goods$color[goods$direction==1 & goods$pval>(-log(level3,10))]=colors[2]
	goods$color[goods$direction==0 & goods$pval>(-log(level3,10))]=colors[1]
	for (i in length(labs):1) {
		ypos=top-step*ii
		ii=ii+1
		if (goods$pval[i]> -log(level3,10)) { 
			text(left,ypos,labs[i],font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family) 
		} else {
			if (goods$pval[i]>-log(level2,10)) { 
				text(left,ypos,labs[i],font=1,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
			} else {
	#			if (goods$pval[i]>cutoff) { 
	#				text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
		#		} else { 
			text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family) 
			#}
			}
		}
	}
	
    par(mar = c(3,1,1,0))
	
	plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
	text(left,top-step*2,paste("p < ",level3,sep=""),font=2,cex=1* txtsize,adj=c(0,0),family=font.family)
	text(left,top-step*3,paste("p < ",level2,sep=""),font=1,cex=0.8* txtsize,adj=c(0,0),family=font.family)
	text(left,top-step*4,paste("p < ",10^(-cutoff),sep=""),font=3,col="grey50",cex=0.8* txtsize,adj=c(0,0),family=font.family)
	
	message("GO terms dispayed: ",length(goods.names))
	message("\"Good genes\" accounted for:  ", ngenes," out of ",totSum, " ( ",round(100*ngenes/totSum,0), "% )")
	par(old.par)	
	goods$pval=10^(-1*goods$pval)
	return(list(goods,cl.goods))
}

#ggtree version of gomwuPlot
gomwuPlot_gg <- function(
    inFile,
    goAnnotations,
    goDivision,
    level1 = 0.1,
    level2 = 0.05,
    level3 = 0.01,
    absValue = 0.001,
    adjusted = TRUE,
    txtsize = 3,
    font.family = "sans",
    treeHeight = 0.5,  # not really needed in ggtree but kept for compatibility
    data_type = "LFC", # specify data type: "LFC" for colors, "kME" for grayscale
    colors = NULL
) {
  # 1) Locate and read input files
  in.mwu    <- paste0("MWU_", goDivision, "_", inFile)
  in.dissim <- paste0("dissim_", goDivision, "_", inFile, "_", goAnnotations)
  
  # 2) Read in MWU results and raw data
  pv <- read.table(in.mwu, header = TRUE)
  rownames(pv) <- pv$term
  
  in.raw <- paste0(goDivision, "_", inFile)
  rsq    <- read.table(in.raw, sep = "\t", header = TRUE)
  rsq$term <- as.factor(rsq$term)
  
  # 3) Decide which p-values to use
  pvals <- if (adjusted) pv$p.adj else pv$pval
  
  # Create a small data frame of key info
  heat <- data.frame(
    rawp      = pvals,
    pval      = -log10(pvals + 1e-15),
    direction = ifelse(pv$delta.rank > 0, 1, 0),
    row.names = pv$term
  )
  
  # 4) Choose cutoff by p-value OR by absValue threshold
  cutoff <- -log10(level1)
  if (cutoff > 0) {
    goods <- subset(heat, pval >= cutoff)
  } else {
    # if cutoff <= 0, use absValue logic
    goods.names <- unique(rsq$term[ abs(rsq$value) >= absValue ])
    goods       <- heat[ rownames(heat) %in% goods.names, ]
  }
  
  # 5) Set color scheme manually based on data_type
  if (data_type == "kME") {
    grayscale <- TRUE
    if (is.null(colors) || length(colors) < 4) {
      colors <- c("black", "black", "grey50", "grey50")
    }
  } else if (data_type == "LFC") {
    grayscale <- FALSE
    if (is.null(colors) || length(colors) < 4) {
      colors <- c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral")
    }
  } else {
    stop("Unknown data_type. Must be either 'LFC' or 'kME'.")
  }
  
  # 6) Subset the dissimilarity matrix to these “good” GO terms
  diss <- read.table(in.dissim, sep = "\t", header = TRUE, check.names = FALSE)
  rownames(diss) <- colnames(diss)
  
  goods.names <- rownames(goods)
  diss.goods  <- diss[goods.names, goods.names, drop = FALSE]
  
  # 7) Count how many genes are covered by these “good” terms
  sel_list   <- list()  
  good.len   <- numeric(0)
  good.genes <- character(0)
  
  for (g in goods.names) {
    # Subset the data for the current term
    current_sel <- subset(rsq, term == g)
    
    # Filter the rows based on the absolute value criteria
    pass <- abs(current_sel$value) >= absValue
    current_sel <- current_sel[pass, ]
    
    # Store the filtered sel object in the list using the term as the key
    sel_list[[g]] <- current_sel
    
    # Append genes and count the number of rows for the current term
    good.genes <- c(good.genes, as.character(current_sel$seq))
    good.len   <- c(good.len, nrow(current_sel))
  }
  flattened_sel <- bind_rows(sel_list, .id = "term")
  ngenes <- length(unique(good.genes))
  totSum <- length(unique(rsq$seq[abs(rsq$value) >= absValue]))
  
  # ------------------ CHECK FOR ZERO OR ONE TERM BEFORE RENAMING ------------------
  n_goods <- nrow(goods)
  
  # Construct a plot_df "shell" to return even if empty
  plot_df <- data.frame(
    pval      = numeric(0),
    direction = integer(0),
    label     = character(0),
    original_term = character(0)
  )
  
  if (n_goods == 0) {
    # NO TERMS: Just print a message, return an empty result
    message("No significant GO terms above the threshold for: ", inFile)
    return(invisible(list(goods = plot_df, tree_plot = NULL)))
    
  } else if (n_goods == 1) {
    # ONE TERM: finalize rename, build minimal plot
    renameVec <- paste0(
      good.len, "/", 
      pv[match(goods.names, rownames(pv)), "nseqs"], " ",
      pv[match(goods.names, rownames(pv)), "name"]
    )
    
    rownames(goods)      <- renameVec
    rownames(diss.goods) <- renameVec
    
    # Build the single-row plot_df
    plot_df <- data.frame(
      rawp      = goods$rawp,
      pval      = goods$pval,
      direction = goods$direction,
      label     = rownames(goods),  # after rename
      original_term = goods.names,
      row.names = rownames(goods)
    )
    
    # Approx. p value and direction
    direction01     <- goods$direction[1]  # 0 => down, 1 => up
    direction_label <- ifelse(direction01 == 1, "Up-regulated", "Down-regulated")
    approx_p        <- signif(10^(-goods$pval[1]), 3)
    
    message("\nOnly one significant GO term found: ", rownames(goods)[1])
    message("Direction: ", direction_label)
    message("Approx. p-value: ", approx_p)
    
    # Minimal dummy plot
    library(ggplot2)
    dummy_plot <- ggplot() +
      annotate("text", x=0, y=0, label=paste(
        "Term:", rownames(goods)[1],
        "\n", direction_label,
        "\np=", approx_p
      ), size=5) +
      theme_void()
    
    print(dummy_plot)
    return(invisible(list(tree_plot=dummy_plot, goods = plot_df, passing_genes = flattened_sel)))
  }
  
  # 8) Now we safely rename rownames for final labeling if we have >= 2 terms
  renameVec <- paste0(
    good.len, "/", 
    pv[match(goods.names, rownames(pv)), "nseqs"], " ",
    pv[match(goods.names, rownames(pv)), "name"]
  )
  rownames(goods)      <- renameVec
  rownames(diss.goods) <- renameVec
  
  # Create a small data frame for plotting (>=2 terms)
  plot_df <- goods %>%
    dplyr::mutate(
      label         = rownames(.),
      original_term = goods.names
    )
  
  # Print summary info (like original)
  message("GO terms displayed: ", nrow(plot_df))
  message(
    "\"Good genes\" accounted for: ", ngenes, " out of ", totSum,
    " (", round(100*ngenes/totSum, 0), "%)"
  )
  
  # 9) Hierarchical clustering
  go_dist  <- as.dist(diss.goods)
  cl_goods <- hclust(go_dist, method = "average")
  
  # Convert to phylo
  phy_goods <- ape::as.phylo(cl_goods)
  
  # 10) Reorder 'goods' according to the hclust order so tips match
  plot_df <- plot_df[match(cl_goods$labels[cl_goods$order], plot_df$label), ]
  
  # 11) Determine color / fontface based on direction & p-value thresholds
  plot_df <- plot_df %>%
    dplyr::mutate(
      cutoff_l3 = -log10(level3),
      cutoff_l2 = -log10(level2),
      cutoff_l1 = -log10(level1),
      
      # "Raw" significance category
      sig_cat_raw = dplyr::case_when(
        pval > cutoff_l3 ~ "HIGH",
        pval > cutoff_l2 ~ "MED",
        pval > cutoff_l1 ~ "LOW",
        TRUE             ~ "NS"
      ),
      # Human-readable label for legend
      sig_cat = dplyr::case_when(
        pval > cutoff_l3 ~ paste0("p < ", level3),
        pval > cutoff_l2 ~ paste0("p < ", level2),
        pval > cutoff_l1 ~ paste0("p < ", level1),
        TRUE             ~ "NS"
      ),
      # direction => "up"/"down"
      direction_factor = factor(
        direction, levels = c(0,1),
        labels = c("down","up")
      ),
      # fontface
      fontface_col = dplyr::case_when(
        sig_cat_raw == "HIGH" ~ "bold",
        sig_cat_raw == "MED"  ~ "plain",
        TRUE                  ~ "italic"
      )
    )
  
  # Factor order
  plot_df$sig_cat <- factor(
    plot_df$sig_cat,
    levels = c(
      paste0("p < ", level3),
      paste0("p < ", level2),
      paste0("p < ", level1),
      "NS"
    )
  )
  
  # 12) Build the ggtree
  library(ggtree)
  library(dplyr)
  
  p <- ggtree(phy_goods, branch.length = "none") +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(5, 200, 5, 5))
  
  # Join the data
  p$data <- p$data %>%
    left_join(plot_df, by = "label")
  
  # 13) Plot with mapped color and size based on data_type
  if (grayscale) {
    p <- p +
      geom_tiplab(
        aes(
          label    = label,
          color    = direction_factor,
          size     = sig_cat,
          fontface = fontface_col
        ),
        show.legend = TRUE,
        family      = font.family
      ) +
      scale_color_manual(
        values = c("down" = colors[1], "up" = colors[2]),
        guide  = "none"  # Hide the legend for direction
      )
  } else {
    p <- p +
      geom_tiplab(
        aes(
          label    = label,
          color    = direction_factor,
          size     = sig_cat,
          fontface = fontface_col
        ),
        show.legend = TRUE,
        family      = font.family
      ) +
      scale_color_manual(
        values = c("down" = colors[1], "up" = colors[2]),
        name   = "Direction"
      )
  }
  
  p <- p +
    scale_size_manual(
      values = c(
        structure(3, names = paste0("p < ", level3)),
        structure(2, names = paste0("p < ", level2)),
        structure(1, names = paste0("p < ", level1)),
        structure(0.5, names = "NS")
      ),
      name = "Significance"
    ) +
    theme(
      axis.line  = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    ) +
    theme(
      legend.position = "left",
      legend.title    = element_text(face = "bold")
    ) +
    labs(
      title    = paste("GO MWU Tree -", inFile),
      subtitle = paste0("Filtered at p < ", level1, " or absValue >= ", absValue)
    )
  
  # 14) Print the final plot, return
  print(p)
  return(invisible(list(
    goods = plot_df,      # Equivalent to results[[1]]
    hcl = cl_goods,       # Equivalent to results[[2]]
    tree_plot = p,        # The plotted tree (if needed)
    passing_genes = flattened_sel   
  )))
  
}

# Extract GOs based on clustering of tree
extract_representative_GOs <- function(results, goDivision, input, gomwu_dir, 
                                       pcut = 1e-2, hcut = 0.9, plot_tree = FALSE) {
  message("\n--- Extracting Representative GO Terms ---\n")
  
  # Set working directory to the GO-MWU directory
  if (!dir.exists(gomwu_dir)) {
    stop("Error: The specified GO-MWU directory does not exist.")
  }
  setwd(gomwu_dir)
  message("Working directory changed to: ", getwd())
  
  # Check that results contains the required 'goods' component and that it has a 'rawp' column
  if (!("goods" %in% names(results))) {
    stop("Error: The results object must contain a 'goods' element.")
  }
  if (!("rawp" %in% colnames(results$goods))) {
    stop("Error: The 'goods' element must contain a 'rawp' column.")
  }
  
  # Initialize annotation storage
  annots <- c()
  num_terms <- nrow(results$goods)
  
  if (num_terms > 1) {
    if ("hcl" %in% names(results)) {
      # Check if the hcl object has enough leaves to form a dendrogram
      if (length(results$hcl$order) < 3) {
        message("Too few terms in hcl for clustering/plotting; selecting all terms passing significance filter.")
        significant_indices <- which(results$goods$rawp <= pcut)
        if (length(significant_indices) > 0) {
          annots <- sub("\\d+/\\d+ ", "", rownames(results$goods)[significant_indices])
        }
      } else {
        # Valid hcl available: cut the tree
        ct <- tryCatch({
          cutree(results$hcl, h = hcut)
        }, error = function(e) {
          message("Error cutting dendrogram for ", input, ": ", e$message)
          return(NULL)
        })
        
        if (is.null(ct)) {
          message("Dendrogram issues for ", input, "; selecting all terms passing significance filter.")
          significant_indices <- which(results$goods$rawp <= pcut)
          if (length(significant_indices) > 0) {
            annots <- sub("\\d+/\\d+ ", "", rownames(results$goods)[significant_indices])
          }
        } else {
          # Optionally plot the dendrogram if requested
          if (plot_tree) {
            plot_title <- gsub("\\.csv$", "", input)
            if (length(results$hcl$order) >= 3) {
              tryCatch({
                plot(results$hcl, cex = 0.6, main = plot_title)
                abline(h = hcut, col = "red")
              }, error = function(e) {
                message("Error plotting dendrogram for ", input, ": ", e$message)
              })
            } else {
              message("Not enough elements to plot dendrogram for ", input, "; skipping plot.")
            }
          }
          
          # Process clusters from the cut tree
          for (ci in unique(ct)) {
            message("Processing Cluster: ", ci)
            rn <- names(ct)[ct == ci]  # Terms in this cluster
            obs <- grep("obsolete", rn)  # Exclude obsolete terms
            if (length(obs) > 0) rn <- rn[-obs]
            if (length(rn) == 0) next  # Skip empty clusters
            
            rr <- results$goods[rn, , drop = FALSE]
            bestrr <- rr[which(rr$rawp == min(rr$rawp, na.rm = TRUE)), , drop = FALSE]
            best <- 1
            
            # If there are ties, break them using a frequency-based selection
            if (nrow(bestrr) > 1) {
              nns <- sub(" .+", "", row.names(bestrr))
              fr <- numeric(0)
              for (i in seq_along(nns)) {
                fr <- c(fr, eval(parse(text = nns[i])))
              }
              best <- which(fr == max(fr))[1]
            }
            
            if (!is.na(bestrr$rawp[best]) && bestrr$rawp[best] <= pcut) {
              annots <- c(annots, sub("\\d+/\\d+ ", "", row.names(bestrr)[best]))
            }
          }
        }
      }
    } else {
      # No hcl element present: fall back to selecting overall best GO term
      message("No hcl element present for ", input, "; selecting overall best GO term.")
      best_index <- which.min(results$goods$rawp)
      if (results$goods$rawp[best_index] <= pcut) {
        annots <- sub("\\d+/\\d+ ", "", rownames(results$goods)[best_index])
      }
    }
  } else if (num_terms == 1) {
    # Single-term case: process without clustering
    message("Only one term is differentially regulated; skipping tree clustering for ", input, ".")
    single_rawp <- results$goods[1, "rawp"]
    if (!is.na(single_rawp) && length(single_rawp) > 0 && single_rawp <= pcut) {
      term_name <- rownames(results$goods)[1]
      if (is.null(term_name) || term_name == "") term_name <- "term1"
      annots <- sub("\\d+/\\d+ ", "", term_name)
    }
  } else {
    stop("Error: No terms available in the 'goods' element.")
  }
  
  # Load the MWU file and filter the results based on selected annotations
  mwus_file <- paste("MWU", goDivision, input, sep = "_")
  if (!file.exists(mwus_file)) {
    stop("Error: MWU file not found in the specified directory: ", mwus_file)
  }
  
  mwus <- read.table(mwus_file, header = TRUE)
  bestGOs <- mwus[mwus$name %in% annots, ]
  
  message("\nRepresentative GO extraction complete for ", input, ".\n")
  return(bestGOs)
}


#------------------
# returns non-overlapping GO categories based on dissimilarity table
indepGO=function(dissim.table,min.similarity=1) {
	tt=read.table(dissim.table,sep="\t", header=TRUE)
	tt=as.matrix(tt)
	diag(tt)=1
	for (i in 1:ncol(tt)) {	
		mins=apply(tt,2,min)
		if (min(mins)>=min.similarity) { break }
		sums=apply(tt,2,sum)
		worst=which(sums==min(sums))[1]
		# cat("\n",worsts,"\n")
		# gw=c()
		# for(w in worsts) { gw=append(gw,sum(tt[,w]==1)) }
		# cat(gw)
		# worst=worsts[gw==min(gw)]
#		cat("\n",i," removing",worst,"; sum =",sums[worst])
		tt=tt[-worst,-worst]		
		mins=mins[-worst]		
#		cat(" new min =",min(mins))
	}
	goods=colnames(tt)
	goods=gsub("GO\\.","GO:",goods)
	goods=gsub("\\.GO",";GO",goods)
}

