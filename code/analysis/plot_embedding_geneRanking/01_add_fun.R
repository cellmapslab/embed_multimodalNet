node='asthma'

plot_pca_within <- function(node){
  print(node)
  newnode <- node_name[old==node, old2]
  
  
  pca_sub_dt2 <- merge(pca_data_sub, node_anno, by.x="V1", by.y="name", all.x=T)
  pca_sub_dt3 <- merge(pca_sub_dt2, all[run==i][node1==node][order(-sim)][,.(node2, sim)][1:150], by.x="V1", by.y="node2",all.x=T)
  pca_sub_dt3[V1==node, type:=node]
  pca_sub_dt3[V1==node, sim:=1]
  pca_sub_dt3[!is.na(sim), mysize:=0.5]
  pca_sub_dt3[is.na(sim), mysize:=0.6]
  
  mysim_cutoff <- pca_sub_dt3[order(-sim)][16,sim]
  
  # change name of nodes
  pca_sub_dt3 <- change_names(pca_sub_dt3, node_name, mycol="V1")
  
  # node_colors <- rev(c("#FFAC1C","#C7CEEA", "#F67280", "#C06C84",
  #                      "#6C5B7B", "#355C7D")) 
  
  # node_types <- c("genes", "pheno", "como", "prs", "tissue", newnode)
  # mycolors <- setNames(node_colors, node_types)
  
  mycolors <- as.character(c(colors$node_colors,"#FFAC1C" ))
  names(mycolors) <- as.character(c(colors$node_types, newnode))
  
  pca_sub_dt3[V1==newnode, type:=newnode]
  library(ggpubr)
  library(ggrepel )
  options(ggrepel.max.overlaps = Inf)
  
  x_val=pca_sub_dt3$Dim.1
  y_val=pca_sub_dt3$Dim.2
  
  x_min=ceiling(min(x_val))#-abs(max(x_val) - min(x_val))/2
  x_max=ceiling(max(x_val))
  
  y_min=ceiling(min(y_val))
  y_max=ceiling(max(y_val))#+abs(max(y_val) - min(y_val))/2

  plt1 <- ggplot(pca_sub_dt3[order(sim,na.last = F)], aes(x=Dim.1, y=Dim.2, color=sim, label = V1)) + 
    geom_point(aes(size=mysize)) +theme_tufte() +
    scale_size_continuous(range = c(0.5, 1.2))+
    scale_x_continuous(breaks=seq(x_min,x_max,2)) +
    scale_y_continuous(breaks=seq(y_min,y_max,2)) +
    guides(size = "none") +
    labs(x="PC1", y="PC2") +  
    scale_colour_gradient(low = "#d9a6b5",high = "#603642", na.value = "#f2e1e6") +
    theme(axis.line = element_line(size = 0.5))
    # theme(axis.line=element_blank(),
    #   axis.text.x=element_blank(),
    #   axis.text.y=element_blank(),
    #   axis.ticks=element_blank(),
    #   axis.title.x=element_blank(),
    #   axis.title.y=element_blank(),
    #   panel.background=element_blank(),
    #   panel.border=element_blank(),
    #   panel.grid.major=element_blank(),
    #   panel.grid.minor=element_blank(),
    #   plot.background=element_blank())
  
  x_val=pca_sub_dt3[!is.na(sim)]$Dim.1
  y_val=pca_sub_dt3[!is.na(sim)]$Dim.2
  
  x_min=min(x_val)#-abs(max(x_val) - min(x_val))/2
  x_max= max(x_val)+abs(max(x_val) - min(x_val))/5
  
  y_min=min(y_val)-abs(max(x_val) - min(x_val))/5
  y_max=max(y_val)#+abs(max(y_val) - min(y_val))/2
  
  plt2<- ggplot(pca_sub_dt3[!is.na(sim),], aes(x=Dim.1, y=Dim.2, color=type)) + 
  geom_point(size=1) +theme_tufte() +
    scale_color_manual(values = mycolors)+ labs(x="PC1", y="PC2") +  
    # xlim(x_min, x_max) +
    # ylim(y_min, y_max) +
    # scale_x_continuous(expand=c(0,0), limits=c(x_min,x_max)) +
    # scale_y_continuous(expand=c(0,0), limits=c(y_min,y_max)) +
    scale_x_continuous(limits=c(x_min,x_max)) +
    scale_y_continuous(limits=c(y_min,y_max)) +
    # scale_x_continuous(breaks=seq(x_min,x_max,2)) +
    # scale_y_continuous(breaks=seq(y_min,y_max,2)) +
    geom_label_repel(aes(label=ifelse(sim>mysim_cutoff,as.character(V1),'')),
                     point.padding = 0.5, 
                     box.padding = 1,
                     size = 2,show_guide  = FALSE)+
    theme(axis.line = element_line(size = 0.5))
  l1 <- get_legend(plt1)
  l2 <- get_legend(plt2)
  plt1_1<-plt1+
    theme(#plot.background = element_rect(colour = "grey"),
          legend.position = "none")
  
  
  g_plt1 <- ggplotGrob(plt1_1)
  plt_both <- plt2+ theme(legend.position = "none") +
  annotation_custom(
    grob = g_plt1,
    xmin = x_max-(x_max + abs(x_min))/3,
    xmax = x_max,#+abs(x_max -x_min)/5,
    ymin = y_min,#-abs(x_max - x_min)/5,
    ymax = y_min+(y_max + abs(y_min))/3
  )
  final_plot <- grid.arrange(plt_both,l1,l2,
                                layout_matrix = rbind(c(1,1,1,1,1,3),
                                                      c(1,1,1,1,1,2)))
  return(final_plot)
}


##############

plot_sim_ranking_all <- function(node){
  # get data
  dt = results[node1==node]#[type_node2=="genes"]
  dt <- dt[order(-max_sim)]
  # only top 50
  dt <- dt[1:50]
  # change node name
  dt <- change_names(dt, node_name, mycol="node2")
  node <- node_name[old==node, old2]
  dt[,type_col:=paste0(type_node2, "_", neighbour)]
  dt[type_node2=="gene", type_node2:=type_col]

  # # order the node type so that genes, prs, como,.. will always have the same color
  # n_type <- unique(dt$type_node2)
  # mycols <- colors[match(n_type, colors$node_types),]

  mycolors <- as.character(c(colors$node_colors))
  names(mycolors) <- as.character(c(colors$node_types))

  idorder <- rev(dt$node2)
  dt$node2 <- factor(dt$node2, levels=idorder)
  x_m <- max(dt$max_sim)

  # idorder <- rev(dt$type_node2)
  # dt$type_node2 <- factor(dt$type_node2, levels=mycols$node_types)
  plt1 <-
    ggplot(dt[1:15], aes(x = node2, y = max_sim, fill = type_node2)) +
    geom_bar(stat = "identity", width = 0.9) +
    scale_fill_manual(values = mycolors) +
    theme_tufte() +
    labs(x = "nodes", y = "max sim score", title=node) +
    #theme(legend.position = "bottom") +
    guides(fill=guide_legend(title="node type"))+
    coord_flip(ylim = c(0.6, x_m))
  return(plt1)
  
}

get_plot <- function(node){
  plt1 <- plot_pca_within(node)
  plt2 <- plot_sim_ranking_all(node)
  return(grid.arrange(plt1,plt2,
               layout_matrix = rbind(c(1,1,1,1,2,2,2))))
}

