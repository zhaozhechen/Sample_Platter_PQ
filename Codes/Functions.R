library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# This function is to get the length of all interpolated length
# Input includes:
# varname: name for the variable to investigate
# df: data frame to process
interpolated_interval <- function(varname,df){
  # Get the target variable
  df_tmp <- Site_df %>%
    filter(var == varname)
  # Get the index of interpolated data
  interpolated_idx <- which(df_tmp$ms_interp == 1)
  if(length(interpolated_idx) != 0){
    # Get the length of each interpolated segment length
    # First get breaks of the segments
    group_id <- cumsum(c(1, diff(interpolated_idx) != 1))
    # Then get lengths of each segment
    segment_length <- table(group_id)
  }else{
    segment_length <- NA
  }
  return(segment_length)
}


# Theme for all plots
my_theme <- theme(
  #axis.line=element_line(color="black"),
  panel.background = element_blank(),
  panel.border = element_rect(colour="black",fill=NA),
  legend.key = element_blank(),
  #legend.key.size = unit(6,"cm"),
  #aspect.ratio = 1/1,
  #legend.key.size = unit(0.3,'cm'),
  legend.text = element_text(size=14),
  plot.title = element_text(size=14),
  axis.text = element_text(size=14),
  axis.title = element_text(size=14)
)

# This function is to print pdf and png figure
# Input is the figure g,title,width, and height
print_g <- function(g,title,w,h){
  pdf(paste0(Output_path,"/",title,".pdf"),
      width=w,height=h)
  print(g)
  dev.off()
  png(paste0(Output_path,"/",title,".png"),
      width=w,height=h,units = "in",
      res=600)
  print(g)
  dev.off()
}

# This function is to make TS plot of target variable
# Input include
# varname: name for the variable to investigate
# df: data frame to process
# Site_ID: Site code
TS_plot <- function(varname,df,Site_ID){
  # Keep target variable
  df_tmp <- df %>%
    filter(var == varname) %>%
    mutate(date = as.Date(date))
  if(varname == "discharge"){
    g <- ggplot(data = df_tmp,aes(x = date,y = val))+
      geom_segment(aes(xend=date,y=0,yend=val),color=my_color[3])+
      #geom_line(color=my_color[3])+
      my_theme+
      labs(x = "",y="Q",color="")+
      ggtitle(Site_ID)
  }else if(varname == "precipitation"){
    g <- ggplot(data = df_tmp,aes(x = date,y = val,color=factor(ms_interp)))+
      geom_segment(aes(xend=date,y=0,yend=val))+
      #geom_line()+
      my_theme+
      labs(x = "",y="P",color="")+
      scale_color_manual(values = c("0" = my_color[1],
                                    "1" = my_color[2]),
                         labels = c("0" = "Not interpolated",
                                    "1" = "Interpolated"))+
      theme(legend.position = c(0.85,0.9)) 
  }
  return(g)
}

# This function is to make histogram of target variable
# Input include
# varname: name for the variable to investigate
# df: data frame to process
Hist_plot <- function(varname,df){
  # Keep target variable
  df_tmp <- df %>%
    filter(var == varname) %>%
    mutate(date = as.Date(date))
  if(varname == "discharge"){
    g <- ggplot(data=df_tmp,aes(x=val))+
      geom_histogram(bins = 11,color="black",fill=my_color[3])+
      my_theme+
      labs(x=varname)
  }else if(varname == "precipitation"){
    g <- ggplot(data=df_tmp,aes(x=val,fill=as.factor(ms_interp)))+
      geom_histogram(bins=11,color="black",alpha=0.7)+
      scale_fill_manual(values = c("0" = my_color[1],
                                   "1" = my_color[2]),
                        labels = c("0" = "Not interpolated",
                                   "1" = "Interpolated"))+
      theme(legend.position = c(0.7,0.8))+
      my_theme+
      labs(x=varname,fill="")
  }
  return(g)
}

