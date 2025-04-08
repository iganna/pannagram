#' Heatmap Plot Function
#'
#' This function generates a heatmap plot using ggplot2 based on a table, data frame, or matrix input.
#'
#' @param tbl A table, data frame, or matrix to be used for generating the heatmap.
#' @param cols A vector of colors to use for the heatmap. If NULL, a predefined color palette will be used based on `c.name`.
#' @param c.name A string indicating the name of the predefined color palette to use. 
#'                Options include 'lila', 'mint', 'sea', 'pastel', 'div1', 'div2', 'div3'. Default is 'lila'.
#' @param xlab A string to label the x-axis. Default is NULL (no label).
#' @param ylab A string to label the y-axis. Default is NULL (no label).
#' @param fill.lab A string to label the fill legend. Default is NULL (no label).
#'
#' @return A ggplot object representing the heatmap.
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(25), nrow=5)
#' p <- heatplot(mat, c.name = 'mint')
#' p
#' }
#'
#' @export
heatplot <- function(tbl, 
                     cols = NULL,
                     c.name = 'gray_sea',
                     xlab=NULL,
                     ylab=NULL,
                     fill.lab=NULL,
                     cord.fix = F,
                     add.label = F,
                     to.norm = 'none'
                     ) {
  
  cols.list = list('lila'     = c('#F5EDED', '#CB80AB', '#8967B3', '#624E88'),
                   'mint'     = c('#F6F6F6', '#97DECE', '#439A97', '#2E4F4F'),
                   'sea'      = c('#FFDC7F', '#78B7D0', '#227B94', '#16325B'),
                   'pastel'   = c('#F4EDCC', '#A4CE95', '#6196A6', '#5F5D9C'),
                   'gray'     = c('white', 'grey90', 'grey70', 'grey50'),
                   'gray_sea' = c('#F8FAFC', '#D9EAFD', '#BCCCDC', '#9AA6B2'),
                   'pale_sea' = c('#F8F4FC', '#E0D2F5', '#C6AEE2', '#A989C7'),
                   
                   
                   'div1' = c("#399918", "#88D66C", "#ECFFE6", "#FFAAAA", "#FF7777"),
                   'div2' = c('#439A97', '#91DDCF', '#F7F9F2', '#F19ED2', '#CD6688'),
                   'div3' = c('#439A97', '#96CEB4', '#FFF7D1', '#FFAD60', '#A66E38'),
                   'redblue' = c('#11468F','#CDDEFF', 'white',  '#FFBCBC', '#EB455F'),
                   'dot' = c("#CE1F6A", 'white', '#27374D'))
    
  if(!(to.norm %in% c('none', 'row', 'col'))) stop('Wrong normalisation')
  
  if(is.null(cols)){
    if(c.name %in% names(cols.list)){
      cols = cols.list[[c.name]]
    } else {
      pokazAttention('Avaliable colorpallets:', names(cols.list))
      stop(paste('Colorpallet', c.name, 'doesnt exist.'))
    }
  }
  
  tbl.class = class(tbl)[1]
  if(tbl.class == 'table'){
    tbl = matrix(tbl, ncol = ncol(tbl), 
                 dimnames = list(rownames(tbl), colnames(tbl)))
  } else  if(tbl.class == "data.frame") {
    tbl = as.matrix(tbl)
  } else if(tbl.class != 'matrix'){
    stop('Class should be table/data.frame/matrix.')
  }
  
  df <- reshape2::melt(tbl)
  df$Var1 = factor(df$Var1)
  df$Var2 = factor(df$Var2)
  df$label = df$value
  
  if(to.norm == 'row'){
    for(tmp in unique(df$Var1)){
      tmp.val = df$value[df$Var1 == tmp]
      tmp.val = (tmp.val - min(tmp.val)) / (max(tmp.val) - min(tmp.val))
      df$value[df$Var1 == tmp] = tmp.val
    }
  }
  
  if(to.norm == 'col'){
    for(tmp in unique(df$Var2)){
      tmp.val = df$value[df$Var2 == tmp]
      tmp.val = (tmp.val - min(tmp.val)) / (max(tmp.val) - min(tmp.val))
      df$value[df$Var2 == tmp] = tmp.val
    }
  }
  
  p <- ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    theme_minimal() +
    labs(x = xlab, y = ylab, fill = fill.lab) +
    scale_x_discrete(expand = c(0, 0)) +   
    scale_y_discrete(expand = c(0, 0))
  
  if(add.label){
    p = p + geom_text(aes(label = round(label, 2)), color = "black")
  }
  
  if(c.name == 'dot'){
    wsize = max(abs(df$value))
    p = p + scale_fill_gradient2(low = "#CE1F6A", mid = "white", high = "#27374D",
                                 breaks = c(-wsize, 0, wsize))
  } else {
    p = p + scale_fill_gradientn(colors = cols)
  }
  
  if(cord.fix){
    p = p + coord_fixed()
  }
  
  
  return(p)
}

