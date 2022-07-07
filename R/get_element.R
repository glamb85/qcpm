#' @title get_element
#' 
#' @details
#' Internal function. \code{get_element} is called by \code{get_paths} 
#' @param path_matrix the matrix of path coefficients 
#' @return the path coefficients labels
#' @keywords internal
#' @export
#'
get_element	<-	function(path_matrix)
{
  y = mat.or.vec(1,sum(ifelse(path_matrix==0,0,1)))
  k = 1
  for (i in 1:dim(path_matrix)[1]) 
  {
    for (j in 1:dim(path_matrix)[1])
    {
      if(path_matrix[i,j] != 0)
      {
        v	= paste(colnames(path_matrix)[j],rownames(path_matrix)[i],sep = "->")
        y[k]	= v
        k	= k+1
      }
    }
  }
  as.vector(y)
}

#' @title get_element
#' 
#' @details
#' Internal function. \code{get_element} is called by \code{get_paths} 
#' @param path_matrix the matrix of path coefficients 
#' @return the path coefficients labels
#' @keywords internal
#' @export
#'
get_names_MV	<-	function(x)
{
  y = mat.or.vec(1,sum(ifelse(x==0,0,1)))
  k = 1
  for (i in 1:dim(x)[1]) 
  {
    for (j in 1:dim(x)[2])
    {
      if(x[i,j] != 0)
      {
        v	= paste(rownames(x)[i],colnames(x)[j],sep = "-")
        y[k]	= v
        k	= k+1
      }
    }
  }
  as.vector(y)
}
