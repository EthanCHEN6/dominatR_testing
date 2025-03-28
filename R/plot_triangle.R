

#' Plotting Triangle with center of mass
#'
#' @param data a table with three columns
#' @param entropyrange a vector with 2 elements in the range [0,1] i.e. c(minimum entropy , maximum entropy)
#' only the entries that has their entropy within this range are highlighted
#' @param magnituderange a vector with 2 elements in the range [0,1] i.e. c(minimum magnitude , maximum magnitude)
#' only the entries that has their magnitude within this range are highlighted
#' @param col a vector of 3 elements that contains color names
#' @param output_table a boolean value of T or F indicating whether to return the output table or not
#' @param plotAll a boolean value
#' If T, all the points that are out of entropyrange or magnitude range are plotted in whitesmoke color, else they are not plotted at all.
#' @param cex a numeric value indicating the size of the points
#' @param pch a integer value indicating different shapes of the points
#'
#' @return
#'  Returns a table with the coordinates of center of mass entropy(in the range of 0 and 1) and magnitude(in the range of 0 and 1)
#' @export
#' @import dplyr
#'
#'
plot_triangle=function(data,
                       entropyrange=c(0,1),
                       magnituderange=c(0,1),
                       col=c("red","green","blue"),
                       output_table=T,
                       plotAll=T,
                       cex=1,
                       pch=16)
{

  ### Making a polygon
  n = 3

  t=data.frame(x=sapply(0:(n-1),  function(x) sin(pi*(360/n)*x/180)),
               y=sapply(0:(n-1),  function(x) cos(pi*(360/n)*x/180)))


  #data=log2(data)
  colnames(data)=c("a","b","c")

  # calculate center of mass
  com = centmass(data, t$x, t$y)

  final=cbind(data,com)

  # calculating entropy
  entrop_data = entropy(data)$Entropy

  #Finding domination score
  final=final %>% mutate(entropy=minmax_normalization(entrop_data,0,1),
                         magnitude=minmax_normalization(a^2+b^2+c^2,0,1))

  # setting colors
  final=final %>% mutate(max=pmax(a,b,c)) %>%
    mutate(color=ifelse(a==max,col[1],ifelse(b==max,col[2],col[3])))

  # background colors
  final$color[!(final$entropy>entropyrange[1] &
                  final$entropy<entropyrange[2] &
                  final$magnitude>magnituderange[1] &
                  final$magnitude<magnituderange[2])]="whitesmoke"

  #plotting
  {
    plot(c(1,-1),c(1,-1),
         pch=16,
         xlab = "",
         ylab="",
         frame.plot = F,
         ann = F,
         axes = F,
         col="transparent")

    if(plotAll==T){
      badpoints=final %>%
        filter(color=="whitesmoke")

      points(badpoints$comx,
             badpoints$comy,
             col=badpoints$color,
             pch=pch,
             cex=cex)
      }

    goodpoints=final %>%
      filter(color!="whitesmoke")
    points(goodpoints$comx,
           goodpoints$comy,
           col=goodpoints$color,
           pch=pch,
           cex=cex)

    polygon(t$x,
            t$y,
            pch=16)
  }##Luka##

  if(output_table==T) return(final)
}
