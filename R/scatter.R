## Matt Galloway


#' @title Scatter
#' @description This function simply streamlines the process of creating a scatterplot with ggplot
#'
#' @param data. data frame
#' @param x. x-axis
#' @param y. y-axis
#' @return a scatterplot
#' @export
#' @examples
#' scatter(iris, Sepal.Length, Sepal.Width)



# we define the scatter plot function
scatter = function(data., x., y.) {
    
    # substitute
    x. = substitute(x.)
    y. = substitute(y.)
    
    # create scatterplot
    ggplot(data.) + geom_point(mapping = aes(x = eval(x., 
        data.), y = eval(y., data.))) + ggtitle("Scatterplot") + 
        ylab(y.) + xlab(x.)
    
}
