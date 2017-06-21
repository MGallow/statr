## Matt Galloway


#' @title Multiple Plot
#' @description Taken from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#'
#' @param ... object can be passed in
#' @param cols number of columns in layout
#' @param layout a matrix specify the layout. If present, 'cols' is ignored
#' @return plots
#' @examples
#' multiplot(p1, p2, cols = 1)

multiplot <- function(..., plotlist = NULL, file, cols = 1, 
    layout = NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel ncol: Number of columns of plots nrow:
        # Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), 
            ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots == 1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), 
            ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain
            # this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, 
                layout.pos.col = matchidx$col))
        }
    }
}



##----------------------------------------------------------------------




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
    ggplot(data.) + geom_point(mapping = aes(x = eval(x., data.), 
        y = eval(y., data.))) + ggtitle("Scatterplot") + ylab(y.) + 
        xlab(x.)
    
}



##----------------------------------------------------------------------


#' @title Diagnostic
#' @description This function simply streamlines the process of creating diagnostic plots with ggplot
#'
#' @param data. data frame
#' @param x. x-axis
#' @param y. y-axis
#' @return a residual plot and QQ plot
#' @export
#' @examples
#' diagnostic(iris, Sepal.Length, Sepal.Width)

# we define the scatter plot function
diagnostic = function(data., x., y.) {
    
    # substitute
    x. = substitute(x.)
    y. = substitute(y.)
    
    # create residual plot
    fit = lm(eval(y., data.) ~ eval(x., data.), data.)
    
    residual_plot = ggplot(data., mapping = aes(x = fit$fitted.values, 
        y = fit$residuals)) + geom_abline(intercept = 0, slope = 0, 
        color = "red") + geom_point() + ggtitle("Residual Plot") + 
        ylab("residuals") + xlab("fitted values")
    
    # create qq plot
    y <- quantile(eval(y., data.), c(0.25, 0.75))  # Find the 1st and 3rd quartiles
    x <- qnorm(c(0.25, 0.75))  # Find the matching normal values on the x-axis
    slope <- diff(y)/diff(x)  # Compute the line slope
    int <- y[1] - slope * x[1]  # Compute the line intercept
    
    qqplot = ggplot(data.) + stat_qq(mapping = aes(sample = eval(y., 
        data.))) + geom_abline(intercept = int, slope = slope, 
        color = "red") + ggtitle("QQ Plot")
    
    # output plots
    multiplot(residual_plot, qqplot, cols = 1)
    
}
