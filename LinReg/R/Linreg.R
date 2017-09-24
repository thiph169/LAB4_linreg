#' @title Linear regression
#' @field formula a formula, format y ~ x_1 + x_2 + ... + x_n.
#' @field data a data frame.
#' @description Returns the result of the Linear Regression
#' @examples 
#' data(iris)
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$print()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$plot()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$resid()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$pred()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$coef()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$summary()
#' @export linreg
#' @export
linreg <- setRefClass(
  Class = "linreg",
  fields = list(
    formula = "formula",
    data = "data.frame",
    beta_hat = "matrix",
    y_hat = "matrix",
    e_hat = "matrix",
    df_linreg = "numeric",
    var_sigma_hat = "numeric",
    p_value = "matrix",
    t_value = "matrix",
    var_beta_hat = "matrix",
    data_name = "character"
  ),
  methods = list(
    initialize = function(formula, data){
      formula <<- formula
      data <<- data
      data_name <<- deparse(substitute(data))
      
      
      # Define matrix
      X <- model.matrix(formula, data)
      
      # extract the dependant variable from data set
      dep_name <- all.vars(expr = formula)[1]
      y <- (data[, dep_name])
      
      # Calculate coefficients
      beta_hat <<- solve((t(X) %*% X)) %*% t(X) %*% y
      
      # find the fitted values of y using beta_hat
      y_hat <<- X %*% beta_hat
      
      # find the residual values
      e_hat <<- y - y_hat
      
      # degrees of freedom
      n <- nrow(X) # number of observations
      p <- ncol(X) # number of parameters
      
      df_linreg <<- n - p
      
      # residual variance
      var_sigma_hat <<- as.numeric((t(e_hat) %*% e_hat) / df_linreg)
      var_beta_hat <<- var_sigma_hat * solve((t(X) %*% X))
      
      #The t-values for each coefficient
      t_value <<- beta_hat / sqrt(diag(var_beta_hat))
      p_value <<- pt(abs(t_value), df = df_linreg,lower.tail=FALSE)
    },
    # Build linreg print function
    print = function() {
      # Print with help of a list
      # return(
      #     list(
      #         Formula_call = formula,
      #         Regression_Coefficients = beta_hat
      #     )
      # )
      
      # Print function call
      cat(sep = "\n")
      cat("Call:")
      cat(sep = "\n")
      cat(paste("linreg(", "formula = ", formula[2], " ", formula[1], " ", formula[3], ", ", "data = ", data_name, ")", sep = "" ))
      
      cat(sep = "\n")
      cat(sep = "\n")
      cat("Coefficients:")
      
      cat(sep = "\n")
      
      # opt1.
      # Print coef with named vector, causes error in testthat
      # print_coef <- as.vector(beta_hat)
      # names(print_coef) <- row.names(beta_hat)
      # return(print_coef)
      
      # opt2.
      # Print with cat(), error cause by too few blankspaces between output, ugly output
      cat(row.names(beta_hat))
      cat(sep = "\n")
      cat(beta_hat)
      
    },
    
    # Build linreg plot function
    plot = function() {
      library(ggplot2)
      
      # Build plotting theme
      liu_blue <- "#54D8E0"
      theme_liu <- theme(plot.margin = unit(c(1,1,1,1), "cm"), 
                         panel.background = element_rect(fill="white"),
                         panel.grid.major.y = element_blank(),
                         panel.grid.minor.y = element_blank(),
                         panel.grid.major.x = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         axis.line = element_line(color= "#58585b", size=0.1),
                         axis.text.x = element_text(color="Black", size="10"),
                         axis.text.y = element_text(color="Black", size="10"),
                         axis.title.x = element_text(color="Black", size="10", face="bold"),
                         axis.title.y = element_text(color="Black", size="10", face="bold"),
                         axis.ticks.y = element_blank(),
                         axis.ticks.x = element_line(color = "#58585b", size = 0.3),
                         plot.title = element_text(color="Black", face="bold", size="14"),
                         legend.position="bottom", legend.title = element_blank(),
                         legend.key = element_blank(),
                         legend.text = element_text(color="Black", size="10"))
      
      # Build plot data
      plot_df <- data.frame(
        df_resid = e_hat,
        df_fitted_values = y_hat)
      
      # Build plot1, Residual vs fitted values
      p1 <-
        ggplot(data = plot_df, aes(x = df_fitted_values, y = df_resid)) +
        geom_point(colour = liu_blue) +
        geom_smooth(method = "loess",
                    color = "red",
                    se = FALSE) +
        geom_abline(slope = 0,
                    intercept = 0,
                    linetype = "dotted") +
        ggtitle("Residual vs Fitted") +
        ylab("Residuals") +
        xlab("Fitted Values") +
        theme_liu
      
      # Build plot1, Scale Location
      p2 <- ggplot(data = plot_df,
                   aes(x = df_fitted_values, y = sqrt(abs((df_resid - mean(df_resid)) / sqrt(var_sigma_hat)
                   )))) +
        geom_point(colour = liu_blue) +
        geom_smooth(method = "loess",
                    color = "red",
                    se = FALSE) +
        ggtitle("Scale-Location") +
        ylab(expression(sqrt(abs("Standardized Residuals")))) +
        xlab("Fitted Values") +
        theme_liu
      
      return(list(Residual_vs_Fitted = p1, 
                  Scale_Location = p2))
    },
    # Build linreg resid print function
    resid = function() {
      return((Residuals = round(e_hat, 2)))
    },
    # Build linreg pred print function
    pred = function() {
      return((Fitted_values = round(y_hat, 2)))
    },
    # Build linreg coef print function
    coef = function() {
      vector <- as.vector(beta_hat)
      vect_names <-  rownames(beta_hat)
      names(vector) <-  vect_names
      return(vector)
    },
    # Build linreg summary print function
    summary = function() {
      # coef_mx <- as.matrix(cbind(
      #     round(beta_hat, 2),
      #     round(sqrt(diag(
      #         var_beta_hat
      #     )), 2),
      #     round(t_value, 2),
      #     round(p_value, 4)
      # ))
      # 
      # colnames(coef_mx) <-
      #     c("Estimate", "Sd.Error", "T-value", "P-value")
      
      coef_mx <- data.frame(
        var = rownames(beta_hat),
        estimate = round(beta_hat, 2),
        std.error = round(sqrt(diag(
          var_beta_hat
        )), 2),
        t_value = round(t_value, 2),
        p_value = round(p_value, 4)
      )
      
      coef_mx$var <- as.character(coef_mx$var)
      # str(coef_mx)
      row.names(coef_mx) <- NULL
      
      coef_mx <- rbind(c(" ", "Estimate", "Std. Error", "t value", "Pr(>|t|)"), coef_mx)
      
      for(i in 2:nrow(coef_mx)){
        if(coef_mx$p_value[i] == 0){
          coef_mx$p_value[i] <- "***"
        } else if(coef_mx$p_value[i] > 0 & coef_mx$p_value[i] <= 0.001){
          coef_mx$p_value[i] <- "**"
        } else if(coef_mx$p_value[i] > 0.001 & coef_mx$p_value[i] <= 0.01){
          coef_mx$p_value[i] <- "*"
        } else if(coef_mx$p_value[i] > 0.01 & coef_mx$p_value[i] <= 0.05){
          coef_mx$p_value[i] <- "."
        } else if(coef_mx$p_value[i] > 0.05 & coef_mx$p_value[i] <= 0.1){
          coef_mx$p_value[i] <- " "
        } else if(coef_mx$p_value[i] > 0.1){
          coef_mx$p_value[i] <- " "
        }
      }
      
      
      
      for(c in 1:ncol(coef_mx)){
        wdth <- max(nchar(as.character(coef_mx[, c])), na.rm = TRUE)
        for(r in 1:nrow(coef_mx)){
          coef_mx[r, c] <- format(coef_mx[r, c], width = wdth, justify = c("right"))
        }
      }
      
      # return(
      #     list(
      #         Call = paste("linreg(", "formula = ", formula[2], " ", formula[1], " ", formula[3], ", ", "data = ", data_name, ")", sep = "" ),
      #         Residuals = c(
      #             Min = round(min(e_hat), 2),
      #             round(quantile(e_hat, .25), 2),
      #             Median = round(median(e_hat), 2),
      #             round(quantile(e_hat, .75), 2),
      #             Max = round(max(e_hat), 2)
      #         ),
      #         Coefficients = coef_mx,
      #         Evaluation = c(
      #             paste(
      #                 "Residual standard error: ",
      #                 round(sqrt(var_sigma_hat), 4),
      #                 " on ",
      #                 df,
      #                 " degrees of freedom"
      #             )
      #         )
      #     )
      # )
      
      # Prints function call
      cat(sep = "\n")
      cat("Call:")
      cat(sep = "\n")
      cat(paste("linreg(", "formula = ", formula[2], " ", formula[1], " ", formula[3], ", ", "data = ", data_name, ")", sep = "" ))
      
      
      
      # Print coef matrix with results and names, still causes errors.
      cat(sep = "\n")
      cat(sep = "\n")
      cat("Coefficients:")
      cat(sep = "\n")
      #cat(format(coef_mx, width=max(nchar(coef_mx), nchar(coef_mx)), justify = c("right")))
      cat(sep = "\n")
      for(i in 1:nrow(coef_mx)){
        cat(paste(as.character(coef_mx[i, ]),collapse = " "))
        cat(sep = "\n")
      }
      
      # Print df
      cat(sep = "\n")
      cat("Residual standard error:", round(sqrt(var_sigma_hat), 2), "on", df_linreg, "degrees of freedom")
    }
  )
)
