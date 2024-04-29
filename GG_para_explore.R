library(flexsurv)
library(ggplot2)

#########################################################################
# GG
#########################################################################
# Define a sequence of y values for the log-transformed variable
x_values <- seq(0, 2, length.out = 400)

# Define shape parameters to compare
shape_parameters <- c(0.1, 1, 2)

# Compute the PDF values for the log-transformed generalized gamma distribution
pdf_values <- sapply(shape_parameters, function(shape) {
  dgengamma(x_values, mu = 1, sigma = 1, Q = shape)
})

# Create a data frame for plotting
plot_data <- data.frame(
  x = rep(x_values, times = length(shape_parameters)),
  pdf = as.vector(pdf_values),
  shape = factor(rep(shape_parameters, each = length(x_values)))
)

# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = pdf, color = shape)) +
  geom_line() +
  labs(title = "PDF of Log-Generalized Gamma Distribution with Various Shape Parameters",
       x = "x",
       y = "Probability Density",
       color = "Shape Parameter") +
  theme_minimal()

#########################################################################
# GLG
#########################################################################

# Define a sequence of y values for the log-transformed variable
y_values <- seq(-2, 5, length.out = 400)

# Define shape parameters to compare
shape_parameters <- c(0.1, 1, 2, 10)

# Compute the PDF values for the log-transformed generalized gamma distribution
pdf_values <- sapply(shape_parameters, function(shape) {
  exp(y_values) * dgengamma(exp(y_values), mu = 1, sigma = 0.5, Q = shape)
})

# Create a data frame for plotting
plot_data <- data.frame(
  y = rep(y_values, times = length(shape_parameters)),
  pdf = as.vector(pdf_values),
  shape = factor(rep(shape_parameters, each = length(y_values)))
)

# Plot using ggplot2
ggplot(plot_data, aes(x = y, y = pdf, color = shape)) +
  geom_line() +
  labs(title = "PDF of Log-Generalized Gamma Distribution with Various Shape Parameters",
       x = "Log-transformed x (y = log(x))",
       y = "Probability Density",
       color = "Shape Parameter") +
  theme_minimal()
