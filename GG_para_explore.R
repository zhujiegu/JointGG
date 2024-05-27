library(flexsurv)
library(ggplot2)



x <- seq(0, 30, length.out = 1000)

# Parameters for the first generalized gamma distribution
mu1 <- 1
sigma1 <- 0.4
q1 <- -0.15

# Parameters for the second generalized gamma distribution
mu2 <- 1.5
sigma2 <- 0.34
q2 <- -0.15

# Calculate the density values for each distribution
pdf1 <- dgengamma(x, mu1, sigma1, q1)
pdf2 <- dgengamma(x, mu2, sigma2, q2)

h1 <- hgengamma(x, mu1, sigma1, q1)
h2 <- hgengamma(x, mu2, sigma2, q2)
# Plot the PDFs
plot(x, pdf1, type = "l", col = "blue", lwd = 2, xlab = "x", ylab = "Density", main = "Generalized Gamma PDFs")
lines(x, pdf2, col = "red", lwd = 2)
legend("topright", legend = c("Distribution 1", "Distribution 2"), col = c("blue", "red"), lwd = 2)

plot(x, h1, type = "l", col = "blue", lwd = 2, xlab = "x", ylab = "hazard", main = "Generalized Gamma Hazard")
lines(x, h2, col = "red", lwd = 2)
legend("topright", legend = c("Distribution 1", "Distribution 2"), col = c("blue", "red"), lwd = 2)



#########################################################################
# GG
#########################################################################
# Define a sequence of y values for the log-transformed variable
x_values <- seq(0, 100, length.out = 100000)

# Define shape parameters to compare
shape_parameters <- c(-0.08)

# Compute the PDF values for the log-transformed generalized gamma distribution
pdf_values <- sapply(shape_parameters, function(shape) {
  dgengamma(x_values, mu = 1.7, sigma = 0.09, Q = shape)
  # pgengamma(x_values, mu = 2.5, sigma = 0.4, Q = shape, lower.tail = F)
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
  labs(title = "PDF of Generalized Gamma Distribution (mu=2, sigma=3, q=1)",
       x = "Survival time",
       y = "Probability Density") +
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
