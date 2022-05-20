library(tidyverse)

# some sample script to see how the carrying capacity
# equation from M'Gonigle et al 2012 actually looks like

# environmental gradient
gradient_x <- seq(0,1,0.01)
gradient_y <- seq(0,0.5,0.01)

# calculate the carrying capacity
b <- 1.0
k0 <- 1.0

sigma_k <- 0.1

# the carrying capacity function
k <- function(x, y) {

    k.out <- 0.0

    for (i.val in c(0,1))
    {
        for (j.val in c(0,1))
        {
            k.out <- k.out + exp(-((x - (0.25 + i.val/2.0))^2 + (y - (0.25 + j.val/2.0))^2)/(2 * sigma_k^2)) * k0
        }
    }

    k.out <- k.out + b * k0
}

# plotting the thing

# making all combinations of x and y coordinates
the.data <- as_tibble(expand.grid(x=gradient_x, y= gradient_y))

# apply the function to the tibble
the.data <- mutate(the.data,
       kval = k(x,y))

# plot the result
ggplot(data=the.data
       ,mapping=aes(x=x,y=y)) +
    geom_tile(mapping=aes(fill=kval))

ggsave("contour_carrying_capacity.pdf")


k0 <- 0.5
b <- 1.0

k.uni <- function(x) {

    k.out <- 0.0

    k.out <- k.out + exp(-((x - 0.5)^2) / (2 * sigma_k^2)) * k0

    k.out <- k.out + b * k0

}

the.data <- tibble(x=seq(0,1,0.01))

the.data <- mutate(the.data,
                   kval = k.uni(x))

ggplot(data=the.data
       ,mapping=aes(x=x,y=kval)) +
    geom_line()

ggsave("unidimensional_carrying_capacity.pdf")
