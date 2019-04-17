f##########################################################################
#########################Problem Solving in R#############################
##########################################################################


##################Exploring Code########################
#figure out what each line of code does
set.seed(2212)

df <- data.frame(replicate(8, round(rnorm(50, 0, 4), 2)))

names(df) <- paste0('run', 1:8)

x <- apply(unname(df[,3:6]), 1, rep, times=4)



##################Types of Errors########################

#######Syntax Errors######
#find and correct all of the syntax errors in these for lines of code

dfunif <- data.frame(replicate(100, runif(20, 0, 10)))

x <- lapply(1:110, function(i), mean(dfunif[,i]) )

shapiro.test(x) #a test of normality
shapiro.text(runif(100, 0, 10, FALSE))


######Runtime Errors######
#find and correct the runtime errors in this section
vector1 <- 1:10
vector2 <- c('cat','whale','horse','owl')

addfive <- function(v) {
    v2 <- v+5
    return(v3)
}

addfive(vector1)
addfive(vector2)





######Logic Errors########

##Problem 1: extract the numbers in random vector that are between 25 and 75
    #save those numbers as a variable named middle

set.seed(32)
randomvector <- sample(1:100, 50)

#####################
#first try
middle <- if (randomvector > 25 & randomvector > 75) {
    randomvector
} else {

}

#####################
#second try
indexes <- which(randomvector>25 & randomvector>75) #identify which elements the ones we want
middle <- randomvector[indexes] #select those elements


###############################################################
##Problem 2: Find the mean of v. Then add 10 to every element of v and find the
            #mean again.

#run the following 2 lines only once
set.seed(292929)
v <- sample(1:100, 20)
############

mean(v) #should be 47.95

add10 <- sapply(1:10, function(i) {
    v[i]+10
})

mean(add10) #should be 57.95





