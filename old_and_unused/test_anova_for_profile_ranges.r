library("reshape")

num.profiles <- 15
profile.len <- 2000


#-------------- Create the data --------------

centers <- runif(num.profiles, min=-5, max=5)
centers <- centers - mean(centers)

data <- t( t( matrix( rnorm(num.profiles * profile.len, mean=0.0, sd=0.7), nrow=num.profiles) + centers ) + 1:profile.len - 20 )
stopifnot( all( dim(data) == c(num.profiles, profile.len)) )

print(rowMeans(data))
print(colMeans(data))

#-------------- Perform the regression --------------

stopifnot( all( dim(data) == c(num.profiles, profile.len)) )

meanProfile <- colMeans(data)
data.normalized <- t( t(data) - meanProfile)
print(data.normalized)
stopifnot( all( dim(data.normalized) == c(num.profiles, profile.len)) )

df <- melt(data)
colnames(df) <- c("Profile", "Pos", "Value")
df$Profile <- as.factor(df$Profile)

df.normalized <- melt(data.normalized)
colnames(df.normalized) <- c("Profile", "Pos", "NormalizedValue")
df.normalized$Profile <- as.factor(df.normalized$Profile)

#print(df.normalized)

print(rowMeans(data.normalized))
print(colMeans(data.normalized))

fit <- lm( NormalizedValue ~ Profile + 0, df.normalized )
print( coef(fit) )


print( cor( coef(fit), centers, method="pearson") )
print( coef(fit) - centers )


