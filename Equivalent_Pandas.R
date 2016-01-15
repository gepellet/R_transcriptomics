# Open the iris data
class(iris)

# Example of numerical indexing (rows)
iris_data = subset(iris,as.numeric(rownames(iris))%%10 == 0) 
iris_data = iris[seq(0,nrow(iris),10),]

# Example of semantic indexing (columns)
sepal_lengths = iris_data$Sepal.Length       # or you can use the attach function
sepal_lengths

# Example of an inline-query, featuring boolean indexing
# Extract all records where sepal_length exceeds the third quartile
boolean_mask = sepal_lengths > quantile(sepal_lengths, .75)
boolean_mask

long_sepals = iris_data[boolean_mask,]

# Create new data_frame
petal_data = data.frame(iris_data$Species,iris_data$Petal.Length,iris_data$Petal.Width)

petal_data = subset(iris_data, select=c('Species','Petal.Length','Petal.Width'))
sepal_data = subset(iris_data, select=c('Species','Sepal.Length','Sepal.Width'))

# example of grouping and group-based operations
# compute the per-species medians of remaining columns
species_medians = aggregate(. ~ Species, iris_data, FUN=median)


petal_medians = subset(species_medians, select=c('Species','Petal.Length','Petal.Width'))
result = merge(sepal_data,petal_medians, by='Species')
