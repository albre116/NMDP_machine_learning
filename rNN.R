library('RSNNS')
data("snnsData")
patterns <- snnsData$art2_tetra_med.pat
inputs <- patterns[, inputColumns(patterns)]
targets <- patterns[, outputColumns(patterns)]
patterns <- splitForTrainingAndTest(patterns, targets, ratio = 0.15)


model <- elman(patterns$inputsTrain, patterns$targetsTrain,
               size = c(8, 8), learnFuncParams = c(0.1), maxit = 500,
               inputsTest = patterns$inputsTest, targetsTest = patterns$targetsTest,
               linOut = FALSE)

plot(inputs, type = "l")
plot(targets[1:100], type = "l")
lines(model$fitted.values[1:100], col = "green")
plotIterativeError(model)
plotRegressionError(patterns$targetsTrain, model$fitted.values)
plotRegressionError(patterns$targetsTest, model$fittedTestValues)
hist(model$fitted.values - patterns$targetsTrain)
