#' @title Compare different machine learning models
#' @description
#' Trains a machine learning model on a single or different
#' train data frames using caret train with given methods,
#' tests the model on a test data frames,
#' and returns a data frame with names of passed parameters
#' and results of model training
#' (Some methods may not be supported)
#'
#'
#' @param train_data a list containing one or multiple data frames
#' @param test_data a single data frame
#' @param caret_methods a string with one or vector with multiple caret's train methods
#'
#' @export
compare = function(train_data, test_data, caret_method) {
  if(!(is.list(train_data) && !is.data.frame(train_data))) stop("Train data must be in list")
  if(is.vector(test_data)) stop("Test data must be passed as a single variable")

  df_len = length(train_data) * length(caret_method)
  df = data.frame(train_data = NA,
                  test_data = NA,
                  method = NA,
                  accuracy = NA,
                  precision = NA,
                  recall = NA,
                  f1 = NA,
                  tmp = c(1:df_len))[,-8]

  count = 0
  train_names = unlist(stringr::str_split(stringr::str_remove(stringr::str_remove(
    stringr::str_remove(trimws(deparse(substitute(train_data))), "list"), "[)]"), "[(]"), ", "))
  train_names = sapply(train_names, function (i) trimws(i), USE.NAMES = F)
  train_names = sapply(train_names, function (i) stringr::str_remove(i, "[,]"), USE.NAMES = F)

  for (i in 1:length(train_data)) {
    for (j in 1:length(caret_method)) {
      count = count + 1
      print(paste0("Now testing on ", train_names[i], " using ", caret_method[j], "..."))
      if(caret_method[j] == "rpart") {
        model = rpart::rpart(sentiment ~ ., data = train_data[[i]], method = "class")
        matr = caret::confusionMatrix(predict(model, na.omit(test_data), type="class"), na.omit(test_data)$sentiment)
        df$train_data[count] = train_names[i]
        df$method[count] = caret_method[j]
        df$accuracy[count] = matr$overall[1]
        df$precision[count] = matr$byClass[5]
        df$recall[count] = matr$byClass[6]
        df$f1[count] = matr$byClass[7]

        rpart.plot::rpart.plot(model, box.palette = "PuPu", shadow.col = "#EDB3C8", type = 5)
        title(main = paste0("train: ", train_names[i], ", \ntest: ", deparse(substitute(test_data))),
              cex.main=1.4, adj=0, col.main="#b33d70")
      }
      else if (caret_method[j] == "glm") {
        model = glm(sentiment ~ ., data=train_data[[i]], family="binomial")
        matr = as.matrix(table("predicted" = ifelse(predict(model, na.omit(test_data), type="response") > 0.5, "positive", "negative"),
                               "actual" = na.omit(test_data)$sentiment))
        df$train_data[count] = train_names[i]
        df$method[count] = caret_method[j]
        df$accuracy[count] = sum(diag(matr))/sum(matr)
        df$precision[count] = matr[1, 1]/sum(matr[1,])
        df$recall[count] = matr[1, 1]/sum(matr[, 1])
        df$f1[count] = (2*df$precision[count]*df$recall[count])/(df$precision[count] + df$recall[count])
      }
      else {
        model = caret::train(sentiment ~ ., data = train_data[[i]], method = caret_method[j])
        matr = caret::confusionMatrix(predict(model, na.omit(test_data)), na.omit(test_data)$sentiment)
        df$train_data[count] = train_names[i]
        df$method[count] = caret_method[j]
        df$accuracy[count] = matr$overall[1]
        df$precision[count] = matr$byClass[5]
        df$recall[count] = matr$byClass[6]
        df$f1[count] = matr$byClass[7]
      }
      print("Done")
    }
  }
  df$test_data = deparse(substitute(test_data))
  return(df)
}
