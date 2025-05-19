
# 读取整个文件为一个字符串
text <- readLines("D:/8701/May/Poi-C60opt.txt", warn = FALSE)
text_combined <- paste(text, collapse = " ")

# 使用正则表达式提取所有 ":1.0:" 后面的数字
weights <- unlist(regmatches(text_combined, gregexpr("(?<=:1\\.0:)\\d+\\.\\d+", text_combined, perl = TRUE)))

# 转换为数值向量
weights_numeric <- as.numeric(weights)

# 输出为 R 的向量形式
cat("c(", paste(weights_numeric, collapse = ", "), ")\n", sep = "")
