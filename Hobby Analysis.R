library(ggplot2)
library(dplyr)
library(plyr)
library(gridExtra)
library(RColorBrewer)
library(ggthemes)
library(gridExtra)
library(qgraph)
library(mnormt)
library(igraph)
library(bootnet) 
library(dplyr)
library(NetworkComparisonTest)
library(mgm)
library(VIM)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(cluster)
library(randomForest)
library(data.table)


Data <- read.csv("./R/new_data.csv")
Data2 <- Data[,2:33] # Select only variables

corMat <- cor_auto(Data2) # Correlate data
names<-c("1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
         "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32")

Graph_pcor <- qgraph( corMat, graph = "pcor", layout = "spring",
                      tuning = 0.25,sampleSize = nrow(Data2), 
                      legend.cex = 0.35, vsize = 5,esize = 10,
                      posCol = "#003399", negCol = "#FF9933",vTrans = 200,
                      nodeNames = colnames(corMat), labels = names)

Layout <- averageLayout(Graph_pcor)

Graph_pcor <- qgraph(corMat, graph = "pcor", layout = Layout, threshold = "bonferroni",
                     sampleSize = nrow(Data2), alpha = 0.05, 
                     legend.cex = 0.35, vsize = 5,esize = 10,
                     posCol = "#003399", negCol = "#FF9933",vTrans = 200,
                     nodeNames = colnames(corMat), labels = names)

iGraph_pcor <- as.igraph(Graph_pcor, attributes = TRUE)


Graph_lasso <- qgraph(corMat, graph = "glasso", layout = Layout, tuning = 0.25,
                      sampleSize = nrow(Data2), 
                      legend.cex = 0.35, vsize = 5,esize = 10,
                      posCol = "#003399", negCol = "#FF9933",vTrans = 200,
                      nodeNames = colnames(corMat), labels = names)


############## Centrality analysis
centRes <- centrality(Graph_lasso)
qplot(centRes$OutDegree, geom = "histogram") + geom_histogram(bins=10) + theme_minimal(12) +
  labs(title="Node Strength") + xlab("Strength")

qplot(centRes$Closeness, geom = "histogram") + geom_histogram(bins=20) + theme_minimal(12) +
  labs(title="Closeness Centrality") + xlab("Closeness")

qplot(centRes$Betweenness, geom = "histogram") + geom_histogram(bins=20) + theme_minimal(12) +
  labs(title="Betweenness Centrality") + xlab("Betweenness")

centralityPlot(Graph_lasso)
centralityPlot(GGM = list(unregularized = Graph_pcor, regularized = Graph_lasso))

# Layout <- averageLayout(Graph_pcor,Graph_lasso)
layout(t(1:2))
qgraph(corMat, graph = "pcor", layout = Layout, threshold = "bonferroni",
       sampleSize = nrow(Data2), minimum = 0,
       cut = 0.15, maximum = 1, details = TRUE,
       title = "Partial correlations", 
       legend.cex = 0.35, vsize = 5,esize = 10,
       posCol = "#003399", negCol = "#FF9933",vTrans = 200, nodeNames = colnames(corMat), labels = names)

qgraph(corMat, graph = "glasso", layout = Layout, tuning = 0.25,
       sampleSize = nrow(Data2), minimum = 0,
       cut = 0.15, maximum = 1, details = TRUE,
       title = "LASSO regularization", 
       legend.cex = 0.35, vsize = 5,esize = 10,
       posCol = "#003399", negCol = "#FF9933",vTrans = 200, nodeNames = colnames(corMat), labels = names)


layout(t(1:1))
g = as.igraph(Graph_lasso, attributes=TRUE)
sgc <- spinglass.community(g)
sgc$membership

group.worktrap<- list(c(1,2,9,10,11,13,14), c(9,10,14), c(19,20,23,24,25,26,27,30,31,32),
                       c(2,11,17,18,21,22,29), c(4,5,6,7,16,28))
Graph_lasso <- qgraph(corMat, graph = "glasso", layout = "spring", tuning = 0.25,
                      sampleSize = nrow(Data2), 
                      legend.cex = 0.35, vsize = 5,esize = 10,
                      posCol = "#003399", negCol = "#FF9933",vTrans = 200,groups=group.spinglass,
                      color=c("red", "orange", "white", "blue", "green"),
                      nodeNames = colnames(corMat), labels = names)



data1 <- slice(Data2, c(1:505))    # n = 505
data2 <- slice(Data2, c(506:1010))   # n = 506

### similarity: visual
network1 <- estimateNetwork(data1, default="EBICglasso")
network2 <- estimateNetwork(data2, default="EBICglasso")
layout(t(1:2))
graph1 <- plot(network1, cut=0)
graph2 <- plot(network2, cut=0)

### similarity: statistical
cor(network1$graph[lower.tri(network1$graph)], network2$graph[lower.tri(network1$graph)], method="spearman")

cor(centrality(network1)$InDegree, centrality(network2)$InDegree)


layout(t(1:1))
Data <- read.csv("./R/data2.csv")

Data2 <- subset(Data, Data$id!=1)[,1:32]
Data2_ <- subset(Data, Data$id==1)[,1:32]

corMat2 <- corMat
corMat2[1,2] <- corMat2[2,1] <- corMat[1,2] + 0.0001 # we change the correlation matrix just a tiny little bit
L1 <- averageLayout(corMat, corMat2)

#making model
fit1 <- mgm(Data2, type=rep('g',32), level=rep(1,32), k=2, lambdaSel = 'EBIC', 
           lambdaGam = 0.25)

Data2 <- subset(Data, Data$id!=2)[,1:32]
Data2_ <- subset(Data, Data$id==2)[,1:32]

corMat2 <- corMat
corMat2[1,2] <- corMat2[2,1] <- corMat[1,2] + 0.0001 # we change the correlation matrix just a tiny little bit
L1 <- averageLayout(corMat, corMat2)

#making model
fit2 <- mgm(Data2, type=rep('g',32), level=rep(1,32), k=2, lambdaSel = 'EBIC', 
           lambdaGam = 0.25)

Data2 <- subset(Data, Data$id!=3)[,1:32]
Data2_ <- subset(Data, Data$id==3)[,1:32]

corMat2 <- corMat
corMat2[1,2] <- corMat2[2,1] <- corMat[1,2] + 0.0001 # we change the correlation matrix just a tiny little bit
L1 <- averageLayout(corMat, corMat2)

#making model
fit3 <- mgm(Data2, type=rep('g',32), level=rep(1,32), k=2, lambdaSel = 'EBIC', 
           lambdaGam = 0.25)

Data2 <- subset(Data, Data$id!=4)[,1:32]
Data2_ <- subset(Data, Data$id==4)[,1:32]

corMat2 <- corMat
corMat2[1,2] <- corMat2[2,1] <- corMat[1,2] + 0.0001 # we change the correlation matrix just a tiny little bit
L1 <- averageLayout(corMat, corMat2)

#making model
fit4 <- mgm(Data2, type=rep('g',32), level=rep(1,32), k=2, lambdaSel = 'EBIC', 
           lambdaGam = 0.25)

#plotting
# fit_plot <- qgraph(fit$pairwise$wadj, cut=0, legend.cex = 0.35,
#                    layout = L1, edge.color = fit$pairwise$edgecolor,
#                    nodeNames = colnames(corMat), labels = names )

#predict
pred_obj1 <- predict(object = fit1, 
                    data = Data2, 
                    errorCon = 'R2')
pred_obj_1 <- predict(object = fit1, 
                    data = Data2_, 
                    errorCon = 'R2')
pred_obj2 <- predict(object = fit2, 
                    data = Data2, 
                    errorCon = 'R2')
pred_obj_2 <- predict(object = fit2, 
                     data = Data2_, 
                     errorCon = 'R2')
pred_obj3 <- predict(object = fit3, 
                    data = Data2, 
                    errorCon = 'R2')
pred_obj_3 <- predict(object = fit3, 
                     data = Data2_, 
                     errorCon = 'R2')
pred_obj4 <- predict(object = fit4, 
                    data = Data2, 
                    errorCon = 'R2')
pred_obj_4 <- predict(object = fit4, 
                     data = Data2_, 
                     errorCon = 'R2')
#error_train
error1 = join(x = pred_obj1$errors, y = pred_obj2$errors,
               by = "Variable", type = "left")
error2 = join(x = error1, y = as.data.frame(pred_obj3$errors), 
               by = "Variable", type = "left")
error3 = join(x = error2, y = as.data.frame(pred_obj4$errors),
               by = "Variable", type = "left")

#error
error_1 = join(x = pred_obj_1$errors, y = pred_obj_2$errors,
              by = "Variable", type = "left")
error_2 = join(x = error_1, y = as.data.frame(pred_obj_3$errors), 
              by = "Variable", type = "left")
error_3 = join(x = error_2, y = as.data.frame(pred_obj_4$errors),
              by = "Variable", type = "left")

row_mean <- apply(error3[,2:5], 1, mean)
row_mean_ <- apply(error_3[,2:5], 1, mean)
write.csv(row_mean_, "R/temp.csv")


#plotting with error in pie chart
fit_plot2 <- qgraph(fit$pairwise$wadj, cut=0, pie = pred_obj$error[,2],legend.cex = .35,
                    layout = L1, edge.color = fit$pairwise$edgecolor, nodeNames = colnames(corMat),
                    labels = names)


#missing_value

data <- read.csv("responses.csv")
data_using <- data[,32:63]

png("missing_.png",width = 3000,height = 2000,res = 400)
missing_data = aggr(data_using, prop = T, number = T)
dev.off()

row <- dim(data_using)[1]
colomn <- dim(data_using)[2]
for (i in 1:row)
  if (sum(is.na(data_using[i,])) < colomn * 0.15){
    for (j in 1:colomn) {
      data_using[,j][is.na(data_using[,j])] = median(data_using[,j], na.rm = T)
    } 
  }else{
    data_using[-i,]
  }
write.csv(data_using, file = "new_data.csv")

#visual_analysis

data <- read.csv("new_data.csv")
data <- data[,2:33]

PC <- as.data.frame(table(data$PC))
Mathematics <- as.data.frame(table(data$Mathematics))
Reading <- as.data.frame(table(data$Reading))
counts <- rbind(PC[,2],Mathematics[,2],Reading[,2])

png("p1.png",width = 3000,height = 2000,res = 400)
barplot(counts, beside = TRUE,
        ylim = c(0,500),
        xlab = "Preferences", ylab = "Amount of People",
        legend.text = c("PC","Mathematics","Reading"), 
        col = brewer.pal(12,"Set3")[1:3],
        names.arg = c(1,2,3,4,5))
dev.off()

data_1 <- data[which(data[,4]=="1"),]
data_2 <- data[which(data[,4]=="2"),]
data_3 <- data[which(data[,4]=="3"),]
data_4 <- data[which(data[,4]=="4"),]
data_5 <- data[which(data[,4]=="5"),]
PC_1 <- as.data.frame(table(data_1$PC))
PC_2 <- as.data.frame(table(data_2$PC))
PC_3 <- as.data.frame(table(data_3$PC))
PC_4 <- as.data.frame(table(data_4$PC))
PC_5 <- as.data.frame(table(data_5$PC))
count <- rbind(PC_1[,2],PC_2[,2],PC_3[,2],PC_4[,2],PC_5[,2])

png("p2.png",width = 3000,height = 2000,res = 400)
par(mai=c(1,1,1,1.5))
barplot(count, 
        xlab = "Preferences of Mathematics", ylab = "Preferences of PC",
        legend.text = c("1","2","3","4","5"), 
        col = brewer.pal(12,"Set3")[1:5],
        names.arg = c(1,2,3,4,5),
        args.legend = c(x=9,y=200),xpd = TRUE)
dev.off()

#correlation_analysis

data <- read.csv("new_data.csv")

png("corr.png",width = 3000,height = 2000,res = 400)
data_corr <- cor(data)
corrplot(data_corr, type = "upper")
dev.off()

corrplot.mixed(data_corr)

#PCA

data <- read.csv("new_data.csv")
data <- data[,2:33]

data.pc <- PCA(data, graph = FALSE)

png("screeplot.png",width = 3000,height = 2000,res = 400)
fviz_screeplot(data.pc, addlabels=TRUE, ncp=14)
dev.off()

var <- get_pca_var(data.pc)
var
eig.val <- get_eigenvalue(data.pc)
eig.val

png("dim1.png",width = 3000,height = 2000,res = 400)
fviz_contrib(data.pc, choice = "var", axes = 1, top = 10)
dev.off()

png("dim2.png",width = 3000,height = 2000,res = 400)
fviz_contrib(data.pc, choice = "var", axes = 2, top = 10)
dev.off()

png("cos2.png",width = 3000,height = 2000,res = 400)
fviz_pca_var(data.pc, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE)
dev.off()


x <- read.table("new_data.csv",sep=",",header=T)
x1 <- x[,2:33]
y <- t(x1)
#write.table(y,"data1.csv",row.names=FALSE,col.names=TRUE,sep=",")
km <- kmeans(y, 5)
print(km)

z <- data.frame(x1)
k = 4
z$id <- sample(1:k, nrow(z), replace = TRUE)
list <- 1:k
z
#write.table(z,"data2.csv",row.names=FALSE,col.names=TRUE,sep=",")
# 每次迭代的预测用数据框，测试用数据框
# the folds
prediction <- data.table()
prediction1 <- data.table()
testsetCopy <- data.table()
testsetCopy1 <- data.table()
z <- z[,c(4,1:3,5:33)]
z
#k层的函数
for(i in 1:k){
  # 删除id为i的行，创建训练集
  # 选id为i的行，创建训练集
  trainingset <- subset(z, id %in% list[-i])
  testset <- subset(z, id %in% c(i))
  #运行一个随机森林模型
  mymodel <- randomForest(trainingset$Mathematics ~ ., data = trainingset, ntree = 100, importance = TRUE)
  #去掉回应列1, Sepal.Length
  temp <- as.data.frame(predict(mymodel, testset[,-1]))
  temp1 <- as.data.frame(predict(mymodel, trainingset[,-1]))
  # 将迭代出的预测结果添加到预测数据框的末尾
  prediction <- rbind(prediction, temp)
  prediction1 <- rbind(prediction1, temp1)
  # 将迭代出的测试集结果添加到测试集数据框的末尾
  # 只保留Sepal Length一列
  testsetCopy <- rbind(testsetCopy, as.data.frame(testset[,1]))
  testsetCopy1 <- rbind(testsetCopy1, as.data.frame(trainingset[,1]))
}
importance1 = importance(x = mymodel, type = 1)
importance2 = importance(x = mymodel, type = 2)
importance1
importance2
varImpPlot(mymodel) #自变量重要性排名

# 将预测和实际值放在一起
result <- cbind(prediction, testsetCopy[, 1])
result1 <- cbind(prediction1, testsetCopy1[, 1])
names(result) <- c("Predicted", "Actual")
names(result1) <- c("Predicted", "Actual")
result$Difference <- (result$Actual - result$Predicted)^2 
result1$Difference <- (result1$Actual - result1$Predicted)^2 
result
result1
summary(result$Difference) #测试集均方误差
summary(result1$Difference) #训练集均方误差



