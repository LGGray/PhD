library(ggplot2)

studies <- c('AD_GSE147424', 'MS_GSE193770', 'pSS_GSE157278', 'SLE_SDY997', 'UC_GSE125527', 'UC_GSE182270')
cells <- c('Regulatory.T.cells', 'Tem.Trm.cytotoxic.T.cells', 'Tcm.Naive.helper.T.cells', 'Tem.Effector.helper.T.cells')
colours <- c('#FBE799', '#E6CA89', '#F7D6A0', '#E3B187', '#F3BF87', '#EEAF92')

LR <- lapply(cells, function(x){
    metrics <- lapply(studies, function(y){
        if(file.exists(paste0(y, '/exp.matrix/metrics/logit_metrics_', x, '.common.csv'))){
            df <- read.csv(paste0(y, '/exp.matrix/metrics/logit_metrics_', x, '.common.csv'))
            nFeatures <- length(read.delim(paste0(y, '/ML.models/features/logit_model_', x, '.common.txt'))$Features)
            df$nFeatures <- nFeatures
            return(df)
        } else{data.frame(Accuracy=0, Precision=0, Recall=0, F1=0, AUC=0, nFeatures=0)}
    })
    names(metrics) <- studies
    return(metrics)
})

LR.df <- lapply(LR, function(x){
    df <- dplyr::bind_rows(x, .id='study')
    return(df)
})
names(LR.df) <- cells

z=1
pdf(paste0('../ML.figures/LR.', names(LR.df)[[z]], '.common.pdf'))
ggplot(LR.df[[z]], aes(x=study, y=F1, fill=study)) + 
  geom_col() +
  scale_fill_manual(values=colours) +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = 'none') +
  geom_text(aes(label=nFeatures), vjust=-1) +
  ggtitle(paste0('LR: ', gsub('\\.', ' ', names(LR.df)[[z]]))) +
  xlab('') + ylab('F1 score')
dev.off()
#--------------------------------------------

RF <- lapply(cells, function(x){
    metrics <- lapply(studies, function(y){
        if(file.exists(paste0(y, '/exp.matrix/metrics/RF_metrics_', x, '.common.csv'))){
            df <- read.csv(paste0(y, '/exp.matrix/metrics/RF_metrics_', x, '.common.csv'))
            nFeatures <- length(read.delim(paste0(y, '/ML.models/features/RF_model_', x, '.common.txt'))$Features)
            df$nFeatures <- nFeatures
            return(df)
        } else{data.frame(Accuracy=0, Precision=0, Recall=0, F1=0, AUC=0, nFeatures=0)}
    })
    names(metrics) <- studies
    return(metrics)
})

RF.df <- lapply(RF, function(x){
    df <- dplyr::bind_rows(x, .id='study')
    return(df)
})
names(RF.df) <- cells

z=1
pdf(paste0('../ML.figures/RF.', names(RF.df)[[z]], '.common.pdf'))
ggplot(RF.df[[z]], aes(x=study, y=F1, fill=study)) + 
  geom_col() +
  scale_fill_manual(values=colours) +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = 'none') +
  geom_text(aes(label=nFeatures), vjust=-1) +
  ggtitle(paste0('RF: ', gsub('\\.', ' ', names(RF.df)[[z]]))) +
  xlab('') + ylab('F1 score')
dev.off()
#--------------------------------------------

SVM <- lapply(cells, function(x){
    metrics <- lapply(studies, function(y){
        if(file.exists(paste0(y, '/exp.matrix/metrics/SVM_metrics_', x, '.common.csv'))){
            df <- read.csv(paste0(y, '/exp.matrix/metrics/SVM_metrics_', x, '.common.csv'))
            nFeatures <- length(read.delim(paste0(y, '/ML.models/features/SVM_model_', x, '.common.txt'))$Features)
            df$nFeatures <- nFeatures
            return(df)
        } else{data.frame(Accuracy=0, Precision=0, Recall=0, F1=0, AUC=0, nFeatures=0)}
    })
    names(metrics) <- studies
    return(metrics)
})

SVM.df <- lapply(SVM, function(x){
    df <- dplyr::bind_rows(x, .id='study')
    return(df)
})
names(SVM.df) <- cells

z=1
pdf(paste0('../ML.figures/SVM.', names(SVM.df)[[z]], '.common.pdf'))
ggplot(SVM.df[[z]], aes(x=study, y=F1, fill=study)) + 
  geom_col() +
  scale_fill_manual(values=colours) +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = 'none') +
  geom_text(aes(label=nFeatures), vjust=-1) +
  ggtitle(paste0('SVM: ', gsub('\\.', ' ', names(SVM.df)[[z]]))) +
  xlab('') + ylab('F1 score')
dev.off()
#--------------------------------------------

GBM <- lapply(cells, function(x){
    metrics <- lapply(studies, function(y){
        if(file.exists(paste0(y, '/exp.matrix/metrics/GBM_metrics_', x, '.common.csv'))){
            df <- read.csv(paste0(y, '/exp.matrix/metrics/GBM_metrics_', x, '.common.csv'))
            nFeatures <- length(read.delim(paste0(y, '/ML.models/features/GBM_model_', x, '.common.txt'))$Features)
            df$nFeatures <- nFeatures
            return(df)
        } else{data.frame(Accuracy=0, Precision=0, Recall=0, F1=0, AUC=0, nFeatures=0)}
    })
    names(metrics) <- studies
    return(metrics)
})

GBM.df <- lapply(GBM, function(x){
    df <- dplyr::bind_rows(x, .id='study')
    return(df)
})
names(GBM.df) <- cells

z=4
pdf(paste0('../ML.figures/GBM.', names(GBM.df)[[z]], '.common.pdf'))
ggplot(GBM.df[[z]], aes(x=study, y=F1, fill=study)) + 
  geom_col() +
  scale_fill_manual(values=colours) +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = 'none') +
  geom_text(aes(label=nFeatures), vjust=-1) +
  ggtitle(paste0('GBM: ', gsub('\\.', ' ', names(GBM.df)[[z]]))) +
  xlab('') + ylab('F1 score')
dev.off()