
tt <-  genefilter::rowttests(X, fac = factor(cl), tstatOnly = TRUE)
tt2 <- rowTstat(X, cl)

tt$statistic - tt2

system.time(rowTstatNull(expmat = X, classes = cl, nperm = 1000))


library(Category)
system.time(ttperm(X, fac = factor(cl), B = 1000))

A <- read_matrix('/Volumes/HCM/Human_LCM/UCSC_hg19/Counts_TPM/All-TPM-UCSChg19-Symbol.txt')

tt <- tpm(counts = A)

A <- A[,sample(1:ncol(A), 30)]

a1 <- basic_rank_norm(counts = A)

a2 <- basicRankNorm(A)
dim(a2)

boxplot(a1)
boxplot(a2)
diag(cor(a1, a2))

At <- tpm(counts = A, log = FALSE)

A1 <- t(t(apply(At, 2, rank, na.last="keep"))/(colSums(!is.na(A))+1))
A2 <- t(apply(A1, 1, rank, na.last="keep"))/(rowSums(!is.na(A1))+1)
summary(A2[,1])

# function to change gene names

load("data/geneInfo.RData")
counts <- readr::read_tsv("../UCSC_hg19/Counts_TPM/AllCounts-UCSC-hg19-knownGene.txt")

test <- any2symbol(as.character(counts$EntrezID))
geneids <- any2entrez(gene_ids = test)

tt <- c('weird', 'weirder', 'weirdest')

idx <- Position(function(i) any(test %in% i), geneInfo)
if(is.na(idx)) stop("Sorry, could not find any of the provided gene ids among the gene symbol
                    and Ensembl reference, respectively!")

entrez <- geneInfo[["ENTREZID"]]
names(entrez) <- geneInfo[[idx]]
trez[test])
sum(test %in% names(entrez))

data("geneInfo")
sum(test %in% geneInfo$SYMBOL)

sum(is.na(geneids))

geneids <- letters[1:25]

# GSEA --------------------------------------------------------------------
library(binilib)
reg <- readRDS("../Analyses/1_CUMC_LCM/Networks/Epi/ARACNE-AP/UCSChg19-Epi-TPMraw-ARACNEap-TF.rds")
gsigs <- read_matrix('../Analyses/External_Data_Sets/IPMN_Precursor/Log2FC-Matrix_Hiraoka.tsv')

gsig <- gsigs[,1]
head(gsig)
library(fgsea)

toi <- reg[['AEBP1']]

gs <- gsea_regulon(gsig, regulon = toi)
plot(gs)

tt <- viper::aREA(gsig, regulon = reg)$nes
tt[which.min(tt)]
head(tt)
ps <- names(toi$tfmode)[toi$tfmode >= 0]
ng <- names(toi$tfmode)[toi$tfmode < 0]

fgsea(pathways = list(PLAG1pos = ps), stats = signature, nperm = 10000)

tt <- gsea1T(signature, gS = ps)
tt$ES
tt <- gsea_null(tt, perm =10000)
tt$pval
plot(tt)


tt <- gsea_regulon(signature, toi)
tt <- gsea_null(tt)
plot(tt, plotSignature = TRUE, signatureNames = c('Precursor', 'PDA'),
     signatureType = 'Z-score')

str(tt)
tfm <- toi[['tfmode']]
pos_targ <- names(tfm)[tfm>=0]
neg_targ <- names(tfm)[tfm<0]

position_pos <- match(pos_targ, names(sort(signature, decreasing = TRUE))) # 1st decreasing for pos. targets
position_neg <- match(neg_targ, names(sort(signature, decreasing = FALSE))) # 2nd increasing for neg. targets
names(position_pos) <- pos_targ
names(position_neg) <- neg_targ

rlistr <- rank(signature)
rlistr <- rlistr[!(names(signature) %in% names(tfm))] # ranks of non-hits
rlistr <- sort(c(rlistr, position_pos, position_neg))
rlist <- signature[match(names(rlistr), names(signature))]
x <- which(names(rlist) %in% names(tfm))
nr <- sum(abs(rlist[x])^1)
nh <- length(rlist)-length(x)
es <- rep(-(1/nh),length(rlist))
es[x] <- abs(rlist[x])^1/nr
rs <- cumsum(es)
max(rs) + min(rs)

null <- gsea1T_null(signature, pos_targ)
sum(null$null_es > max(rs))
max(rs)/null$pos_null_es
plot(res)

tt <- gsea_null(res)

tt$NES_neg

ps <-

res <- gsea1T(signature, pos_targ)
class(res)
res <- gsea_null(res)
res$ES
res$nes

res
pdf('test.pdf', width = 5, height = 4, useDingbats = FALSE)
plot(res, plotSignature = TRUE, main = 'KLF5 targets',
     signatureNames = c('DMSO', 'PTC596'), signatureType = 'Z-score')
dev.off()

res <- gsea1T(ptc596_signature, names(ps), sorting = 'increasing')
pdf('test.pdf', width = 6, height = 4, useDingbats = FALSE)
plot(res, plotSignature = TRUE, main = 'Positive KLF5 targets',
     signatureNames = c('DMSO', 'PTC596'),
     signatureType = 'Z-score')
dev.off()


data("homologyInfo")
class(homologyInfo$Entrez)
class(names(homologyInfo$Entrez))

nms <- names(homologyInfo$Entrez)
human <- as.character(homologyInfo$Entrez)
names(human) <- nms
head(human)
class(nms)
class(human)

homologyInfo$Entrez <- as.character(homologyInfo$Entrez)


iris %>% mutate(bb = sqrt(Sepal.Length))
iris %<>% mutate(bb = sqrt(Sepal.Length))

library(binilib)
marinaLH <- readxl::read_excel('../AG_Rad/Data/Human-Signatures.xlsx', sheet = 3)
set1 <- readRDS('set1.rds')
set2 <- readRDS('set2.rds')

lhsig <- marinaLH$nes
names(lhsig) <- marinaLH$gene
gs_lc <- gsea2T(lhsig, set1, set2)
plot(gs_lc, legend = FALSE,
     plotSignature = TRUE,
     signatureNames = c('High-grade','Low-grade'),
     signatureType = 'Z-score')


xs <- rnorm(800)
s1 <- gl(n = 4, k = 200, length = 800, labels = paste('very long factor name', sample(letters[1:26], size = 4)))
s2 <- gl(n = 2, k = 100, length = 800, labels = LETTERS[25:26])

xs <- do.call('c', lapply(1:8, function(i) rnorm(n = 100, mean = i)))

split_violin(x = xs, s1 = s1, s2 = s2, rug = TRUE, rotate_xlabels = 60, main = 'Whats up')


# Gene sets ---------------------------------------------------------------
library(binilib)
fs <- list.files('~/Downloads/msigdb/')
nms <- str_extract(fs, '^\\w\\d{0,1}\\.\\w{2,3}')
path <- '~/Downloads/msigdb/'

l <- lapply(seq(length(fs)), function(i){
    hm <- strsplit(readLines(paste0(path, fs[i])), split = "\t")
    nm <- purrr::map_chr(hm, 1)
    hm <- purrr::map(hm, ~ .[3:length(.)])
    names(hm) <- nm
    hm
})
names(l) <- nms
msigdb <- l
save(msigdb, file = 'data/msigdb.RData')


# Mouse
data("homologyTable")

l2 <- lapply(l, function(i){

    lapply(i, function(j){

        tmp <- tibble(human = j) %>% inner_join(homologyTable$symbol, by = 'human')
        unique(return(tmp$mouse)) # in case the same mouse gene is returned multiple times

    })

})

# Check that names are still the same
identical(names(l), names(l2))
identical(unlist(lapply(l, function(i) names(i)), use.names = FALSE), unlist(lapply(l2, function(i) names(i)), use.names = FALSE))
# check lengths
hu <- unlist(lapply(l, function(i){
    vapply(i, function(j) length(j), FUN.VALUE = numeric(1))
}), use.names = FALSE)

mo <- unlist(lapply(l2, function(i){
    vapply(i, function(j) length(j), FUN.VALUE = numeric(1))
}), use.names = FALSE)

plot(hu, mo, pch = 16, cex = 0.7) # looks good ; numbers are quite different per se

summary(hu - mo) # 2 genes median difference
which(hu - mo == 357)
hu[145] - mo[145]

# check hallmark
l2$h.all$HALLMARK_TNFA_SIGNALING_VIA_NFKB
l$h.all$HALLMARK_TNFA_SIGNALING_VIA_NFKB

msigdbMouse <- l2
save(msigdbMouse, file = 'data/msigdbMouse.RData')

# Build hallmark founder set information ----------------------------------
fs <- list.files('~/Desktop/ADVOCATE_NewPaper/ADVOCATE-Paper/External_cohorts/Screening-Pathways/Hallmark_Founder_GeneSets/')

load("~/Desktop/ADVOCATE_NewPaper/ADVOCATE-Paper/External_cohorts/Screening-Pathways/Hallmark-Founder-GeneSets-Symbol.rda")



l <- map(hallmark_founder, ~ names(.))
hallmark_founder <- l
save(hallmark_founder, file = 'data/hallmarkFounder.RData')



sigs <- readRDS('../Laise-PDA-Activity/Data/bulkSubtypes-SignatureMatrix.rds')
sigs[1:5, 1:3]


tt <- specific_n(sigs)



# fgseaResults --------------------------------------------------------------------------------
library(DESeq2)
library(binilib)
dset <- readRDS('/home/carlo/Documents/Analyses/AG_Keller/Habringer_Stefan/Richard_2020-06-24/Data/Lewis-Eset-Prefilter-Counts.rds')
dset <- dset[,dset$Tissue == 'lymphoma']

dds <- DESeqDataSetFromMatrix(apply(exprs(dset), 2, round), pData(dset), ~ Genotype)
dds <- DESeq(dds)
resultsNames(dds)

res <- results(object = dds, contrast = c('Genotype', 'EMU_TCL1_WHIM', 'EMU_TCL1'), tidy = TRUE) %>%
    filter(!is.na(padj)) %>%
    tbl_df() %>% arrange(desc(stat))
res

gsig <- p2z(res$pvalue, res$stat)
names(gsig) <- res$row
fres_hm <- fgsea::fgsea(pathways = msigdbMouse$h.all,
                        stats = gsig, nperm = 10000)[order(NES, decreasing = TRUE), ]
fres_hm <- subset(fres_hm, padj < 0.1)

gs <- msigdbMouse$h.all[fres_hm$pathway]

tmp <- lapply(gs, function(i){
    gsea1T(signature = gsig, gS = i)
})

!any(names(gs) %in% fres_hm$pathway)

absrange <- cut(abs(gsig), breaks = quantile(abs(gsig)), include.lowest = TRUE, labels = FALSE)
length(absrange)
head(absrange)
table(absrange)

cols <- t(apply(col2rgb(c('firebrick2', 'royalblue')), 1, function(i){
    seq(i[1], i[2], length.out = length(gsig))
}))
dim(cols)

cols <- rbind(cols, als)

tmp <- cols[,1]
tmp
rgb(238.0000, 44, 44, maxColorValue = 255, alpha = 244.5688)

cols[1,]
als[1]

plot_fgseaRes(ges = gsig,
              geneSets = msigdbMouse$h.all[fres_hm$pathway],
              signatureNames = c('Emu TCL1 WHIM', 'Emu TCL1'),
              fgseaRes = fres_hm)

sigmat <- binilib::read_matrix('../Analyses/External_Data_Sets/Somerville_CellRep_2018/Data/log2FCMatrix-Comparisons.tsv')
reg <- readRDS('../Analyses/1_CUMC_LCM/Salmon_Gencode_hg38/Networks_Salmon/Epi/CUMC-Salmonhg38-DEtransform-Epi-TFPlusCoTF-AnimalTFDb.rds')


vpres <- viper::aREA(sigmat, reg)$nes
fdrs <- apply(vpres, 2, function(i) p.adjust(p = pnorm(abs(i), lower.tail = FALSE)*2, 'fdr'))

fdrs['TP63',]
plot_OneReg_MultSig(sigmat = sigmat,
                    tf = 'TP63', regulon = reg)

tt <- abs(seq(-255, 255, length.out = 15e3))

# Keller Lab Favorite Pathways ----------------------------------------------------------------
vst <- getVarianceStabilizedData(dds)
boxplot(split(vst['Dmtn', ], f = dds$Genotype))

reg <- readRDS('../../KellerLab_Pathways/KellerLab-FavoritePathway-Regulons.rds')
reg <- reg[-7]

vpres <- lapply(reg, function(i){
    viper(vst, i, method = 'ttest')
})
vpres <- do.call(rbind, vpres)

dd <- model.matrix(~ dds$Genotype)

library(limma)
colnames(dd) <- levels(dds$Genotype)
fit <- eBayes(lmFit(vpres, dd))
res <- topTable(fit, coef = 'EMU_TCL1_WHIM',
                number = nrow(fit),
                genelist = rownames(fit)) %>%
    tbl_df() %>% arrange(desc(t))
res
boxplot(split(vpres['HALLMARK_WNT_BETA_CATENIN_SIGNALING', ], f = dds$Genotype))





# Import --------------------------------------------------------------------------------------
dset <- readRDS('../Analyses/AG_Reichert/Data/RNAseq_Organoids_LogTransform-22584x54.rds')
colaps <- gsub('_\\d$', '', dset$sample)
vmat <- t(apply(Biobase::exprs(dset), 1, function(i) tapply(i, colaps, mean)))

purIST(emat = vmat, raw = TRUE)


# DEtransform -------------------------------------------------------------
gg <- readRDS('../Analyses/External_Data_Sets/Cell_line_expression/Klijn_Cells/Data/140625_Klijn_counts_coding.rds')
gg <- readRDS('../Analyses/Califano_Lab/Counts-Califano-NET.rds')

tt <- readRDS('../Analyses/External_Data_Sets/ICGC-Canada/EGA_Downloads/RNA_Seq/Counts_TPM/ICGC-PACA-CA-EsetCounts.rds')

exprs(gg)[1:5, 1:5]

mode(exprs(gg)[,1])



emat <- gg@assayData$exprs

vv <- DEtransform(emat = emat)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = emat,
                                      colData = data.frame(condition = rep(1, ncol(emat)),
                                                           row.names = colnames(emat)),
                                      design = ~ 1)
dds

