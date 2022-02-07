library(Biostrings)
library(ggplot2)
library(cowplot)
library(dplyr)
setwd("~/cfo155/whole_genome_lrTAPS")
id <- "bc1001_cmb"
# Read Length ----
dat <- read.table(paste0("fastq/", id, ".read_all.txt"), header=FALSE)
n50 <- N50(dat$V1)
dat <- read.table(paste0("fastq/", id, ".read_len.txt"), header=TRUE)
colnames(dat) <- c("len","nread")
dat<- dat[order(dat$len),]
p <- ggplot(dat, aes(x=len,y=nread)) +
  geom_bar(stat="identity") +
  annotate("text", y = 200000, x = 1000, label = paste0("N50: ", n50, " bp")) +
  xlab("read length") +
  theme_light() + xlim(0,7500) +
  theme(legend.position = "bottom")
ggsave(paste0("plots/",id,".read_len.pdf"), p, width = 6, height = 4)


# QC ----
dat<- read.table(paste0("sta/", id, ".all.sta"), header=TRUE, stringsAsFactors = FALSE)
dat$conversion <- gsub(";.*","",dat$kb4_meth) %>% as.character() %>% as.numeric() 
dat$fp <-  gsub(";.*","",dat$kb4_unmeth) %>% as.character() %>% as.numeric() 


p1 <- ggplot(dat, aes(x = id, y = conversion*100)) + geom_bar(stat = "identity", fill = "orangered3", col = "orangered3") + 
  geom_text(aes(label=round(conversion*100, 2)), vjust=-0.4) +
  theme_light() + xlab("") + ylab("mCG conversion on 4kb (%)") + theme(axis.text = element_text(color = "black")) + ylim(0,100)
p2 <- ggplot(dat, aes(x = id, y = fp*100)) + geom_bar(stat = "identity", fill = "steelblue", col = "steelblue") + 
  geom_text(aes(label=round(fp*100, 2)), vjust=-0.4) +
  theme_light() + xlab("") + ylab("false positive on unmodified cg(4kb) (%)") + theme(axis.text = element_text(color = "black"))  + ylim(0,0.5)
p <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
ggsave(paste0("plots/",id,"_conversion_fp.pdf"),p,width=3,height = 4)

dat$genome_meth <- gsub(";.*","",dat$genome_meth) %>% as.character() %>% as.numeric() 
p1 <- ggplot(dat, aes(x = id, y = genome_meth*100, fill=id)) + geom_bar(stat = "identity") + 
  geom_text(aes(label=round(genome_meth*100, 2)), vjust=-0.4) +
  theme_light() + xlab("") + ylab("genome methylation (%)") + theme(axis.text = element_text(color = "black"), legend.position = "none") +
  scale_fill_manual(values=brewer.pal(2,"Dark2"))

dat$pmap <- (dat$nq20map_bam)/(dat$nfwd_read + dat$nrev_read)
dat$pq1map <- (dat$nq1map_bam)/(dat$nfwd_read + dat$nrev_read)
p2 <- ggplot(dat, aes(x = id, y = pmap*100)) + geom_bar(stat = "identity", fill = "steelblue", col = "steelblue") + 
  geom_text(aes(label=round(pmap*100, 2)), vjust=-0.4) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  theme_light() + xlab("") + ylab(" Q20 mapping rate (%)") + theme(axis.text = element_text(color = "black"))

p <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
ggsave(paste0("plots/",id,"genome_sta.pdf"),p,width=3,height = 4)

dat$prm_map <- dat$nrm_bam / dat$nmap_bam
p1 <- ggplot(dat,aes(x=id,y=(1-prm_map)*100))+ 
  geom_bar(stat="identity")+
  geom_text(aes(label=round((1-prm_map)*100, 2)), vjust=-0.1) +
  ylab("PCR duplicates%") +
  xlab("") +
  theme(legend.position = "bottom")
p1
ggsave(paste0("plots/",id,"_pcr_rm.pdf"),p1,width = 1.5, height = 4)




# comparison with short-read taps ----
taps <- read_delim("shortread_taps/taps_pub.average.methratio.txt", delim = "\t") %>% as.data.frame()
lrtaps <- read_delim(paste0("meth/", id,".uniflag.md_CpG.merge.bedGraph"), delim = "\t") %>% as.data.frame()
all_taps <- merge(taps,lrtaps, by=c("chr","start"))
colnames(all_taps) <- c("chr","start","taps_mC","taps_aC","taps_meth","lrtaps_mC","lrtaps_aC","lrtaps_meth")
e14_snp <- read_delim("resource/e14.cgsnp.bed",delim="\t", col_names=FALSE)
colnames(e14_snp) <- c("chr","start","end","snp")
all_taps <- merge(all_taps,e14_snp,by=c("chr","start"),all.x=TRUE)
all_taps <- all_taps[!complete.cases(all_taps$end),-c(ncol(all_taps)-1,ncol(all_taps))]
cg_states <- read_delim("resource/mm9_genome.cg.chrHMM.info", delim = "\t", col_names =FALSE) %>% as.data.frame()
all_taps <- merge(all_taps, cg_states, by.x=c("chr","start"), by.y=c("X1","X2") )

# methylation correlation ----
for(depth in c(5, 20)){
  seltaps <- all_taps[all_taps$taps_aC>=depth & all_taps$lrtaps_aC>=depth & complete.cases(all_taps) & all_taps$X3!=".",]
  cors <- cor(seltaps %>% dplyr::select(taps_meth, lrtaps_meth))
  pdf(paste0("plots/",id,"_depth.",depth,".scatter.pdf"), width = 3.2, height = 3.5)
  smoothScatter(x=seltaps$taps_meth,
                y=seltaps$lrtaps_meth,
                xlab=paste0("TAPS (cor:", round(cors[2,1],2),")"),
                ylab= "lrTAPS",
                xlim=c(0,1),
                ylim=c(0,1))
  dev.off()
  print(dim(seltaps))
  print(cors)
  
  cors1 <- c ("all", cors[2,1])%>% t()%>% as.data.frame(); cors1$V2 <- as.numeric(as.character(cors1$V2))
  cors2 <- as.table(by(seltaps[,c(5,8)], seltaps$X3, function(x) {cor(x$taps_meth, x$lrtaps_meth)})) %>% as.data.frame()
  colnames(cors1) <- c("feature","cor")
  colnames(cors2) <- c("feature","cor")
  allcors <- rbind(cors1,cors2)
  write.table(allcors,paste0("plots/", id,"_depth.",depth,".cor.txt"),sep="\t", col.names =TRUE, row.names = FALSE, quote = FALSE)
  
  pdf(paste0("plots/",id,"_depth.",depth,".feature.scatter.pdf"), width = 3.2, height = 3.5)
  for( i in sort(unique(seltaps$X3))[c(1,8:15,2:7)]){
    selcpg <- seltaps[seltaps$X3==i, ]
    smoothScatter(x=selcpg$taps_meth,
                  y=selcpg$lrtaps_meth,
                  xlab=paste0("TAPS (cor:", round(cor(selcpg$taps_meth, selcpg$lrtaps_meth),2),")\n",i),
                  ylab= "lrTAPS",
                  xlim=c(0,1),
                  ylim=c(0,1))
  }
  dev.off()
}





# covered sites ----

all_taps$taps <- ifelse(all_taps$taps_aC >= 5, 1, 0); all_taps$taps[is.na(all_taps$taps)] <- 0
all_taps$lrtaps <- ifelse(all_taps$lrtaps_aC >= 5, 1, 0); all_taps$lrtaps[is.na(all_taps$lrtaps)] <- 0
options(scipen = 100)
write.table(all_taps,paste0("meth/", id,"_shortread_taps.info.txt"),sep="\t", col.names =TRUE, row.names = FALSE, quote = FALSE)
mean(all_taps$lrtaps_aC[!is.na(all_taps$lrtaps_aC)])

pdf(paste0("plots/", id,"_all.coveredC.venn.pdf"), width = 5, height = 4)
num_taps <- sum(all_taps$taps=="1",na.rm=TRUE)
num_lrtaps <- sum(all_taps$lrtaps=="1",na.rm=TRUE)
num_taps_lrtaps <- sum(all_taps$taps=="1"&all_taps$lrtaps=="1",na.rm=TRUE)
draw.pairwise.venn(
  area1     = num_taps,
  area2     = num_lrtaps,
  cross.area       = num_taps_lrtaps,
  category  = c('taps', 'lrtaps'),
  fill      = c("#446db4","#23b177"),
  cat.col   = c("#446db4","#23b177"),
  euler.d = TRUE,
  scaled    = TRUE
)
dev.off()

all_taps[all_taps$taps!="1"&all_taps$lrtaps=="1",] %>% head()

depth_sta <- all_taps[all_taps$X3!=".",] %>% dplyr::select(taps_aC,lrtaps_aC, X3) %>% melt(id.vars="X3") 
p <-ggplot(depth_sta, aes(value, fill=variable)) +
  geom_bar(stat="count", position=position_dodge() ) + theme(legend.position = "bottom") +
  xlab("depth") +
  facet_wrap(~X3, scales = "free", nrow=3) +
  xlim(0,40)
  
ggsave(paste0("plots/", id,"_all.depth_sta.pdf"),p,width =12, height = 8)


# covered repeats ----
sta <- data.frame("non_repeats"=c(6497927, 47295, 6154800, 12954564), "repeats"=c(3954621, 251304, 2596746, 8076613))
sta$group <- c("overlap","lrTAPS_only","TAPS_only","all")
sta$non_repeat_pct <- sta$non_repeats / (sta$non_repeats + sta$repeats)
sta$repeat_pct <- sta$repeats / (sta$non_repeats + sta$repeats)
sta_melt <- melt(sta %>% dplyr::select("group","non_repeat_pct","repeat_pct"), by=c("group"))
sta_melt$group <- factor(sta_melt$group, levels=c("all","TAPS_only","overlap","lrTAPS_only"))
p <- ggplot(sta_melt, aes(x=group, y=value, fill=variable)) +  geom_bar(stat="identity") + theme_light() + scale_fill_manual( values = c("#446db4","#23b177")) + ylab("percentage") + theme(legend.position = "bottom")
ggsave(paste0("plots/", id,"repeats_sta.pdf"),p,width =4, height = 4)




# Methylation on chrHMM ----
dat1 <- read.table("meth/bc1001_reseq.chrHMM.meth.sta", header=FALSE)
dat2 <- read.table("shortread_taps/mESC_cStates_HMM.taps.shortread.sta",header=FALSE)
dat <- merge(dat1,dat2,by=c("V1"))
colnames(dat) <- c("chrHMM","lrtaps_meth","lrtaps_mC","taps_aC","taps_meth","srtaps_mC","taps_aC")
p <- ggplot(dat[,c(1,grep("meth",colnames(dat)))] %>% melt(id.vars=c("chrHMM")),aes(x=chrHMM,y=value*100,fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") +
  xlab("chrHMM state") +
  ylab("methylation (%)") +
  scale_fill_manual(values=brewer.pal(3,"Dark2"))

ggsave("plots/bc1001_reseq_mESC_chrHMM.meth.pdf",p,width = 6, height = 4)




# Aligned Length ----





# Coverage on chrHMM ----
bc1001_cov <- read.table("align/bc1001_reseq.chrHMM.cov.sta", header=FALSE)
short_cov <- read.table("shortread_taps/merged_bam.md.cov.sta",header=FALSE)
len <- read.table("resource/mESC_cStates_HMM.len.sta",header=FALSE)

dat <- merge(bc1001_cov,short_cov,by=c("V1")) %>% merge(len,by=c("V1"))
colnames(dat) <- c("chrHMM","lrtaps_depth","lrtaps_pcov","lrtaps_cov","taps_depth","taps_pcov","taps_cov","len")
dat$lrtaps_p <- dat$lrtaps_cov/dat$len
dat$taps_p <- dat$taps_cov/dat$len

p1 <- ggplot(dat[,c(1,grep("depth",colnames(dat)))] %>% melt(id.vars=c("chrHMM")),aes(x=chrHMM,y=value*100,fill=variable))+
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7) +
  ylab("# total read") +
  xlab("chr HMM") +
  facet_grid(rows = vars(variable), scales="free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")

p2 <- ggplot(dat[,c(1,grep("p$",colnames(dat)))] %>% melt(id.vars=c("chrHMM")),aes(x=chrHMM,y=value*100,fill=variable))+
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7) +
  ylab("% cov") +
  xlab("chr HMM") +
  facet_grid(rows = vars(variable)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
p <- plot_grid(p1, p2)
p
ggsave("plots/bc1001_reseq.mESC_chrHMM.coverage.cmb.pdf",p,width = 12, height = 6)
ggsave("plots/bc1001_reseq.mESC_chrHMM.coverage.cmb.png",p,width = 12, height = 6)


# covered sites ---
grid.newpage()
pdf("plots/covered_sites.pdf",width = 4,height = 4)
bc1001 <- 1425646958
bc1002 <- 1484131961
shortr <- 2542284168  
bc1001_bc1002 <- 910008142 
bc1001_shortr <- 1423644222
bc1002_shortr <- 1481521517  
bc1001_bc1002_shortr <-  908875523
draw.triple.venn(
  area1     = bc1001,
  area2     = bc1002,
  area3     = shortr,
  n12       = bc1001_bc1002,
  n23       = bc1002_shortr,
  n13       = bc1001_shortr,
  n123      = bc1001_bc1002_shortr,
  category  = c('bc1001', 'bc1002', 'short_read'),
  fill      = c("#446db4","#23b177","#fac133"),
  cat.col   = c("#446db4","#23b177","#fac133"),
  euler.d = TRUE,
  scaled    = TRUE
)
dev.off()

dat <- read.table("align/bc1002.chrHMM.len.txt")
p <- ggplot(dat,aes(x=V8,y=V4)) +  geom_violin() + geom_boxplot(outlier.shape = NA, width=0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
ggsave("plots/frag_len.chrhmm.bc1002.pdf",p,width = 6, height = 4)
dat <- read.table("align/bc1001.chrHMM.len.txt")
p <- ggplot(dat,aes(x=V8,y=V4)) +  geom_violin() + geom_boxplot(outlier.shape = NA, width=0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
ggsave("plots/frag_len.chrhmm.bc1001.pdf",p,width = 6, height = 4)

# covered sites ----
bc1001_cov <- read.table("align/bc1001.lc_extrap.txt", header=TRUE)
bc1002_cov <- read.table("align/bc1002.lc_extrap.txt", header=TRUE)
bc1003_cov <- read.table("align/bc1003.lc_extrap.txt", header=TRUE)


dat <- merge(bc1001_cov,bc1002_cov,by="TOTAL_READS") %>% merge(bc1003_cov,by="TOTAL_READS")
colnames(dat) <- c("total","bc1001_distinct","bc1001_lower_ci","bc1001_upper_ci","bc1002_distinct","bc1002_lower_ci","bc1002_upper_ci","subbc1002_distinct","subbc1002_lower_ci","subbc1002_higher_ci")
dat$bc1001_p <- dat$bc1001_distinct/dat$total
dat$bc1002_p <- dat$bc1002_distinct/dat$total
dat$subbc1002_p <- dat$subbc1002_distinct/dat$total


p1 <- ggplot(dat[dat$total<10000000 & dat$total>0,c(grep("_p|total",colnames(dat)))] %>% melt(id.vars=c("total")),aes(x=total,y=value*100,fill=variable))+
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7) +
  ylab("distinct reads%") +
  xlab("reads sequenced") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
p1
ggsave("plots/lc_extrap.pdf",p1,width = 6, height = 4)


# Meth context ----
bc1001_meth <- read.table("meth/bc1001.uniflag_CpG.context.sta", header=FALSE)
bc1002_meth <- read.table("meth/bc1002.uniflag_CpG.context.sta", header=FALSE)

dat <- merge(bc1001_meth[,c(1,4)], bc1002_meth[,c(1,4)],by=c("V1"))
colnames(dat) <- c("context","bc1001","bc1002")

p1 <- ggplot(melt(dat),aes(x=context,y=value,fill=variable))+ 
  geom_bar(stat="identity", position = "dodge")+
  ylab("Methylation%") +
  xlab("smp") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
p1
ggsave("plots/meth_context.pdf",p1,width = 6, height = 4)

# Compare SV ----
pdf(paste0("plots/", id,"_SV.cmp.pdf"), width = 4.5, height = 4)
grid.newpage()
draw.pairwise.venn(
  area1     = 641+1170,
  area2     = 27299+1170,
  cross.area       = 1170,
  category  = c('ins_taps', 'ins_lrtaps'),
  fill      = c("#446db4","#23b177"),
  cat.col   = c("#446db4","#23b177"),
  euler.d = TRUE,
  scaled    = TRUE
)
grid.newpage()
draw.pairwise.venn(
  area1     = 6625+11418,
  area2     = 24723+11418,
  cross.area       = 11418,
  category  = c('del_taps', 'del_lrtaps'),
  fill      = c("#446db4","#23b177"),
  cat.col   = c("#446db4","#23b177"),
  euler.d = TRUE,
  scaled    = TRUE
)

dat <- read.table("plots/all_insertion_sta.txt", header=TRUE, stringsAsFactors = FALSE)
dat$svlen <- as.numeric(dat$svlen)
insertion_num <- spread(dat[,-c(4)], svlen, num)
colnames(insertion_num) <- c("type","0-50","50-100", "100+")
insertion_num[is.na(insertion_num)] <- 0
insertion_num <- melt(insertion_num,id.vars=("type"), variable.name = "svlen")
p1 <- ggplot(insertion_num, aes(x=svlen, y=value, fill=type)) + 
  geom_bar(stat="identity", position = "dodge")  + 
  ylab("num") + theme_light() + xlab("insertion len") + 
  scale_fill_manual(values=brewer.pal(3,"Dark2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

insertion_pct <- spread(dat[,-c(3)], svlen, pct)
insertion_pct[is.na(insertion_pct)] <- 0
colnames(insertion_pct) <- c("type","0-50","50-100", "100+")
insertion_pct <- melt(insertion_pct,id.vars=("type"), variable.name = "svlen")
p2 <- ggplot(insertion_pct, aes(x=svlen, y=value, fill=type)) + geom_bar(stat="identity", position = "dodge")  + ylab("pct") + theme_light()+ xlab("insertion len") + scale_fill_manual(values=brewer.pal(3,"Dark2"))

plot_grid(p1, p2, labels = c('A', 'B'), ncol=1)

dat <- read.table("plots/all_deletion_sta.txt", header=TRUE, stringsAsFactors = FALSE)
dat$svlen <- as.numeric(dat$svlen)
deletion_num <- spread(dat[,-c(4)], svlen, num)
colnames(deletion_num) <- c("type","0-50","50-100", "100+")
deletion_num[is.na(deletion_num)] <- 0
deletion_num <- melt(deletion_num,id.vars=("type"), variable.name = "svlen")
p3 <- ggplot(deletion_num, aes(x=svlen, y=value, fill=type)) + geom_bar(stat="identity", position = "dodge")  + ylab("num") + theme_light()+ xlab("deletion len") + 
  scale_fill_manual(values=brewer.pal(3,"Dark2")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

deletion_pct <- spread(dat[,-c(3)], svlen, pct)
deletion_pct[is.na(deletion_pct)] <- 0
colnames(deletion_pct) <-c("type","0-50","50-100", "100+")
deletion_pct <- melt(deletion_pct,id.vars=("type"), variable.name = "svlen")
p4 <- ggplot(deletion_pct, aes(x=svlen, y=value, fill=type)) + geom_bar(stat="identity", position = "dodge")  + ylab("pct") + theme_light()+ xlab("deletion len") + scale_fill_manual(values=brewer.pal(3,"Dark2"))

plot_grid(p3, p4, labels = c('C', 'D'),ncol=1)
dev.off()
