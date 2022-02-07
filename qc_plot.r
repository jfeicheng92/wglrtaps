library(Biostrings)
# QCs ----
setwd("~/cfo155/whole_genome_lrTAPS")
dat1 <- read.table("sta/bc1001.uniflag.sta", header=TRUE)
dat2 <- read.table("sta/bc1002.uniflag.sta", header=TRUE)
dat3 <- read.table("sta/bc1003.uniflag.sta", header=TRUE)
dat <- rbind(dat1,dat2)
rownames(dat) <- c("bc1001","bc1002")
dat$id <- c("bc1001","bc1002")
dat$conversion <- (gsub(".*;\\(|/.*","",dat$kb4_meth_c) %>% as.character() %>% as.numeric() +
                  gsub(".*;\\(|/.*","",dat$kb4_meth_g) %>% as.character() %>% as.numeric()) /
                  (gsub(".*/|\\)","",dat$kb4_meth_c) %>% as.character() %>% as.numeric() +
                  gsub(".*/|\\)","",dat$kb4_meth_g) %>% as.character() %>% as.numeric()) 

dat$fp <- (gsub(".*;\\(|/.*","",dat$kb4_unmeth_c) %>% as.character() %>% as.numeric() +
             gsub(".*;\\(|/.*","",dat$kb4_unmeth_g) %>% as.character() %>% as.numeric()) /
  (gsub(".*/|\\)","",dat$kb4_unmeth_c) %>% as.character() %>% as.numeric() +
     gsub(".*/|\\)","",dat$kb4_unmeth_g) %>% as.character() %>% as.numeric()) 


p <- ggplot(dat, aes(x = id, y = conversion*100)) + geom_bar(stat = "identity", fill = "orangered3", col = "orangered3") + 
  geom_text(aes(label=round(conversion*100, 2)), vjust=-0.4) +
  theme_light() + xlab("") + ylab("mCG conversion on 4kb (%)") + theme(axis.text = element_text(color = "black")) + ylim(0,100)
ggsave("plots/mC_conversion.pdf",p,width=2,height = 4)
p <- ggplot(dat, aes(x = id, y = fp*100)) + geom_bar(stat = "identity", fill = "steelblue", col = "steelblue") + 
  geom_text(aes(label=round(fp*100, 2)), vjust=-0.4) +
  theme_light() + xlab("") + ylab("false positive on unmodified cg(4kb) (%)") + theme(axis.text = element_text(color = "black"))  + ylim(0,0.5)
ggsave("plots/mC_fp.pdf",p,width=2,height = 4)



dat$genome_meth <- (gsub(".*;\\(|/.*","",dat$genome_meth_c) %>% as.character() %>% as.numeric() +
                     gsub(".*;\\(|/.*","",dat$genome_meth_g) %>% as.character() %>% as.numeric()) /
  (gsub(".*/|\\)","",dat$genome_meth_c) %>% as.character() %>% as.numeric() +
     gsub(".*/|\\)","",dat$genome_meth_g) %>% as.character() %>% as.numeric()) 
p <- ggplot(dat, aes(x = id, y = genome_meth*100, fill=id)) + geom_bar(stat = "identity") + 
  geom_text(aes(label=round(genome_meth*100, 2)), vjust=-0.4) +
  theme_light() + xlab("") + ylab("genome methylation (%)") + theme(axis.text = element_text(color = "black"), legend.position = "none") +
  scale_fill_manual(values=brewer.pal(2,"Dark2"))
ggsave("plots/genome_meth.pdf",p,width=2,height = 4)


dat$pmap <- (dat$nfwd_bam + dat$nrev_bam)/dat$nraw
p <- ggplot(dat, aes(x = id, y = pmap*100)) + geom_bar(stat = "identity", fill = "steelblue", col = "steelblue") + 
  geom_text(aes(label=round(pmap*100, 2)), vjust=-0.4) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  theme_light() + xlab("") + ylab(" Q20 mapping rate (%)") + theme(axis.text = element_text(color = "black"))
  
ggsave("plots/map_rate.pdf",p,width=2,height = 4)



# Methylation on chrHMM ----
dat1 <- read.table("meth/bc1001.uniflag_CpG.chrHMM.sta", header=FALSE)
dat2 <- read.table("meth/bc1002.uniflag_CpG.chrHMM.sta", header=FALSE)
dat3 <- read.table("resource/mESC_cStates_HMM.taps.shortread.sta",header=FALSE)
dat <- merge(dat1,dat2,by=c("V1")) %>% merge(dat3,by=c("V1"))
colnames(dat) <- c("chrHMM","bc1001_meth","bc1001_mC","bc_1001_aC","bc1002_meth","bc1002_mC","bc_1002_aC","taps_meth","taps_mC","taps_aC")
p <- ggplot(dat[,c(1,grep("meth",colnames(dat)))] %>% melt(id.vars=c("chrHMM")),aes(x=chrHMM,y=value*100,fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") +
  xlab("chrHMM state") +
  ylab("methylation (%)") +
  scale_fill_manual(values=brewer.pal(3,"Dark2"))

ggsave("plots/mESC_chrHMM.meth.pdf",p,width = 6, height = 4)


# Depth vs. CpG density ----
bc1001_tg_cov <- read.table("align/bc1001.tg.cov.bedGraph", header=FALSE)[-c(1:2),1:5]
bc1001_ca_cov <- read.table("align/bc1001.ca.cov.bedGraph", header=FALSE)[-c(1:2),1:5]
bc1002_tg_cov <- read.table("align/bc1002.tg.cov.bedGraph", header=FALSE)[-c(1:2),1:5]
bc1002_ca_cov <- read.table("align/bc1002.ca.cov.bedGraph", header=FALSE)[-c(1:2),1:5]

cov <- merge(bc1001_tg_cov, bc1001_ca_cov,by=c("V1","V2","V3")) %>% 
  merge(bc1002_tg_cov,by=c("V1","V2","V3")) %>% 
  merge(bc1002_ca_cov,by=c("V1","V2","V3"))
colnames(cov) <- c("chr","start","end","bc1001_tg_depth","bc1001_tg_cov","bc1001_ca_depth","bc1001_ca_cov","bc1002_tg_depth","bc1002_tg_cov","bc1002_ca_depth","bc1002_ca_cov")
cov$bc1001_tg_depth <- cov$bc1001_tg_depth/7.45
cov$bc1001_ca_depth <- cov$bc1001_ca_depth/7.45
cov$bc1002_tg_depth <- cov$bc1002_tg_depth/8.08
cov$bc1002_ca_depth <- cov$bc1002_ca_depth/8.08

taps <- read_delim("/gpfs2/well/ludwig/users/cfo155/ncbi/public/sra/taps_pub.average.methratio.txt",delim="\t",col_names=TRUE)
taps$taps <- (taps$taps_meth>0.1 & (taps$mC + taps$uC >3)) %>% as.numeric()
taps$idx <- floor(taps$start/50000)*50000
taps_dens <- aggregate(taps ~ chr+idx, taps, sum)
cov_dens <- merge(cov,taps_dens, by.x=c("chr","start"), by.y=c("chr","idx"),all.x=TRUE)
cov_dens$taps_lvl <- cut(cov_dens$taps, breaks=quantile(cov_dens$taps,na.rm=TRUE, probs=seq(0,1,1/10)), include.lowest = TRUE, labels = seq(1,10))


p <- ggplot(na.omit(cov_dens[,c(1:3,grep("depth|lvl",colnames(cov_dens)))]) %>% melt(id.vars=c("chr","start","end","taps_lvl")),aes(x=taps_lvl,y=value)) +
  geom_boxplot(outlier.shape = NA) + 
  ylim(0,4)+
  facet_wrap(factor(variable) ~ ., scales="free", nrow=2) +
  xlab("mCG density") +
  ylab("depth in 50k bin") 

ggsave("plots/cpg_dens.pdf",p,width = 6, height = 5)


# Read Length ----
setwd("~/cfo155/whole_genome_lrTAPS")
dat1 <- read.table("fastq/bc1001.read_len.txt", header=TRUE)
dat2 <- read.table("fastq/bc1002.read_len.txt", header=TRUE)
dat1[,2] <- dat1[,2]/7.45
dat2[,2] <- dat2[,2]/8.08

dat <- merge(dat1,dat2,by=c("len"), all.x = TRUE, all.y=TRUE) %>% melt(id.vars="len")
dat$variable <- gsub("fastq.|.read_len.txt","",dat$variable )
p <- ggplot(dat, aes(x=len,y=value,fill=variable)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7) +
  ylab("normalized reads number") +
  xlab("read length") +
  theme_light() + xlim(0,10000) +
  theme(legend.position = "bottom")

ggsave("plots/read_lens.pdf",p,width = 6, height = 4)
ggsave("plots/read_lens.png",p,width = 6, height = 4)
# Aligned Length ----
dat1 <- read.table("align/bc1001.ca.sta.txt", header=TRUE)
dat2 <- read.table("align/bc1001.tg.sta.txt", header=TRUE)
dat3 <- read.table("align/bc1002.ca.sta.txt", header=TRUE)
dat4 <- read.table("align/bc1002.tg.sta.txt", header=TRUE)
dat1[,2] <- dat1[,2]/7.45; dat2[,2] <- dat2[,2]/7.45
dat3[,2] <- dat3[,2]/8.08; dat4[,2] <- dat4[,2]/7.45


dat <- merge(dat1,dat2,by=c("len"), all.x = TRUE, all.y=TRUE) %>% 
  merge(dat3,by=c("len"), all.x = TRUE, all.y=TRUE )%>% 
  merge(dat4, by=c("len"), all.x = TRUE, all.y=TRUE)  %>% melt(id.vars="len")

dat$variable <- gsub("align.|.len.txt","",dat$variable )
p <- ggplot(dat, aes(x=len,y=value,fill=variable)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7) +
  ylab("normalized reads number") +
  xlab("align length") +
  theme_light() + xlim(0,8000) +
  theme(legend.position = "bottom")

ggsave("plots/align_lens.pdf",p,width = 6, height = 4)
ggsave("plots/align_lens.png",p,width = 6, height = 4)


dat1 <- read.table("align/bc1001.ca.len.txt", header=FALSE)
dat2 <- read.table("align/bc1001.tg.len.txt", header=FALSE)
N50(c(dat1$V1,dat2$V1))
dat3 <- read.table("align/bc1002.ca.len.txt", header=FALSE)
dat4 <- read.table("align/bc1002.tg.len.txt", header=FALSE)
N50(c(dat3$V1,dat4$V1))

# Coverage on chrHMM ----
bc1001_tg_cov <- read.table("align/bc1001.tg.chrHMM.cov.sta", header=FALSE)
bc1001_ca_cov <- read.table("align/bc1001.ca.chrHMM.cov.sta", header=FALSE)
bc1002_tg_cov <- read.table("align/bc1002.tg.chrHMM.cov.sta", header=FALSE)
bc1002_ca_cov <- read.table("align/bc1002.ca.chrHMM.cov.sta", header=FALSE)
bc1003_tg_cov <- read.table("align/bc1003.tg.chrHMM.cov.sta", header=FALSE)
bc1003_ca_cov <- read.table("align/bc1003.ca.chrHMM.cov.sta", header=FALSE)
len <- read.table("resource/mESC_cStates_HMM.len.sta",header=FALSE)


dat <- merge(bc1001_tg_cov,bc1001_ca_cov,by=c("V1")) %>% merge(bc1002_tg_cov,by=c("V1")) %>% merge(bc1002_ca_cov,by=c("V1")) %>% merge(bc1003_tg_cov,by=c("V1")) %>% merge(bc1003_ca_cov,by=c("V1")) %>% merge(len,by=c("V1"))
colnames(dat) <- c("chrHMM","bc1001_tg_depth","bc1001_tg_pcov","bc1001_tg_cov", "bc1001_ca_depth","bc1001_ca_pcov","bc1001_ca_cov","bc1002_tg_depth","bc1002_tg_pcov","bc1002_tg_cov", "bc1002_ca_depth","bc1002_ca_pcov","bc1002_ca_cov", "subbc1002_tg_depth","subbc1002_tg_pcov","subbc1002_tg_cov", "subbc1002_ca_depth","subbc1002_ca_pcov","subbc1002_ca_cov","len")
dat$bc1001_tg_p <- dat$bc1001_tg_cov/dat$len
dat$bc1001_ca_p <- dat$bc1001_ca_cov/dat$len
dat$bc1002_tg_p <- dat$bc1002_tg_cov/dat$len
dat$bc1002_ca_p <- dat$bc1002_ca_cov/dat$len
dat$subbc1002_tg_p <- dat$subbc1002_tg_cov/dat$len
dat$subbc1002_ca_p <- dat$subbc1002_ca_cov/dat$len


p1 <- ggplot(dat[,c(1,grep("depth",colnames(dat)))] %>% melt(id.vars=c("chrHMM")),aes(x=chrHMM,y=value*100,fill=variable))+
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7) +
  ylab("# total read") +
  xlab("chr HMM") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")

p2 <- ggplot(dat[,c(1,grep("p$",colnames(dat)))] %>% melt(id.vars=c("chrHMM")),aes(x=chrHMM,y=value*100,fill=variable))+
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7) +
  ylab("% cov") +
  xlab("chr HMM") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
p <- plot_grid(p1, p2)
ggsave("plots/mESC_chrHMM.coverage.pdf",p,width = 15, height = 6)
ggsave("plots/mESC_chrHMM.coverage.png",p,width = 15, height = 6)

## combine
bc1001_cov <- read.table("align/bc1001.chrHMM.cov.sta", header=FALSE)
bc1002_cov <- read.table("align/bc1002.chrHMM.cov.sta", header=FALSE)
bc1003_cov <- read.table("align/bc1003.chrHMM.cov.sta", header=FALSE)
short_cov <- read.table("shortread_taps/merged_bam.md.cov.sta",header=FALSE)
len <- read.table("resource/mESC_cStates_HMM.len.sta",header=FALSE)


dat <- merge(bc1001_cov,bc1002_cov,by=c("V1")) %>% merge(bc1003_cov,by=c("V1"))  %>% merge(short_cov,by=c("V1")) %>% merge(len,by=c("V1"))
colnames(dat) <- c("chrHMM","bc1001_depth","bc1001_pcov","bc1001_cov","bc1002_depth","bc1002_pcov","bc1002_cov","subbc1002_depth","subbc1002_pcov","subbc1002_cov","short_depth","short_pcov","short_cov","len")
dat$bc1001_p <- dat$bc1001_cov/dat$len
dat$bc1002_p <- dat$bc1002_cov/dat$len
dat$subbc1002_p <- dat$subbc1002_cov/dat$len
dat$short_p <- dat$short_cov/dat$len

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
ggsave("plots/mESC_chrHMM.coverage.cmb.pdf",p,width = 12, height = 6)
ggsave("plots/mESC_chrHMM.coverage.cmb.png",p,width = 12, height = 6)
# CGI cov ----
bc1001_cov <- read.table("align/bc1001.cgi.cov.sta", header=FALSE)
bc1002_cov <- read.table("align/bc1002.cgi.cov.sta", header=FALSE)
bc1003_cov <- read.table("align/bc1003.cgi.cov.sta", header=FALSE)
len <- read.table("resource/cpgIslandExt.txt",header=FALSE)


dat <- cbind(bc1001_cov,bc1002_cov) %>% cbind(bc1003_cov)  %>% cbind(sum(len$V4-len$V3))
colnames(dat) <- c("bc1001_depth","bc1001_pcov","bc1001_cov","bc1002_depth","bc1002_pcov","bc1002_cov","subbc1002_depth","subbc1002_pcov","subbc1002_cov","len")
dat$bc1001_p <- dat$bc1001_cov/dat$len
dat$bc1002_p <- dat$bc1002_cov/dat$len
dat$subbc1002_p <- dat$subbc1002_cov/dat$len
dat$feature <- "cgi"

p1 <- ggplot(dat[,c(grep("depth|feature",colnames(dat)))] %>% melt(id.vars=c("feature")),aes(x=feature,y=value*100,fill=variable))+
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7) +
  ylab("# total read") +
  xlab("feature") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")

p2 <- ggplot(dat[,c(grep("p$|feature",colnames(dat)))] %>% melt(id.vars=c("feature")),aes(x=feature,y=value*100,fill=variable))+
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7) +
  ylab("% cov") +
  xlab("feature") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
p <- plot_grid(p1, p2)
ggsave("plots/mESC_cgi.coverage.cmb.pdf",p,width = 4, height = 6)
ggsave("plots/mESC_cgi.coverage.cmb.png",p,width = 4, height = 6)

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

# PCR duplicates ----
bc1001_dup <- read.table("align/bc1001.rm.sta", header=FALSE)
bc1002_dup <- read.table("align/bc1002.rm.sta", header=FALSE)


dat <- rbind(bc1001_dup, bc1002_dup)
dat$smp <- c("bc1001","bc1002")

p1 <- ggplot(dat,aes(x=smp,y=(1-V3)*100))+ 
  geom_bar(stat="identity")+
  geom_text(aes(label=round((1-V3)*100, 2)), vjust=-0.4) +
  ylab("PCR duplicates%") +
  xlab("smp") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
p1
ggsave("plots/pcr_rm.pdf",p1,width = 2, height = 4)

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



