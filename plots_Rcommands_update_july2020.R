#Code to reproduce all of the plots

#Tab_new: this is the expression tab with the RNAseq data mapped to TEs
#Annotation_tab_new: this is the information from the repeatmasker annotations of the different TEs
load("Annotation_Tab_new.Rdata")
library(clinfun)
#required for the jonckheere tests for statistics

#figure 1 and supplemental figure 1
load("TE_genebodies_FPKM.Rdata")

 sums_all_new_tab<-apply(Tab_new,2,sum)
 
 GROUPS_texpsn<-c(rep(0,length=20),rep(1,length=10),rep(2,length=5))
 GROUPS_texpsn_bp<-GROUPS_texpsn+1
 COLS_texpsn<-GROUPS_texpsn
 COLS_texpsn[COLS_texpsn==0]<-"red"
 COLS_texpsn[COLS_texpsn=="1"]<-"green"
 COLS_texpsn[COLS_texpsn=="2"]<-"blue"
 
 #figure 1B
  boxplot(sums_all_new_tab[1:20],sums_all_new_tab[21:30],sums_all_new_tab[31:35], names=c(0,1,2), boxwex=0.5,ylab="TE fpkm", xlab="log10 population size")
points(GROUPS_texpsn_bp,sums_all_new_tab[1:35], col=COLS_texpsn, pch=18,cex=1.5)
 abline(h=sums_all_new_tab[36],col="gray")
 fit_lm_allte<-lm(sums_all_new_tab[1:35]~GROUPS_texpsn_bp)
 abline(fit_lm_allte)
 cor_all_lm_te<-cor.test(sums_all_new_tab[1:35],GROUPS_texpsn_bp)
 
 text(1.5,fit_lm_allte$coefficients[2]*1.5+fit_lm_allte$coefficients[1]+0.2*(fit_lm_allte$coefficients[2]*1.5+fit_lm_allte$coefficients[1]), paste("p=", round(cor_all_lm_te$p.value,digits=3)))
 text(1.5,fit_lm_allte$coefficients[2]*1.5+fit_lm_allte$coefficients[1]+0.3*(fit_lm_allte$coefficients[2]*1.5+fit_lm_allte$coefficients[1]), paste("r=", round(cor_all_lm_te$estimate,digits=3)))
 dev.copy(pdf, "figure1B_sept.pdf")
 dev.off()
 
 
 #supplemental figure 1A
 
 Tab_new_norm<-(Tab_new[,1:35]+1)/(Tab_new[,36]+1)
 Tab_new_norm_means<-apply(Tab_new_norm,2,mean)
 boxplot(log2(Tab_new_norm_means[1:20]),log2(Tab_new_norm_means[21:30]),log2(Tab_new_norm_means[31:35]),ylab="log2(TE/pma", names=c(0,1,2), xlab="log10 population size")
 points(GROUPS_texpsn_bp, log2(Tab_new_norm_means),col=COLS_texpsn,pch=18,cex=1.5, ylab="log2(TE/pma)",xlab="log10 population size")
 
 fit_lm_normte<-lm(log2(Tab_new_norm_means)~GROUPS_texpsn_bp)
 abline(fit_lm_normte)
 jh_test_normte<-jonckheere.test(Tab_new_norm_means, GROUPS_texpsn_bp)
 cor_norm_te<-cor.test(log2(Tab_new_norm_means),GROUPS_texpsn_bp)
 legend("topright", paste("jh p=", round(jh_test_normte$p.value, digits=3), sep=""),bty="n")
legend("topright", paste("\n\nLM p=", round(cor_norm_te$p.value, digits=3),sep=""),bty="n")
 dev.copy(pdf, "Supplemental_figure1A.pdf")
 dev.off()
 
 #supplemental figure 1B
 norm_list<-list()
 for(i in 1:35){norm_list[[i]]<-log2(Tab_new_norm[,i])
 }
 
 meanE<-sapply(norm_list,mean)
 ordered_list_norm<-list()
 cols_order<-c()
 for(i in 1:35){ordered_list_norm[[i]]<-norm_list[[order(meanE)[i]]]
 cols_order<-c(cols_order, COLS_texpsn[order(meanE)[i]])
 
 }
 
 library(vioplot)
 
 vioplot(ordered_list_norm, col=cols_order,pchMed="", ylim=c(-2,5), ylab="log2 TE/pma")
 abline(h=0, col="gray")
 
 legend("topleft", c("N=1","N=10","N=100"),fill=c("red","blue","green"))
 dev.copy(pdf, "Supplemental_figure_1B.pdf")
 
 
 #figure 1C
 wilcoxon_table_xpsn<-matrix(0, ncol=3, nrow=nrow(Tab_new))
 for(i in 1:nrow(Tab_new)){
 M<-wilcox.test(Tab_new[i,1:20],Tab_new[i,21:36])$p.value
 wilcoxon_table_xpsn[i,1]<-M
 }
 wilcoxon_table_xpsn[,2]<-apply(Tab_new[,1:20],1,median)
 wilcoxon_table_xpsn[,3]<-apply(Tab_new[,21:36],1,median)
 colnames(wilcoxon_table_xpsn)<-c("pval","median_Ne1","median_Ne>1")
 row.names(wilcoxon_table_xpsn)<-row.names(Tab_new)
 #plot(wilcoxon_table_xpsn[,2]-wilcoxon_table_xpsn[,3],(-1)*log10(wilcoxon_table_xpsn[,1]),pch=18,cex=1.5,ylab="log10 p val",xlab="median Ne1-median Ne>1")
 
 
  load("Annotations_Tab_new.Rdata")
 TE_cats<-list("TC1","TC4","TURMOIL2",c("LINE","VINGI"),"CER")
 names(TE_cats)<-c("TC1","TC4","TURMOIL2","non-LTR","LTR")
 TE_cols<-c("purple","hotpink","red","green","cyan")
 for(i in 1:length(TE_cats)){A<-c();for(j in 1:length(TE_cats[[i]])){A<-c(A,grep(TE_cats[[i]][j], Annotations_Tab_new[,2]))}
 
 points(wilcoxon_table_xpsn[A,2]-wilcoxon_table_xpsn[A,3],(-1)*log10(wilcoxon_table_xpsn[A,1]),pch=18,cex=1.5,col=TE_cols[i])}
 
 
 legend("topright", names(TE_cats), col=TE_cols,cex=1.5,pch=18)
 
abline(h=(-1)*log10(0.05),lty=2, col="orange")
text(1,1.5,"p=0.05",col="orange")
dev.copy(pdf, "figure1C.pdf")
dev.off()



#supplemental figure 2


jh_tests_all<-matrix(0, ncol=3,nrow=nrow(Tab_new))

for(i in 1:nrow(Tab_new)){Test_jh<-jonckheere.test(Tab_new[i,1:35],g=GROUPS_texpsn,nperm=1000);
Test_fit<-lm(Tab_new[i,1:35]~GROUPS_texpsn)
jh_tests_all[i,1]<-Test_jh$p.value
jh_tests_all[i,2]<-summary(Test_fit)$coefficients[2,4]
jh_tests_all[i,3]<-Test_fit$coefficients[2]
}
colnames(jh_tests_all)<-c("jonckeere pval","linear model pval","gradient expression vs popsize")

plot((-1)*jh_tests_all[,3],(-1)*log10(jh_tests_all[,2]), pch=18,cex=1.5,ylab="log10 pval",xlab="gradient with decreasing popsize")
for(i in 1:length(TE_cats)){A<-c();for(j in 1:length(TE_cats[[i]])){A<-c(A,grep(TE_cats[[i]][j], Annotations_Tab_new[,2]))}
 
 points((-1)*jh_tests_all[A,3],(-1)*log10(jh_tests_all[A,2]),pch=18,cex=1.5,col=TE_cols[i])}
legend("topright",names(TE_cats), col=TE_cols, cex=1.5,pch=18)
abline(h=(-1)*log10(0.05),col="orange",lty=2)
dev.copy(pdf, "Supplemental_figure2A.pdf")
dev.off()

q<-which(jh_tests_all[,1]<0.05)
subset<-Tab_new[q,]
row.names(subset)<-paste(Annotations_Tab_new[q,2],Annotations_Tab_new[q,1], sep="|")
pdf("Supplemental_figure2B.pdf")
heatmap.2(log2(subset+1)-log2(subset[,36]+1), col=colorRampPalette(c("blue","white","orange"))(51), trace="none", Colv=F,ColSideColors=c(COLS_texpsn,"gray"),cexRow=0.5)
dev.off()


#figure 1D
Turmoil2_genes<-Annotations_Tab_new[grep("TURMOIL2", Annotations_Tab_new[,2]),1]
Turmoil2_total_expression<-apply(Tab_new[which(row.names(Tab_new)%in%Turmoil2_genes==T),],2,sum)
boxplot(Turmoil2_total_expression[1:20],Turmoil2_total_expression[21:30],Turmoil2_total_expression[31:35], boxwex=0.5, ylab="Turmoil2_FPKM", xlab="log10 population size")
points(GROUPS_texpsn_bp,Turmoil2_total_expression[1:35],pch=18,col=COLS_texpsn,cex=1.5,ylab="Turmoil2 FPKM",xlab="log10 population size")
abline(h=Turmoil2_total_expression[36], col="gray")
LM_turmoil2<-lm(Turmoil2_total_expression[1:35]~GROUPS_texpsn_bp)
abline(LM_turmoil2)
text(1.5,15, paste("p=",round(summary(LM_turmoil2)$coefficients[2,4],digits=3),sep=""))
jh_test_turmoil2<-jonckheere.test(Turmoil2_total_expression[1:35],GROUPS_texpsn_bp,nperm=10000)
text(1.5,30,paste("jh p=",round(jh_test_turmoil2$p.value,digits=3),sep=""),col="blue")
dev.copy(pdf, "Figure1D.pdf")
dev.off()


#figure 1E
Tc1_genes<-Annotations_Tab_new[grep("TC1", Annotations_Tab_new[,2]),1]

Tc1_total_expression<-apply(Tab_new[which(row.names(Tab_new)%in%Tc1_genes==T),],2,sum)
boxplot(Tc1_total_expression[1:20],Tc1_total_expression[21:30],Tc1_total_expression[31:35], boxwex=0.5, xlab="log10 population size", ylab="Tc1 FPKM")
points(GROUPS_texpsn_bp,Tc1_total_expression[1:35],pch=18,col=COLS_texpsn,cex=1.5,ylab="Tc1 FPKM",xlab="log10 population size")
abline(h=Tc1_total_expression[36],col="gray")
jh_test_tc1<-jonckheere.test(Tc1_total_expression[1:35],GROUPS_texpsn_bp)
dev.copy(pdf, "Figure1E.pdf")
dev.off()


#Figure 1F,G,H
load("Counts_100g_fpkm_norm.Rdata")
load("Counts_25g_fpkm_norm.Rdata")

counts_all_25_norm_pma<-(counts_all_25_fpkm+1)/(counts_all_25_fpkm[,12]+1)
counts_all_100_norm_pma<-(counts_all_100_fpkm+1)/(counts_all_25_fpkm[,12]+1)

TE_counts_25<-counts_all_25_norm_pma[which(row.names(counts_all_25_norm_pma)%in%row.names(Tab_new)==T),]
TE_counts_100<-counts_all_100_norm_pma[which(row.names(counts_all_100_norm_pma)%in%row.names(Tab_new)==T),]


TE_means_25g<-apply(TE_counts_25[,1:11],2,mean)
TE_means_100g<-apply(TE_counts_100[,1:11],2,mean)

TE_means_400g<-Tab_new_norm_means[1:20]



VALS_timecourse<-c(TE_means_25g,TE_means_100g,TE_means_400g)
TIMES_timecourse<-c(rep(25,length=11),rep(100,length=11),rep(400, length=20))
cols_timecourse<-c(rep("darkorange2",length=11),rep("tomato",length=11),rep("red",length=20))

plot(TIMES_timecourse,(VALS_timecourse), col=cols_timecourse, ylab="mean TE expression/pma", xlab="generations", log="x",pch=18,cex=1.5,xaxt="n")
axis(side=1, at=c(25,100,400), labels=c(25,100,400))

lm_generations<-lm((VALS_timecourse)~log10(TIMES_timecourse))
abline(lm_generations)
text(250,1.2, paste("LM p=", format(summary(lm_generations)$coefficients[2,4],digits=3,scientific=T)))
library(clinfun)
jh_generations<-jonckheere.test(VALS_timecourse, TIMES_timecourse)
text(250, 1.5, paste("JH p=", format(jh_generations$p.value, digits=3, scientific=T)))
dev.copy(pdf, "Figure1F.pdf")
dev.off()

TE_sums_25g_turmoil<-apply(counts_all_25_fpkm[which(row.names(counts_all_25_fpkm)%in%Turmoil2_genes==T),],2,sum)
TE_sums_100g_turmoil<-apply(counts_all_100_fpkm[which(row.names(counts_all_100_fpkm)%in%Turmoil2_genes==T),],2,sum)

VALS_turmoil_generations<-c(TE_sums_25g_turmoil[1:11]/TE_sums_25g_turmoil[12],TE_sums_100g_turmoil/TE_sums_25g_turmoil[12],Turmoil2_total_expression[1:20]/Turmoil2_total_expression[36])

plot(TIMES_timecourse,(VALS_turmoil_generations),col=cols_timecourse,ylab="total turmoil2 expression/pma", xlab="generations", log="x", pch=18, cex=1.5,xaxt="n")
axis(side=1, at=c(25,100,400), labels=c(25,100,400))

lm_turmoil_generations<-lm((VALS_turmoil_generations)~log10(TIMES_timecourse))
abline(lm_turmoil_generations)
text(250,3,paste("p=", format(summary(lm_turmoil_generations)$coefficients[2,4], digits=3, scientific=T)))

jh_turmoil_generations<-jonckheere.test(VALS_turmoil_generations, TIMES_timecourse)
text(250,6,paste("jh p=", format(jh_turmoil_generations$p.value, digits=3, scientific=T)))
dev.copy(pdf, "Figure1G.pdf")
dev.off()

TE_sums_25g_tc1<-apply(counts_all_25_fpkm[which(row.names(counts_all_25_fpkm)%in%Tc1_genes==T),],2,sum)
TE_sums_100g_tc1<-apply(counts_all_100_fpkm[which(row.names(counts_all_100_fpkm)%in%Tc1_genes==T),],2,sum)

VALS_tc1_generations<-c(TE_sums_25g_tc1[1:11]/TE_sums_25g_tc1[12], TE_sums_100g_tc1/TE_sums_25g_tc1[12], Tc1_total_expression[1:20]/Tc1_total_expression[36])

plot(TIMES_timecourse,(VALS_tc1_generations), col=cols_timecourse, ylab="total tc1 expression/pma", xlab="generations", log="x", pch=18, cex=1.5,xaxt="n")
axis(side=1, at=c(25,100,400), labels=c(25,100,400))

lm_tc1_generations<-lm((VALS_tc1_generations)~log10(TIMES_timecourse))

abline(lm_tc1_generations)
text(250,2,paste("p=", format(summary(lm_tc1_generations)$coefficients[2,4], digits=3, scientific=F)))
jh_tc1_generations<-jonckheere.test(VALS_tc1_generations, TIMES_timecourse)
text(250, 300, paste("p=", format(jh_tc1_generations$p.value, digits=3, scientific=T)))
dev.copy(pdf, "Figure1H.pdf")
dev.off()


#Supplemental Figure 1 

#gene expression fpkm 

load("fpkm_genes.Rdata")
#nb tes have been removed from this set of genes











Vm_tes<-cbind(apply(Tab_new[,1:20],1,var),apply(Tab_new[,21:30],1,var), apply(Tab_new[,31:35],1,var))

#concern that this is unduly influenced by the number of samples
#thus we will use a downsampling approach in order to estimate the sd for each
Samples_vm_size1<-combn(1:20,5)
Samples_vm_size10<-combn(21:30,5)
Samples_vm_size1_1000<-sample(1:ncol(Samples_vm_size1),1000,replace=F)

All_size1_ff<-c()
All_size10_ff<-c()
for(i in 1:1000){
temp<-Tab_new[,Samples_vm_size1[,Samples_vm_size1_1000[i]]]	
	All_size1_ff<-cbind(All_size1_ff, apply(temp,1,var)/apply(temp,1,mean))
	
				if(i<253){
				All_size10_ff<-cbind(All_size10_ff, apply(Tab_new[,Samples_vm_size10[,i]],1,var)/apply(Tab_new[,Samples_vm_size10[,i]],1,mean))
				}
				}
				
est_ff<-cbind(apply(All_size1_ff,1,mean),apply(All_size10_ff,1,mean),apply(Tab_new[,31:35],1,var)/apply(Tab_new[,31:35],1,mean))
boxplot(est_ff, outline=F, names=c("N1","N10","N100"),col=c("red","green","blue"),ylab="fano factors")





est_sd<-cbind(apply(All_size1_vm,1,mean),apply(All_size10_vm,1,mean),apply(Tab_new[,31:35],1,var))
boxplot(log10(est_sd+1),outline=F, names=c("N1","N10","N100"),col=c("red","green","blue"),ylab="log10 sd TEfpkm")
print(wilcox.test(apply(All_size1_vm,1,mean),apply(All_size10_vm,1,mean),paired=T))
data:  apply(All_size1_vm, 1, mean) and apply(All_size10_vm, 1, mean)
V = 2605, p-value = 0.003473
alternative hypothesis: true location shift is not equal to 0
print(wilcox.test(apply(All_size1_vm,1,mean),apply(Tab_new[,31:35],1,sd),paired=T))
data:  apply(All_size1_vm, 1, mean) and apply(Tab_new[, 31:35], 1, sd)
V = 2576, p-value = 0.001047
alternative hypothesis: true location shift is not equal to 0
print(wilcox.test(apply(All_size10_vm,1,mean),apply(Tab_new[,31:35],1,sd),paired=T))
data:  apply(All_size10_vm, 1, mean) and apply(Tab_new[, 31:35], 1, sd)
V = 1912, p-value = 0.7128
alternative hypothesis: true location shift is not equal to 0

#fano factors are a better measure of noise for low expression.
fano_factors<-cbind(apply(Tab_new[,1:20],1,var)/apply(Tab_new[,1:20],1,mean),apply(Tab_new[,21:30],1,var)/apply(Tab_new[,21:30],1,mean),apply(Tab_new[,31:35],1,var)/apply(Tab_new[,31:35],1,mean))

boxplot(fano_factors, outline=F, names=c("N1","N10","N100"),col=c("red","green","blue"),ylab="Fano factors")
dev.copy(pdf, "FigureS2A.pdf")
dev.off()

wilcox.test(fano_factors[,1], fano_factors[,2], paired=T)
#p=0.01544|0.03088 adj

wilcox.test(fano_factors[,1],fano_factors[,3], paired=T)
#p=0.004253|0.01359 adj
wilcox.test(fano_factors[,2], fano_factors[,3], paired=T)

#p=0.7423

Vm_genes<-cbind(apply(fpkm_val[,1:20],1,var), apply(fpkm_val[,21:30],1,var), apply(fpkm_val[,31:35],1,var))

fano_factors_genes<-cbind(Vm_genes[,1]/apply(fpkm_val[,1:20],1,mean),Vm_genes[,2]/apply(fpkm_val[,21:30],1,mean),Vm_genes[,3]/apply(fpkm_val[,31:35],1,mean))

boxplot(fano_factors_genes, outline=F, names=c("N1","N10","N100"), col=c("red","green","blue"),ylab="Fano factors_genes")
dev.copy(pdf, "FigureS2B.pdf")
dev.off()
wilcox.test(fano_factors_genes[,1], fano_factors_genes[,2], paired=T)$p.value
#p<1e-16




Var_per_strain_tes<-apply(log2(Tab_new[,1:35]+1)-log2(Tab_new[,36]+1),2,var)
Var_per_strain_genes<-apply(log2(fpkm_val[,1:35]+1)-log2(fpkm_val[,36]+1),2,var)


Var_per_strain_tes2<-apply(log2(Tab_new[,1:35]+1)-log2(Tab_new[,36]+1),2,var)
Var_per_strain_genes2<-apply(log2(fpkm_val[,1:35]+1)-log2(fpkm_val[,36]+1),2,var)

plot(GROUPS_texpsn, Var_per_strain_tes[1:35], pch=18, col=COLS_texpsn, cex=1.5, ylab="var texpsn",xlab="log10 pop size")
lm_var_per_strain<-lm(Var_per_strain_tes[1:35]~GROUPS_texpsn)
abline(lm_var_per_strain)
text(1.5,0.2, paste("p=", round(summary(lm_var_per_strain)$coefficients[2,4], digits=3), sep=""))
jh_test_var_te<-jonckheere.test(Var_per_strain_tes[1:35],GROUPS_texpsn)
text(1.5, 0.35, paste("jh p=", round(jh_test_var_te$p.value,digits=3), sep=""))


dev.copy(pdf, "FigureS2C.pdf")
dev.off()



plot(GROUPS_texpsn, Var_per_strain_genes[1:35], pch=18, col=COLS_texpsn, cex=1.5, ylab="var_genexpsn", xlab="log10 pop size")
lm_var_per_strain_gene<-lm(Var_per_strain_genes[1:35]~GROUPS_texpsn)
abline(lm_var_per_strain_gene)
text(1.5,0.2, paste("p=", round(summary(lm_var_per_strain)$coefficients[2,4],digits=3), sep=""))
jh_test_var_gene<-jonckheere.test(Var_per_strain_genes[1:35],GROUPS_texpsn)
text(1.5, 0.35, paste("jh p=", round(jh_test_var_gene$p.value, digits=3), sep=""))




dev.copy(pdf, "FigureS2D.pdf")
dev.off()




plot(Var_per_strain_genes,Var_per_strain_tes, col=COLS_texpsn, pch=18,cex=2)

legend("topright", c("size1","size10","size100"),col=c("red","green","blue"),pch=18, cex=2)
dev.copy(pdf, "FigureS2E.pdf")
dev.off()



#figure 2
te_types<-table(Annotations_Tab_new[,2])

output_matrix_te_sums<-matrix(0, nrow=length(te_types), ncol=ncol(Tab_new))

for(i in 1:nrow(te_types)){individual_tes<-Annotations_Tab_new[which(Annotations_Tab_new[,2]==names(te_types)[i]),1]
A<-which(row.names(Tab_new)%in%individual_tes==T)
if(length(A)>1){
output_matrix_te_sums[i,]<-apply(Tab_new[A,],2,sum)}
if(length(A)==1){output_matrix_te_sums[i,]<-Tab_new[A,]}
}
row.names(output_matrix_te_sums)<-names(te_types)

output_matrix_te_sums<-output_matrix_te_sums[which(apply(output_matrix_te_sums,1,max)>0),]
output_matrix_correlations<-matrix(1,ncol=nrow(output_matrix_te_sums),nrow=nrow(output_matrix_te_sums))

for(i in 1:nrow(output_matrix_te_sums)){
for(j in 1:nrow(output_matrix_te_sums)){
if(!i==j){
A<-cor.test(output_matrix_te_sums[i,],output_matrix_te_sums[j,],method="spearman")
output_matrix_correlations[i,j]<-(-1*log10(A$p.value)*sign(A$estimate))

}else{output_matrix_correlations[i,j]<-10}

}
}

row.names(output_matrix_correlations)<-row.names(output_matrix_te_sums)
colnames(output_matrix_correlations)<-row.names(output_matrix_te_sums)

library(gplots)
heatmap.2(output_matrix_correlations, col=colorRampPalette(c("blue","white","orange"))(20), trace="none")

 dev.copy(pdf, "Figure2.pdf")
 dev.off()
 

#figure 3

load("cnv_pe_sorted.Rdata")
cnv_pe_2<-cnv_pe_2[1:181,]
#so that the cnv only contains TEs
cnv_mean<-apply(cnv_pe_2[,29:34],1,mean)

cnv_norm_mean<-cnv_pe_2/cnv_mean

cnv_total<-apply(cnv_pe_2,2,sum,na.rm=T)

groups_cnv_bp<-c(rep(0, 19),rep(1,10),rep(2,5))+1

cols_cnv<-COLS_texpsn[-1]
boxplot(log2(cnv_total[1:19]),log2(cnv_total[20:29]),log2(cnv_total[30:34]),outline=F,xlab="log10 population size", ylab="log2 te cnv", names=c(0,1,2),boxwex=0.5)
points(groups_cnv_bp,log2(cnv_total), col=cols_cnv, pch=18, cex=1.5)
 
lm_cnv<-lm(log2(cnv_total)~groups_cnv_bp)
jh_cnv<-jonckheere.test(cnv_total, groups_cnv_bp)

abline(lm_cnv)
text(2.5,18, paste("p=", format(summary(lm_cnv)$coefficients[2,4],digits=2,scientific=TRUE)))
text(3,16 , paste("p=",format(jh_cnv$p.value, digits=2, scientific=TRUE)))
dev.copy(pdf, "Figure3A_new2020.pdf")
dev.off()

cnv_to_sum_te<-match(names(sums_all_new_tab), names(cnv_total))


plot(log2(sums_all_new_tab),log2(cnv_total[cnv_to_sum_te]),xlab="log2 te expression", ylab="log2 te cnv", col=COLS_texpsn,cex=1.5,pch=18)

fit_xpsn_to_cnv<-lm(log2(cnv_total[cnv_to_sum_te])~log2(sums_all_new_tab))

abline(fit_xpsn_to_cnv)

text(5.2,0,paste("p=",round(summary(fit_xpsn_to_cnv)$coefficients[2,4],digits=2), sep=""))


legend("topleft", c("N.1","N.10","N.100"),col=c("red","green","blue"),pch=18,bty="n")
dev.copy(pdf, "Figure3B.pdf")
dev.off()
turmoil2_cnv<-cnv_norm_mean[grep("TURMOIL2", row.names(cnv_norm_mean)),]
a<-match(names(Turmoil2_total_expression), names(turmoil2_cnv))
plot(log2(Turmoil2_total_expression), log2(turmoil2_cnv)[a], col=cols_cnv, pch=18, cex=1.5, ylab="log2 Turmoilcnv", xlab="log2 turmoil expression")
legend("topright", c("N.1","N.10","N.100"), col=c("red","green","blue"), pch=18, cex=1.5, bty="n")
dev.copy(pdf, "Figure3C.pdf")
dev.off()



 
 #Figure4
 load("sRNA_clean_normalized_miRNA.Rdata")
 load("piRNA_targeted_genes_all.Rdata")
 
 piRNA_targeted_tes<-Tab_new[which(row.names(Tab_new)%in%piRNA_targets[,5]==T|row.names(Tab_new)%in%piRNA_targets[,13]==T),]
 non_piRNA_targeted_tes<-Tab_new[which(row.names(Tab_new)%in%piRNA_targets[,5]==F&row.names(Tab_new)%in%piRNA_targets[,13]==F),]
 sum_piRNA_targeted_tes_xpsn<-apply(piRNA_targeted_tes,2,sum)
 sum_nonpiRNA_targeted_tes<-apply(non_piRNA_targeted_tes,2,sum)
 boxplot(sum_piRNA_targeted_tes_xpsn[1:20],sum_piRNA_targeted_tes_xpsn[21:30],sum_piRNA_targeted_tes_xpsn[31:35],ylab="fpkm piRNA targets", xlab="log10 population size",outline=F,box.wex=0.5)
 points(GROUPS_texpsn_bp,sum_piRNA_targeted_tes_xpsn[1:35], col=COLS_texpsn, pch=18, cex=1.5, ylab="fpkm piRNA targets", xlab="log10 population size")
 
 lm_pitargs<-lm(sum_piRNA_targeted_tes_xpsn[1:35]~GROUPS_texpsn_bp)
 jh_test_pitargs<-jonckheere.test(sum_piRNA_targeted_tes_xpsn[1:35], GROUPS_texpsn)
 
 abline(lm_pitargs)
abline(h=sum_piRNA_targeted_tes_xpsn[36], col="gray")
 text(1.5,40, paste("p=", round(summary(lm_pitargs)$coefficients[2,4], digits=2), sep=""))
 text(1.5, 60, paste("jh p=", round(jh_test_pitargs$p.value, digits=2), sep=""))
 
 dev.copy(pdf, "Figure4A.pdf")
 dev.off()
 boxplot(sum_nonpiRNA_targeted_tes[1:20],sum_nonpiRNA_targeted_tes[21:30],sum_nonpiRNA_targeted_tes[31:35], box.wex=0.5, ylab="fpkm nonpiRNA targets", xlab="log10 population size")
  points(GROUPS_texpsn_bp,sum_nonpiRNA_targeted_tes[1:35], col=COLS_texpsn, pch=18, cex=1.5, ylab="fpkm piRNA targets", xlab="log10 population size")
lm_nontargs<-lm(sum_nonpiRNA_targeted_tes[1:35]~GROUPS_texpsn_bp)
jh_test_nontargs<-jonckheere.test(sum_nonpiRNA_targeted_tes[1:35], GROUPS_texpsn)
abline(lm_nontargs)
abline(h=sum_nonpiRNA_targeted_tes[36], col="gray")

text(1.5,15, paste("p=", round(summary(lm_nontargs)$coefficients[2,4], digits=2), sep=""))
text(1.5,30, paste("jh p=", round(jh_test_nontargs$p.value, digits=2), sep=""))

dev.copy(pdf, "Figure4B.pdf")
dev.off()


#supp fig expression hrde targets
hrde_targeted_tes<-Tab_new[which(row.names(Tab_new)%in%hrdeDown==T|row.names(Tab_new)%in%hrdeDown==T),]
 non_hrde_targeted_tes<-Tab_new[which(row.names(Tab_new)%in%hrdeDown==F&row.names(Tab_new)%in%hrdeDown==F),]
 sum_hrdeDown<-apply(hrde_targeted_tes,2,sum)
 sum_non_hrdeDown-apply(non_hrde_targeted_tes,2,sum)
 boxplot(sum_hrdeDown[1:20],sum_hrdeDown[21:30],sum_hrdeDown[31:35],ylab="fpkm hrde targets", xlab="log10 population size",outline=F,box.wex=0.5,ylim=c(0,2))
 points(GROUPS_texpsn_bp,sum_hrdeDown[1:35], col=COLS_texpsn, pch=18, cex=1.5, ylab="fpkm hrde targets", xlab="log10 population size")
 
 lm_hrdeTarg<-lm(sum_hrdeDown[1:35]~GROUPS_texpsn_bp)
 jh_hrdeTarg<-jonckheere.test(sum_hrdeDown[1:35], GROUPS_texpsn)
 
 abline(lm_hrdeTarg)
abline(h=sum_hrdeDown[36], col="gray")
 text(1.5,1.2, paste("p=", round(summary(lm_hrdeTarg)$coefficients[2,4], digits=5), sep=""))
 text(1.5, 1.6, paste("jh p=", round(jh_hrdeTarg$p.value, digits=7), sep=""))
 
 dev.copy(pdf, "Hrde_targets.pdf")
 dev.off()




piRNAs_sequences<-read.table("piRNA_sequences_all.fa", stringsAsFactors=F, fill=T)

load("piRNA_analysis_norm.Rdata")
piRNAs_targeting_tes<-piRNA_targets[which(piRNA_targets[,5]%in%row.names(piRNA_targeted_tes)==T),2]

piRNA_seqs_mapping_to_tes<-piRNAs_sequences[which(piRNAs_sequences[,1]%in%piRNAs_targeting_tes==T)-1,1]

piRNA_seqs_mapping_to_tes<-substr(piRNA_seqs_mapping_to_tes,2,nchar(piRNA_seqs_mapping_to_tes))

piRNA_reads_te_mapping<-piRNAs_norm[which(row.names(piRNAs_norm)%in%piRNA_seqs_mapping_to_tes==T),]
piRNA_sums_te_mapping<-apply(piRNA_reads_te_mapping,2,sum)
Groups_piRNA_clean<-c(rep(0, length=16), rep(1, length=8), rep(2,length=5))
Groups_piRNA_clean_bp<-Groups_piRNA_clean+1
cols_piRNA_clean<-Groups_piRNA_clean
cols_piRNA_clean[which(cols_piRNA_clean==0)]<-"red"
cols_piRNA_clean[which(cols_piRNA_clean==1)]<-"green"
cols_piRNA_clean[which(cols_piRNA_clean==2)]<-"blue"
boxplot(piRNA_sums_te_mapping[1:16],piRNA_sums_te_mapping[17:24], piRNA_sums_te_mapping[25:29], box.wex=0.5, outline=F, ylab="piRNA_reads_DEseq_normalized", xlab="log10 population size")
points(Groups_piRNA_clean_bp, piRNA_sums_te_mapping[1:29], col=cols_piRNA_clean, cex=1.5, pch=18, ylab="piRNA reads DEseq normalized", xlab="log10 population size")
abline(h=piRNA_sums_te_mapping[30], col="gray")
lm_piRNAs_all_tes<-lm(piRNA_sums_te_mapping[1:29]~Groups_piRNA_clean_bp)

abline(lm_piRNAs_all_tes)
text(1.5, 1700, paste("p=", round(summary(lm_piRNAs_all_tes)$coefficients[2,4],digits=3)))

jh_pis_all<-jonckheere.test(piRNA_sums_te_mapping[1:29], Groups_piRNA_clean)

text(1.5, 2500, paste("jh p=", round(jh_pis_all$p.value, digits=3), sep=""))

dev.copy(pdf, "Figure_4C.pdf")
dev.off()
#sup fig_total_piRNAs
piRNA_sums_all<-apply(piRNAs_norm,2,sum)
> boxplot(piRNA_sums_all[1:16],piRNA_sums_all[17:24],piRNA_sums_all[25:29])
> boxplot(piRNA_sums_all[1:16],piRNA_sums_all[17:24],piRNA_sums_all[25:29], ylim=c(0,10000), ylab="total_piRNA_seqs_DEseq_norm",  xlab="log10 population size")
> boxplot(piRNA_sums_all[1:16],piRNA_sums_all[17:24],piRNA_sums_all[25:29], ylim=c(0,100000), ylab="total_piRNA_seqs_DEseq_norm",  xlab="log10 population size")
> boxplot(piRNA_sums_all[1:16],piRNA_sums_all[17:24],piRNA_sums_all[25:29], ylim=c(0,100000), ylab="total_piRNA_seqs_DEseq_norm",  xlab="log10 population size",outline=F)
> boxplot(piRNA_sums_all[1:16],piRNA_sums_all[17:24],piRNA_sums_all[25:29], ylab="total_piRNA_seqs_DEseq_norm",  xlab="log10 population size",outline=F)
> points(Groups_piRNA_clean_bp, piRNA_sums_all[1:29], col=cols_piRNA_clean, cex=1.5,pch=18)
> abline(h=piRNA_sums_all[30], col="gray")










turmoil_piRNA_seqs<-piRNA_targets[which(piRNA_targets[,5]%in%Turmoil2_genes==T),2]

piRNA_seqs_mapping_to_turmoil<-piRNAs_sequences[which(piRNAs_sequences[,1]%in%turmoil_piRNA_seqs==T)-1,1]

piRNA_seqs_mapping_to_turmoil<-substr(piRNA_seqs_mapping_to_turmoil,2,nchar(piRNA_seqs_mapping_to_turmoil))

piRNA_reads_turmoil_mapping<-piRNAs_norm[which(row.names(piRNAs_norm)%in%piRNA_seqs_mapping_to_turmoil==T),]

piRNA_turmoil_mapping_sums<-apply(piRNA_reads_turmoil_mapping,2,sum)

plot(Groups_piRNA_clean,piRNA_turmoil_mapping_sums[1:29], col=cols_piRNA_clean,cex=1.5,pch=18,ylab="turmoil piRNA reads DEseq normalized", xlab="log10 population size")
abline(piRNA_turmoil_mapping_sums[30], col="gray")

lm_turmoil_only_pis<-lm(piRNA_turmoil_mapping_sums[1:29]~Groups_piRNA_clean)

abline(lm_turmoil_only_pis)

text(1.5,400, paste("p=", round(summary(lm_turmoil_only_pis)$coefficients[2,4],digits=3)))

jh_turmoil_only_pis<-jonckheere.test(piRNA_turmoil_mapping_sums[1:29], Groups_piRNA_clean)

text(1.5,600, paste("jh p=", round(jh_turmoil_only_pis$p.value, digits=3)))

dev.copy(pdf, "Figure4D.pdf")
dev.off()




sRNA_clean_TEs<-sRNA_clean_normalized_miRNA[which(row.names(sRNA_clean_normalized_miRNA)%in%Annotations_Tab_new[,1]==T),]










sRNA_clean_te_sums<-apply(sRNA_clean_TEs,2,sum)
Groups_sRNA_clean<-c(rep(0, length=16), rep(1, length=8), rep(2,length=5)) 
Groups_sRNA_clean_bp<-Groups_sRNA_clean+1
cols_sRNA_clean<-Groups_sRNA_clean
cols_sRNA_clean[which(cols_sRNA_clean==0)]<-"red"
cols_sRNA_clean[which(cols_sRNA_clean==1)]<-"green"
cols_sRNA_clean[which(cols_sRNA_clean==2)]<-"blue"
boxplot(sRNA_clean_te_sums[1:16], sRNA_clean_te_sums[17:24], sRNA_clean_te_sums[25:29], box.wex=0.5, ylab="total 22G-RNA reads normalized", xlab="log10 population size")
points(Groups_sRNA_clean_bp,sRNA_clean_te_sums[1:29], col=cols_sRNA_clean, pch=18, cex=1.5, ylab="total 22G-RNA reads normalized", xlab="log10 population size")

abline(h=sRNA_clean_te_sums[30], col="gray")
lm_all_tes<-lm(sRNA_clean_te_sums[1:29]~Groups_sRNA_clean_bp)
abline(lm_all_tes)
text(1.5,20000, paste("p=", round(summary(lm_all_tes)$coefficients[2,4],digits=2),  sep=""))
jh_test_all_tes<-jonckheere.test(sRNA_clean_te_sums[1:29], Groups_sRNA_clean)

text(1.5, 40000, paste("jh p=", round(jh_test_all_tes$p.value, digits=2),sep=""))
dev.copy(pdf, "Figure4E.pdf")
dev.off()

jh_tests_individual_tes<-matrix(0, ncol=3, nrow=nrow(sRNA_clean_TEs))

for(i in 1:nrow(sRNA_clean_TEs)){
										jh_t<-jonckheere.test(sRNA_clean_TEs[i,1:29], Groups_sRNA_clean)
										lm_t<-lm(sRNA_clean_TEs[i,1:29]~Groups_sRNA_clean)
										
										jh_tests_individual_tes[i,1]<-jh_t$p.value
										jh_tests_individual_tes[i,2]<-summary(lm_t)$coefficients[2,4]
										jh_tests_individual_tes[i,3]<-lm_t$coefficients[2]
										
										}

row.names(jh_tests_individual_tes)<-row.names(sRNA_clean_TEs)

colnames(jh_tests_individual_tes)<-c("jh test", "lm p", "lm gradient")

plot(jh_tests_individual_tes[,3], -1*log10(jh_tests_individual_tes[,2]), cex=1.5,pch=18, ylab="-log10 p value", xlab="gradient")


abline(h=(-1)*log10(0.05),lty=2)

TE_cats<-list("TC1","TC4",c("MARINER","TC3","TC2","TC5"),"TURMOIL1",c("LINE","VINGI","NESL","RTE"),"CER")
 names(TE_cats)<-c("TC1","TC4","TC/MARINER(other)","TURMOIL1","non-LTR","LTR")
 TE_cols<-c("purple","hotpink","red","green","cyan")
 for(i in 1:length(TE_cats)){A<-c();for(j in 1:length(TE_cats[[i]])){A<-c(A,grep(TE_cats[[i]][j], Annotations_Tab_new[,2]))}
 A<-which(row.names(jh_tests_individual_tes)%in%Annotations_Tab_new[A,1]==T)
 
 points(jh_tests_individual_tes[A,3],(-1)*log10(jh_tests_individual_tes[A,2]),pch=18,cex=1.5,col=TE_cols[i])}
 
 
 legend("bottomleft", names(TE_cats), col=TE_cols,cex=1.2,pch=18)
 dev.copy(pdf, "Figure4F.pdf")
dev.off()


#supplementary figure sRNAs against hrde1 22G-depleted TEs
hrde1Down_22Gs<-apply(sRNA_clean_normalized_miRNA[which(row.names(sRNA_clean_normalized_miRNA)%in%hrdeDown==T),],2,sum)

boxplot(hrde1Down_22Gs[1:16], hrde1Down_22Gs[17:24], hrde1Down_22Gs[25:29], box.wex=0.5, ylab="hrde1_down 22G-RNA readsnormalized", xlab="log10 population size")
points(Groups_sRNA_clean_bp,hrde1Down_22Gs[1:29], col=cols_sRNA_clean, pch=18, cex=1.5, ylab="hrde1_ddown 22G-RNA reads normalized", xlab="log10 population size")

lm_hrde<-lm(hrde1Down_22Gs[1:29]~Groups_sRNA_clean_bp)
jh_hrde<-jonckheere.test(hrde1Down_22Gs[1:29],Groups_sRNA_clean_bp)

abline(lm_hrde)

text(3,8000, "jh p=0.044")
text(3,6000, "LM p=0.021")
dev.copy(pdf, "hrde1_down_22Gs.pdf")
dev.off()

#supplementary figure sRNAs against piRNA targeted TEs

piRNA_targeted_TEs_22Gs<-apply(sRNA_clean_normalized_miRNA[which(row.names(sRNA_clean_normalized_miRNA)%in%row.names(piRNA_targeted_tes)==T),],2,sum)
boxplot(piRNA_targeted_TEs_22Gs[1:16],piRNA_targeted_TEs_22Gs[17:24],piRNA_targeted_TEs_22Gs[25:29],box.wex=0.5, ylab="piRNA_targ_TEs_norm22Gs",xlab="log10_pop_size",names=c(0,1,2),outline=F)
points(Groups_sRNA_clean_bp, piRNA_targeted_TEs_22Gs[1:29],col=cols_sRNA_clean,pch=18,cex=1.5)
lm_pi_22G<-lm(piRNA_targeted_TEs_22Gs[1:29]~Groups_sRNA_clean_bp)
text(3,2000, "lm p=0.074")
text(3,20000,"jh p=0.055")


#Extra figure (not included in final figure set)
Line_elements<-grep("LINE", Annotations_Tab_new[,2])
Line_elements<-c(Line_elements, grep("VINGI", Annotations_Tab_new[,2]))
Line_elements<-c(Line_elements, grep("NESL-1", Annotations_Tab_new[,2]))
Line_elements<-Annotations_Tab_new[Line_elements,1]
Sums_sRNA_lines<-apply(sRNA_clean_TEs[which(row.names(sRNA_clean_TEs)%in%Line_elements==T),],2,sum)
#Supplementary Figure 3
plot(Groups_sRNA_clean, Sums_sRNA_lines[1:29], pch=18,col=cols_sRNA_clean, cex=1.5, ylab="total line 22G-RNA reads normalized" , xlab="log10 population size")
abline(h=Sums_sRNA_lines[30], col="gray")
lm_all_lines_sRNA<-lm(Sums_sRNA_lines[1:29]~Groups_sRNA_clean)
abline(lm_all_lines_sRNA)
text(1.5,2200, paste("p=", round(summary(lm_all_lines_sRNA)$coefficients[2,4], digits=3), sep=""))

jh_test_all_lines<-jonckheere.test(Sums_sRNA_lines[1:29], Groups_sRNA_clean)
text(1.5,7000, paste("jh p=", round(jh_test_all_lines$p.value, digits=3), sep=""))

dev.copy(pdf, "Extra_figure.pdf")
dev.off()

#Figure 4G


Sums_sRNA_turmoil2<-apply(sRNA_clean_TEs[which(row.names(sRNA_clean_TEs)%in%Turmoil2_genes==T),],2,sum)

q<-match(names(Sums_sRNA_turmoil2), names(Turmoil2_total_expression))

plot(log2(Sums_sRNA_turmoil2),log2(Turmoil2_total_expression)[q], col=cols_sRNA_clean,pch=18, cex=1.5, ylab="log2 Turmoil2 RNAseq FPKM", xlab="log2 Turmoil2 22G-RNA" )
legend("bottomleft", c("N.1","N.10","N.100"), col=c("red","green","blue"), pch=18, cex=1.5,bty="n")
lm_fit_sRNAturmoil<-lm(log2(Turmoil2_total_expression)[q]~log2(Sums_sRNA_turmoil2))
abline(lm_fit_sRNAturmoil)
text(7,4, paste("LR p=", format(summary(lm_fit_sRNAturmoil)$coefficients[2,4], scientific=T, digits=2)))

dev.copy(pdf, "Figure4G.pdf")
dev.off()

#supplemental figure 3

setwd("piRNA_mismatch_analysis")

Mismatch.reads<-function(file_in){A<-read.table(file_in, sep="_", stringsAsFactors=F, fill=T)

A<-A[grep(">", A[,1]),]

A_exact<-A[which(A[,3]=="-1"),]

A_mm<-A[-which(A[,3]=="-1"),]

A_real_mm<-A_mm[which(A_mm[,3]%in%A_exact[,5]==F&A_mm[,4]>0),]

if(dim(A_real_mm)[1]>0){Tot_mm<-nrow(A_real_mm)

Tot_reads_mm<-unlist(strsplit(A_real_mm[,1], "-"))

Tot_reads_mm<-sum(as.numeric(Tot_reads_mm[-grep(">", Tot_reads_mm)]))}


else{Tot_reads_mm<-0

Tot_mm<-0}


return(c(nrow(A_exact),Tot_mm,Tot_reads_mm))

}

mm_files<-list.files()[grep(".fa", list.files())]
Output_all<-c()
for(i in 1:length(mm_files)){Output_all<-rbind(Output_all, Mismatch.reads(mm_files[i]))}
row.names(Output_all)<-mm_files

mm_reads_per_loci<-Output_all[,3]/Output_all[,2]

ne1<-grep("1[A-Z]", names(mm_reads_per_loci))
ne10<-grep("10[A-Z]", names(mm_reads_per_loci))
ne100<-grep("100[A-Z]", names(mm_reads_per_loci))

mm_reads_per_loci_ordered<-mm_reads_per_loci[c(ne1,ne10,ne100)]
groups_reads_per_loci<-c(rep(0, length=length(ne1)), rep(1,length=length(ne10)),rep(2,length=length(ne100)))
cols_mm<-groups_reads_per_loci
cols_mm[cols_mm==0]<-"red"
cols_mm[cols_mm==1]<-"green"
cols_mm[cols_mm==2]<-"blue"

plot(groups_reads_per_loci, mm_reads_per_loci_ordered, pch=18, col=cols_mm, cex=1.5, xlab="log10 population size", ylab="piRNA mismatch reads/mismatch only locus detected")

abline(h=mm_reads_per_loci[grep("PMA", names(mm_reads_per_loci))], col="gray")

lm_mm<-lm(mm_reads_per_loci_ordered~groups_reads_per_loci)
abline(lm_mm)
text(1.5,1.4, paste("p=", round(summary(lm_mm)$coefficients[2,4], digits=3)))

jh_test_mm<-jonckheere.test(mm_reads_per_loci_ordered,groups_reads_per_loci)
text(1.5, 2.2, paste("jh p=", round(jh_test_mm$p.value, digits=3)))

dev.copy(pdf, "FigureS3.pdf")
dev.off()





#Figure 5
load("Regulated_domain_TEs.Rdata")
Reg_TEs<-Reg_TEs[,9]
load("Active_domain_TEs.Rdata")
load("Het_tes.Rdata")

Reg_xpsn<-apply(Tab_new[which(row.names(Tab_new)%in%Reg_TEs==T&row.names(Tab_new)%in%piRNA_targets[,13]==T),],2,sum)
boxplot(Rex_xpsn[1:20],Reg_xpsn[21:30],Reg_xpsn[31:35], names=c(0,1,2), ylab="Regulated domain TE FPKM",xlab="log10 population size",box.wex=0.5)
 points(GROUPS_texpsn_bp,Reg_xpsn[1:35], col=COLS_texpsn, pch=18,cex=1.5, ylab="Regulated domain TEs fpkm", xlab="log10 population size")

lm_reg<-lm(Reg_xpsn[1:35]~GROUPS_texpsn_bp)
jh_reg<-jonckheere.test(Reg_xpsn[1:35], GROUPS_texpsn)

abline(h=Reg_xpsn[36], col="gray")
abline(lm_reg)

text(1.5,15, paste("p=",round(summary(lm_reg)$coefficients[2,4], digits=2), sep=""))
text(1.5,40, paste("jh p=", round(jh_reg$p.value,digits=2), sep=""))
dev.copy(pdf, "Figure5A.pdf")
dev.off()

Het_xpsn<-apply(Tab_new[which(row.names(Tab_new)%in%Heterochromatin_tes==T&row.names(Tab_new)%in%piRNA_targets[,13]==T),],2,sum)
boxplot(Het_xpsn[1:20],Het_xpsn[21:30], Het_xpsn[31:35],names=c(0,1,2), ylab="Heterochromatin TE FPKM", xlab="log10 population size")
points(GROUPS_texpsn_bp, Het_xpsn[1:35], col=COLS_texpsn, pch=18, cex=1.5, ylab="Heterochromatin TEs fpkm", xlab="log10 population size",ylim=c(0,40))
abline(h=Het_xpsn[36],col="gray")

lm_het<-lm(Het_xpsn[1:35]~GROUPS_texpsn_bp)
jh_het<-jonckheere.test(Het_xpsn[1:35],GROUPS_texpsn)

abline(lm_het)
text(1.5, 5, paste("p=", round(summary(lm_het)$coefficients[2,4], digits=2), sep=""))
text(1.5, 32, paste("jh p=", round(jh_het$p.value, digits=2), sep=""))

dev.copy(pdf, "Figure5B.pdf")
dev.off()

Active_xpsn<-apply(Tab_new[which(row.names(Tab_new)%in%Active_TEs==T&row.names(Tab_new)%in%piRNA_targets[,13]==T),],2,sum)

plot(GROUPS_texpsn, Active_xpsn[1:35], col=COLS_texpsn, pch=18, cex=1.5, ylab="Active domain TEs fpkm", xlab="log10 population size")

lm_active<-lm(Active_xpsn[1:35]~GROUPS_texpsn)

jh_active<-jonckheere.test(Active_xpsn[1:35], GROUPS_texpsn)

abline(h=Active_xpsn[36], col="gray")
abline(lm_active)

text(1.5, 15, paste("p=", round(summary(lm_active)$coefficients[2,4], digits=3)))

text(1.5, 50, paste("jh p=", round(jh_active$p.value, digits=3)))

dev.copy(pdf, "Figure5C.pdf")
dev.off()

sRNA_clean_TEs<-sRNA_clean_normalized_miRNA[which(row.names(sRNA_clean_normalized_miRNA)%in%row.names(Tab_new)==T),]

sRNAs_tes_reg_sums<-apply(sRNA_clean_TEs[which(row.names(sRNA_clean_TEs)%in%Reg_TEs==T),],2,sum)
plot(Groups_sRNA_clean, sRNAs_tes_reg_sums[1:29], pch=18,col=cols_sRNA_clean, cex=1.5, ylab="total Regulated TE piRNA targets 22G-RNA reads normalized" , xlab="log10 population size")


lm_reg_sRNA<-lm(sRNAs_tes_reg_sums[1:29]~Groups_sRNA_clean)

abline(lm_reg_sRNA)
text(1.5, 7500, paste("p=", round(summary(lm_reg_sRNA)$coefficients[2,4], digits=3)))




abline(h=sRNAs_tes_reg_sums[30], col="gray")
jh_test_sRNA<-jonckheere.test(sRNAs_tes_reg_sums[1:29], Groups_sRNA_clean)

text(1.5, 15000, paste("jh p=", round(jh_test_sRNA$p.value, digits=3)))

dev.copy(pdf, "Figure5D.pdf")
dev.off()

sRNAs_test_het_sums<-apply(sRNA_clean_TEs[which(row.names(sRNA_clean_TEs)%in%Heterochromatin_tes==T&row.names(sRNA_clean_TEs)%in%piRNA_targets[,13]==T),],2,sum)

plot(Groups_sRNA_clean, sRNAs_test_het_sums[1:29], pch=18, col=cols_sRNA_clean, cex=1.5, ylab="total Heterochromatin TE piRNA targets 22G-RNAs normalized", xlab="log10 population size")
abline(h=sRNAs_test_het_sums[30], col="gray")

lm_het_sRNA<-lm(sRNAs_test_het_sums[1:29]~Groups_sRNA_clean)
abline(lm_het_sRNA)

text(1.5, 3200, paste("p=", round(summary(lm_het_sRNA)$coefficients[2,4], digits=3)))
jh_het_sRNA<-jonckheere.test(sRNAs_test_het_sums[1:29], Groups_sRNA_clean)

text(1.5, 7500, paste("jh p=", round(jh_het_sRNA$p.value, digits=3)))
dev.copy(pdf, "Figure5E.pdf")
dev.off()

sRNAs_tes_active_sums<-apply(sRNA_clean_TEs[which(row.names(sRNA_clean_TEs)%in%Active_TEs==T&row.names(sRNA_clean_TEs)%in%piRNA_targets[,13]==T),],2,sum)

plot(Groups_sRNA_clean, sRNAs_tes_active_sums[1:29], pch=18, col=cols_sRNA_clean, cex=1.5, ylab="total Active TE piRNA targets 22G-RNAs normalized", xlab="log10 population size")
lm_act_sRNA<-lm(sRNAs_tes_active_sums[1:29]~Groups_sRNA_clean)
jh_act_sRNA<-jonckheere.test(sRNAs_tes_active_sums[1:29], Groups_sRNA_clean)
abline(lm_act_sRNA)
abline(h=sRNAs_tes_active_sums[36], col="gray")

text(1.5, 4000,paste("p=", round(summary(lm_act_sRNA)$coefficients[2,4], digits=3)))

text(1.5, 8000, paste("jh p=", round(jh_act_sRNA$p.value,digits=3)))
dev.copy(pdf, "Figure5F.pdf")
dev.off()


#Figure 6 PATC data

PATC<-read.table("TE_PATC_scores.txt", sep="\t", stringsAsFactors=F)

PATC_covered_xpsn<-PATC[which(PATC[,1]%in%row.names(Tab_new)==T),]

#sort into bins

PATC_covered_xpsn<-PATC_covered_xpsn[order(PATC_covered_xpsn[,4], decreasing=T),]

patc_bin<-c()

n<-1
while(length(patc_bin)<=nrow(PATC_covered_xpsn)){for(i in 1:21){patc_bin<-c(patc_bin, n)};n<-n+1}

patc_bin<-patc_bin[which(patc_bin<5)]
length(patc_bin)
patc_bin<-c(patc_bin,4,4,4)
PATC_covered_xpsn<-cbind(PATC_covered_xpsn,patc_bin)
patc_bins_te<-list()
for(i in 1:4){A<-which(row.names(Tab_new)%in%PATC_covered_xpsn[patc_bin==i,1]==T);patc_bins_te[[i]]<-apply(Tab_new[A,],2,sum)}


#expression bins

plotfunc<-function(data_in,groups){tabtemp<-cbind(data_in[1:length(groups)], groups);plot(tabtemp[,2],tabtemp[,1], col=c(rep("red", length=length(which(groups==0))),rep("green", length=length(which(groups==1))),rep("blue", length=length(which(groups==2)))),pch=18,cex=1.5,xlab="log_popsize", ylab="fpkm",ylim=c(min(data_in)-1,max(data_in)+1));abline(h=data_in[36],col="grey");fit<-lm(tabtemp[,1]~tabtemp[,2]);abline(fit);return(summary(fit))}


FITS_patc_bins<-list();for(i in 1:4){FITS_patc_bins[[i]]<-plotfunc(patc_bins_te[[i]],patc_bin_plot_1[,2]);dev.copy(pdf, paste("TE_expression_patc_bins",i,".pdf", sep=""));dev.off()}
##These figures then assembled into Figure 6A

jonck_tests_for_patc<-c()

for(i in 1:4){jonck_tests_for_patc<-c(jonck_tests_for_patc, jonckheere.test(patc_bins_te[[i]][1:35], g=GROUPS_texpsn)$p.value)}


All_data_normalized_miRNA_orig<-All_data_normalized_miRNA
All_data_normalized_miRNA<-sRNA_clean_TEs


for(i in 1:4){

				A<-which(row.names(sRNA_clean_TEs)%in%PATC_covered_xpsn[which(PATC_covered_xpsn[,5]==i),1])
				sRNA_data_temp<-apply(sRNA_clean_TEs[A,],2,sum)
				
				plotfunc(sRNA_data_temp,Groups_sRNA_clean)
				jh_t<-jonckheere.test(sRNA_data_temp[1:29],Groups_sRNA_clean)$p.value
				lm_t<-lm(sRNA_data_temp[1:29]~Groups_sRNA_clean)
				abline(h=sRNA_data_temp[30], col="gray")
				text(1.5,max(sRNA_data_temp)-200, paste(paste("jh p=", round(jh_t, digits=3),sep=""), paste("lm p=",round(summary(lm_t)$coefficients[2,4], digits=3),sep=""),sep="\n"))
				
				
				dev.copy(pdf, paste("sRNAs_TEs_PATC_bin_",i, ".pdf", sep=""))
				dev.off()
				}
##these figures then assembled into Figure 6B


#Figures 6C and D assembled from the following pdfs
domain_list<-list(Reg_TEs, Heterochromatin_tes, Active_TEs)
names(domain_list)<-c("Reg", "Het", "Act")
Top_patc_sRNA<-sRNA_clean_TEs[which(row.names(sRNA_clean_TEs)%in%PATC_covered_xpsn[which(PATC_covered_xpsn[,5]==1),1]==T),]

Lowest_patc_sRNA<-sRNA_clean_TEs[which(row.names(sRNA_clean_TEs)%in%PATC_covered_xpsn[which(PATC_covered_xpsn[,5]==4),1]==T),]

for(i in 1:3){

A<-which(row.names(Top_patc_sRNA)%in%domain_list[[i]]==T)
B<-which(row.names(Lowest_patc_sRNA)%in%domain_list[[i]]==T)

sum_top_temp<-apply(Top_patc_sRNA[A,],2,sum)
sum_lowest_temp<-apply(Lowest_patc_sRNA[B,],2,sum)


plotfunc(sum_top_temp,Groups_sRNA_clean)
jh_t_a<-jonckheere.test(sum_top_temp[1:29], Groups_sRNA_clean)$p.value
lm_t_a<-lm(sum_top_temp[1:29]~Groups_sRNA_clean)

abline(h=sum_top_temp[29], col="gray")

text(1.5, max(sum_top_temp)-200, paste(paste("jh p=", round(jh_t_a,digits=3), sep=""), paste("lm p=", round(summary(lm_t_a)$coefficients[2,4],digits=3), sep=""), sep="\n"))

dev.copy(pdf, paste("top_patc", names(domain_list)[i], "sRNAs.pdf", sep=""))
dev.off()


plotfunc(sum_lowest_temp, Groups_sRNA_clean)

jh_t_b<-jonckheere.test(sum_lowest_temp[1:29], Groups_sRNA_clean)$p.value

lm_t_b<-lm(sum_lowest_temp[1:29]~Groups_sRNA_clean)

abline(h=sum_lowest_temp[30], col="gray")
text(1.5, max(sum_lowest_temp)-200, paste(paste("jh p=", round(jh_t_b,digits=3), sep=""), paste("lm p=", round(summary(lm_t_b)$coefficients[2,4], digits=3)), sep="\n"))

dev.copy(pdf, paste("low_patc", names(domain_list)[i], "sRNA.pdf", sep=""))
dev.off()

}


#Supplemental Figure 6B

Top_patc_xpsn<-Tab_new[which(row.names(Tab_new)%in%PATC_covered_xpsn[patc_bin==1,1]==T),]
Low_patc_xpsn<-Tab_new[which(row.names(Tab_new)%in%PATC_covered_xpsn[patc_bin==4,1]==T),]

for(i in 1:3){

A<-which(row.names(Top_patc_xpsn)%in%domain_list[[i]]==T)
B<-which(row.names(Low_patc_xpsn)%in%domain_list[[i]]==T)

sum_top_temp<-apply(Top_patc_xpsn[A,],2,sum)
sum_lowest_temp<-apply(Low_patc_xpsn[B,],2,sum)


plotfunc(sum_top_temp,GROUPS_te_xpsn)
jh_t_a<-jonckheere.test(sum_top_temp[1:35], GROUPS_te_xpsn)$p.value
lm_t_a<-lm(sum_top_temp[1:35]~GROUPS_te_xpsn)

abline(h=sum_top_temp[36], col="gray")

text(1.5, max(sum_top_temp)-0.5, paste(paste("jh p=", round(jh_t_a,digits=3), sep=""), paste("lm p=", round(summary(lm_t_a)$coefficients[2,4],digits=3), sep=""), sep="\n"))

dev.copy(pdf, paste("top_patc", names(domain_list)[i], "xpsn.pdf", sep=""))
dev.off()



