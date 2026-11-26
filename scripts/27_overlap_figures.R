# Run from the project root after scripts/26_overlap_reanalysis.R.
library(ggplot2)
library(data.table)
library(patchwork)
r <- file.path("results", "overlap_reanalysis")
out <- file.path(r, "figures")
dir.create(out, showWarnings=FALSE)
read <- function(n) fread(file.path(r, paste0(n,".csv")))
ids <- c("IFX_ARJS_5","OSM_OSMR","IL22_NEUTROPHIL_6","VDZ_NR_3","VDZ_4","ETRO_2","VDZ_R_5")
labels <- c("Arijs anti-TNF","OSM-OSMR","IL-22 chemokines","VDZ non-response","VDZ 4-gene",
            "Etrolizumab","VDZ response")
prepare <- function(d) d[, signature := factor(signature_id, levels=rev(ids),labels=rev(labels))]
ink <- "#183D50"
accent <- "#B45322"
style <- theme_classic(base_size=12) +
  theme(axis.title=element_text(size=11), axis.text=element_text(color="#263C46"),
        plot.title=element_text(face="bold",size=12),
        panel.grid.major.x=element_line(color="#e2e9ea",linewidth=.3),
        strip.background=element_blank(), strip.text=element_text(face="bold"),
        legend.position="bottom",legend.title=element_blank(),
        plot.margin=margin(8,10,8,10))
forest <- function(d,x_label="Odds ratio per 1 SD",subtitle=NULL) {
  d <- prepare(copy(d))
  ggplot(d, aes(x=or,y=signature)) +
    geom_vline(xintercept=1,linetype=2,color="#728A94") +
    geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=.15,
                  color=ink,na.rm=TRUE,linewidth=.7) +
    geom_point(color=ink,size=2.7,na.rm=TRUE) +
    scale_y_discrete(drop=FALSE) + scale_x_log10() +
    labs(x=x_label,y=NULL,subtitle=subtitle) + style
}
write_plot <- function(p, file, w, h) {
  for (ext in c("pdf","tiff","png")) {
    args <- list(filename=file.path(out,paste0(file,".",ext)),plot=p,width=w,height=h,
                 units="in",bg="white",dpi=if(ext=="png") 170 else 300)
    if (ext=="pdf") args$device <- grDevices::pdf
    if (ext=="tiff") args$compression <- "lzw"
    do.call(ggsave,args)
  }
}
active <- read("active_meta")[scenario=="primary"]
external <- read("source_excluded_meta")[scenario=="primary"]
common_scale <- function(breaks) scale_x_log10(limits=range(breaks),breaks=breaks)
make_panel <- function(d,title,breaks) {
  p <- forest(d)
  p$scales$scales <- Filter(function(s)!any("x" %in% s$aesthetics),p$scales$scales)
  p + common_scale(breaks) + ggtitle(title)
}
ticks <- c(.5,1,2,4)
p1 <- make_panel(active,"A   Retained active arms",ticks)
p2 <- make_panel(external,"B   Publication sources excluded",ticks) +
  geom_text(data=prepare(copy(external[cohorts==1])),
            aes(x=1.38,y=signature,label="1 cohort only"),inherit.aes=FALSE,
            size=3.3,color="#5C7379",hjust=0)
write_plot(p1+p2,"Figure_2",9.0,4.0)

placebo <- read("placebo_meta")
interaction <- read("interaction_precision")
p3 <- make_panel(placebo,"A   Placebo arms",ticks)
p4 <- make_panel(interaction,"B   Treatment-by-signature interaction",ticks) +
  labs(x="Interaction odds ratio")
write_plot(p3+p4,"Figure_3",9.1,4.0)

stability <- prepare(read("stability")[analysis=="all_active"])
predict <- ggplot(stability,aes(or,signature)) +
  geom_vline(xintercept=1,linetype=2,color="#728A94") +
  geom_errorbar(aes(xmin=prediction_low,xmax=prediction_high,color="Prediction interval"),
                orientation="y",width=.18,linewidth=.6) +
  geom_errorbar(aes(xmin=ci_low,xmax=ci_high,color="Confidence interval"),
                orientation="y",width=0,linewidth=1.65) +
  geom_point(size=2.5,color=ink) +
  scale_color_manual(values=c("Prediction interval"="#93A9B0","Confidence interval"=ink)) +
  scale_x_log10(breaks=c(.5,1,2,4)) + labs(x="Odds ratio per 1 SD", y=NULL) + style
write_plot(predict,"Figure_S1",6.8,4.2)
limits <- ggplot(prepare(copy(interaction)),aes(or,signature)) +
  annotate("rect",xmin=2/3,xmax=1.5,ymin=-Inf,ymax=Inf,fill="#E0EEE9") +
  geom_vline(xintercept=1,linetype="dashed",color="#728A94") +
  geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=.14,color=ink,linewidth=.7) +
  geom_point(color=ink,size=2.6) + scale_x_log10(breaks=ticks) +
  labs(x="Interaction odds ratio",y=NULL,subtitle="Shading shows the exploratory interval from 0.67 to 1.50.") + style
write_plot(limits,"Figure_S3",6.8,4.1)
ep <- read("endpoint_stratified")
p5 <- make_panel(ep[stratum=="objective_or_composite_mucosal"],"A   Mucosal or composite",c(.3,.5,1,2,4,8))
p6 <- make_panel(ep[stratum=="clinical_only"],"B   Clinical response or remission",c(.3,.5,1,2,4,8))
write_plot(p5+p6,"Figure_S4",9.0,4)
sc <- prepare(read("active_meta"))
sc[, scheme:=factor(scenario,levels=c("primary","switch_ACT","switch_Leuven","switch_both"),
                    labels=c("Primary","Switch ACT 1","Switch Leuven IFX","Switch both"))]
swaps <- ggplot(sc,aes(or,signature)) +
  geom_vline(xintercept=1,linetype=2,color="#728A94") +
  geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=.14,color=ink) +
  geom_point(color=ink,size=2.3) +
  facet_wrap(~scheme,ncol=2) + scale_x_log10(breaks=c(.5,1,2,4)) +
  labs(y=NULL,x="Odds ratio per 1 SD") + style
write_plot(swaps,"Figure_S6",7.9,6.4)
