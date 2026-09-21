#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(susieR); library(coloc)})
root <- normalizePath(getwd())
od <- file.path(root, "results/phase2BR")
fr <- fread(file.path(root, "config/phase2B_PILOT_LOCUS_FREEZE.tsv"), na.strings=c("", "NA"))
fr <- fr[pilot_id %in% sprintf("P2BP%02d", 2:8)]
qc <- fread(file.path(od, "PANUKBB_VARIANT_MATCH_QC.tsv"), na.strings=c("", "NA"))
meta <- fread(file.path(root, "config/phase1D_dataset_freeze.tsv"), na.strings=c("", "NA"))
old <- fread(file.path(root, "results/phase2B/SUSIE_LD_MISMATCH_DIAGNOSTICS.tsv"), na.strings=c("", "NA"))
old_pip <- fread(file.path(root, "results/phase2B/SUSIE_VARIANT_PIP.tsv"), na.strings=c("", "NA"))
old_fits <- readRDS(file.path(root, "results/phase2B/SUSIE_FITS.rds"))
access <- fread(file.path(od, "PANUKBB_ACCESS_AUDIT.tsv"), na.strings=c("", "NA"))
B <- as.numeric(access[check == "n_samples", value][1])
stopifnot(is.finite(B), B > 0, nrow(fr) == 7L)

cs <- function(f, ids) {
  if (inherits(f, "fit_error") || is.null(f$sets) || is.null(f$sets$cs)) return(list())
  lapply(f$sets$cs, function(x) ids[as.integer(x)])
}
pur <- function(f, R) tryCatch(as.data.frame(susie_get_cs(f, Xcorr=R)$purity), error=function(e) NULL)
fit <- function(b, se, R, n, rf, mm) tryCatch(
  susie_rss(bhat=b, shat=se, R=R, n=n, L=10, R_finite=rf, R_mismatch=mm,
            estimate_residual_variance=FALSE, max_iter=1000),
  error=function(e) structure(list(error=conditionMessage(e)), class="fit_error"))
fit_summary <- function(f, p, pair, trait, model, ids, R, N, Neff) {
  if (inherits(f, "fit_error")) return(data.table(pilot_id=p,pair=pair,trait=trait,model=model,
    n_variants=length(ids),N_total=N,N_effective=Neff,n_supplied=Neff,B=B,converged=FALSE,
    n_credible_sets=NA,effective_rank=NA,r_over_B=NA,Q_art=NA,R_sensitivity_flag=NA,
    R_reliability_flag=NA,B_corrected=NA,lambda_bias=NA,artifact_flag=NA,mode_label="ERROR",
    per_variable_penalty_median=NA,per_variable_penalty_95th=NA,per_variable_penalty_max=NA,
    purity_pass=FALSE,reliability_class="ERROR",error=f$error))
  d <- f$R_finite_diagnostics; pu <- pur(f,R)
  ppass <- !is.null(pu) && nrow(pu)>0 && all(is.finite(pu$min.abs.corr) & pu$min.abs.corr>=0.5)
  rf <- if (is.null(d$R_reliability_flag)) FALSE else isTRUE(d$R_reliability_flag)
  pen <- if (is.null(d$per_variable_penalty)) numeric() else as.numeric(d$per_variable_penalty)
  reliable <- isTRUE(f$converged) && length(cs(f,ids))>0 && ppass && !rf
  data.table(pilot_id=p,pair=pair,trait=trait,model=model,n_variants=length(ids),N_total=N,
    N_effective=Neff,n_supplied=Neff,B=B,converged=isTRUE(f$converged),n_credible_sets=length(cs(f,ids)),
    effective_rank=d$effective_rank,r_over_B=d$r_over_B,Q_art=d$Q_art,
    R_sensitivity_flag=d$R_sensitivity_flag,R_reliability_flag=rf,B_corrected=d$B_corrected,
    lambda_bias=d$lambda_bias,artifact_flag=d$artifact_flag,mode_label=d$mode_label,
    per_variable_penalty_median=if(length(pen))median(pen,na.rm=TRUE) else NA,
    per_variable_penalty_95th=if(length(pen))as.numeric(quantile(pen,.95,na.rm=TRUE,names=FALSE)) else NA,
    per_variable_penalty_max=if(length(pen))max(pen,na.rm=TRUE) else NA,purity_pass=ppass,
    reliability_class=if(rf)"LD_MISMATCH_LIMITED" else if(reliable)"RELIABLE" else "INDETERMINATE",error="")
}
cred <- function(f,p,pair,trait,model,ids,R,reliable) {
  z <- cs(f,ids); if(!length(z)) return(data.table())
  pu <- pur(f,R); ii <- if(!is.null(f$sets$cs_index)) f$sets$cs_index else seq_along(z)
  rbindlist(lapply(seq_along(z),function(i){x<-z[[i]]; pp<-f$pip[match(x,names(f$pip))]
    data.table(pilot_id=p,pair=pair,trait=trait,model=model,credible_set=paste0("L",i),
      lead_SNP=x[which.max(pp)],PIP=max(pp),credible_set_size=length(x),
      min_abs_LD=if(!is.null(pu))pu$min.abs.corr[i] else NA,
      mean_abs_LD=if(!is.null(pu))pu$mean.abs.corr[i] else NA,
      median_abs_LD=if(!is.null(pu))pu$median.abs.corr[i] else NA,
      logBF=if(length(f$lbf)>=ii[i])f$lbf[ii[i]] else NA,primary_reliable=reliable,
      cs_variants=paste(x,collapse=";"))}),fill=TRUE)
}
ser <- function(f,p,pair,trait,ids) {
  d <- f$R_finite_diagnostics
  if(is.null(d) || !isTRUE(d$R_reliability_flag) || is.null(d$ser_model)) return(NULL)
  s <- d$ser_model; z <- if(!is.null(s$sets) && !is.null(s$sets$cs))lapply(s$sets$cs,function(x)ids[as.integer(x)]) else list()
  k <- which.max(s$pip)
  data.table(pilot_id=p,pair=pair,trait=trait,model="SER_FALLBACK_DIAGNOSTIC",top_SER_PIP=max(s$pip),
    lead_variant=names(s$pip)[k],n_SER_credible_sets=length(z),SER_credible_set_width=if(length(z))length(z[[1]]) else NA,
    SER_credible_set=if(length(z))paste(z[[1]],collapse=";") else NA,diagnostic_only=TRUE)
}
old_top <- function(p,trait) {
  trait_name <- trait
  z <- old_pip[pilot_id==p & get("trait")==trait_name]
  if(!nrow(z)) return(NA_character_)
  z[which.max(PIP),variant_id]
}
old_fit_for <- function(p,trait) {
  x <- fr[pilot_id==p][1]; tr0 <- strsplit(x$pair,"-",fixed=TRUE)[[1]]
  old_fits[[p]][[if(tr0[1]==trait)"trait1" else "trait2"]][["mismatch"]]
}
old_fit_summary <- function(p,trait) {
  f <- old_fit_for(p,trait); ids0 <- old_fits[[p]]$variant_id; z <- cs(f,ids0)
  data.table(old_top_PIP=if(length(f$pip))max(f$pip,na.rm=TRUE) else NA_real_,
    old_credible_set_sizes=if(length(z))paste(lengths(z),collapse=";") else NA_character_)
}
all_diag <- list(); all_cred <- list(); all_pip <- list(); all_ser <- list(); fits <- list()
comparisons <- list(); coloc_out <- list(); priors <- list(); direction <- list()

for(ii in seq_len(nrow(fr))) {
  x <- fr[ii]; p <- x$pilot_id; q <- qc[pilot_id==p]
  stopifnot(nrow(q)==1L, q$coverage_gate=="PASS")
  ss <- fread(file.path(od,p,"SUMMARY_STATS_PANUKBB_ALIGNED.tsv"),na.strings=c("","NA")); ids <- ss$variant_id
  R <- as.matrix(fread(file.path(od,p,"PANUKBB_LD.tsv"),header=FALSE)); storage.mode(R)<-"double"; diag(R)<-1
  dimnames(R)<-list(ids,ids)
  stopifnot(nrow(R)==nrow(ss),ncol(R)==nrow(ss),all(is.finite(R)),isTRUE(all.equal(unname(R),unname(t(R)),tolerance=1e-10)))
  tr <- strsplit(x$pair,"-",fixed=TRUE)[[1]]; m1<-meta[trait==tr[1]][1]; m2<-meta[trait==tr[2]][1]
  N1<-as.numeric(m1$N);N2<-as.numeric(m2$N);E1<-as.numeric(m1$effective_N);E2<-as.numeric(m2$effective_N);s1<-as.numeric(m1$cases)/N1;s2<-as.numeric(m2$cases)/N2
  d1<-list(beta=ss$BETA_1,varbeta=ss$SE_1^2,snp=ids,position=ss$BP,type="cc",N=N1,s=s1,LD=R,trait=tr[1])
  d2<-list(beta=ss$BETA_2,varbeta=ss$SE_2^2,snp=ids,position=ss$BP,type="cc",N=N2,s=s2,LD=R,trait=tr[2])
  ck1<-tryCatch({check_dataset(d1,req="LD");"PASS"},error=function(e)paste0("FAIL: ",conditionMessage(e)))
  ck2<-tryCatch({check_dataset(d2,req="LD");"PASS"},error=function(e)paste0("FAIL: ",conditionMessage(e)))
  f1<-list(pan_r_finite=fit(ss$BETA_1,ss$SE_1,R,E1,B,"none"),pan_mismatch=fit(ss$BETA_1,ss$SE_1,R,E1,B,"eb"))
  f2<-list(pan_r_finite=fit(ss$BETA_2,ss$SE_2,R,E2,B,"none"),pan_mismatch=fit(ss$BETA_2,ss$SE_2,R,E2,B,"eb"))
  fits[[p]]<-list(trait1=f1,trait2=f2,variant_id=ids,summary=ss)
  for(j in 1:2) {
    fs<-if(j==1)f1 else f2; trait<-tr[j]; N<-if(j==1)N1 else N2; Ne<-if(j==1)E1 else E2
    pan_fit<-fs$pan_mismatch
    dr<-lapply(names(fs),function(md)fit_summary(fs[[md]],p,x$pair,trait,md,ids,R,N,Ne)); all_diag[[length(all_diag)+1]]<-rbindlist(dr,fill=TRUE)
    all_cred[[length(all_cred)+1]]<-rbindlist(lapply(seq_along(fs),function(k)cred(fs[[k]],p,x$pair,trait,names(fs)[k],ids,R,dr[[k]]$reliability_class=="RELIABLE")),fill=TRUE)
    z<-ser(pan_fit,p,x$pair,trait,ids); if(!is.null(z))all_ser[[length(all_ser)+1]]<-z
    new<-dr[[2]]; trait_name<-trait; oldz<-old[pilot_id==p & model=="mismatch" & get("trait")==trait_name]
    stopifnot(nrow(oldz)==1L)
    oldflag<-isTRUE(as.logical(oldz$R_reliability_flag[1])); newflag<-isTRUE(new$R_reliability_flag[1])
    qimp<-is.finite(oldz$Q_art[1])&&is.finite(new$Q_art[1])&&new$Q_art[1]<=.5*oldz$Q_art[1]
    odist<-if(is.finite(oldz$B_corrected[1])&&oldz$B_corrected[1]>0)abs(log(oldz$B_corrected[1]/as.numeric(oldz$B_reference[1]))) else Inf
    ndist<-if(is.finite(new$B_corrected[1])&&new$B_corrected[1]>0)abs(log(new$B_corrected[1]/B)) else Inf
    bimp<-is.finite(odist)&&is.finite(ndist)&&ndist<=.5*odist
    rc<-if(oldflag&&!newflag)"RESCUED" else if(oldflag&&newflag&&(qimp||bimp))"PARTIALLY_RESCUED" else "NOT_RESCUED"
    os<-old_fit_summary(p,trait); nz<-cs(pan_fit,ids)
    comparisons[[length(comparisons)+1]]<-data.table(pilot_id=p,pair=x$pair,trait=trait,old_model="1000G_mismatch",new_model="PanUKBB_mismatch",old_B=oldz$B_reference[1],new_B=B,
      old_effective_rank=oldz$effective_rank[1],new_effective_rank=new$effective_rank[1],old_r_over_B=oldz$r_over_B[1],new_r_over_B=new$r_over_B[1],old_Q_art=oldz$Q_art[1],new_Q_art=new$Q_art[1],old_B_corrected=oldz$B_corrected[1],new_B_corrected=new$B_corrected[1],old_R_sensitivity_flag=oldz$R_sensitivity_flag[1],new_R_sensitivity_flag=new$R_sensitivity_flag[1],old_R_reliability_flag=oldz$R_reliability_flag[1],new_R_reliability_flag=new$R_reliability_flag[1],old_n_credible_sets=oldz$n_credible_sets[1],new_n_credible_sets=new$n_credible_sets[1],old_top_PIP_SNP=old_top(p,trait),new_top_PIP_SNP=names(pan_fit$pip)[which.max(pan_fit$pip)],rescue_class=rc,Q_art_improved_50pct=qimp,B_ratio_improved_50pct=bimp)
    comparisons[[length(comparisons)]]$old_top_PIP<-os$old_top_PIP
    comparisons[[length(comparisons)]]$new_top_PIP<-if(length(pan_fit$pip))max(pan_fit$pip,na.rm=TRUE) else NA_real_
    comparisons[[length(comparisons)]]$old_credible_set_sizes<-os$old_credible_set_sizes
    comparisons[[length(comparisons)]]$new_credible_set_sizes<-if(length(nz))paste(lengths(nz),collapse=";") else NA_character_
    all_pip[[length(all_pip)+1]]<-data.table(pilot_id=p,pair=x$pair,trait=trait,variant_id=ids,PIP=pan_fit$pip,beta=ss[[paste0("BETA_",j)]],SE=ss[[paste0("SE_",j)]])
  }
  dr1<-tail(all_diag,2)[[1]][model=="pan_mismatch"]; dr2<-tail(all_diag,1)[[1]][model=="pan_mismatch"]
  r1<-nrow(dr1)==1L&&dr1$reliability_class=="RELIABLE"; r2<-nrow(dr2)==1L&&dr2$reliability_class=="RELIABLE"
  if(!r1||!r2||ck1!="PASS"||ck2!="PASS") {
    coloc_out[[length(coloc_out)+1]]<-data.table(pilot_id=p,pair=x$pair,signal1=NA,signal2=NA,PP.H3=NA,PP.H4=NA,H4_H3=NA,classification="LD_LIMITED",status=if(ck1!="PASS"||ck2!="PASS")"DATASET_QC_FAIL" else "LD_LIMITED",upgrade_label="NO_UPGRADE",direction_status="NOT_RUN")
    for(pr in c(1e-6,5e-6,1e-5))priors[[length(priors)+1]]<-data.table(pilot_id=p,pair=x$pair,signal1=NA,signal2=NA,p1=1e-4,p2=1e-4,p12=pr,PP.H3=NA,PP.H4=NA,H4_H3=NA,status="NOT_RUN_LD_LIMITED")
    next
  }
  z<-tryCatch(coloc.susie(f1$pan_mismatch,f2$pan_mismatch,p1=1e-4,p2=1e-4,p12=5e-6),error=function(e)NULL)
  if(is.null(z)||is.null(z$summary)){coloc_out[[length(coloc_out)+1]]<-data.table(pilot_id=p,pair=x$pair,signal1=NA,signal2=NA,PP.H3=NA,PP.H4=NA,H4_H3=NA,classification="INDETERMINATE",status="NO_RESULT",upgrade_label="NO_UPGRADE",direction_status="NOT_RUN");next}
  su<-as.data.table(z$summary)
  for(k in seq_len(nrow(su))) {
    h3<-if("PP.H3"%in%names(su))su$PP.H3[k] else if("PP.H3.abf"%in%names(su))su$PP.H3.abf[k] else NA_real_; h4<-if("PP.H4"%in%names(su))su$PP.H4[k] else if("PP.H4.abf"%in%names(su))su$PP.H4.abf[k] else NA_real_; h43<-if(is.finite(h3)&&h3>0)h4/h3 else Inf
    s1<-if("idx1"%in%names(su))su$idx1[k] else k;s2<-if("idx2"%in%names(su))su$idx2[k] else k
    for(pr in c(1e-6,5e-6,1e-5)) {
      zs<-tryCatch(coloc.susie(f1$pan_mismatch,f2$pan_mismatch,p1=1e-4,p2=1e-4,p12=pr),error=function(e)NULL)
      if(is.null(zs)||is.null(zs$summary)) {
        priors[[length(priors)+1]]<-data.table(pilot_id=p,pair=x$pair,signal1=s1,signal2=s2,p1=1e-4,p2=1e-4,p12=pr,PP.H3=NA,PP.H4=NA,H4_H3=NA,status="NO_RESULT")
      } else {
        ss<-as.data.table(zs$summary)
        rr<-if("idx1"%in%names(ss))ss[idx1==s1 & idx2==s2] else ss[1]
        if(!nrow(rr))rr<-ss[1]
        hh3<-if("PP.H3"%in%names(rr))rr$PP.H3[1] else if("PP.H3.abf"%in%names(rr))rr$PP.H3.abf[1] else NA_real_;hh4<-if("PP.H4"%in%names(rr))rr$PP.H4[1] else if("PP.H4.abf"%in%names(rr))rr$PP.H4.abf[1] else NA_real_
        priors[[length(priors)+1]]<-data.table(pilot_id=p,pair=x$pair,signal1=s1,signal2=s2,p1=1e-4,p2=1e-4,p12=pr,PP.H3=hh3,PP.H4=hh4,H4_H3=if(is.finite(hh3)&&hh3>0)hh4/hh3 else Inf,status="PASS")
      }
    }
    cl<-if(isTRUE(h4>=.8&&h43>=3))"STRONG_SHARED_SIGNAL" else if(isTRUE(h4>=.5))"MODERATE_SHARED_SIGNAL" else if(isTRUE(h3>=.8))"DISTINCT_SIGNALS" else "INDETERMINATE"
    coloc_out[[length(coloc_out)+1]]<-data.table(pilot_id=p,pair=x$pair,signal1=s1,signal2=s2,PP.H3=h3,PP.H4=h4,H4_H3=h43,classification=cl,status="PASS",upgrade_label="NO_UPGRADE",direction_status="NOT_RUN")
  }
}
bind<-function(x)if(length(x))rbindlist(x,fill=TRUE)else data.table()
fwrite(bind(all_diag),file.path(od,"PANUKBB_SUSIE_DIAGNOSTICS.tsv"),sep="\t",na="NA")
fwrite(bind(all_cred),file.path(od,"PANUKBB_CREDIBLE_SETS.tsv"),sep="\t",na="NA")
fwrite(bind(all_pip),file.path(od,"PANUKBB_VARIANT_PIP.tsv"),sep="\t",na="NA")
fwrite(bind(all_ser),file.path(od,"SER_FALLBACK_DIAGNOSTIC.tsv"),sep="\t",na="NA")
fwrite(bind(comparisons),file.path(od,"REFERENCE_COMPARISON.tsv"),sep="\t",na="NA")
fwrite(bind(coloc_out),file.path(od,"COLOC_SUSIE_PANUKBB.tsv"),sep="\t",na="NA")
fwrite(bind(priors),file.path(od,"COLOC_PRIOR_SENSITIVITY.tsv"),sep="\t",na="NA")
fwrite(bind(direction),file.path(od,"COLOC_DIRECTION_AUDIT.tsv"),sep="\t",na="NA")
fwrite(data.table(pilot_id=character(),trait=character(),check_dataset=character(),N_total=numeric(),N_effective=numeric(),s=numeric(),B=numeric(),L=integer()),file.path(od,"PANUKBB_DATASET_QC.tsv"),sep="\t",na="NA")
saveRDS(fits,file.path(od,"PANUKBB_SUSIE_FITS.rds"))
cat("Pan-UKBB rescue core complete; B=",B,"\n",sep="")
