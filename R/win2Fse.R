#'Compute power for Simple Effects in Two Factor Within Subjects ANOVA with up to two by four levels.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alterative values may be entered by user
#'@param m1.1 Mean of first level factor 1, 1st level factor two
#'@param m2.1 Mean of second level factor 1, 1st level factor two
#'@param m3.1 Mean of third level factor 1, 1st level factor two
#'@param m4.1 Mean of fourth level factor 1, 1st level factor two
#'@param m1.2 Mean of first level factor 1, 2nd level factor two
#'@param m2.2 Mean of second level factor 1, 2nd level factor two
#'@param m3.2 Mean of third level factor 1, 2nd level factor two
#'@param m4.2 Mean of fourth level factor 1, 2nd level factor two
#'@param s1.1 Standard deviation of first level factor 1, 1st level factor two
#'@param s2.1 Standard deviation of second level factor 1, 1st level factor two
#'@param s3.1 Standard deviation of third level factor 1, 1st level factor two
#'@param s4.1 Standard deviation of forth level factor 1, 1st level factor two
#'@param s1.2 Standard deviation of first level factor 1, 2nd level factor two
#'@param s2.2 Standard deviation of second level factor 1, 2nd level factor two
#'@param s3.2 Standard deviation of third level factor 1, 2nd level factor two
#'@param s4.2 Standard deviation of forth level factor 1, 2nd level factor two
#'@param r12 correlation Factor 1, Level 1 and Factor 1, Level 2
#'@param r13 correlation Factor 1, Level 1 and Factor 1, Level 3
#'@param r14 correlation Factor 1, Level 1 and Factor 1, Level 4
#'@param r15 correlation Factor 1, Level 1 and Factor 2, Level 1
#'@param r16 correlation Factor 1, Level 1 and Factor 2, Level 2
#'@param r17 correlation Factor 1, Level 1 and Factor 2, Level 3
#'@param r18 correlation Factor 1, Level 1 and Factor 2, Level 4
#'@param r23 correlation Factor 1, Level 2 and Factor 1, Level 3
#'@param r24 correlation Factor 1, Level 2 and Factor 1, Level 4
#'@param r25 correlation Factor 1, Level 2 and Factor 2, Level 1
#'@param r26 correlation Factor 1, Level 2 and Factor 2, Level 2
#'@param r27 correlation Factor 1, Level 2 and Factor 2, Level 3
#'@param r28 correlation Factor 1, Level 2 and Factor 2, Level 4
#'@param r34 correlation Factor 1, Level 3 and Factor 1, Level 4
#'@param r35 correlation Factor 1, Level 3 and Factor 2, Level 1
#'@param r36 correlation Factor 1, Level 3 and Factor 2, Level 2
#'@param r37 correlation Factor 1, Level 3 and Factor 2, Level 3
#'@param r38 correlation Factor 1, Level 3 and Factor 2, Level 4
#'@param r45 correlation Factor 1, Level 4 and Factor 2, Level 1
#'@param r46 correlation Factor 1, Level 4 and Factor 2, Level 2
#'@param r47 correlation Factor 1, Level 4 and Factor 2, Level 3
#'@param r48 correlation Factor 1, Level 4 and Factor 2, Level 4
#'@param r56 correlation Factor 2, Level 1 and Factor 2, Level 2
#'@param r57 correlation Factor 2, Level 1 and Factor 2, Level 3
#'@param r58 correlation Factor 2, Level 1 and Factor 2, Level 4
#'@param r67 correlation Factor 2, Level 2 and Factor 2, Level 3
#'@param r68 correlation Factor 2, Level 2 and Factor 2, Level 4
#'@param r78 correlation Factor 2, Level 3 and Factor 2, Level 4
#'@param r sets same correlations between DVs on all factor levels (seriously, just use this)
#'@param s sets same standard deviation for factor levels (see comment above)
#'@param n Sample size for first group
#'@param alpha Type I error (default is .05)
#'@return Power for Simple Effects for Two Factor Within Subjects ANOVA
#'@export

win2Fse<-function(m1.1,m2.1,m3.1=NA,m4.1=NA,m1.2,m2.2,m3.2=NA,m4.2=NA,
                  s1.1=NA,s2.1=NA,s3.1=NA,s4.1=NA,s1.2=NA,s2.2=NA,s3.2=NA,s4.2=NA,
                  r12=NULL, r13=NULL, r14=NULL, r15=NULL, r16=NULL, r17=NULL, r18=NULL,
                  r23=NULL, r24=NULL, r25=NULL, r26=NULL, r27=NULL, r28=NULL,
                  r34=NULL, r35=NULL, r36=NULL, r37=NULL, r38=NULL,
                  r45=NULL, r46=NULL, r47=NULL, r48=NULL,
                  r56=NULL, r57=NULL, r58=NULL,
                  r67=NULL, r68=NULL,
                  r78=NULL, r=NULL, s = NULL, n, alpha=.05)

{

  levels<-NA
  levels[is.na(m4.1) & is.na(m3.1)]<-2
  levels[is.na(m4.1) & !is.na(m3.1)]<-3
  levels[!is.na(m4.1)]<-4

  if (levels=="4"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s3.1<-s;s4.1<-s;s1.2<-s;s2.2<-s;s3.2<-s;s4.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2;var5<-s^2;var6<-s^2;var7<-s^2;var8<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s4.1^2;var5<-s1.2^2;var6<-s2.2^2;var7<-s3.2^2;var8<-s4.2^2}
    r12<-r;r13<-r;r14<-r;r15<-r;r16<-r;r17<-r;r18<-r;r23<-r;r24<-r;r25<-r;r26<-r;r27<-r;r28<-r
    r34<-r;r35<-r;r36<-r;r37<-r;r38<-r;r45<-r;r46<-r;r47<-r;r48<-r;r56<-r;r57<-r;r58<-r
    r67<-r;r68<-r;r78<-r
    cov12<-r12*s1.1*s2.1;cov13<-r13*s1.1*s3.1;cov14<-r14*s1.1*s4.1;cov15<-r15*s1.1*s1.2;cov16<-r16*s1.1*s2.2;cov17<-r17*s1.1*s3.2;cov18<-r18*s1.1*s4.2
    cov23<-r23*s2.1*s3.1;cov24<-r24*s2.1*s4.1;cov25<-r25*s2.1*s1.2;cov26<-r26*s2.1*s2.2;cov27<-r27*s2.1*s3.2;cov28<-r28*s2.1*s4.2
    cov34<-r34*s3.1*s4.1;cov35<-r35*s3.1*s1.2;cov36<-r36*s3.1*s2.2;cov37<-r37*s3.1*s3.2;cov38<-r38*s3.1*s4.2
    cov45<-r45*s4.1*s1.2;cov46<-r46*s4.1*s2.2;cov47<-r47*s4.1*s3.2;cov48<-r48*s4.1*s4.2
    cov56<-r56*s1.2*s2.2;cov57<-r57*s1.2*s3.2;cov58<-r58*s1.2*s4.2
    cov67<-r67*s2.2*s3.2;cov68<-r68*s2.2*s4.2
    cov78<-r78*s3.2*s4.2
    out <- MASS::mvrnorm(n, mu = c(m1.1,m2.1,m3.1,m4.1,m1.2,m2.2,m3.2,m4.2),
                   Sigma = matrix(c(var1,cov12,cov13, cov14, cov15, cov16, cov17, cov18,
                                    cov12,var2,cov23, cov24, cov25, cov26, cov27, cov28,
                                    cov13, cov23,var3,cov34, cov35, cov36, cov37, cov38,
                                    cov14, cov24, cov34, var4, cov45, cov46, cov47, cov48,
                                    cov15, cov25, cov35, cov45, var5, cov56, cov57, cov58,
                                    cov16, cov26, cov36, cov46, cov56, var6, cov67, cov68,
                                    cov17, cov27, cov37, cov47, cov57, cov67, var7, cov78,
                                    cov18, cov28, cov38, cov48, cov58, cov68, cov78, var8), ncol = 8),
                   empirical = TRUE)

    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3 = V3, y4 = V4, y5 = V5, y6 = V6, y7 = V7, y8 = V8)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="time",value="dv",-id)
    out$time<-as.factor(out$time)
    out$time<-as.numeric(out$time)
    out$iv1<-NA
    out$iv1[out$time==1|out$time==5]<-1
    out$iv1[out$time==2|out$time==6]<-2
    out$iv1[out$time==3|out$time==7]<-3
    out$iv1[out$time==4|out$time==8]<-4
    out$iv2<-NA
    out$iv2[out$time==1|out$time==2|out$time==3|out$time==4]<-1
    out$iv2[out$time==5|out$time==6|out$time==7|out$time==8]<-2
    out$iv1<-as.ordered(out$iv1)
    out$iv2<-as.ordered(out$iv2)
    options(contrasts=c("contr.helmert", "contr.poly"))

    #split stuff here...
    data.ab1<-subset(out, iv2==1)
    modelab1<-ez::ezANOVA(data=data.ab1, dv=.(dv), wid=.(id), within = .(iv1), type=3, detailed=TRUE)
    dfab1<-modelab1$ANOVA$DFn[2]
    dfWab1<-modelab1$ANOVA$DFd[2]
    SSab1<-modelab1$ANOVA$SSn[2]
    SSWab1<-modelab1$ANOVA$SSd[2]
    eta2ab1<-SSab1/(SSab1+SSWab1)
    f2ab1<-eta2ab1/(1-eta2ab1)
    lambdaab1<-f2ab1*dfWab1
    minusalpha<-1-alpha
    Ftab1<-qf(minusalpha, dfab1, dfWab1)
    powerab1<-round(1-pf(Ftab1, dfab1,dfWab1,lambdaab1),3)
    ggeab1<-round(modelab1$`Sphericity Corrections`$GGe[1],3)
    hfeab1<-round(modelab1$`Sphericity Corrections`$HFe[1],3)
    hfeab1[hfeab1>1]<-1
    ggdfab1<-ggeab1*dfab1
    ggdfWab1<-ggeab1*dfWab1
    hfdfab1<-hfeab1*dfab1
    hfdfWab1<-hfeab1*dfWab1
    lambdaggab1<-f2ab1*ggdfWab1
    lambdahfab1<-f2ab1*hfdfWab1
    Ftggab1<-qf(minusalpha, ggdfab1, ggdfWab1)
    Fthfab1<-qf(minusalpha, hfdfab1, hfdfWab1)
    powerggab1<-round(1-pf(Ftggab1, ggdfab1,ggdfWab1,lambdaggab1),3)
    powerhfab1<-round(1-pf(Fthfab1, hfdfab1,hfdfWab1,lambdahfab1),3)

    data.ab2<-subset(out, iv2==2)
    modelab2<-ez::ezANOVA(data=data.ab2, dv=.(dv), wid=.(id), within = .(iv1), type=3, detailed=TRUE)
    dfab2<-modelab2$ANOVA$DFn[2]
    dfWab2<-modelab2$ANOVA$DFd[2]
    SSab2<-modelab2$ANOVA$SSn[2]
    SSWab2<-modelab2$ANOVA$SSd[2]
    eta2ab2<-SSab2/(SSab2+SSWab2)
    f2ab2<-eta2ab2/(1-eta2ab2)
    lambdaab2<-f2ab2*dfWab2
    minusalpha<-1-alpha
    Ftab2<-qf(minusalpha, dfab2, dfWab2)
    powerab2<-round(1-pf(Ftab2, dfab2,dfWab2,lambdaab2),3)
    ggeab2<-round(modelab2$`Sphericity Corrections`$GGe[1],3)
    hfeab2<-round(modelab2$`Sphericity Corrections`$HFe[1],3)
    hfeab2[hfeab2>1]<-1
    ggdfab2<-ggeab2*dfab2
    ggdfWab2<-ggeab2*dfWab2
    hfdfab2<-hfeab2*dfab2
    hfdfWab2<-hfeab2*dfWab2
    lambdaggab2<-f2ab2*ggdfWab2
    lambdahfab2<-f2ab2*hfdfWab2
    Ftggab2<-qf(minusalpha, ggdfab2, ggdfWab2)
    Fthfab2<-qf(minusalpha, hfdfab2, hfdfWab2)
    powerggab2<-round(1-pf(Ftggab2, ggdfab2,ggdfWab2,lambdaggab2),3)
    powerhfab2<-round(1-pf(Fthfab2, hfdfab2,hfdfWab2,lambdahfab2),3)

    data.ba1<-subset(out, iv1==1)
    modelba1<-ez::ezANOVA(data=data.ba1, dv=.(dv), wid=.(id), within = .(iv2), type=3, detailed=TRUE)
    dfba1<-modelba1$ANOVA$DFn[2]
    dfWba1<-modelba1$ANOVA$DFd[2]
    SSba1<-modelba1$ANOVA$SSn[2]
    SSWba1<-modelba1$ANOVA$SSd[2]
    eta2ba1<-SSba1/(SSba1+SSWba1)
    f2ba1<-eta2ba1/(1-eta2ba1)
    lambdaba1<-f2ba1*dfWba1
    minusalpha<-1-alpha
    Ftba1<-qf(minusalpha, dfba1, dfWba1)
    powerba1<-round(1-pf(Ftba1, dfba1,dfWba1,lambdaba1),3)

    data.ba2<-subset(out, iv1==2)
    modelba2<-ez::ezANOVA(data=data.ba2, dv=.(dv), wid=.(id), within = .(iv2), type=3, detailed=TRUE)
    dfba2<-modelba2$ANOVA$DFn[2]
    dfWba2<-modelba2$ANOVA$DFd[2]
    SSba2<-modelba2$ANOVA$SSn[2]
    SSWba2<-modelba2$ANOVA$SSd[2]
    eta2ba2<-SSba2/(SSba2+SSWba2)
    f2ba2<-eta2ba2/(1-eta2ba2)
    lambdaba2<-f2ba2*dfWba2
    minusalpha<-1-alpha
    Ftba2<-qf(minusalpha, dfba2, dfWba2)
    powerba2<-round(1-pf(Ftba2, dfba2,dfWba2,lambdaba2),3)

    data.ba3<-subset(out, iv1==3)
    modelba3<-ez::ezANOVA(data=data.ba3, dv=.(dv), wid=.(id), within = .(iv2), type=3, detailed=TRUE)
    dfba3<-modelba3$ANOVA$DFn[2]
    dfWba3<-modelba3$ANOVA$DFd[2]
    SSba3<-modelba3$ANOVA$SSn[2]
    SSWba3<-modelba3$ANOVA$SSd[2]
    eta2ba3<-SSba3/(SSba3+SSWba3)
    f2ba3<-eta2ba3/(1-eta2ba3)
    lambdaba3<-f2ba3*dfWba3
    minusalpha<-1-alpha
    Ftba3<-qf(minusalpha, dfba3, dfWba3)
    powerba3<-round(1-pf(Ftba3, dfba3,dfWba3,lambdaba3),3)

    data.ba4<-subset(out, iv1==4)
    modelba4<-ez::ezANOVA(data=data.ba4, dv=.(dv), wid=.(id), within = .(iv2), type=3, detailed=TRUE)
    dfba4<-modelba4$ANOVA$DFn[2]
    dfWba4<-modelba4$ANOVA$DFd[2]
    SSba4<-modelba4$ANOVA$SSn[2]
    SSWba4<-modelba4$ANOVA$SSd[2]
    eta2ba4<-SSba4/(SSba4+SSWba4)
    f2ba4<-eta2ba4/(1-eta2ba4)
    lambdaba4<-f2ba4*dfWba4
    minusalpha<-1-alpha
    Ftba4<-qf(minusalpha, dfba4, dfWba4)
    powerba4<-round(1-pf(Ftba4, dfba4,dfWba4,lambdaba4),3)



    {print(paste("Power Factor A at B1 (Unadjusted) for n =",n,"=", powerab1))}
    {print(paste("Power Factor A at B1 H-F Adjusted (Epsilon = ",hfeab1 ,") for n =",n, "=", powerhfab1))}
    {print(paste("Power Factor A at B1 G-G Adjusted (Epsilon = ", ggeab1,") for n =",n, "=", powerggab1))}
    {print(paste("Power Factor A at B2 (Unadjusted) for n =",n,"=", powerab2))}
    {print(paste("Power Factor A at B2 H-F Adjusted (Epsilon = ",hfeab2 ,") for n =",n, "=", powerhfab2))}
    {print(paste("Power Factor A at B2 G-G Adjusted (Epsilon = ", ggeab2,") for n =",n, "=", powerggab2))}
    {print(paste("Power Factor B at A1 for n =",n, "=", powerba1))}
    {print(paste("Power Factor B at A2 for n =",n, "=", powerba2))}
    {print(paste("Power Factor B at A3 for n =",n, "=", powerba3))}
    {print(paste("Power Factor B at A4 for n =",n, "=", powerba4))}

  }}

