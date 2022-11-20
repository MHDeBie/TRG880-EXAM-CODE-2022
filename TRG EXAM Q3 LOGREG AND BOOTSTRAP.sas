*TRG880 Q3 LOGISTIC REGRESSION AND BOOTSTRAP;
quit;
proc import datafile="C:/Users/Matthew/Documents/1 UNI STUFF/TRG880/TRG_EXAM ACTUAL/QUESTION 3/q3.csv"
dbms=csv out=log_dat replace;
run;
/*proc print data=log_dat;*/
/*run;*/

data rlog_dat;
	set log_dat;
n_glucose = input(glucose, comma9.);
run;
proc print data=rlog_dat (obs=10);
run;
******NORMAL LOG REG PROCEDURE;
ods output clear;
ods output parameterestimates=est1; 
proc logistic data=rlog_dat plots(only)=Effect outest=logdat_pe;
   model TenYearCHD(Event='1') = male age education currentSmoker cigsPerDay BPMeds
						prevalentStroke prevalentHyp diabetes totChol sysBP
						diaBP BMI  heartRate n_glucose / clparm=wald ;
run;
ods output close;

proc print data=est1;
run;
data _null_;
set est1;
call execute ('data non_signit;set est1(keep=variable Estimate ProbChiSq);where ProbChiSq>=0.05;run;');
call execute ('data signit;set est1(keep=variable Estimate ProbChiSq);where ProbChiSq<0.05;run;');
call execute('data varnames;set signit (keep=Variable); where Variable^="Intercept";run;');
call execute('proc transpose data=varnames out=wide prefix=x; var Variable;run;');
call execute('data var_all;set wide;length all_pred $500.;all_pred=catx("", of x:);run;');
call execute('proc print data=var_all;run;');
run;
proc print data=non_signit;
run;
proc print data=signit;
run;

proc print data=var_all;
run;

ods output clear;
ods output parameterestimates=est2;
data _null_;
set var_all;
/*where all_pred =: 'x';*/
call execute ('proc logistic data=rlog_dat plots(only)=Effect; model TenYearCHD(Event="1") = '|| all_pred||'; run;quit;');
run;
ods output close;

data clean_logdat;
 set rlog_dat;
 if cmiss(of _all_) then delete;
run;
proc iml;
*NEWTON RAPHSON;
use clean_logdat;
read all into xy;
kk=ncol(xy);
n=nrow(xy);
y=xy[,(kk-1)];
x=xy[,1:(kk-2)]||xy[,(kk)];
/*print (xy[1:5,]);*/
/*print (x[1:5,]) (y[1:5,]);*/
x_org=x;
y_org=y;

xm=J(n,1,1)||x_org;
bhi = j(ncol(xm),1,0);
bho = bhi;
diff = 11111;
i=1;
x=x_org;
namr={"INT" "male" "age" "education" "currentSmoker" "cigsPerDay" "BPMeds"
						"prevalentStroke" "prevalentHyp" "diabetes" "totChol" "sysBP"
						"diaBP" "BMI"  "heartRate" "n_glucose"};
do i=1 to 15 while (diff>0.00001);
	y=y_org;
	x=xm;
	lo = x*bho; 
	o = exp(lo); 
	p = o/(1+o);
	w = diag(p#(1-p));
	bhn = bho + inv(x`*w*x)*x`*(y-p); /*Newton-Raphson algorithm*/
	diff = max(abs(bhn-bho));
	bho=bhn;
end;
print bhn[r=namr];

yhat=xm*bhn; *vector of predicted response;
call sort(yhat);
print (yhat[1:10,]);

print diff;
/*estimate the variance*/
cov=inv(x`*w*x);
variance=vecdiag(cov);
print cov, variance;
bh_print=bhn`;


quit;


*PLOT THE ROC CURVE;
proc logistic data=rlog_dat plots(only)=ROC ;
   ods select ROCcurve;
   model TenYearCHD(Event='1') = male age education currentSmoker cigsPerDay BPMeds
						prevalentStroke prevalentHyp diabetes totChol sysBP
						diaBP BMI  heartRate n_glucose/ outroc=roc  ;
  output out=b p=pred;
run;

*LOG REG ON ONLY SIGNIT;
/*proc logistic data=rlog_dat plots(only)=ROC ;*/
/*   ods select ROCcurve;*/
/*   model TenYearCHD(Event='1') = male age cigsPerDay totChol sysBP n_glucose/ outroc=roc  ;*/
/*  output out=b p=pred;*/
/*run;*/

proc print data=roc;
   var _1mSpec_ _Sensit_;
run;


proc gplot data=roc;
title "ROC for Logistic Regression Model";
 plot ( _1mSpec_ _Sensit_)*_1mSpec_ ;
run ;

proc iml;
use roc;
read all var {"_1MSPEC_"} into x;
read all var {"_SENSIT_"} into y;



xy=x||y;
nm={"x" "y"};
maxy= round((xy[,2])[<>]+0.1);
print maxy;

a=0;
b=1;
npoint=nrow(xy);
ds=J(npoint,1,0);
*these are coordinate initialisers;
*xc are the x_coords of the points;
*yc are the ycoordinates;
xc=ranuni(ds)*(b-a)+a;
yc=ranuni(ds)*maxy;

call sort(xc);
yy=y;
/*print yy yc xc;*/
cprop=(yc<yy)[+]/npoint;
grp=J(nrow(yy),1,.);
grp_calc=J(nrow(yy),1,.);
do ii=1 to nrow(yy);
	if yy[ii]>yc[ii] then grp[ii]=1;
	if yy[ii]<yc[ii] then grp[ii]=2;

	if yy[ii]>yc[ii] then grp_calc[ii]=1;
	if yy[ii]<yc[ii] then grp_calc[ii]=0;
end;
mc_area=grp_calc[+,]/npoint;
print mc_area;
auc=mc_area;
*the idea is to saturate the area with so many points;
area=cprop*((b-a)*maxy);



npoints=nrow(yc);
gini=(2*area)-1;
print npoints gini area;
datall=(xy||yc||xc||grp);
nm1={"x" "y" "yc" "xc" "gr"};

create dall from datall[colname=nm1];
append from datall;
close dall;

create areadat from area[colname="area"];
append from area;

quit;


symbol1 interpol=none width=1
color=blue
value=dot
height=2;

symbol2 interpol=none width=1
color=red
value=dot
height=2;


proc sgplot data=dall;
title "Monte Carlo Integration";
series y=y x=xc/ lineattrs=(color=green thickness=4);
scatter y=yc x=xc/ group=gr;
run;

proc gplot data=dall;
plot y*x=gr;
run;



********************************
*COMPUTE THE ROUGH ROC CURVE HERE;
proc means data=b min p10 p20  p30 p40 p50 p60 p70 p80 max ;
 var pred ;
 output out=cc min= p10= p20=  p30= p40= p50= p60= p70= p80= max= /autoname;
run ;
proc print data=cc;
run;

data dd (drop = _TYPE_ _FREQ_);
 set b ;
 if _N_ =1 then set cc ; 
 if pred <= pred_max then pred_c = "p80 - max" ;
 if pred <= pred_p80 then pred_c = "p60 - p80" ;
 if pred <= pred_p60 then pred_c = "p40 - p60" ;
 if pred <= pred_p40 then pred_c = "p20 - p40" ;
 if pred <= pred_p20 then pred_c = "min - p20" ;
run ;


proc freq data=dd;
 tables pred_c*TenYearCHD ;
run ;

proc freq data=dd;
 tables pred_c*TenYearCHD / out=ee1;
 where TenYearCHD=0 ;
run ;
data ee1 (keep = pred_c freq0 perc0 cumper0);
 set ee1 ;
 freq0 = count ;
 perc0 = percent ;
 retain cumper0;
 cumper0+perc0;
run ;
proc freq data=dd;
 tables pred_c*TenYearCHD / out=ee2;
 where TenYearCHD=1 ;
run ;
data ee2 (keep = pred_c freq1 perc1 cumper1);
 set ee2 ;
 freq1 = count ;
 perc1 = percent ;
 retain cumper1;
 cumper1+perc1;
run ;

data ff ;
 pred_c = "min      " ;
 cumper0 = 0 ;
 cumper1 = 0 ;
run ;

data gg ;
 merge ff ee1 ee2 ;
 by pred_c ;
run ;


data gg ; 
 set gg ;
 ks = abs(cumper1-cumper0);
 ksr = ks/100 ;
run ;

proc print data=gg ;
run ;
proc means data=gg max;
 var ks ksr;
run ;

data b ; 
set b ;
chd_pTenYearCHD_p = 0 ;
if pred > 0.5 then chd_pTenYearCHD_p= 1 ;
run ;

proc freq data=b ;
 tables TenYearCHD*TenYearCHD_p ;
run ;

symbol1 interpol=join width=4 color=blue value=dot  height=1;
symbol2 interpol=join width=3 color=red  value=dot  height=1;
data nm_gg;
	set gg;
TPR=cumper0;
FPR=cumper1;
run;

proc gplot data=nm_gg;
title "ROC for Logistic Regression Model";
 plot ( TPR FPR)*FPR / overlay href=(17.857 25.466 39.596
       62.422) chref=purple ;
run ;

proc gplot data=gg;
 plot (cumper1 cumper0 )*cumper1 / overlay href=(17.857 25.466 39.596
       62.422) chref=purple areas=3;
run ;
*END OF COMPUTATION OF ROUGH ROC CURVE;
**************************************;



*************************************
BOOTSTRAPPING GINI;
proc iml;
use roc;
read all var {"_1MSPEC_"} into x;
read all var {"_SENSIT_"} into y;

x_org=x;
y_org=y;
xy=x||y;
xy_org=xy;
nm={"x" "y"};
n=nrow(xy); ss_n=n;
ss_m=J(ss_n,1,0);
do i=1 to ss_n;
	ss_m[i,]=i;
end;
/*print ss_m;*/
sam_size=ss_n;
NSAM=1000; *Number of Bootstraps;
AUC_res=J(NSAM,1,.);
GIN_res=J(NSAM,1,.);

do jj=1 to NSAM;

	x=x_org;
	y=y_org;
	ss=sample((ss_m),sam_size, "Rep");
	ss=ss`;
	call sort(ss);
	xy_sam=xy[ss,];
	x=xy_sam[,1];
	y=xy_sam[,2];


	maxy= round((xy[,2])[<>]+0.1);
/*	print maxy;*/

	a=0;
	b=1;
	npoint=nrow(xy_sam);
	ds=J(npoint,1,0);
	xc=ranuni(ds)*(b-a)+a;
	yc=ranuni(ds)*maxy;

	call sort(xc);
	yy=y;
	/*print yy yc xc;*/
	cprop=(yc<yy)[+]/npoint;
	grp=J(nrow(yy),1,.);
	grp_calc=J(nrow(yy),1,.);
	do ii=1 to nrow(yy);
		if yy[ii]>yc[ii] then grp[ii]=1;
		if yy[ii]<yc[ii] then grp[ii]=2;
	end;
	mc_area=cprop;
	auc=mc_area;
	*the idea is to saturate the area with so many points;
	area=cprop*((b-a)*maxy);
	npoints=nrow(yc);
	gini=(2*area)-1;
/*	print jj npoints gini area;*/

AUC_res[jj]=AUC;
GIN_res[jj]=gini;
end;

print AUC_res GIN_res;

*AUC CONFIDENCE;
call sort(AUC_res);
	a=1-0.05;
	p_vec={0.025,0.975};
    call qntl(CI_AUC,AUC_res,p_vec);
	print "Confidence Interval:" a ;
	print p_vec CI_AUC;

*GINI CONFIDENCE;
call sort(GIN_res);
	a=1-0.05;
	p_vec={0.025,0.975};
    call qntl(CI_GIN,GIN_res,p_vec);
	print "Confidence Interval:" a ;
	print p_vec CI_GIN;


create dall from datall[colname=nm1];
append from datall;
close dall;

create areadat from area[colname="area"];
append from area;

quit;

*************************************
BOOTSTRAPPING REGRESSION;
*NEWTON RAPHSON;
proc iml;
use clean_logdat;
read all into xy;
kk=ncol(xy);
n=nrow(xy);
y=xy[,(kk-1)];
x=xy[,1:(kk-2)]||xy[,(kk)];
x_org=x;
y_org=y;
n=nrow(xy); ss_n=n;
ss_m=J(ss_n,1,0);
do i=1 to ss_n;
	ss_m[i,]=i;
end;
sam_size=ss_n;
NSAM=30;
yhat_res=J(n,NSAM,.);
beta_res=J(ncol(xy),NSAM,.);
print ss_n;

print (x_org[1:5,]);

do jj=1 to NSAM;
/*print JJ;*/
	x=x_org;
	y=y_org;
	ss=sample((ss_m),sam_size, "Rep");
	ss=ss`;
	x_sam=x[ss,];
	y_sam=y[ss,];
	
	do qq=1 to ncol(x_sam);
		xcol=x_sam[,qq];
		if xcol[+,]=0 then x_sam[,qq]=x_org[,(qq)];
		else x_sam[,qq]=xcol;
	end;
	x=x_sam;
	y=y_sam;
		n=nrow(x);
	xm=J(n,1,1)||x;
/*	det_x=det(xm);*/
/*	if det_x^=0 then do;*/
	bhi = j(ncol(xm),1,0);
	bho = bhi;
	diff = 11111;
	i=1;
/*	print (xm[1:3,]);*/

	do i=1 to 15 while (diff>0.00001);
		lo = xm*bho; 
		o = exp(lo); 
		p = o/(1+o);
		w = diag(p#(1-p));
		bhn = bho + inv(xm`*w*xm)*xm`*(y-p); /*Newton-Raphson algorithm*/
		diff = max(abs(bhn-bho));
		bho=bhn;
	end;
	yhat=xm*bhn; *vector of predicted response;
	prob=1/(1+exp(-yhat));
	call sort(prob);
	yhat_res[,jj]=prob;
	beta_res[,jj]=bhn;
end;

*95% CI FOR REGRESSION COEFFICIENTS;
a=1-0.05;
p_vec={0.025,0.975};
beta_res95CI=J(nrow(beta_res),2,.);
do ii=1 to nrow(beta_res);
	beta_int=beta_res[ii,]; *ALL YHAT ii's;
	beta_T=beta_int`;
	call sort(beta_T);
	    call qntl(CI_beta,beta_T,p_vec);
	CI_betaT=CI_beta`;
	beta_res95CI[ii,]=CI_betaT;
end;
namr={"INT" "male" "age" "education" "currentSmoker" "cigsPerDay" "BPMeds"
						"prevalentStroke" "prevalentHyp" "diabetes" "totChol" "sysBP"
						"diaBP" "BMI"  "heartRate" "n_glucose"};
cpvec=char(p_vec);
print beta_res95CI[r=namr c=cpvec l="Regression Coefficients 95% Confidence Interval"];





print (yhat_res[1:5,1:5]);

a=1-0.05;
p_vec={0.025,0.975};

CI95_yhat=J(nrow(yhat_res),2,.);
AVG_yhat=J(nrow(yhat_res),1,.);
do ii=1 to nrow(yhat_res);
	yhat_int=yhat_res[ii,]; *ALL YHAT ii's;
	yhat_T=yhat_int`;
	ayhat=mean(yhat_T);
	AVG_yhat[ii]=ayhat;
	call sort(yhat_T);
	    call qntl(CI_yhat,yhat_T,p_vec);
	CI_yhatT=CI_yhat`;
	CI95_yhat[ii,]=CI_yhatT;
end;

/*print CI95_yhat;*/
nmcy={"LCI" "UCI" "AV_PROB" "PVEC"};
call sort(x_org[,2]);
namr={"male" "age" "education" "currentSmoker" "cigsPerDay" "BPMeds"
						"prevalentStroke" "prevalentHyp" "diabetes" "totChol" "sysBP"
						"diaBP" "BMI"  "heartRate" "n_glucose"};
nmy={"CHD"};
nmcy=nmcy||namr||nmy;
jump=1/nrow(CI95_yhat);
a=0 ;b=1;
prob_vec=J(nrow(CI95_yhat),1,.);

prob_vec=do((a+jump),(b),jump)`;
ntest1=nrow(CI95_yhat);
ntest2=nrow(prob_vec);
print ntest1 ntest2;

outdat=CI95_yhat||AVG_yhat||prob_vec||x_org||y_org;
/*print (outdat[1:5,]);*/


create dout from outdat[colname=nmcy];
append from outdat;
close dout;
quit;

proc gplot data=dout ;
plot LCI*PVEC;
run ;

proc sgplot data=dout;
title "Regression Confidence Interval";
/*scatter y=CHD x=age;*/
series y=LCI x=PVEC;
series y=UCI x=PVEC;
series y=AV_PROB x=PVEC;
run;

*************************************
BOOTSTRAPPING REGRESSION STANDARD ERRORS;
*NEWTON RAPHSON;
proc iml;
use clean_logdat;
read all into xy;
kk=ncol(xy);
n=nrow(xy);
y=xy[,(kk-1)];
x=xy[,1:(kk-2)]||xy[,(kk)];
x_org=x;
y_org=y;
n=nrow(xy); ss_n=n;
ss_m=J(ss_n,1,0);
do i=1 to ss_n;
	ss_m[i,]=i;
end;
sam_size=ss_n;
NSAM=50;
print ss_n;

print (x_org[1:5,]);
namr={"INT" "male" "age" "education" "currentSmoker" "cigsPerDay" "BPMeds"
						"prevalentStroke" "prevalentHyp" "diabetes" "totChol" "sysBP"
						"diaBP" "BMI"  "heartRate" "n_glucose"};

Res_std=J((kk),nsam,.);
do jj=1 to NSAM;
/*print JJ;*/
	x=x_org;
	y=y_org;
	ss=sample((ss_m),sam_size, "Rep");
	ss=ss`;
	x_sam=x[ss,];
	y_sam=y[ss,];
	
	do qq=1 to ncol(x_sam);
		xcol=x_sam[,qq];
		if xcol[+,]=0 then x_sam[,qq]=x_org[,(qq)];
		else x_sam[,qq]=xcol;
	end;
	x=x_sam;
	y=y_sam;
		n=nrow(x);
	xm=J(n,1,1)||x;
/*	det_x=det(xm);*/
/*	if det_x^=0 then do;*/
	bhi = j(ncol(xm),1,0);
	bho = bhi;
	diff = 11111;
	i=1;
/*	print (xm[1:3,]);*/

	do i=1 to 15 while (diff>0.00001);
		lo = xm*bho; 
		o = exp(lo); 
		p = o/(1+o);
		w = diag(p#(1-p));
		bhn = bho + inv(xm`*w*xm)*xm`*(y-p); /*Newton-Raphson algorithm*/
		diff = max(abs(bhn-bho));
		bho=bhn;
	end;
	yhat=xm*bhn; *vector of predicted response;
	prob=1/(1+exp(-yhat));
	cov=inv(xm`*w*xm);
	variance=vecdiag(cov);
	std_err=sqrt(variance);
	Res_std[,jj]=std_err;
end;

print (Res_std[1:kk,1:5]);

a=1-0.05;
p_vec={0.025,0.975};

std95=J(nrow(Res_std),2,.);
stdAV=J(nrow(Res_std),1,.);
do ii=1 to nrow(Res_std);
	std_int=Res_std[ii,]; *ALL YHAT ii's;
	std_T=std_int`;
	a_std=mean(std_T);
	call sort(std_T);
	    call qntl(CI_std,std_T,p_vec);
	CI_stdT=CI_std`;
	std95[ii,]=CI_stdT;
	stdAV[ii,]=a_std;
end;

std_all=std95||stdAV;
cn={"LCI" "UCI" "Average"};
print std_all[r=namr c=cn];

quit;



