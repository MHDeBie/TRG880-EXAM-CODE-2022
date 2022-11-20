*******************************
TRG 880 EXAM u18119543
Gaussian Mixture Regression;

quit;
proc import datafile="C:/Users/Matthew/Documents/1 UNI STUFF/TRG880/TRG_EXAM ACTUAL/QUESTION 2/q2.csv"
dbms=csv out=gmr_dat replace;
run;

proc print data=gmr_dat (obs=10);
run;

proc contents data=gmr_dat noprint out=varnames (keep=name);
run;

proc print data=varnames;
title "VARNAMES";
run;
******
GENERATE SCATTERPLOTS FOR ALL VARIABLES;

data _null_;
set varnames;
where name =: 'x';
call execute ('proc sgplot data=gmr_dat; title "Scatter Plot of Y and Predictor:'||name||'";
				scatter y=y x='||name||'; run;');
run;

******
INITIAL FEATURE SELECT INVESTIGATIONS;

data _null_;
call execute('data pred;set gmr_dat (KEEP=x:);run;');
call execute('proc contents data=pred noprint out=varnames (keep=name);run;');
call execute('proc transpose data=varnames out=wide prefix=x;var NAME; run;');
call execute('data var_all;set wide;length all_pred $500.; all_pred=catx("", of x:);run;');
run;

ods output clear;
ods output parameterestimates=est1; 
title "Initial Feature Selection Investigations";
data _null_;
set var_all;
where all_pred =: 'x';
call execute ('proc reg data=gmr_dat; model y= '|| all_pred||'/vif tol collin; run;quit;');
call execute ('proc print data=est1;run;');
call execute ('data non_signit;set est1(keep=variable ProbChiSq);where ProbChiSq>=0.05;run;');
call execute ('data signit;set est1(keep=variable ProbChiSq);where ProbChiSq<0.05;run;');
call execute('data varnames;set signit (keep=Variable);run;');
call execute('proc transpose data=varnames out=wide prefix=x; var Variable;run;');
call execute('data var_all;set wide;length all_pred $500.;all_pred=catx("", of x:);run;');
call execute('proc print data=var_all;run;');
run;
ods output close;


********
GAUSSIAN MIXTURE REGRESSION MODEL ESTIMATION;

/*proc fastclus data=gmr_dat out=Clust maxclusters=2 maxiter=3 random=321;*/
/*   var x;*/
/*run;*/
data _null_;
set var_all;
where all_pred =: 'x';
call execute ('proc fastclus data=gmr_dat out=Clust maxclusters=2 maxiter=3 random=321; var '|| all_pred||'; run;');
run;
proc print data=Clust (obs=5);
run;

proc sgplot data=Clust; 
scatter x=x1 y=y/group=CLUSTER;
run;


quit;
proc iml;
use gmr_dat;
read all into xy;
pp=ncol(xy);
print pp;
x=xy[,1:(pp-1)];
y=xy[,pp];
n=nrow(x);
og_x=x;
og_y=y;
u=og_x;
*****************NUMBER OF COMPONENTS;
k=2;
****************;

use Clust;
read all into clusdat;
group=clusdat[,(pp+1)];
*QUICK N DIRTY ASSIGNMENT OF INITIAL CLUSTERS;

count=J(n,1,.);
do ii=1 to n;
	count[ii]=ii;
end;

start LinReg(x,y);
	bh=inv(x`*x)*x`*y;
	yh=x*bh;
    return(bh);
finish;

start MVNormalPDF(y, mu, sigma);
	n=nrow(y);
	p = nrow(sigma);
	phi=j(n,1,.);
	do i=1 to n;
	   yi=y[i,];
	   mui=mu[i,];
	   const = (2*constant("PI"))##(p/2) * sqrt(det(sigma));
	   d = Mahalanobis(yi, mui, sigma);
	   phi[i,] = exp( -0.5*d#d ) / const;
	end;
	 return( phi );
finish;

x=og_x;
ni=j(1,k,.);
b_old=j(pp, k, .);
s2=j(1,k,.);
xvec=J(n,1,1)||og_x;
yhat=J(n,k,.);
do jj=1 to k;

   xvec=J(n,1,1)||x; *define the design matrix;
   id=loc(group=jj);
   xi=xvec[id,];
   yi=y[id,];
   s2[,jj]=var(yi);
   ni[,jj]=nrow(xi);
   BETA=LinReg(xi, yi);
   b_old[,jj]=LinReg(xi, yi);
   yhat[,jj]=xvec*BETA;
end;

mp=ni/n;
var=sum(s2);
p=ncol(x);



var=sum(var);
init_beta=b_old;
init_mp=mp;
init_var=var;
init_group=group;

idat=yhat||og_y||og_x||group;
/*print (idat[1:10,]);*/

*INITIAL RESULTS;
*************************;
xnam="X1":"X20";
bnam="B0":"B20";
compnam="K1":"K2";
nm= {"yhat1" "yhat2" "y"};
gm={"group"};
nmx=nm||xnam||gm;
create initdat from idat[colname=nmx];
append from idat;
close initdat;
***********************;



bet_old=init_beta;
mp_old=init_mp;
var_old=init_var;
print init_beta[r=bnam c={"BK1" "BK2"} l="Initial:"] init_mp[c={"Pi1" "Pi2"}] init_var[c="Var"] ;

x=xvec;
y=og_y;
p=ncol(x);
gamma=J(n,k,.);
maxiter=100;
conv=99999999;
Log_pdf=j(maxiter,1,.);

*****LOOP GOES HERE;
do qq=1 to maxiter while(conv>0.0000001);
/*print "Iter" qq;*/

	do tt=1 to k;
		   mu_k=x*bet_old[,tt];
		   gamma_k=mp_old[tt] # MVNormalPDF(y, mu_k, var_old);
		   gamma[,tt]=gamma_k;
	end;
	gamma=gamma/gamma[,+];
	group=gamma[,<:>];
	mp_new=(gamma[+,] / n);*mixprob;
	bet_new=j(p, k, .);


	var_new=J(n, k, .);
	yhat=J(n,k,.);
	MMG=j(n,k,.);
	do jj=1 to k;
		   wk=diag(gamma[,jj]); *BETA;
		   beta=inv(x`*wk*x)*x`*wk*y;
		   bet_new[,jj]=beta;*b_old*;
		   sigma2=gamma[,jj]#((y-x*beta)##2); *VARIANCE;
		   var_new[,jj]=sigma2;
		   yhat_k=x*beta;
		   yhat[,jj]=yhat_k;
	end;
		var_new=var_new[+,+]/n;*s_old*;

	do jj=1 to k;
		mu_knew=x*bet_new[,jj];
		MMG[,jj]=mp_new[jj]#MVNormalPdf(y, mu_knew, var_new);
	end;
		sum_MMG=MMG[,+];
		Log_pdf[qq]=(log(sum_MMG))[+,];

	diff1=abs(bet_old-bet_new);
	diff2=abs(var_old-var_new);
	diff3=abs(mp_old-mp_new);
	diff=abs(diff1+diff2+diff3);
	conv=diff;

	mp_old=mp_new;
	bet_old=bet_new;
	var_old=var_new;
/*	print bet_old mp_old var_old;*/
	group_end=group;
end;
print "Number of Iter:" qq;
fin_log_pdf=log_pdf[(qq-1)];
num_parm=k*3;
AIC=2*num_parm - 2*fin_log_pdf;
BIC=num_parm*log(n) -2*fin_log_pdf;

bet_final=bet_old;
mp_final=mp_old;
var_final=var_old;

print bet_final[r=bnam c={"BK1" "BK2"} l="Final:"]"" mp_final[c={"Pi1" "Pi2"} l=""]""  var_final[c="Var" l=""] ;
act_y=og_y;
rss_pp=J(1,k,.);
******************************;
*CALCULATE COMPONENT WISE RSS;
******************************;
/*print group_end;*/
do hh=1 to ncol(yhat);
/*			if hh=2 then hh=3;*/
			id_yh=loc(group_end=hh)`;
			act_y=og_y[id_yh,];
			cmp_y=yhat[id_yh,hh];
		
		rss_k=(act_y-cmp_y)`*(act_y-cmp_y);
		rss_pp[,hh]=rss_k;
end;
print AIC BIC rss_pp  qq;


fdat=yhat||og_y||og_x||group;



nm= {"yhat1" "yhat2" "y" "x" "group" };
/*print (fdat[1:10,])[c=nmx];*/
create gmrdat from fdat[colname=nmx];
append from fdat;
close gmrdat;


quit;

*WE JUST USE THIS TO CHECK AND MAKE SURE ALL IS WORKING AS INTENDED;
proc sgplot data=initdat;
title "INITIAL GROUP ASSIGNMENT";
	scatter x=x1 y=y/ group=group;
run;

proc contents data=gmrdat noprint out=varnames (keep=name);
run;

data _null_;
set varnames;
where name =: 'X';
call execute ('proc sgplot data=gmrdat; title "Mixture Regression Results of Y and:'||name||'";
				scatter x='||name||' y=y/ group=group markeroutlineattrs=(thickness=2)
   								  markerattrs=(symbol=circlefilled size=10 );
				scatter x='||name||' y=yhat1/  filledoutlinedmarkers 
   markerfillattrs=(color=green) 
   markeroutlineattrs=(color=green thickness=2)
   markerattrs=(symbol=circlefilled size=3);
				scatter x='||name||' y=yhat2/ filledoutlinedmarkers 
   markerfillattrs=(color=purple) 
   markeroutlineattrs=(color=purple thickness=2)
   markerattrs=(symbol=circlefilled size=3);;
				run;');
run;

*******************
INTRODUCING RIDGE REGRESSION;

proc contents data=gmr_dat noprint out=varnames (keep=name);
run;

proc print data=varnames;
title "VARNAMES";
run;
data _null_;
set var_all;
where all_pred =: 'x';
call execute ('proc fastclus data=gmr_dat out=Clust maxclusters=2 maxiter=3 random=321; var '|| all_pred||'; run;');
run;



quit;
proc iml;
use gmr_dat;
read all into xy;
pp=ncol(xy);
print pp;
x=xy[,1:(pp-1)];
y=xy[,pp];
n=nrow(x);
ogg_x=x;
ogg_y=y;
u=og_x;
*****************NUMBER OF COMPONENTS;
k=2;
****************;

use Clust;
read all into clusdat;
group=clusdat[,(pp+1)];
*QUICK N DIRTY ASSIGNMENT OF INITIAL CLUSTERS;

count=J(n,1,.);
do ii=1 to n;
	count[ii]=ii;
end;

start stand(x); ****FOR STANDARDISING;
	xBar = mean(X);           
	X = X - xBar;                           
	std=std(X);
	X_std = X / std;  
	return (X_std);
finish;

x=stand(x);
y=stand(y);
og_x=x;
og_y=y;



start MVNormalPDF(y, mu, sigma);
	n=nrow(y);
	p = nrow(sigma);
	phi=j(n,1,.);
	do i=1 to n;
	   yi=y[i,];
	   mui=mu[i,];
	   const = (2*constant("PI"))##(p/2) * sqrt(det(sigma));
	   d = Mahalanobis(yi, mui, sigma);
	   phi[i,] = exp( -0.5*d#d ) / const;
	end;
	 return( phi );
finish;


start RidgeReg(x,y,lambda);
	kn=ncol(x);
	Iden=I(kn);
	bh=inv(x`*x+lambda*Iden)*x`*y;
	yh=x*bh;
    return(bh);
finish;


*************************;
*DEFINE LAMBDA;
lambda_all=(do(0.0,10,0.5))`;
/*lambda_all={0.3,0.5,0.75, 0.8};*/
n_lam=nrow(lambda_all);
Iden=I(ncol(x));
print n_lam;
mm=ncol(og_x)+1;
rss_res=J(n_lam,k,.);
beta_res1=J(mm,n_lam,.); flg_res1=J(mm,n_lam,.);
beta_res2=J(mm,n_lam,.); flg_res2=J(mm,n_lam,.);
grp_res=J(n,n_lam,.);
lam_res=J(n_lam,1,.);
itr_res=J(n_lam,1,.);
stp_res=J(n_lam,1,.);

do lam_it=1 to n_lam;
lambda=lambda_all[lam_it];
/*lambda=0.5;*/
print "Possible Lambda Value:" lam_it lambda;
lam_res[lam_it]=lambda;
stp_res[lam_it]=lam_it;
*************************;

x=og_x;
ni=j(1,k,.);
b_old=j(pp, k, .);
s2=j(1,k,.);
xvec=J(n,1,1)||og_x;
yhat=J(n,k,.);


do jj=1 to k;

   xvec=J(n,1,1)||x; *define the design matrix;
   id=loc(group=jj);
   xi=xvec[id,];
   yi=y[id,];
   s2[,jj]=var(yi);
   ni[,jj]=nrow(xi);
   BETA=RidgeReg(xi, yi,lambda);
   b_old[,jj]=RidgeReg(xi, yi,lambda);
   yhat[,jj]=xvec*BETA;
end;

mp=ni/n;
var=sum(s2);
p=ncol(x);



var=sum(var);
init_beta=b_old;
init_mp=mp;
init_var=var;
init_group=group;

idat=yhat||og_y||og_x||group;


*INITIAL RESULTS;
*************************;
xnam="X1":"X20";
bnam="B0":"B20";
compnam="K1":"K2";
nm= {"yhat1" "yhat2" "y"};
gm={"group"};
nmx=nm||xnam||gm;
/*print (idat[1:10,])[c=nmx];*/
/*create initdat from idat[colname=nmx];*/
/*append from idat;*/
/*close initdat;*/
***********************;



bet_old=init_beta;
mp_old=init_mp;
var_old=init_var;
/*print init_beta[r=bnam c={"BK1" "BK2"} l="Initial:"] init_mp[c={"Pi1" "Pi2"}] init_var[c="Var"] ;*/

x=xvec;
y=og_y;
p=ncol(x);
gamma=J(n,k,.);
maxiter=100;
conv=99999999;
Log_pdf=j(maxiter,1,.);


do qq=1 to maxiter while(conv>0.0000001);


	do tt=1 to k;
		   mu_k=x*bet_old[,tt];
		   gamma_k=mp_old[tt] # MVNormalPDF(y, mu_k, var_old);
		   gamma[,tt]=gamma_k;
	end;
	gamma=gamma/gamma[,+];
	group=gamma[,<:>];
	mp_new=(gamma[+,] / n);*mixprob;
	bet_new=j(p, k, .);


	var_new=J(n, k, .);
	yhat=J(n,k,.);
	MMG=j(n,k,.);
	do jj=1 to k;
		   wk=diag(gamma[,jj]); *BETA;
		   xvec=J(n,1,1)||og_x;
		   x=xvec;
		   Iden=ncol(xvec);
		   beta=inv(x`*wk*x+lambda*Iden)*x`*wk*y; *NOW WE HAVE TO INTRODUCE THE RIDGE PARM INTO THE CALCULATION;
		   bet_new[,jj]=beta;*b_old*;
		   sigma2=gamma[,jj]#((y-x*beta)##2); *VARIANCE;
		   var_new[,jj]=sigma2;
		   yhat_k=x*beta;
		   yhat[,jj]=yhat_k;
	end;
		var_new=var_new[+,+]/n;*s_old*;

	do jj=1 to k;
		mu_knew=x*bet_new[,jj];
		MMG[,jj]=mp_new[jj]#MVNormalPdf(y, mu_knew, var_new);
	end;
		sum_MMG=MMG[,+];
		Log_pdf[qq]=(log(sum_MMG))[+,];

	diff1=abs(bet_old-bet_new);
	diff2=abs(var_old-var_new);
	diff3=abs(mp_old-mp_new);
	diff=abs(diff1+diff2+diff3);
	conv=diff;

	mp_old=mp_new;
	bet_old=bet_new;
	var_old=var_new;
/*	print bet_old mp_old var_old;*/
	group_end=group;
end;
/*print "Number of Iter:" qq;*/
fin_log_pdf=log_pdf[(qq-1)];
num_parm=k*3;
AIC=2*num_parm - 2*fin_log_pdf;
BIC=num_parm*log(n) -2*fin_log_pdf;

bet_final=bet_old;
mp_final=mp_old;
var_final=var_old;

*****
BASIC THRESH HOLD ANALYSIS;
thr=0.05;
flg_k1=(abs(bet_final[,1])>thr)#bet_final[,1];
flg_k2=(abs(bet_final[,2])>thr)#bet_final[,2];
rss_pp=J(1,k,.);
******************************;
*CALCULATE COMPONENT WISE RSS;
******************************;
do hh=1 to ncol(yhat);
			id_yh=loc(group_end=hh)`;
			act_y=og_y[id_yh,];
			cmp_y=yhat[id_yh,hh];
		
		rss_k=(act_y-cmp_y)`*(act_y-cmp_y)+lambda*(bet_final[+,hh])**2;
		rss_pp[,hh]=rss_k;
end;
rss_res[lam_it,]=rss_pp;
beta_res1[,lam_it]=bet_final[,1]; flg_res1[,lam_it]=flg_k1;
beta_res2[,lam_it]=bet_final[,2]; flg_res2[,lam_it]=flg_k2;
grp_res[,lam_it]=group_end;
itr_res[lam_it,]=qq;
/*lam_it*/
************************;
end;*END OF BIG LOOP;
************************;
lam_nas=char(lam_res);
total_rss=rss_res[,+];
print rss_res[c=compnam L="RSS"] total_rss itr_res lam_res[l="Lambda"];
print flg_res1[l="Parameters Not at Risk of Elimination" r=bnam c=lam_nas];
print flg_res2[l="Parameters Not at Risk of Elimination" r=bnam c=lam_nas];
print beta_res1[r=bnam];

parm_dat=lam_res||(beta_res1)`;
rnam="Lambda"||bnam;
fdat=yhat||og_y||og_x||group;



nm= {"yhat1" "yhat2" "y" "x" "group" };
/*print (fdat[1:10,])[c=nmx];*/


create parmdat from parm_dat[c=rnam];
append from parm_dat;
close parmdat;


quit;
proc sgplot data=parmdat;
title "Parameter Trajectories component k=1";
	series x=Lambda y=B1;
	series x=Lambda y=B2;
	series x=Lambda y=B3;
	series x=Lambda y=B4;
	series x=Lambda y=B5;
	series x=Lambda y=B6;
	series x=Lambda y=B7;
	series x=Lambda y=B8;
	series x=Lambda y=B9;
	series x=Lambda y=B10;
	series x=Lambda y=B11;
	series x=Lambda y=B12;
	series x=Lambda y=B13;
	series x=Lambda y=B14;
	series x=Lambda y=B15;
	series x=Lambda y=B16;
	series x=Lambda y=B17;
	series x=Lambda y=B18;
	series x=Lambda y=B19;
	series x=Lambda y=B20;
run;
proc print data=initdat (obs=10);
run;

proc sgplot data=initdat;
title "INITIAL GROUP ASSIGNMENT";
	scatter x=x1 y=y/ group=group;
	scatter x=x1 y=yhat1;
	scatter x=x1 y=yhat2;
run;

*******************ORIGINAL FEATURE SET;
*************************RIDGE REGRESSION WITH ONLY NON-ZERO FEATURES;

quit;
proc fastclus data=gmr_dat out=Clust maxclusters=2 maxiter=3 random=321; var x1 x2 x3 x4 x8; run;
proc iml;
use gmr_dat;
read all into xy;
pp=ncol(xy);
print pp;
read all var {"x1" "x2" "x3" "x4" "x8"} into x;
y=xy[,pp];
n=nrow(x);
ogg_x=x;
ogg_y=y;

*****************NUMBER OF COMPONENTS;
k=2;
****************;

use Clust;
read all into clusdat;
group=clusdat[,(pp+1)];
print (group[1:10,]);
*QUICK N DIRTY ASSIGNMENT OF INITIAL CLUSTERS;

count=J(n,1,.);
do ii=1 to n;
	count[ii]=ii;
end;

start stand(x); ****FOR STANDARDISING;
	xBar = mean(X);           
	X = X - xBar;                           
	std=std(X);
	X_std = X / std;  
	return (X_std);
finish;

x=stand(x);
y=stand(y);
og_x=x;
og_y=y;
pp=ncol(x)+1;


start MVNormalPDF(y, mu, sigma);
	n=nrow(y);
	p = nrow(sigma);
	phi=j(n,1,.);
	do i=1 to n;
	   yi=y[i,];
	   mui=mu[i,];
	   const = (2*constant("PI"))##(p/2) * sqrt(det(sigma));
	   d = Mahalanobis(yi, mui, sigma);
	   phi[i,] = exp( -0.5*d#d ) / const;
	end;
	 return( phi );
finish;


start RidgeReg(x,y,lambda);
	kn=ncol(x);
	Iden=I(kn);
	bh=inv(x`*x+lambda*Iden)*x`*y;
	yh=x*bh;
    return(bh);
finish;


*************************;
*DEFINE LAMBDA;
lambda_all=(do(0.0,10,0.5))`;
/*lambda_all={0.3,0.5,0.75, 0.8};*/
n_lam=nrow(lambda_all);
Iden=I(ncol(x));
print n_lam;
mm=ncol(og_x)+1;
rss_res=J(n_lam,k,.);
beta_res1=J(mm,n_lam,.); flg_res1=J(mm,n_lam,.);
beta_res2=J(mm,n_lam,.); flg_res2=J(mm,n_lam,.);
grp_res=J(n,n_lam,.);
lam_res=J(n_lam,1,.);
itr_res=J(n_lam,1,.);
stp_res=J(n_lam,1,.);
yhat_res1=J(n,n_lam,.);
yhat_res2=J(n,n_lam,.);


do lam_it=1 to n_lam;
lambda=lambda_all[lam_it];
print "Possible Lambda Value:" lam_it lambda;
lam_res[lam_it]=lambda;
stp_res[lam_it]=lam_it;
*************************;

x=og_x;
ni=j(1,k,.);
b_old=j(pp, k, .);
s2=j(1,k,.);
xvec=J(n,1,1)||og_x;
yhat=J(n,k,.);


do jj=1 to k;

   xvec=J(n,1,1)||x; *define the design matrix;
   id=loc(group=jj);
   xi=xvec[id,];
   yi=y[id,];
   s2[,jj]=var(yi);
   ni[,jj]=nrow(xi);
   BETA=RidgeReg(xi, yi,lambda);
   b_old[,jj]=RidgeReg(xi, yi,lambda);
   yhat[,jj]=xvec*BETA;
end;

mp=ni/n;
var=sum(s2);
p=ncol(x);



var=sum(var);
init_beta=b_old;
init_mp=mp;
init_var=var;
init_group=group;

idat=yhat||og_y||og_x||group;


*INITIAL RESULTS;
*************************;
xnam="X1":"X20";
bnam="B0":"B20";
compnam="K1":"K2";
nm= {"yhat1" "yhat2" "y"};
gm={"group"};
nmx=nm||xnam||gm;
***********************;



bet_old=init_beta;
mp_old=init_mp;
var_old=init_var;
/*print init_beta[r=bnam c={"BK1" "BK2"} l="Initial:"] init_mp[c={"Pi1" "Pi2"}] init_var[c="Var"] ;*/

x=xvec;
y=og_y;
p=ncol(x);
gamma=J(n,k,.);
maxiter=100;
conv=99999999;
Log_pdf=j(maxiter,1,.);


do qq=1 to maxiter while(conv>0.0000001);


	do tt=1 to k;
		   mu_k=x*bet_old[,tt];
		   gamma_k=mp_old[tt] # MVNormalPDF(y, mu_k, var_old);
		   gamma[,tt]=gamma_k;
	end;
	gamma=gamma/gamma[,+];
	group=gamma[,<:>];
	mp_new=(gamma[+,] / n);*mixprob;
	bet_new=j(p, k, .);


	var_new=J(n, k, .);
	yhat=J(n,k,.);
	MMG=j(n,k,.);
	do jj=1 to k;
		   wk=diag(gamma[,jj]); *BETA;
		   xvec=J(n,1,1)||og_x;
		   x=xvec;
		   Iden=ncol(xvec);
		   beta=inv(x`*wk*x+lambda*Iden)*x`*wk*y; *NOW WE HAVE TO INTRODUCE THE RIDGE PARM INTO THE CALCULATION;
		   bet_new[,jj]=beta;*b_old*;
		   sigma2=gamma[,jj]#((y-x*beta)##2); *VARIANCE;
		   var_new[,jj]=sigma2;
		   yhat_k=x*beta;
		   yhat[,jj]=yhat_k;
	end;
		var_new=var_new[+,+]/n;*s_old*;

	do jj=1 to k;
		mu_knew=x*bet_new[,jj];
		MMG[,jj]=mp_new[jj]#MVNormalPdf(y, mu_knew, var_new);
	end;
		sum_MMG=MMG[,+];
		Log_pdf[qq]=(log(sum_MMG))[+,];

	diff1=abs(bet_old-bet_new);
	diff2=abs(var_old-var_new);
	diff3=abs(mp_old-mp_new);
	diff=abs(diff1+diff2+diff3);
	conv=diff;

	mp_old=mp_new;
	bet_old=bet_new;
	var_old=var_new;
/*	print bet_old mp_old var_old;*/
	group_end=group;
end;
/*print "Number of Iter:" qq;*/
fin_log_pdf=log_pdf[(qq-1)];
num_parm=k*3;
AIC=2*num_parm - 2*fin_log_pdf;
BIC=num_parm*log(n) -2*fin_log_pdf;

bet_final=bet_old;
mp_final=mp_old;
var_final=var_old;

*****
BASIC THRESH HOLD ANALYSIS;
thr=0.05;
flg_k1=(abs(bet_final[,1])>thr)#bet_final[,1];
flg_k2=(abs(bet_final[,2])>thr)#bet_final[,2];
rss_pp=J(1,k,.);
******************************;
*CALCULATE COMPONENT WISE RSS;
******************************;
do hh=1 to ncol(yhat);
			id_yh=loc(group_end=hh)`;
			act_y=og_y[id_yh,];
			cmp_y=yhat[id_yh,hh];
		
		rss_k=(act_y-cmp_y)`*(act_y-cmp_y)+lambda*(bet_final[+,hh])**2;
		rss_pp[,hh]=rss_k;
end;
rss_res[lam_it,]=rss_pp;
beta_res1[,lam_it]=bet_final[,1]; flg_res1[,lam_it]=flg_k1;
beta_res2[,lam_it]=bet_final[,2]; flg_res2[,lam_it]=flg_k2;
grp_res[,lam_it]=group_end;
itr_res[lam_it,]=qq;


/*lam_it*/
************************;
end;*END OF BIG LOOP;
************************;
lam_nas=char(lam_res);
total_rss=rss_res[,+];
print rss_res[c=compnam L="RSS"] total_rss itr_res lam_res[l="Lambda"];
print flg_res1[l="Parameters Not at Risk of Elimination" r=bnam c=lam_nas];
print flg_res2[l="Parameters Not at Risk of Elimination" r=bnam c=lam_nas];
print beta_res1[r=bnam];

parm_dat=lam_res||(beta_res1)`;
rnam="Lambda"||bnam;


nm= {"yhat1" "yhat2" "y" "x1" "x2" "x3" "x4" "x5" "group" };


create parmdat from parm_dat[c=rnam];
append from parm_dat;
close parmdat;

idx_lamk1=(rss_res[,1])[>:<]; *comp1;
idx_lamk2=(rss_res[,1])[>:<]; 

print idx_lamk1 idx_lamk2;
beta_fin1=beta_res1[,idx_lamk1];
beta_fin2=beta_res2[,idx_lamk2];
beta_finr=beta_fin1||beta_fin2;
grp_fin=grp_res[,idx_lamk1];

yhat_fin=J(n,k,.);
do jj=1 to k;
   xvec=J(n,1,1)||og_x; *define the design matrix;
   Beta=beta_finr[,jj];
   yhat_fin[,jj]=xvec*BETA;
end;


idx_lamk1=(rss_res[,1])[>:<]; *comp1;
idx_lamk2=(rss_res[,1])[>:<]; 

print idx_lamk1 idx_lamk2;
beta_fin1=beta_res1[,idx_lamk1];
beta_fin2=beta_res2[,idx_lamk2];
beta_finr=beta_fin1||beta_fin2;
grp_fin=grp_res[,idx_lamk1];

yhat_fin=J(n,k,.);
do jj=1 to k;
   xvec=J(n,1,1)||og_x; *define the design matrix;
   Beta=beta_finr[,jj];
   yhat_fin[,jj]=xvec*BETA;
end;

print beta_fin1[r=bnam c="Parameters of Component k=1" l=""] beta_fin2[c="Parameters of Component k=2" l=""];


fdat=yhat_fin||og_y||og_x||group;
print (fdat[1:10,])[c=nm];
create gmrdat from fdat[colname=nm];
append from fdat;
close gmrdat;

quit;

proc sgplot data=parmdat;
title "Parameter Trajectories component k=1";
	series x=Lambda y=B1;
	series x=Lambda y=B2;
	series x=Lambda y=B3;
	series x=Lambda y=B4;
	series x=Lambda y=B5;
run;



proc contents data=gmrdat noprint out=varnames (keep=name);
run;

proc print data=varnames;
run;

data _null_;
set varnames;
where name =: 'x';
call execute ('proc sgplot data=gmrdat; title "Mixture Regression Results of Y and:'||name||'";
				scatter x='||name||' y=y/ group=group markeroutlineattrs=(thickness=2)
   								  markerattrs=(symbol=circlefilled size=10 );
				scatter x='||name||' y=yhat1/  filledoutlinedmarkers 
   markerfillattrs=(color=green) 
   markeroutlineattrs=(color=green thickness=2)
   markerattrs=(symbol=circlefilled size=3);
				scatter x='||name||' y=yhat2/ filledoutlinedmarkers 
   markerfillattrs=(color=purple) 
   markeroutlineattrs=(color=purple thickness=2)
   markerattrs=(symbol=circlefilled size=3);;
				run;');
run;








**************************************************;
*FEATURE SELECTION USING EVEN LESS VARIABLES;
**************************************************;




quit;
proc fastclus data=gmr_dat out=Clust maxclusters=2 maxiter=3 random=321; var x1 x2 x3; run;
proc iml;
use gmr_dat;
read all into xy;
pp=ncol(xy);
print pp;
read all var {"x1" "x2" "x3"} into x;
y=xy[,pp];
n=nrow(x);
ogg_x=x;
ogg_y=y;

*****************NUMBER OF COMPONENTS;
k=2;
****************;

use Clust;
read all into clusdat;
group=clusdat[,(pp+1)];
print (group[1:10,]);
*QUICK N DIRTY ASSIGNMENT OF INITIAL CLUSTERS;

count=J(n,1,.);
do ii=1 to n;
	count[ii]=ii;
end;

start stand(x); ****FOR STANDARDISING;
	xBar = mean(X);           
	X = X - xBar;                           
	std=std(X);
	X_std = X / std;  
	return (X_std);
finish;

x=stand(x);
y=stand(y);
og_x=x;
og_y=y;
pp=ncol(x)+1;


start MVNormalPDF(y, mu, sigma);
	n=nrow(y);
	p = nrow(sigma);
	phi=j(n,1,.);
	do i=1 to n;
	   yi=y[i,];
	   mui=mu[i,];
	   const = (2*constant("PI"))##(p/2) * sqrt(det(sigma));
	   d = Mahalanobis(yi, mui, sigma);
	   phi[i,] = exp( -0.5*d#d ) / const;
	end;
	 return( phi );
finish;


start RidgeReg(x,y,lambda);
	kn=ncol(x);
	Iden=I(kn);
	bh=inv(x`*x+lambda*Iden)*x`*y;
	yh=x*bh;
    return(bh);
finish;


*************************;
*DEFINE LAMBDA;
lambda_all=(do(0.0,10,0.5))`;
/*lambda_all={0.3,0.5,0.75, 0.8};*/
n_lam=nrow(lambda_all);
Iden=I(ncol(x));
print n_lam;
mm=ncol(og_x)+1;
rss_res=J(n_lam,k,.);
beta_res1=J(mm,n_lam,.); flg_res1=J(mm,n_lam,.);
beta_res2=J(mm,n_lam,.); flg_res2=J(mm,n_lam,.);
grp_res=J(n,n_lam,.);
lam_res=J(n_lam,1,.);
itr_res=J(n_lam,1,.);
stp_res=J(n_lam,1,.);
yhat_res1=J(n,n_lam,.);
yhat_res2=J(n,n_lam,.);


do lam_it=1 to n_lam;
lambda=lambda_all[lam_it];
print "Possible Lambda Value:" lam_it lambda;
lam_res[lam_it]=lambda;
stp_res[lam_it]=lam_it;
*************************;

x=og_x;
ni=j(1,k,.);
b_old=j(pp, k, .);
s2=j(1,k,.);
xvec=J(n,1,1)||og_x;
yhat=J(n,k,.);


do jj=1 to k;

   xvec=J(n,1,1)||x; *define the design matrix;
   id=loc(group=jj);
   xi=xvec[id,];
   yi=y[id,];
   s2[,jj]=var(yi);
   ni[,jj]=nrow(xi);
   BETA=RidgeReg(xi, yi,lambda);
   b_old[,jj]=RidgeReg(xi, yi,lambda);
   yhat[,jj]=xvec*BETA;
end;

mp=ni/n;
var=sum(s2);
p=ncol(x);



var=sum(var);
init_beta=b_old;
init_mp=mp;
init_var=var;
init_group=group;

idat=yhat||og_y||og_x||group;


*INITIAL RESULTS;
*************************;
xnam="X1":"X20";
bnam="B0":"B20";
compnam="K1":"K2";
nm= {"yhat1" "yhat2" "y"};
gm={"group"};
nmx=nm||xnam||gm;
***********************;



bet_old=init_beta;
mp_old=init_mp;
var_old=init_var;
/*print init_beta[r=bnam c={"BK1" "BK2"} l="Initial:"] init_mp[c={"Pi1" "Pi2"}] init_var[c="Var"] ;*/

x=xvec;
y=og_y;
p=ncol(x);
gamma=J(n,k,.);
maxiter=100;
conv=99999999;
Log_pdf=j(maxiter,1,.);


do qq=1 to maxiter while(conv>0.0000001);


	do tt=1 to k;
		   mu_k=x*bet_old[,tt];
		   gamma_k=mp_old[tt] # MVNormalPDF(y, mu_k, var_old);
		   gamma[,tt]=gamma_k;
	end;
	gamma=gamma/gamma[,+];
	group=gamma[,<:>];
	mp_new=(gamma[+,] / n);*mixprob;
	bet_new=j(p, k, .);


	var_new=J(n, k, .);
	yhat=J(n,k,.);
	MMG=j(n,k,.);
	do jj=1 to k;
		   wk=diag(gamma[,jj]); *BETA;
		   xvec=J(n,1,1)||og_x;
		   x=xvec;
		   Iden=ncol(xvec);
		   beta=inv(x`*wk*x+lambda*Iden)*x`*wk*y; *NOW WE HAVE TO INTRODUCE THE RIDGE PARM INTO THE CALCULATION;
		   bet_new[,jj]=beta;*b_old*;
		   sigma2=gamma[,jj]#((y-x*beta)##2); *VARIANCE;
		   var_new[,jj]=sigma2;
		   yhat_k=x*beta;
		   yhat[,jj]=yhat_k;
	end;
		var_new=var_new[+,+]/n;*s_old*;

	do jj=1 to k;
		mu_knew=x*bet_new[,jj];
		MMG[,jj]=mp_new[jj]#MVNormalPdf(y, mu_knew, var_new);
	end;
		sum_MMG=MMG[,+];
		Log_pdf[qq]=(log(sum_MMG))[+,];

	diff1=abs(bet_old-bet_new);
	diff2=abs(var_old-var_new);
	diff3=abs(mp_old-mp_new);
	diff=abs(diff1+diff2+diff3);
	conv=diff;

	mp_old=mp_new;
	bet_old=bet_new;
	var_old=var_new;
/*	print bet_old mp_old var_old;*/
	group_end=group;
end;
/*print "Number of Iter:" qq;*/
fin_log_pdf=log_pdf[(qq-1)];
num_parm=k*3;
AIC=2*num_parm - 2*fin_log_pdf;
BIC=num_parm*log(n) -2*fin_log_pdf;

bet_final=bet_old;
mp_final=mp_old;
var_final=var_old;

*****
BASIC THRESH HOLD ANALYSIS;
thr=0.05;
flg_k1=(abs(bet_final[,1])>thr)#bet_final[,1];
flg_k2=(abs(bet_final[,2])>thr)#bet_final[,2];
rss_pp=J(1,k,.);
******************************;
*CALCULATE COMPONENT WISE RSS;
******************************;
do hh=1 to ncol(yhat);
			id_yh=loc(group_end=hh)`;
			act_y=og_y[id_yh,];
			cmp_y=yhat[id_yh,hh];
		
		rss_k=(act_y-cmp_y)`*(act_y-cmp_y)+lambda*(bet_final[+,hh])**2;
		rss_pp[,hh]=rss_k;
end;
rss_res[lam_it,]=rss_pp;
beta_res1[,lam_it]=bet_final[,1]; flg_res1[,lam_it]=flg_k1;
beta_res2[,lam_it]=bet_final[,2]; flg_res2[,lam_it]=flg_k2;
grp_res[,lam_it]=group_end;
itr_res[lam_it,]=qq;


/*lam_it*/
************************;
end;*END OF BIG LOOP;
************************;
lam_nas=char(lam_res);
total_rss=rss_res[,+];
print rss_res[c=compnam L="RSS"] total_rss itr_res lam_res[l="Lambda"];
print flg_res1[l="Parameters Not at Risk of Elimination" r=bnam c=lam_nas];
print flg_res2[l="Parameters Not at Risk of Elimination" r=bnam c=lam_nas];
print beta_res1[r=bnam];

parm_dat=lam_res||(beta_res1)`;
rnam="Lambda"||bnam;


nm= {"yhat1" "yhat2" "y" "x1" "x2" "x3" "group" };


create parmdat from parm_dat[c=rnam];
append from parm_dat;
close parmdat;

idx_lamk1=(rss_res[,1])[>:<]; *comp1;
idx_lamk2=(rss_res[,1])[>:<]; 

print idx_lamk1 idx_lamk2;
beta_fin1=beta_res1[,idx_lamk1];
beta_fin2=beta_res2[,idx_lamk2];
beta_finr=beta_fin1||beta_fin2;
grp_fin=grp_res[,idx_lamk1];

yhat_fin=J(n,k,.);
do jj=1 to k;
   xvec=J(n,1,1)||og_x; *define the design matrix;
   Beta=beta_finr[,jj];
   yhat_fin[,jj]=xvec*BETA;
end;

print beta_fin1[r=bnam c="Parameters of Component k=1" l=""] beta_fin2[c="Parameters of Component k=2" l=""];


fdat=yhat_fin||og_y||og_x||group;
print (fdat[1:10,])[c=nm];
create gmrdat from fdat[colname=nm];
append from fdat;
close gmrdat;

quit;

proc sgplot data=parmdat;
title "Parameter Trajectories component k=1";
	series x=Lambda y=B1;
	series x=Lambda y=B2;
	series x=Lambda y=B3;
run;



proc contents data=gmrdat noprint out=varnames (keep=name);
run;

proc print data=varnames;
run;

data _null_;
set varnames;
where name =: 'x';
call execute ('proc sgplot data=gmrdat; title "Mixture Regression Results of Y and:'||name||'";
				scatter x='||name||' y=y/ group=group markeroutlineattrs=(thickness=2)
   								  markerattrs=(symbol=circlefilled size=10 );
				scatter x='||name||' y=yhat1/  filledoutlinedmarkers 
   markerfillattrs=(color=green) 
   markeroutlineattrs=(color=green thickness=2)
   markerattrs=(symbol=circlefilled size=3);
				scatter x='||name||' y=yhat2/ filledoutlinedmarkers 
   markerfillattrs=(color=purple) 
   markeroutlineattrs=(color=purple thickness=2)
   markerattrs=(symbol=circlefilled size=3);;
				run;');
run;














**************************************************;
*FEATURE SELECTION USING ONLY ONE VARIABLE;
**************************************************;




quit;
proc fastclus data=gmr_dat out=Clust maxclusters=2 maxiter=3 random=321; var x1; run;
proc iml;
use gmr_dat;
read all into xy;
pp=ncol(xy);
print pp;
read all var {"x1"} into x;
y=xy[,pp];
n=nrow(x);
ogg_x=x;
ogg_y=y;

*****************NUMBER OF COMPONENTS;
k=2;
****************;

use Clust;
read all into clusdat;
group=clusdat[,(pp+1)];
print (group[1:10,]);
*QUICK N DIRTY ASSIGNMENT OF INITIAL CLUSTERS;

count=J(n,1,.);
do ii=1 to n;
	count[ii]=ii;
end;

start stand(x); ****FOR STANDARDISING;
	xBar = mean(X);           
	X = X - xBar;                           
	std=std(X);
	X_std = X / std;  
	return (X_std);
finish;

x=stand(x);
y=stand(y);
og_x=x;
og_y=y;
pp=ncol(x)+1;


start MVNormalPDF(y, mu, sigma);
	n=nrow(y);
	p = nrow(sigma);
	phi=j(n,1,.);
	do i=1 to n;
	   yi=y[i,];
	   mui=mu[i,];
	   const = (2*constant("PI"))##(p/2) * sqrt(det(sigma));
	   d = Mahalanobis(yi, mui, sigma);
	   phi[i,] = exp( -0.5*d#d ) / const;
	end;
	 return( phi );
finish;


start RidgeReg(x,y,lambda);
	kn=ncol(x);
	Iden=I(kn);
	bh=inv(x`*x+lambda*Iden)*x`*y;
	yh=x*bh;
    return(bh);
finish;


*************************;
*DEFINE LAMBDA;
lambda_all=(do(0.0,10,0.5))`;
/*lambda_all={0.3,0.5,0.75, 0.8};*/
n_lam=nrow(lambda_all);
Iden=I(ncol(x));
print n_lam;
mm=ncol(og_x)+1;
rss_res=J(n_lam,k,.);
beta_res1=J(mm,n_lam,.); flg_res1=J(mm,n_lam,.);
beta_res2=J(mm,n_lam,.); flg_res2=J(mm,n_lam,.);
grp_res=J(n,n_lam,.);
lam_res=J(n_lam,1,.);
itr_res=J(n_lam,1,.);
stp_res=J(n_lam,1,.);
yhat_res1=J(n,n_lam,.);
yhat_res2=J(n,n_lam,.);


do lam_it=1 to n_lam;
lambda=lambda_all[lam_it];
print "Possible Lambda Value:" lam_it lambda;
lam_res[lam_it]=lambda;
stp_res[lam_it]=lam_it;
*************************;

x=og_x;
ni=j(1,k,.);
b_old=j(pp, k, .);
s2=j(1,k,.);
xvec=J(n,1,1)||og_x;
yhat=J(n,k,.);


do jj=1 to k;

   xvec=J(n,1,1)||x; *define the design matrix;
   id=loc(group=jj);
   xi=xvec[id,];
   yi=y[id,];
   s2[,jj]=var(yi);
   ni[,jj]=nrow(xi);
   BETA=RidgeReg(xi, yi,lambda);
   b_old[,jj]=RidgeReg(xi, yi,lambda);
   yhat[,jj]=xvec*BETA;
end;

mp=ni/n;
var=sum(s2);
p=ncol(x);



var=sum(var);
init_beta=b_old;
init_mp=mp;
init_var=var;
init_group=group;

idat=yhat||og_y||og_x||group;


*INITIAL RESULTS;
*************************;
xnam="X1":"X20";
bnam="B0":"B20";
compnam="K1":"K2";
nm= {"yhat1" "yhat2" "y"};
gm={"group"};
nmx=nm||xnam||gm;
***********************;



bet_old=init_beta;
mp_old=init_mp;
var_old=init_var;
/*print init_beta[r=bnam c={"BK1" "BK2"} l="Initial:"] init_mp[c={"Pi1" "Pi2"}] init_var[c="Var"] ;*/

x=xvec;
y=og_y;
p=ncol(x);
gamma=J(n,k,.);
maxiter=100;
conv=99999999;
Log_pdf=j(maxiter,1,.);


do qq=1 to maxiter while(conv>0.0000001);


	do tt=1 to k;
		   mu_k=x*bet_old[,tt];
		   gamma_k=mp_old[tt] # MVNormalPDF(y, mu_k, var_old);
		   gamma[,tt]=gamma_k;
	end;
	gamma=gamma/gamma[,+];
	group=gamma[,<:>];
	mp_new=(gamma[+,] / n);*mixprob;
	bet_new=j(p, k, .);


	var_new=J(n, k, .);
	yhat=J(n,k,.);
	MMG=j(n,k,.);
	do jj=1 to k;
		   wk=diag(gamma[,jj]); *BETA;
		   xvec=J(n,1,1)||og_x;
		   x=xvec;
		   Iden=ncol(xvec);
		   beta=inv(x`*wk*x+lambda*Iden)*x`*wk*y; *NOW WE HAVE TO INTRODUCE THE RIDGE PARM INTO THE CALCULATION;
		   bet_new[,jj]=beta;*b_old*;
		   sigma2=gamma[,jj]#((y-x*beta)##2); *VARIANCE;
		   var_new[,jj]=sigma2;
		   yhat_k=x*beta;
		   yhat[,jj]=yhat_k;
	end;
		var_new=var_new[+,+]/n;*s_old*;

	do jj=1 to k;
		mu_knew=x*bet_new[,jj];
		MMG[,jj]=mp_new[jj]#MVNormalPdf(y, mu_knew, var_new);
	end;
		sum_MMG=MMG[,+];
		Log_pdf[qq]=(log(sum_MMG))[+,];

	diff1=abs(bet_old-bet_new);
	diff2=abs(var_old-var_new);
	diff3=abs(mp_old-mp_new);
	diff=abs(diff1+diff2+diff3);
	conv=diff;

	mp_old=mp_new;
	bet_old=bet_new;
	var_old=var_new;
/*	print bet_old mp_old var_old;*/
	group_end=group;
end;
/*print "Number of Iter:" qq;*/
fin_log_pdf=log_pdf[(qq-1)];
num_parm=k*3;
AIC=2*num_parm - 2*fin_log_pdf;
BIC=num_parm*log(n) -2*fin_log_pdf;

bet_final=bet_old;
mp_final=mp_old;
var_final=var_old;

*****
BASIC THRESH HOLD ANALYSIS;
thr=0.05;
flg_k1=(abs(bet_final[,1])>thr)#bet_final[,1];
flg_k2=(abs(bet_final[,2])>thr)#bet_final[,2];
rss_pp=J(1,k,.);
******************************;
*CALCULATE COMPONENT WISE RSS;
******************************;
do hh=1 to ncol(yhat);
			id_yh=loc(group_end=hh)`;
			act_y=og_y[id_yh,];
			cmp_y=yhat[id_yh,hh];
		
		rss_k=(act_y-cmp_y)`*(act_y-cmp_y)+lambda*(bet_final[+,hh])**2;
		rss_pp[,hh]=rss_k;
end;
rss_res[lam_it,]=rss_pp;
beta_res1[,lam_it]=bet_final[,1]; flg_res1[,lam_it]=flg_k1;
beta_res2[,lam_it]=bet_final[,2]; flg_res2[,lam_it]=flg_k2;
grp_res[,lam_it]=group_end;
itr_res[lam_it,]=qq;


/*lam_it*/
************************;
end;*END OF BIG LOOP;
************************;
lam_nas=char(lam_res);
total_rss=rss_res[,+];
print rss_res[c=compnam L="RSS"] total_rss itr_res lam_res[l="Lambda"];
print flg_res1[l="Parameters Not at Risk of Elimination" r=bnam c=lam_nas];
print flg_res2[l="Parameters Not at Risk of Elimination" r=bnam c=lam_nas];
print beta_res1[r=bnam];

parm_dat=lam_res||(beta_res1)`;
rnam="Lambda"||bnam;


nm= {"yhat1" "yhat2" "y" "x1" "group" };


create parmdat from parm_dat[c=rnam];
append from parm_dat;
close parmdat;

idx_lamk1=(rss_res[,1])[>:<]; *comp1;
idx_lamk2=(rss_res[,1])[>:<]; 

print idx_lamk1 idx_lamk2;
beta_fin1=beta_res1[,idx_lamk1];
beta_fin2=beta_res2[,idx_lamk2];
beta_finr=beta_fin1||beta_fin2;
grp_fin=grp_res[,idx_lamk1];

yhat_fin=J(n,k,.);
do jj=1 to k;
   xvec=J(n,1,1)||og_x; *define the design matrix;
   Beta=beta_finr[,jj];
   yhat_fin[,jj]=xvec*BETA;
end;

print beta_fin1[r=bnam c="Parameters of Component k=1" l=""] beta_fin2[c="Parameters of Component k=2" l=""];


fdat=yhat_fin||og_y||og_x||group;
print (fdat[1:10,])[c=nm];
create gmrdat from fdat[colname=nm];
append from fdat;
close gmrdat;

quit;

proc sgplot data=parmdat;
title "Parameter Trajectories component k=1";
	series x=Lambda y=B1;
run;



proc contents data=gmrdat noprint out=varnames (keep=name);
run;

proc print data=varnames;
run;

data _null_;
set varnames;
where name =: 'x';
call execute ('proc sgplot data=gmrdat; title "Mixture Regression Results of Y and:'||name||'";
				scatter x='||name||' y=y/ group=group markeroutlineattrs=(thickness=2)
   								  markerattrs=(symbol=circlefilled size=10 );
				scatter x='||name||' y=yhat1/  filledoutlinedmarkers 
   markerfillattrs=(color=green) 
   markeroutlineattrs=(color=green thickness=2)
   markerattrs=(symbol=circlefilled size=3);
				scatter x='||name||' y=yhat2/ filledoutlinedmarkers 
   markerfillattrs=(color=purple) 
   markeroutlineattrs=(color=purple thickness=2)
   markerattrs=(symbol=circlefilled size=3);;
				run;');
run;










