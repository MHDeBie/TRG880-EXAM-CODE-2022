proc import datafile="C:/Users/Matthew/Documents/1 UNI STUFF/TRG880/TRG_EXAM ACTUAL/QUESTION 4/q4_corr.csv"
dbms=csv out=quant_reg replace;
run;
dm 'odsresults;clear';

proc print data=quant_reg (obs=10);
run;

data q_dat;
	set quant_reg;
		x=x;
		XSQ=x*x; *create the additional variable;
run;

proc print data=q_dat (obs=10);
run;


proc sgplot data=q_dat;
	title "Quantile Regression Data XSQ and Y";
	scatter y=y x=XSQ;
run;


ods graphics on;
proc quantreg data=q_dat 
   plots=(rdplot ddplot reshistogram) OUTEST=Quant_parm;
   model y =  x x*x /  quantile = 0.05 0.1 0.5 0.9 0.95 seed=5432;
   output out=ModelFit p=Predicted;
run;

proc print data=Quant_parm;
run;

data parms;
	set quant_parm(keep= Intercept x xx );
run;
proc print data=parms;
run;

proc reg data=q_dat;
model y = x xsq/ vif tol collin;
title 'OLS Regression for Proc Reg';
run;
************HETEROSCEDASTICITY TEST;
proc autoreg data=q_dat;
   model y =  x x*x / archtest;
   output out=r r=yresid;
run;

proc iml;
	use q_dat;
	read all var "y" into y;
	read all var {"x" "XSQ"} into x;
	n=nrow(x);
	print n;
	y_org=y;
	x_org=x;

use parms;
	read all into beta;
	beta_q=beta`;
	
	kk=ncol(beta_q);
/*	a=0.3; b=0.9;*/
/*	jump=0.3;*/
	quants={0.05 0.1 0.5 0.90 0.95};
	print quants;
	

	yhat=J(n,kk,.);
	yhat_alt={};
	xh=J(n,1,1)||x;
	do ii=1 to kk;
		beta_int=beta_q[,ii];
		yhat_ii=xh*beta_int;
/*		call sort(yhat_ii);*/
		yhat[,ii]=yhat_ii;

		q_int=quants[,ii];
		cq_int=char(q_int);
		nm_q=rowcatc("Q"||cq_int);
		cn1=cn1||nm_q;
	end;
print cn1;

cnb={"B0" "B1" "B2"};
beta_print=beta_q`;
print beta_print[r=cn1 c=cnB L="Beta for Each Quantile"];


qreg_dat=yhat||x||y;
cn2={"x" "XSQ" "y"};
cn=cn1||cn2;
print (qreg_dat[1:10,])[c=cn];

create qreg from qreg_dat[c=cn];
append from qreg_dat;
close qreg;



start LR( x , y);
    n=nrow(x);
	p=ncol(x);
    xm=j(n,1,1)||x;
	bh=inv(xm`*xm)*xm`*y;
	yh=xm*bh;
	sse=(y-yh)`*(y-yh);
	res=y-yh;
	return (res);
finish;

x=x_org;
y=y_org;
resid=LR(x,y);
print (resid[1:10,]);
cn={"y" "Residual"};
resid_dat=y||resid;

create res_dat from resid_dat[c=cn];
append from resid_dat;
close res_dat;

quit;

proc print data=qreg (obs=10);
run;

proc sgplot data=qreg;
title "Different Quantile Regression Lines with Regards to X";
	scatter y=Q0_05 x=x;
	scatter y=Q0_1 x=x;
	scatter y=Q0_5 x=x;
	scatter y=Q0_9 x=x;
	scatter y=Q0_95 x=x;
run;

proc sgplot data=qreg;
title "Different Quantile Regression Lines with Regards to XSQ";
	scatter y=Q0_05 x=XSQ;
	scatter y=Q0_1 x=XSQ;
	scatter y=Q0_5 x=XSQ;
	scatter y=Q0_9 x=XSQ;
	scatter y=Q0_95 x=XSQ;
run;

proc sgplot data=res_dat;
title "Residual Plot";
	scatter y=y x=Residual;
run;
*************************
*****QUANT REG FOR JUST X;
ods graphics on;
proc quantreg data=q_dat OUTEST=Quant_parm;
   model y =  x  /  quantile = 0.05 0.1 0.5 0.9 seed=5432;
   output out=ModelFit p=Predicted;
run;

data parms;
	set quant_parm(keep= Intercept x);
run;

proc iml;
*FOR X ONLY;
	use q_dat;
	read all var "y" into y;
	read all var {"x" } into x;
	n=nrow(x);
	print n;
	y_org=y;
	x_org=x;

use parms;
	read all into beta;
	beta_q=beta`;
	
	kk=ncol(beta_q);
	quants={0.05 0.1 0.5 0.90};
/*	print quants;*/
	

	yhat=J(n,kk,.);
	yhat_alt={};
	xh=J(n,1,1)||x;
	do ii=1 to kk;
		beta_int=beta_q[,ii];
		yhat_ii=xh*beta_int;
		yhat[,ii]=yhat_ii;

		q_int=quants[,ii];
		cq_int=char(q_int);
		nm_q=rowcatc("Q"||cq_int);
		cn1=cn1||nm_q;
	end;
print cn1;

cnb={"B0" "B1" "B2"};
beta_print=beta_q`;
print beta_print[r=cn1 c=cnB L="Beta for Each Quantile"];


qreg_dat=yhat||x||y;
cn2={"x" "y"};
cn=cn1||cn2;
print (qreg_dat[1:10,])[c=cn];

create qreg from qreg_dat[c=cn];
append from qreg_dat;
close qreg;



start LR( x , y);
    n=nrow(x);
	p=ncol(x);
    xm=j(n,1,1)||x;
	bh=inv(xm`*xm)*xm`*y;
	yh=xm*bh;
	sse=(y-yh)`*(y-yh);
	res=y-yh;
	return (res);
finish;

x=x_org;
y=y_org;
resid=LR(x,y);
cn={"y" "Residual"};
resid_dat=y||resid;

create res_dat from resid_dat[c=cn];
append from resid_dat;
close res_dat;

quit;


proc sgplot data=qreg;
title "Different Quantile Regression Lines with Regards to X ONLY";
	scatter y=Q0_05 x=X;
	scatter y=Q0_1 x=X;
	scatter y=Q0_5 x=X;
	scatter y=Q0_9 x=X;
run;

proc sgplot data=res_dat;
title "Residual Plot for X ONLY";
	scatter y=y x=Residual;
run;


