function[Kappa]=mnl_FindCohensKappa(Tp,Fp,Fn,Tn)
%Tp - Number of True Positives (A)
%Fp - Number of False Positives (C)
%Fn - Number of False Negatives (B)
%Tn - Number of True Negatives (D)

TotalNum=Tp+Fp+Fn+Tn;
Ppos=(Tp+Fp)/TotalNum;
Pneg=(Fn+Tn)/TotalNum;
%Calculate Prob Observed
Po=(Tp+Tn)/TotalNum; %The probability that it will score the trace correctly
%Find Prob expected if random
Pcor=((Tp+Fn)/TotalNum)*Ppos;
Pincor=((Fp+Tn)/TotalNum)*Pneg;
Pexp=Pcor+Pincor;
%Now we can calculate Cohen's Kappa
Kappa=(Po-Pexp)/(1-Pexp);