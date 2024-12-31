%**************************************************************************
%   "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2015)
    
%   In this program, we calculate critical values for the generalized sup 
%   ADF statistic.

%  Note: This function is edited by J. Ma to provide more details of 
%  critical values
% *************************************************************************
 
function [cv_gsadf,cv_bsadf]=MonteCarlo_SIM(T,swindow0,quantiles,SIM_ROUNDS)

m=SIM_ROUNDS;
dim=T-swindow0+1;

%% %%%% DATA GENERATING PROCESS %%%%%%
rng(0)
e=randn(T,m);
a=T^(-1);
y=cumsum(e+a);

%% %% THE GENERALIZED SUP ADF TEST %%%%%

gsadf=ones(m,1);  
bsadfs=cell(m,1); 

for j=1:m
    % if(mod(j,50)==0)
    %     disp(j)
    % end
    bsadfs_temp = zeros(1, dim);
    for r2=swindow0:1:T
        dim0=r2-swindow0+1;
        rwadft=zeros(dim0,1);
        for r1=1:1:dim0
            rwadft(r1)= ADF_Test(y(r1:r2,j),0,1,1);  % two tail 5% significant level
        end
        bsadfs_temp(r2-swindow0+1)=max(rwadft);
    end
    bsadfs{j} = bsadfs_temp;
    gsadf(j)=max(bsadfs_temp);
end
bsadfs = cell2mat(bsadfs);
cv_gsadf=quantile(gsadf,quantiles);
cv_bsadf=quantile(bsadfs,quantiles);
end
