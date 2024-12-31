% -----------------------------------------------------------------------
%  "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2011)

%    In this program, we calculate the ADF statistic with a fixed
%    lag order.

%  if mflag=1 use ADF with constant & without trend
%  if mflag=2 use ADF with constant & trend
%  if mflag=3 use ADF without constant & trend
%  
%  Note: This code have been adjusted by J. Ma 
%  to properly account for lags.
%-------------------------------------------------------------------------

function estm = ADF_test(y,IC,adflag,mflag)

t1=size(y,1)-1;
const=ones(t1,1);
trend=1:1:t1;trend=trend';

y1  = y(size(y,1)-t1:size(y,1)-1);
dy  = y(2:size(y,1)) - y(1:size(y,1)-1);
dy0 = dy(size(dy,1)-t1+1:size(dy,1));
x=y1;

if mflag==1;  x=[x const]; end
if mflag==2;  x=[x const trend]; end

x1=x;

t2=t1-adflag;
dy01=dy0(size(dy0,1)-t2+1:size(dy0,1));      % from k+1 to the end (including dy0)-@

if IC>0
    ICC=zeros(adflag+1,1);
    estm=zeros(adflag+1,1);
    for k = 0:1:adflag
        x2=x1(size(x1,1)-t2+1:size(x1,1),:);         % from k+1 to the end (including y1 and x)-@
        adflag_k = k;
        if adflag_k>0
            j=1;
            while j<=adflag_k
                x2=[x2 dy(size(dy,1)-t2+1-j:size(dy,1)-j)];     % including k lag variables of dy in x2-@
                j=j+1;
            end
        end

        beta = (x2'*x2)\(x2'*dy01);                       % model A-@
        eps  = dy01 - x2*beta;

        if mflag==1 
            sig = sqrt(diag(eps'*eps/(t1-adflag_k-2)*(x2'*x2)^(-1))); 
        elseif mflag==2 
            sig = sqrt(diag(eps'*eps/(t1-adflag_k-3)*(x2'*x2)^(-1)));
        else
             sig = sqrt(diag(eps'*eps/(t1-adflag_k-1)*(x2'*x2)^(-1)));
        end

        npdf=sum(-1/2*log(2*pi)-1/2*(eps.^2));
        if IC==1               %@ AIC @
            ICC(k+1)=-2*npdf/t2+2*size(beta,1)/t2;
        elseif IC==2           %@ BIC @
            ICC(k+1)=-2*npdf/t2+size(beta,1)*log(t2)/t2;
        end

        estm(k+1) = beta(1)./sig(1);
    end
    [~, lag0]=min(ICC);
    estm=estm(lag0);
else
    x2=x1(size(x1,1)-t2+1:size(x1,1),:);
    if adflag>0
        j=1;
        while j<=adflag
            x2=[x2 dy(size(dy,1)-t2+1-j:size(dy,1)-j)];     %including k lag variables of dy in x2-@
            j=j+1;
        end
    end

    beta = (x2'*x2)\(x2'*dy01);                       % model A-@
    eps  = dy01 - x2*beta;

    if mflag==1
        sig = sqrt(diag(eps'*eps/(t1-adflag-2)*(x2'*x2)^(-1)));
    elseif mflag==2
        sig = sqrt(diag(eps'*eps/(t1-adflag-3)*(x2'*x2)^(-1)));
    else
        sig = sqrt(diag(eps'*eps/(t1-adflag-1)*(x2'*x2)^(-1)));
    end
    estm = beta(1)./sig(1);

end

end
