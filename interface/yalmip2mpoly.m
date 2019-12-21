function P=yalmip2mpoly(p,x)
    %yalmip2mpoly
    %% DESCRIPTION:
    %  convert a yalmip polynomial matrix to mpoly polynomial matrix
    %
    %% SYNTAX:
    %  P = yalmip2mpoly(p,x);
    %
    %% INPUTS:
    %  p: yalmip polynomial matrix
    %  x: sdpvar for polynomial variables
    %
    %% OUTPUTS:
    %  P: mpoly polynomial matrix
    %
    %% EXAMPLE:
    %   x = sdpvar(3,1);
    %   p = [x(1)+x(2)*x(3)-2*x(3)^3; sum(x)];
    %   P = yalmip2mpoly(p,x);
    %   mpoly_sdisplay(P);
    %
    %%
    %  See also
    %
    %% COPYRIGHT:
    %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/08/13    Initial Coding
    
    n=length(x);
    P = MPOLY.zeros(n,size(p,1),size(p,2));
    for i=1:numel(p)
        [P(i).coef,monos]=coefficients(p(i),x);
        k=length(P(i).coef); % number of monomials
        P(i).pow=zeros(k,n);
        P(i).k=k;
        for j=1:k
            P(i).pow(j,:)=degree(monos(j),x);
        end
    end
end