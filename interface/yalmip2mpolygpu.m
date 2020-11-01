function P=yalmip2mpolygpu(p,x)
    %yalmip2mpolygpu
    %% DESCRIPTION:
    %  convert a yalmip polynomial matrix to mpoly_gpu polynomial matrix
    %
    %% SYNTAX:
    %  P = yalmip2mpolygpu(p,x);
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
    %   P = yalmip2mpolygpu(p,x);
    %
    %%
    %  See also
    %
    %% COPYRIGHT:
    %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2020/10/31    Initial Coding
    
    n=length(x);
    P = MPOLY_GPU.zeros(n,size(p,1),size(p,2));
    for i=1:numel(p)
        [coef,monos]=coefficients(p(i),x);
        P(i).coef = gpuArray(full(coef));
        k=length(P(i).coef); % number of monomials
        P(i).pow=zeros(k,n,'gpuArray');
        P(i).k=k;
        for j=1:k
            P(i).pow(j,:)=gpuArray(degree(monos(j),x));
        end
    end
end