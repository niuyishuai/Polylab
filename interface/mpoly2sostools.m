function p=mpoly2sostools(P,x)
    %MPOLY2SOSTOOLS
    %% DESCRIPTION:
    %  Convert a polynomial matrix to SOSTOOLS object
    %
    %% SYNTAX:
    %  p = mpoly2sostools(P,x);
    %
    %% INPUTS:
    %  P: matrix of polynomial structure, see mpoly for more details on P.
    %  x: polynomial variables, must be a n x 1 mpvar object
    %
    %% OUTPUTS:
    %  p: yalmip polynomial matrix
    %
    %% EXAMPLE:
    %   n = 3;
    %   P = mpoly(n,[1;-2],[1 0 0;1 2 1]);
    %   x = mpvar('x',[n,1]);
    %   p = mpoly2sostools(P,x);
    %   p
    %
    %%
    %  See also MPOLY, mpvar
    %
    %% COPYRIGHT:
    %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/08/13    Initial Coding
    
    if prod(size(x))~=P(1).n
        error('dimensions of P and x must be matched!');
    end
    p=mpvar('p',size(P));
    for i=1:numel(P)
        p(i)=prod((x(:)*ones(1,P(i).k)).^(P(i).pow'))*P(i).coef(:);
    end
end