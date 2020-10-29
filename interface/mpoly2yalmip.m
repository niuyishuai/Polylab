function p=mpoly2yalmip(P,x)
    %MPOLY2YALMIP
    %% DESCRIPTION:
    %  Convert a polynomial matrix to yalmip object
    %
    %% SYNTAX:
    %  p = mpoly2yalmip(P,x);
    %
    %% INPUTS:
    %  P: matrix of polynomial structure, see mpoly for more details on P.
    %  x: polynomial variables, must be a n x 1 Yalmip sdpvar object
    %
    %% OUTPUTS:
    %  p: yalmip polynomial matrix
    %
    %% EXAMPLE:
    %   n = 3;
    %   P = mpoly(n,[1;-2],[1 0 0;1 2 1]);
    %   x = sdpvar(n,1);
    %   p = mpoly2yalmip(P,x);
    %   sdisplay(p);
    %
    %%
    %  See also MPOLY, MPOLY.zeros, sdpvar
    %
    %% COPYRIGHT:
    %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/08/13    Initial Coding
    
    if numel(x)~=P(1).n
        error('dimensions of P and x must be matched!');
    end
    p=sdpvar(size(P,1),size(P,2));
    for i=1:numel(P)
        p(i)=prod((x(:)*ones(1,P(i).k)).^(P(i).pow'))*P(i).coef(:);
    end
end