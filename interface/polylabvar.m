function x=polylabvar(n,type)
    %polylabvar
    %% DESCRIPTION:
    %  Create a polylab variable of dimension n x 1
    %
    %% SYNTAX:
    %  x = polylabvar(n,type);
    %
    %% INPUTS:
    %  n: number of variables (scalar only).
    %  type: 'cpu' or 'gpu'
    %
    %% OUTPUTS:
    %  x: polylab matrix
    %
    %% EXAMPLE:
    %   x = polylabvar(3);
    %   x.disp;
    %
    %%
    %  See also MPOLY, MPOLY.mpolyvars, MPOLY.disp
    %
    %% COPYRIGHT:
    %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/08/13    Initial Coding
    %  2020/10/31    Support GPU
    if nargin<2
        type='cpu';
    end
    switch type
        case 'cpu'
          x = MPOLY.mpolyvars(n);
        case 'gpu'
          x = MPOLY_GPU.mpolyvars(n);
    end
end