function x=polylabvar(n)
    %polylabvar
    %% DESCRIPTION:
    %  Create a polylab variable of dimension n x 1
    %
    %% SYNTAX:
    %  x = polylabvar(n);
    %
    %% INPUTS:
    %  n: number of variables (scalar only).
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
    
    x = MPOLY.mpolyvars(n);
end