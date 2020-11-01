classdef MPOLY_GPU
    %MPOLY_GPU
    %% DESCRIPTION:
    %  Class of multivariate polynomial
    %
    %% SYNTAX:
    %   P = MPOLY_GPU(n)    return a zero constant polynomial with n variables
    %   P = MPOLY_GPU(n,coef,pow)   return a polynomial with n variables,
    %     coefficients coef, and power pow,
    %
    %% INPUTS:
    %  n: number of variable
    %  coef: list of coefficients k x 1 vector, where k is the number
    %    of monomials
    %  pow: matrix of monomial powers k x n matrix
    %
    %% OUTPUTS:
    %  P: MPOLY polynomial object
    %   -P.n: number of variables
    %   -P.k: number of monomials
    %   -P.coef: list of coefficients
    %   -P.pow: matrix of powers for monomials
    %
    %% EXAMPLE:
    %       P = MPOLY_GPU(3)
    %       P = MPOLY_GPU(3,1,zeros(1,3))
    %
    %% COPYRIGHT:
    %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2020/10/30   Initial Coding (modified from MPOLY)
    
    properties
        n % number of variables
        k % number of monomials
        coef % coefficients
        pow % powers k*n matrix
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = MPOLY_GPU(n,coef,pow)
            switch nargin
                case 0
                    % this case is only for adapting parallel computing e.g., parfor
                    obj.n=0;
                    obj.k=1;
                    obj.coef=0;
                    obj.pow=0;
                case 1
                    if n<0 || norm(n-round(n),1)~=0
                        error('The number of variables n must be nonnegative integer!');
                    end
                    obj.n=n;
                    obj.k=1;
                    obj.coef=gpuArray(0);
                    obj.pow=zeros(1,n,'gpuArray');
                case 2
                    error('coef and pow must be both provided!');
                case 3
                    if n<0
                        error('The number of variables n must be nonnegative integer!');
                    end
                    obj.n=n;
                    obj.k=length(coef);
                    if size(coef,2)>1
                        coef=coef(:);
                    end
                    if isa(coef,'gpuArrays')
                        obj.coef=coef;
                    elseif issparse(coef)
                        obj.coef=gpuArray(full(coef));
                    else
                        obj.coef=gpuArray(coef);
                    end
                    if isa(pow,'gpuArrays')
                        obj.pow=pow;
                    elseif issparse(pow)
                        obj.pow=gpuArray(full(pow));
                    else
                        obj.pow=gpuArray(pow);
                    end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Print functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function disp(obj)
            %MPOLY_GPU.disp
            %% DESCRIPTION:
            %  Display brief information of a polynomial matrix. For displaying
            %  polynomial expression, please use 'sdisp'
            %% SYNTAX:
            %   obj.disp
            %   disp(obj)
            %% INPUTS:
            %  obj: polynomial matrix
            %% OUTPUTS:
            %  matrix information
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,[1;2],[1 0 0;2 1 1]);
            %       p.disp;
            %
            %%
            %  See also MPOLY_GPU, disp
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            [mm,nn]=size(obj);
            fprintf('--------------------------------------------------\n')
            if mm*nn==1
                fprintf(' Polynomial with %d variables and %d monomials.\n',obj.n,obj.k);
            else
                fprintf(' Polynomial matrix of size %d x %d with %d variables.\n',mm,nn,obj(1).n);
            end
            fprintf('--------------------------------------------------\n')
        end
        
        function sdisp(obj,format)
            %MPOLY_GPU.sdisp
            %% DESCRIPTION:
            %  Display the expression of polynomial matrix with given precision
            %% SYNTAX:
            %   obj.sdisp           display polynomial with default format '%.6e'
            %   obj.sdisp(format)   display polynomial with given format
            %% INPUTS:
            %  obj: polynomial matrix
            %  format: string for display format, e.g., '%.3f', '%8.5e'. Default '%.6e'
            %% OUTPUTS:
            %  display polynomial expression
            %% EXAMPLE:
            %       x = MPOLY_GPU.mpolyvars(3);
            %       x.sdisp;
            %       x.sdisp('%.2f');
            %%
            %  See also mpoly, disp
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            if nargin<2
                format='%.6e';
            end
            [nrows,ncols]=size(obj);
            if nrows*ncols==1
                fprintf(' %s\n',func2str(obj,format));
            else
                fprintf(' Matrix of size %d x %d : \n',nrows,ncols);
                for j = 1:ncols
                    for i = 1:nrows
                        fprintf(' (%d, %d): %s\n',i,j,func2str(obj(i,j),format));
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Elements operation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function b = eq(obj1,obj2)
            %MPOLY_GPU.eq
            %% DESCRIPTION:
            %  Check equality of two polynomial matrix
            %% SYNTAX:
            %   obj1 == obj2
            %% INPUTS:
            %  obj1: any
            %  obj2: any
            %% OUTPUTS:
            %  b: true or false
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,[1,2],[1 2 3;4 5 6]);
            %       q = p;
            %       p==q
            %%
            %  See also MPOLY_GPU, MPOLY_GPU.ne
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            b=false;
            if isa(obj1,'MPOLY_GPU') && isa(obj2,'MPOLY_GPU')
                if size(obj1)==size(obj2)
                    if obj1(1).n == obj2(1).n
                        newobj = obj1-obj2;
                        if MPOLY_GPU.iszero(newobj)
                            b=true;
                            return;
                        end
                    end
                end
            end
        end
        
        function b = ne(obj1,obj2)
            %MPOLY_GPU.ne
            %% DESCRIPTION:
            %  Check inequality of two polynomial matrix
            %% SYNTAX:
            %   obj1 ~= obj2
            %% INPUTS:
            %  obj1: any
            %  obj2: any
            %% OUTPUTS:
            %  b: true or false
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,[1,2],[1 2 3;4 5 6]);
            %       q = p+1;
            %       p~=q
            %%
            %  See also MPOLY_GPU, MPOLY_GPU.eq
            b=~eq(obj1,obj2);
        end
        
        function v=degree(obj)
            %MPOLY_GPU.degree
            %% DESCRIPTION:
            %  Get degree of polynomial matrix
            %% SYNTAX:
            %   v = obj.degree
            %% INPUTS:
            %  obj: polynomial matrix
            %% OUTPUTS:
            %  v: scalar for degree
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,[1,2],[1 2 3;4 5 6]);
            %       p.degree
            %%
            %  See also mpoly
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            if numel(obj)==1 % for polynomial scalar
                v = norm(obj(1).pow,inf);
            else % for polynomial matrix
                i = 1:numel(obj);
                v = norm(vertcat(obj(i).pow),inf);
            end
        end
        
        function [coef,monolst] = coefficients(obj)
            %coefficients
            %% DESCRIPTION:
            %  Get list of coefficients and monomials of a polynomial (not
            %  for polynomial matrix)
            %% SYNTAX:
            %   coef = obj.coefficients             get list of coefficients
            %   [coef,monolst] = obj.coefficients   get list of coefficients and monomials.
            %% INPUTS:
            %  obj: polynomial (not polynomial matrix)
            %% OUTPUTS:
            %  coef: list of coefficients monolst: list of monomials (with
            %  1 for leading coefficients)
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,[1,2],[1 2 3;4 5 6]);
            %       coef = p.coefficients
            %       [coef,monolst] = p.coefficients;
            %       monolst.sdisp;
            %%
            %  See also mpoly, mono
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            if numel(obj) > 1
                error('Only for scalar polynomial!');
            end
            if nargout == 1
                coef = obj.coef;
                return;
            else
                coef = obj.coef;
                monolst = obj.mono(1:obj.k);
            end
        end
        
        function newobj = mono(obj,idx)
            %mono
            %% DESCRIPTION:
            %   Get list of monomials for a polynomial (not for matrix)
            %   with given idx list
            %% SYNTAX:
            %   newobj = obj.mono(idx)
            %% INPUTS:
            %  obj: polynomial
            %% OUTPUTS:
            %  newobj: monomial list (with 1 for leading coefficients)
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,[1,2],[1 2 3;4 5 6]);
            %       lst = p.mono(1:2);
            %       lst.sdisp;
            %
            %%
            %  See also mpoly, coefficients
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            if numel(obj)>1
                error('mono is not for matrix');
            end
            if nargin == 1
                idx=1:obj.k;
            end
            nn=length(idx);
            newobj=MPOLY_GPU.zeros(obj.n,nn,1);
            for i=1:nn
                if idx(i) > obj.k || idx(i) <=0
                    error('Monomial index must be positive and not exceed the length of monomials in polynomial');
                end
                newobj(i).coef=gpuArray(1);
                newobj(i).pow=obj.pow(idx(i),:);
            end
        end
        
        function v=eval(obj,x)
            method  = 1; % 1 arrayfun, 2 cellfun, 3 loops
            x=gpuArray(x);
            switch method
                case 1
                    % first method arrayfun
                    i=1:1:numel(obj);
                    f = @(i) obj(i).coef'*prod(x'.^(obj(i).pow),2);
                    v = arrayfun(f,i);
                    v = reshape(v,size(obj));
                case 2
                    % second method cellfun
                    f=@(coef,pow) coef'*prod(x'.^pow,2);
                    v = cellfun(f,{obj.coef},{obj.pow});
                    v = reshape(v,size(obj));
                case 3
                    % third method loops
                    v = nan(size(obj));
                    for i=1:numel(v)
                        %v(i)=obj(i).coef'*prod((ones(obj(i).k,1)*x').^(obj(i).pow),2);
                        v(i)=obj(i).coef'*prod(x'.^(obj(i).pow),2);
                    end
            end
        end
        
        function newobj=subs(obj,x)
            %MPOLY_GPU.subs
            %% DESCRIPTION:
            %  Substitute the default variables of the polynomial/polynomial matrix 
            %  by new expression x.
            %% SYNTAX:
            %   newobj = obj.subs(x)
            %% INPUTS:
            %  x: the expression to replace the default variable
            %% OUTPUTS:
            %  newobj: polynomial/polynomial matrix/double scalar/double matrix
            %% EXAMPLE:
            %       x = polylabvar(2);
            %       f = x(1) + x(1)*x(2);
            %       y = polylabvar(2);
            %       g = f.subs((1+y));
            %       g.sdisp;
            %%
            %  See also MPOLY_GPU, sdisp, plus, mtimes, polylabvar
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2020/08/21    Initial Coding
            method = 1;
            switch method
                case 1
                    % first method (using cellfun)
                    X=repmat(x',obj.k,1);
                    f=@(coef,pow) coef'*prod(X.^pow,2);
                    newobj = cellfun(f,{obj.coef},{obj.pow});
                    newobj = reshape(newobj,size(obj));
                case 2
                    % second method (using arrayfun)
                    i=1:1:numel(obj);
                    X=repmat(x',obj.k,1);
                    f = @(i) obj(i).coef'*prod(X.^(obj(i).pow),2);
                    newobj = arrayfun(f,i);
                    newobj = reshape(newobj,size(obj));
            end
            newobj=newobj.simplify;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Polynomial operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function newobj = plus(obj1,obj2)
            %MPOLY_GPU.plus
            %% DESCRIPTION:
            %  Addition of two polynomial matrices, double scalar with a
            %  polynomial matrix, or double matrix with a polynomial matrix
            %% SYNTAX:
            %   newobj = obj1 + obj2
            %% INPUTS:
            %  obj1: polynomial matrix/double scalar/double matrix
            %  obj2: polynomial matrix/double scalar/double matrix
            %% OUTPUTS:
            %  newobj: polynomial matrix/double scalar/double matrix
            %% EXAMPLE:
            %       m1 = [MPOLY_GPU(3,[1;2],[1 0 0; 0 1 0]); ...
            %               MPOLY_GPU(3,[-1;3],[2 1 2; 1 0 1])];
            %       m2 = m1;
            %       m3 = m1 + m2;
            %       m1.sdisp;
            %       m2.sdisp;
            %       m3.sdisp;
            %%
            %  See also MPOLY_GPU, sdisp, uplus, minus, uminus, sum
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            % convert all type into polynomial matrices
            [op1,op2] = sizeunified(obj1,obj2);
            % polynomial matrix + polynomial matrix
            if op1(1).n ~= op2(1).n || norm(size(op1)-size(op2),1)~=0
                error('dimension not match!');
            else
                newobj = MPOLY_GPU.zeros(op1(1).n,size(op1,1),size(op1,2));
                for i=1:numel(op1)
                    newobj(i).pow=[op1(i).pow; op2(i).pow];
                    newobj(i).coef=[op1(i).coef; op2(i).coef];
                    newobj(i).k=op1(i).k + op2(i).k;
                end
                %newobj=newobj.simplify; % comment for speed
            end
        end
        
        function newobj = uplus(obj)
            %MPOLY_GPU.uplus
            %% DESCRIPTION:
            %  Unary plus for polynomial matrix
            %% SYNTAX:
            %   newobj = + obj
            %   newobj = uplus(obj)
            %% INPUTS:
            %  obj: polynomial matrix
            %% OUTPUTS:
            %  newobj: polynomial matrix
            %% EXAMPLE:
            %       m1 = MPOLY_GPU(3,[1;2],[1 0 0; 0 1 0]);
            %       m2 = + m1;
            %       m1.sdisp;
            %       m2.sdisp;
            %%
            %  See also MPOLY_GPU, sdisp, plus, minus, uminus, sum
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            newobj = obj;
        end
        
        function newobj = sum(obj,d)
            %MPOLY_GPU.sum
            %% DESCRIPTION:
            %  Addition of polynomials.
            %% SYNTAX:
            %   newobj = obj.sum
            %% INPUTS:
            %  obj: polynomial matrix
            %  d: sum directions (1: by column (default), 2: by row)
            %% OUTPUTS:
            %  newobj: polynomial / polynomial matrix
            %% EXAMPLE:
            %%
            %  See also mpoly, simplify, plus
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            %  2020/08/21    Sum by direction
            
            if nargin<2
                d=1;
            end
            [~,nn]=size(obj); % get polynomial matrix size
            
            if isvector(obj) % sum one row or column vector
                newobj=MPOLY_GPU(obj(1).n,vertcat(obj.coef),vertcat(obj.pow));
            elseif d==1 % sum a matrix by column
                newobj=MPOLY_GPU.zeros(obj(1).n,1,nn);
                for j=1:nn
                    newobj(j)=obj(:,j).sum;
                end
            else % sum a matrix by row
                newobj=sum((obj'),1);
                newobj=newobj';
            end
            %newobj=newobj.simplify; % comment for speed
        end
        
        function newobj = prod(obj,d)
            %MPOLY_GPU.prod
            %% DESCRIPTION:
            %  Multiply of a polynomial matrix.
            %% SYNTAX:
            %   newobj = obj.prod
            %% INPUTS:
            %  obj: polynomial matrix
            %  d: product direction (1: by column (default), 2: by row)
            %% OUTPUTS:
            %  newobj: polynomial / polynomial matrix
            %% EXAMPLE:
            %%
            %  See also mpoly, simplify, mtimes
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2020/06/26    Initial Coding
            %  2020/08/21    Prod by direction
            
            if nargin<2
                d=1;
            end
            if min(size(obj))==1 %isvector(obj) % multiply one row or column vector
                newobj=obj(1);
                for i=2:numel(obj)
                    newobj = newobj*obj(i);
                end
            elseif d==1 % multiply a matrix by column
                [mm,nn]=size(obj); % get polynomial matrix size
                newobj=obj(1,:); %MPOLY_GPU.zeros(obj(1).n,1,nn);
                for j=1:nn
                    for i=2:mm
                        newobj(j) = newobj(j)*obj(i,j);
                    end
                    %newobj(j)=obj(:,j).prod;
                end
            else % multiply a matrix by row
                newobj=prod((obj'),1);
                newobj=newobj';
            end
            %newobj=newobj.simplify;
        end
        
        function newobj = minus(obj1,obj2)
            %minus
            %% DESCRIPTION:
            %  Subtraction of two polynomial matrices, double scalar with
            %  a polynomial matrix, or double matrix with a polynomial matrix
            %% SYNTAX:
            %   newobj = obj1 - obj2
            %% INPUTS:
            %  obj1: polynomial matrix/double scalar/double matrix
            %  obj2: polynomial matrix/double scalar/double matrix
            %% OUTPUTS:
            %  newobj: polynomial matrix
            %% EXAMPLE:
            %       x = MPOLY_GPU.mpolyvars(3);
            %       p1 = x(1)+2*x(2);
            %       p2 = x(2);
            %       p3 = p1 - p2;
            %       p1.sdisp;
            %       p2.sdisp;
            %       p3.sdisp;
            %%
            %  See also mpoly, MPOLY_GPU.zeros, MPOLY_GPU.plus, MPOLY_GPU.sdisp,
            %  MPOLY_GPU.uminus
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            newobj = plus(obj1,-obj2);
        end
        
        function newobj = uminus(obj)
            %MPOLY_GPU.uminus
            %% DESCRIPTION:
            %  Unary minus. Return -obj
            %% SYNTAX:
            %   newobj = -obj
            %   newobj = uminus(obj)
            %% INPUTS:
            %  obj: polynomial matrix
            %% OUTPUTS:
            %  newobj: negates the elements of obj
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,1,[1 1 1]);
            %       q = - p;
            %       p.sdisp;
            %       q.sdisp;
            %%
            %  See also mpoly, minus, plus, uplus, sdisp
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            newobj = obj;
            for idx=1:numel(obj)
                newobj(idx).coef = -1*newobj(idx).coef;
            end
        end
        
        function newobj = mtimes(obj1,obj2)
            %MPOLY_GPU.mtimes
            %% DESCRIPTION:
            %  Polynomial matrix multiplication (not for Matlab Coder)
            %
            %% SYNTAX:
            %   mtimes(obj1,obj2)
            %% INPUTS:
            %  obj1: polynomial /double matrix/ double scalar
            %  obj2: polynomial /double matrix/ double scalar
            %
            %% OUTPUTS:
            %  newobj: polynomial matrix
            %
            %% EXAMPLE:
            %       m1 = [x(1),x(2);x(3),x(1)];
            %       m2 = -m1;
            %       m3 = m1*m2;
            %       m1.sdisp;
            %       m2.sdisp;
            %       m3.sdisp;
            %%
            %  See also MPOLY_GPU, sdisp, times
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            %  2019/08/23    Improve for speed
            %  2020/10/21    Improve using indexing and vectorization
            len1 = numel(obj1);
            len2 = numel(obj2);
            % case 1: polynomial scalar * polynomial scalar
            if len1==1 && len2 ==1 && isa(obj1,'MPOLY_GPU') && isa(obj2,'MPOLY_GPU')
                if obj1.n ~= obj2.n
                    error('Dimension not match!');
                else
                    newobj = MPOLY_GPU(obj1.n);
                    newobj.k=obj1.k*obj2.k;
                    % using indexing and vectorization
                    isvectorize=true;
                    if isvectorize
                        ridx=repelem(1:obj1.k,obj2.k)';
                        cidx=repmat(1:obj2.k,1,obj1.k)';
                        newobj.coef=obj1.coef(ridx).*obj2.coef(cidx);
                        newobj.pow=obj1.pow(ridx,:)+obj2.pow(cidx,:);
                    else
                        % using loops
                        newobj.coef= reshape(obj2.coef*obj1.coef',newobj.k,1); % kron(obj1.coef,obj2.coef);
                        %newobj.pow=repmat(obj2.pow,obj1.k,1);
                        cc = cell(obj1.k,1);
                        for i=1:obj1.k
                            cc{i} = obj1.pow(i,:)+obj2.pow;
                        end
                        newobj.pow = cell2mat(cc);
                    end
                end
                return;
            end
            
            % case 2: double scalar * polynomial matrix (no need simplify)
            if isnumeric(obj1) && len1 == 1
                newobj = obj2;
                for i=1:numel(obj2)
                    newobj(i).coef=obj1*newobj(i).coef;
                end
                return;
            elseif isnumeric(obj2) && len2 == 1
                newobj = obj2*obj1;
                return;
            end
            
            % case 3: polynomial scalar * polynomial matrix
            if isa(obj1,'MPOLY_GPU') && len1 == 1 && len2 > 1
                newobj = obj2;
                for i=1:numel(obj2)
                    newobj(i)=obj1*newobj(i);
                end
                return;
            elseif isa(obj2,'MPOLY_GPU') && len2 == 1 && len1 > 1
                newobj = obj2*obj1;
                return;
            end
            
            % case 4: double matrix * polynomial matrix
            if isnumeric(obj1) && len1 > 1 % double matrix * polynomial matrix
                if size(obj1,2) ~= size(obj2,1)
                    error('Dimension not match!');
                elseif isvector(obj1) && isvector(obj2) % double row vector * polynomial column vector
                    aux = obj2;
                    for ii=1:size(obj2,1)
                        aux(ii).coef = obj1(ii)*obj2(ii).coef;
                    end
                    newobj = aux.sum;
                    return;
                else % double matrix * polynomial matrix
                    newobj = MPOLY_GPU.zeros(obj2(1).n,size(obj1,1),size(obj2,2));
                    for i=1:size(obj1,1)
                        for j=1:size(obj2,2)
                            newobj(i,j) = obj1(i,:)*obj2(:,j);
                        end
                    end
                    return;
                end
            elseif isnumeric(obj2) && len2 > 1 % polynomial matrix * double matrix
                newobj = (obj2'*obj1')';
                return;
            end
            
            % case 5: polynomial matrix * polynomial matrix
            if size(obj1,2) ~= size(obj2,1) || obj1(1).n ~= obj2(1).n
                error('Dimension not match!');
            end
            if size(obj1,1)==1 && size(obj2,2)==1 % row vector * col vector
                aux = MPOLY_GPU.zeros(obj1(1).n,1,size(obj2,1));
                for ii=1:size(obj2,1) % we can not vectorize here since it is elementary to pass to monomial multiplications
                    
                    aux(ii).k=obj1(ii).k*obj2(ii).k;
                    aux(ii).coef=reshape(obj2(ii).coef*obj1(ii).coef',aux(ii).k,1);
                    cc = cell(obj1(ii).k,1);
                    for i=1:obj1(ii).k
                        cc{i} = obj1(ii).pow(i,:)+obj2(ii).pow;
                    end
                    aux(ii).pow = cell2mat(cc);
                    %aux(ii).coef=kron(op1(ii).coef,op2(ii).coef); %sparse(newobj.k,1);
                    %aux(ii).pow=repmat(op2(ii).pow,op1(ii).k,1);
                    
                end
                newobj = aux.sum;
            else % matrix * matrix
                newobj = MPOLY_GPU.zeros(obj1(1).n,size(obj1,1),size(obj2,2));
                % vectorize for speed
                for i = 1:size(newobj,1)
                    for j = 1:size(newobj,2)
                        newobj(i,j) = obj1(i,:)*obj2(:,j);
                    end
                end
            end
        end
        
        function newobj = times(obj1,obj2)
            %MPOLY_GPU.times
            %% DESCRIPTION:
            %  Elementwise polynomial matrix multiplication
            %
            %% SYNTAX:
            %   newobj = obj1.*obj2
            %% INPUTS:
            %  obj1: polynomial /double matrix/ double scalar
            %  obj2: polynomial /double matrix/ double scalar
            %
            %% OUTPUTS:
            %  newobj: polynomial matrix
            %
            %% EXAMPLE:
            %       m2 = mpolymat(3,2,1);
            %       m2(1).coef = [1;2];
            %       m2(1).pow = [1 0 0; 0 1 0];
            %       m2(2).coef = [-1;3];
            %       m2(2).pow = [2 1 2; 1 0 1];
            %       m1 = m2;
            %       m3 = mpoly_times(m1,m2);
            %       mpoly_sdisplay(m1);
            %       mpoly_sdisplay(m2);
            %       mpoly_sdisplay(m3);
            %%
            %  See also mpoly, sdisp
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            [obj1,obj2] = sizeunified(obj1,obj2);
            % polynomial matrix .* polynomial matrix
            if obj1(1).n ~= obj2(1).n || norm(size(obj1)-size(obj2),1)~=0
                error('dimension not match!');
            else
                newobj = obj1;
                for i=1:numel(obj1)
                    newobj(i) = obj1(i)*obj2(i);
                end
            end
        end
        
        
        function newobj = mpower(obj,a)
            %MPOLY_GPU.mpower
            %% DESCRIPTION:
            %  Polynomial matrix power
            %% SYNTAX:
            %   newobj = obj^a    compute obj to the a power
            %% INPUTS:
            %  obj: square polynomial matrix
            %  a: exponent (scalar)
            %% OUTPUTS:
            %  newobj: sqaure polynomial matrix
            %% EXAMPLE:
            %       x = MPOLY_GPU.mpolyvars(4);
            %       m = [x(1)^2, x(3)^2; x(2)^2, x(4)^2];
            %       a = 3;
            %       p = m^a;
            %       m.sdisp;
            %       p.sdisp;
            %%
            %  See also mpoly, power
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            %  2020/10/21    Improve fast power
            
            if size(obj,1)~=size(obj,2)
                error('obj must be a square matrix!');
            elseif numel(a)~=1
                error('a must be a double scalar!');
            else
                if a==0
                    if MPOLY_GPU.iszero(obj)
                        % return 0
                        newobj = obj;
                    else
                        % return identity matrix
                        newobj = MPOLY_GPU.identity(obj(1).n,size(obj,1));
                    end
                elseif a==1
                    % return obj
                    newobj = obj;
                elseif a>1 && norm(a-round(a),1)==0
                    % return obj^a
                    method=1; % 1 for big power, 2 for naive power
                    switch method
                        case 1
                            % big power
                            newobj=obj*obj;
                            if a==2
                                % return obj^2
                                %newobj=newobj.simplify;
                                return;
                            end
                            if mod(a,2)==0
                                % return (obj^2)^(floor(a/2))
                                newobj = newobj^(floor(a/2));
                                %newobj=newobj.simplify;
                            else
                                % return obj*(obj^2)^(floor(a/2))
                                newobj = newobj^(floor(a/2));
                                newobj = obj*newobj;
                                %newobj=newobj.simplify;
                            end
                        case 2
                            % naive power
                            newobj = obj;
                            for i=2:a
                                newobj = obj*newobj;
                            end
                    end
                else
                    error('Exponent can not be fractional or negative!');
                end
            end
        end
        
        function newobj = power(obj,a)
            %MPOLY_GPU.power
            %% DESCRIPTION:
            %  Elementwise polynomial matrix power .^
            %% SYNTAX:
            %   newobj = obj.^a    raises each element of obj to the corresponding power in a.
            %% INPUTS:
            %  obj: polynomial matrix
            %  a: exponent matrix
            %% OUTPUTS:
            %  newobj: polynomial matrix
            %% EXAMPLE:
            %       x = MPOLY_GPU.mpolyvars(2);
            %       m = [x(1)^2; x(2)];
            %       a = 3;
            %       p = m.^a;
            %       m.sdisp;
            %       p.sdisp;
            %
            %       a = [2;3];
            %       p = m.^a;
            %       m.sdisp;
            %       p.sdisp;
            %
            %%
            %  See also mpoly, mpower, sdisp, MPOLY_GPU.mpolyvars
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            
            newobj=obj;
            if isnumeric(a) && numel(a) == 1
                for i=1:numel(obj)
                    newobj(i)=newobj(i)^a;
                end
            elseif isnumeric(a) && norm(size(obj)-size(a),1)==0
                for i=1:numel(obj)
                    newobj(i)=newobj(i)^a(i);
                end
            else
                error('Dimension not match!');
            end
        end
        
        function newobj = mrdivide(obj,a)
            % /
            if isnumeric(a) && numel(a)==1
                if a==0
                    error('Can not divide zero');
                end
                newobj = (1/a)*obj;
            else
                error('Can only divide double scalar!');
            end
        end
        
        function newobj = rdivide(obj,a)
            % ./
            if isnumeric(a) && norm(size(obj)-size(a),1)==0
                if any(a==0)
                    error('Can not divide zero');
                else
                    newobj = MPOLY_GPU.zeros(obj(1).n,size(obj,1),size(obj,2));
                    for i=1:numel(newobj)
                        newobj(i) = obj(i)/a(i);
                    end
                end
            else
                error('Dimension not match!');
            end
        end
        
        function newobj = diff(obj,idx)
            %diff
            %% DESCRIPTION:
            %  Partial differential for polynomial (not for matrix, using
            %  jacobian for matrix)
            %% SYNTAX:
            %   newobj = obj.diff(idx)
            %% INPUTS:
            %  obj: polynomial idx: a scalar index
            %% OUTPUTS:
            %  newobj: polynomial
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,[1;-2],[1 1 1;2 0 1]);
            %       J = p.diff(1);
            %       p.sdisp; J.sdisp;
            %
            %%
            %  See also MPOLY_GPU, jacobian
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            if numel(obj) >1
                error('Not for matrix! using jacabian for matrix.');
            end
            if idx<=0 || idx>obj.n
                error('Variable index invalid');
            end
            newobj = obj;
            ithpow=obj.pow(:,idx);
            newobj.coef = obj.coef.*ithpow;
            newobj.pow(:,idx) = max(ithpow-1,0);
            newobj = newobj.simplify; % comment this for speed
        end
        
        function newobj = jacobian(obj)
            %jacobian
            %% DESCRIPTION:
            %  Compute jacobian of a polynomial or polynomial matrix
            %% SYNTAX:
            %   newobj = obj.jacobian
            %% INPUTS:
            %  obj: polynomial / polynomial matrix
            %% OUTPUTS:
            %  newobj: polynomial matrix
            %% EXAMPLE:
            %       x = MPOLY_GPU.mpolyvars(3); p = x(1)*x(2) + x(3); J = p.jacobian;
            %       J.sdisp;
            %%
            %  See also MPOLY_GPU, MPOLY_GPU.mpolyvars, MPOLY_GPU.zeros, plus, mtimes, sdisp
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved. 2019/08/13
            %  Initial Coding
            
            newobj=MPOLY_GPU.zeros(obj(1).n,numel(obj),obj(1).n);
            for i=1:size(newobj,1)
                for j=1:size(newobj,2)
                    newobj(i,j) = obj(i).diff(j);
                end
            end
            newobj = newobj.simplify;
        end
        
        function newobj = simplify(obj,zeroprec)
            %MPOLY_GPU.simplify
            %% DESCRIPTION:
            %  simplify a polynomial matrix (by sortrows) and eliminate all monomials with zero
            %  leading coefficients with given precision.
            %% SYNTAX:
            %   newobj = obj.simplify
            %   newobj = obj.simplify(zeroprec)
            %% INPUTS:
            %  zeroprec: precesion of zero, default 1e-8, this is used for eliminate
            %  monomials with zero coefficient, i.e., |coef|<zeroprec
            %% OUTPUTS:
            %  newobj: simplified polynomial matrix
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,[1;2],[1 0 1;1 0 1]);
            %       p.disp;
            %       p.sdisp;
            %       p = p.simplify;
            %       p.disp;
            %       p.sdisp;
            %%
            %  See also MPOLY_GPU, disp, sdisp
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            if nargin<2
                zeroprec = 1e-8;
            end
            nvar=obj(1).n;
            if isscalar(obj) % for scalar polynomial
                tpow = obj.pow;
                tcoef = obj.coef;
                
                % Sort pow
                [tpow,sortidx]= sortrows(tpow);
                % in the sorted power, find non repeated rows
                repidx=[1;any(tpow(2:end,:)~=tpow(1:end-1,:),2)];
                % extract non repeated pow
                tpow=tpow(logical(repidx),:);
                % sum repeated coefs (using sparse matrix)
                tcoef=full(sparse(cumsum(repidx),sortidx,1)*tcoef);
                
                % Eliminate zero coefs
                idx = find(abs(tcoef)>zeroprec);
                tpow=tpow(idx,:);
                tcoef=tcoef(idx);
                
                % if we get empty list, which means zero polynomial
                if isempty(tcoef)
                    newobj = MPOLY_GPU(nvar);
                else
                    % create result polynomial
                    newobj = MPOLY_GPU(nvar,tcoef,tpow);
                end
            else % for matrix polynomial
                newobj=MPOLY_GPU.zeros(nvar,size(obj,1),size(obj,2));
                for i=1:numel(obj)
                    newobj(i) = obj(i).simplify(zeroprec);
                end
            end
        end

        function newobj = simplify_byunique(obj,zeroprec)
            %MPOLY_GPU.simplify
            %% DESCRIPTION:
            %  simplify a polynomial matrix (by unique) and eliminate all monomials 
            %  with zero leading coefficients with given precision.
            %% SYNTAX:
            %   newobj = obj.simplify_byunique
            %   newobj = obj.simplify_byunique(zeroprec)
            %% INPUTS:
            %  zeroprec: precesion of zero, default 1e-8, this is used for eliminate
            %  monomials with zero coefficient, i.e., |coef|<zeroprec
            %% OUTPUTS:
            %  newobj: simplified polynomial matrix
            %% EXAMPLE:
            %       p = MPOLY_GPU(3,[1;2],[1 0 1;1 0 1]);
            %       p.disp;
            %       p.sdisp;
            %       p = p.simplify;
            %       p.disp;
            %       p.sdisp;
            %%
            %  See also MPOLY_GPU, disp, sdisp
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2020/10/21    Initial Coding
            
            if nargin<2
                zeroprec = 1e-8;
            end
            nvar=obj(1).n;
            if isscalar(obj) % for scalar polynomial
                % unique monomials
                [tpow,~,ic]=unique(obj.pow,'rows');
                % sum repeated coefs (using sparse matrix)
                tcoef=full(sparse(ic,1:length(ic),1)*obj.coef);
                
                % Eliminate zero coefs
                idx = find(abs(tcoef(:))>zeroprec);
                tpow=tpow(idx,:);
                tcoef=tcoef(idx);
                
                % if we get empty list, which means zero polynomial
                if isempty(tcoef)
                    newobj = MPOLY_GPU(nvar);
                else
                    % create result polynomial
                    newobj = MPOLY_GPU(nvar);
                    newobj.k=length(tcoef);
                    newobj.coef = tcoef;
                    newobj.pow = tpow;
                end
            else % for matrix polynomial
                newobj=MPOLY_GPU.zeros(nvar,size(obj,1),size(obj,2));
                for i=1:numel(obj)
                    newobj(i) = obj(i).simplify(zeroprec);
                end
            end
        end

        function newobj = diag(obj)
            %MPOLY_GPU.diag
            %% DESCRIPTION:
            %  Create diagonal matrix or get diagonal elements of matrix
            %% SYNTAX:
            %   newobj = obj.diag
            %   newobj = diag(obj)
            %% INPUTS:
            %  obj: polynomial matrix
            %% OUTPUTS:
            %  newobj: polynomial matrix
            %% EXAMPLE:
            %       p = MPOLY_GPU.identity(3,2);
            %       p.sdisp;
            %       p = p.diag;
            %       p.sdisp;
            %%
            %  See also MPOLY_GPU, disp, sdisp, MPOLY_GPU.identity
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/25    Initial Coding
            
            [mm,nn]=size(obj);
            if mm==1 && nn==1
                newobj = obj;
            elseif mm==1 || nn ==1 % row or column vector
                newobj = MPOLY_GPU.zeros(obj(1).n,max([mm,nn]));
                for i=1:max([mm,nn])
                    newobj(i,i) = obj(i);
                end
            else % matrix
                newobj = MPOLY_GPU.zeros(obj(1).n,min([mm,nn]),1);
                for i=1:min([mm,nn])
                    newobj(i) = obj(i,i);
                end
            end
        end
        
        function newobj = trace(obj)
            %MPOLY_GPU.trace
            %% DESCRIPTION:
            %  Sum of diagonal elements
            %% SYNTAX:
            %   newobj = obj.trace
            %   newobj = trace(obj)
            %% INPUTS:
            %  obj: polynomial matrix
            %% OUTPUTS:
            %  newobj: polynomial scalar
            %% EXAMPLE:
            %       p = MPOLY_GPU.identity(3,2);
            %       p.trace
            %       p.sdisp;
            %%
            %  See also MPOLY_GPU, disp, sdisp, MPOLY_GPU.identity, MPOLY_GPU.diag
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/25    Initial Coding
            
            [mm,nn]=size(obj);
            if mm~=nn
                error('Must be square matrix!');
            end
            newobj=obj.diag.sum;
        end
        
        function latexstr = latex(obj,varname,format)
            %MPOLY_GPU.latex
            %% DESCRIPTION:
            %  Convert polynomial to latex
            %% SYNTAX:
            %   latexstr = obj.latex(varname,format)
            %   latexstr = latex(obj,varname,format)
            %% INPUTS:
            %  obj: polynomial matrix
            %  varname: variable name string ('x' by default)
            %  format: format for coefficients ('%.6e' by default)
            %% OUTPUTS:
            %  latexstr: latex code
            %% EXAMPLE:
            %       p = MPOLY.identity(3,2);
            %       p.latex
            %%
            %  See also MPOLY_GPU, MPOLY_GPU.identity
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2020/11/1    Initial Coding
            
            if nargin<2
                varname='x';
            end
            if nargin<3
                format='%.6e'; 
            end
            [nrows,ncols]=size(obj);
            if nrows*ncols==1
                latexstr=func2str(obj,format,varname,'latexstr');
            else
                latexstr='\begin{bmatrix} ';
                for i = 1:nrows
                    for j = 1:ncols
                        latexstr= [latexstr, func2str(obj(i,j),format,varname,'latexstr'), ' & '];
                    end
                    latexstr(end-1:end)='';
                    latexstr = [latexstr,'\\ '];
                end
                latexstr(end:end-1)='';
                latexstr = [latexstr,'\end{bmatrix}'];
            end
        end
        
    end
    
    methods(Static)
        
        function mat = zeros(n,mm,nn)
            %MPOLY_GPU.zeros
            %% DESCRIPTION:
            %  Static function for creating a MPOLY polynomial matrix with
            %  all zero polynomials
            %
            %% SYNTAX:
            %   mat = MPOLY_GPU.zeros(n)        create a zero polynomial scalar
            %   mat = MPOLY_GPU.zeros(n,mm)     create a mm x mm zero
            %     polynomial matrix with n variables
            %   mat = MPOLY_GPU.zeros(n,mm,nn)  create a mm x nn zero polynomial
            %     matrix with n variables
            %
            %% INPUTS:
            %  n: number of variables mm,nn: rows and columns of polynomial
            %    matrix
            %
            %% OUTPUTS:
            %  mat: MPOLY polynomial matrix
            %
            %% EXAMPLE:
            %       mat = MPOLY_GPU.zeros(3);
            %       mat.disp;
            %       mat = MPOLY_GPU.zeros(3,2);
            %       mat.disp;
            %       mat = MPOLY_GPU.zeros(3,2,5);
            %       mat.disp;
            %
            %%
            %  See also MPOLY_GPU, disp, sdisp, MPOLY_GPU.ones, MPOLY_GPU.identity
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            switch nargin
                case 0
                    error('Number of variables n must be given!');
                case 1
                    mat = MPOLY_GPU(n);
                    return;
                case 2
                    nn = mm;
                case 3
            end
            a=MPOLY_GPU(n); % a is a zero polynomial
            mat=repmat(a,mm,nn);
        end
        
        function obj = identity(n,nn)
            %MPOLY_GPU.identity
            %% DESCRIPTION:
            %  Static function for creating an identity square polynomial
            %  matrix
            %% SYNTAX:
            %   obj = MPOLY_GPU.identity(n)     create an one constant
            %     polynomial with n variables
            %   obj = MPOLY_GPU.identity(n,nn)  create an nn x nn identity
            %     polynomial matrix with n variables
            %% INPUTS:
            %  n: number of variables
            %  nn: number of rows and columns
            %% OUTPUTS:
            %  obj: polynomial matrix
            %% EXAMPLE:
            %       p = MPOLY_GPU.identity(3);
            %       p.sdisp;
            %       p = MPOLY_GPU.identity(3,2);
            %       p.sdisp;
            %
            %%
            %  See also MPOLY_GPU, sdisp, MPOLY_GPU.zeros, MPOLY_GPU.ones
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            switch nargin
                case 0
                    error('Number of variables n must be given!');
                case 1
                    nn = 1;
            end
            obj = MPOLY_GPU.zeros(n,nn,nn);
            for i=1:nn
                obj(i,i).coef = gpuArray(1);
            end
        end
        
        function obj = ones(n,mm,nn)
            %MPOLY_GPU.ones
            %% DESCRIPTION:
            %  Static function for creating an all ones polynomial matrix
            %% SYNTAX:
            %   obj = MPOLY_GPU.ones(n)         create a scalar polynomial with
            %     constant one
            %   obj = MPOLY_GPU.ones(n,mm)      create an mm x mm
            %     square polynomial matrix with all ones
            %   obj = MPOLY_GPU.ones(n,mm,nn)   create an mm x nn polynomial
            %     matrix with all ones
            %% INPUTS:
            %  n: number of variables mm: number of rows nn: number of
            %    columns
            %% OUTPUTS:
            %  obj: polynomial matrix
            %% EXAMPLE:
            %       p = MPOLY_GPU.ones(3);
            %       p.sdisp;
            %       p = MPOLY_GPU.ones(3,2);
            %       p.sdisp;
            %       p = MPOLY_GPU.ones(3,2,3);
            %       p.sdisp;
            %%
            %  See also MPOLY_GPU, sdisp, MPOLY_GPU.zeros, MPOLY_GPU.identity
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            switch nargin
                case 0
                    error('Number of variables n must be given!');
                case 1
                    obj = MPOLY_GPU(n);
                    obj.coef = gpuArray(1);
                    return;
                case 2
                    nn = mm;
            end
            a = MPOLY_GPU(n);
            a.coef = gpuArray(1);
            obj=repmat(a,mm,nn);
        end
        
        function lst = monolist(n,d)
            %MPOLY_GPU.monolist
            %% DESCRIPTION:
            %  Static function for creating a list of monomials with n
            %  variables and degree up to d
            %% SYNTAX:
            %   lst = MPOLY_GPU.monolist(n,d)
            %% INPUTS:
            %  n: number of variables d: max degree
            %% OUTPUTS:
            %  lst: list of monomials (columnwise polynomial vector)
            %% EXAMPLE:
            %       lst = MPOLY_GPU.monolist(3,5);
            %       lst.disp;
            %       lst.sdisp;
            %%
            %  See also MPOLY_GPU, nextmonopow, MPOLY_GPU.disp, MPOLY_GPU.sdisp, nchoosek
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            lst=MPOLY_GPU.ones(n,nchoosek(n+d,d),1);
            pow=zeros(1,n,'gpuArray');
            for i=1:nchoosek(n+d,d)
                lst(i).pow = pow;
                pow=nextmonopow(n,pow);
            end
        end
        
        function [P,Ph,B,mncoefs] = transmatconvex(n,d)
            %MPOLY_GPU.transmatconvex
            %% DESCRIPTION:
            %  Generate a transformation matrix for convex basis of 
            %  homogenious polynomials with n variables of degree d.
            %% SYNTAX:
            %   P = MPOLY_GPU.transmatconvex(n,d)
            %% INPUTS:
            %  n: number of variables 
            %  d: max degree
            %% OUTPUTS:
            %  P: transformation matrix
            %  Ph: transformation matrix without multinomial coefs.
            %  B: set of all choices to put d balls into n bins.
            %  mncoefs: multinomial coefficients
            %% EXAMPLE:
            %       [P,Ph,B,mncoefs] = MPOLY_GPU.transmatconvex(3,2);
            %%
            %  See also MPOLY_GPU, nextmonopow, MPOLY_GPU.monolist, nchoosek
            %
            %% COPYRIGHT:
            %  Copyright 2020, Yi-Shuai NIU. All Rights Reserved.
            %  2020/07/18    Initial Coding
            
            s=nchoosek(d+n-1,d);
            P = zeros(s,s); % transformation matrix
            Ph = P; % transformation matrix without multinomial coefs.
            B = zeros(n,n); % set of all choices to put d balls into n bins.
            % create B
            B(1,1) = d; 
            for i=1:s-1
                B(i+1,:)=nextmonopow(n,B(i,:));
            end
            % create Ph
            for i=1:s
                for j=1:s
                    Ph(i,j) = prod(B(i,:).^B(j,:));
                end
            end
            % create P
            mncoefs=factorial(d)./prod(factorial(B'));
            P=mncoefs.*Ph;
        end
        
        function vars = mpolyvars(n)
            %MPOLY_GPU.mpolyvars
            %% DESCRIPTION:
            %  Create mpoly variables
            %
            %% SYNTAX:
            %  vars = MPOLY_GPU.mpolyvars(n)
            %
            %% INPUTS:
            %  n: number of variables
            %
            %% OUTPUTS:
            %  vars: MPOLY vector
            %
            %% EXAMPLE:
            %       vars = MPOLY_GPU.mpolyvars(5);
            %       var.sdisp;
            %
            %%
            %  See also mpoly, MPOLY_GPU.zeros
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            if n< 0 && norm(n-round(n),1)~=0
                error('n must be a nonnegative integer!');
            end
            vars = MPOLY_GPU.zeros(n,n,1);
            for i=1:n
                vars(i).coef = gpuArray(1);
                vars(i).pow(i)=1;
            end
        end
        
        function b = iszero(obj)
            %MPOLY_GPU.iszero
            %% DESCRIPTION:
            %  Check a zero polynomial or zero matrix
            %% SYNTAX:
            %   b = MPOLY_GPU.iszero(obj)
            %% INPUTS:
            %  obj: polynomial matrix/ double matrix
            %% OUTPUTS:
            %  b: true or false
            %% EXAMPLE:
            %       MPOLY_GPU.iszero(0)
            %       MPOLY_GPU.iszero([0 0])
            %       MPOLY_GPU.iszero(1)
            %       MPOLY_GPU.iszero(MPOLY_GPU(3))
            %       MPOLY_GPU.iszero(MPOLY_GPU(3,1,[0 0 0]))
            %       MPOLY_GPU.iszero(MPOLY_GPU.zeros(3,5));
            %%
            %  See also mpoly, MPOLY_GPU.zeros
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            if isnumeric(obj)
                if norm(obj,1)==0
                    b=true;
                else
                    b=false;
                end
            elseif isa(obj,'MPOLY_GPU')
                for i=1:numel(obj)
                    if ~(obj(i).k==1 && obj(i).coef(1)==0)
                        b=false;
                        return;
                    end
                end
                b=true;
            else
                error('unknown object type');
            end
        end
        
        function q = quadprod(Q,b)
            %MPOLY_GPU.quadprod
            %% DESCRIPTION:
            %  Fast computation of b'Qb
            %
            %% SYNTAX:
            %   q = MPOLY_GPU.quadprod(Q,b)
            %% INPUTS:
            %  Q: a double symmetric matrix
            %  b: a list of monomials
            %
            %% OUTPUTS:
            %  q: result polynomial
            %
            %% EXAMPLE:
            %       Q = rand(3,3); Q=Q+Q';
            %       x = polylabvar(5);
            %       b = [x(1);x(2)*x(3);x(4)*x(5)];
            %       q = MPOLY_GPU.quadprod(Q,b);
            %       q.sdisp;
            %%
            %  See also MPOLY_GPU, sdisp, times, polylabvar
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2020/10/11    Initial Coding
            
            % get the number of variables
            n=b(1).n;
            % get the lenght of monomials
            lenb=length(b);
            % the size of Q is lenb*lenb
            % the number of monomials in q (without simplification) is lenb*(lenb+1)/2
            lenq=lenb*(lenb+1)/2;
            powq=zeros(lenq,n); % powq is the pow list of quadratic form b'*Q*b
            coefq=ones(lenq,1); % coefq is the leading coefficients of b'*Q*b
            idx=1;
            for i=1:lenb
                for j=i:lenb
                    powq(idx,:)=b(i).pow+b(j).pow;
                    if (i==j)
                        coefq(idx)=Q(i,i); % for Qii*bi*bi
                    else
                        coefq(idx)=2*Q(i,j); % for 2*Qij*bi*bj with i<j
                    end
                    idx=idx+1;
                end
            end
            q=MPOLY_GPU(n,coefq,powq);
        end
        
        function obj = gpulize(cpuobj)
            %MPOLY_GPU.gpulize
            %% DESCRIPTION:
            %  Convert a MPOLY object to MPOLY_GPU object
            %
            %% SYNTAX:
            %   g = MPOLY_GPU.gpulize(p)
            %% INPUTS:
            %  p: a MPOLY object
            %
            %% OUTPUTS:
            %  g: a MPOLY_GPU object
            %
            %% EXAMPLE:
            %       
            %%
            %  See also MPOLY_GPU, MPOLY
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2020/10/11    Initial Coding
            obj= repmpolygpu(cpuobj(1).n,MPOLY_GPU(cpuobj(1).n),size(cpuobj,1),size(cpuobj,2));
            for i=1:numel(obj)
                obj(i)=MPOLY_GPU(cpuobj(i).n,cpuobj(i).coef,cpuobj(i).pow);
            end
        end
        
        function obj = cpulize(gpuobj)
            %MPOLY_GPU.cpulize
            %% DESCRIPTION:
            %  Convert a MPOLY_GPU object to MPOLY_CPU object
            %
            %% SYNTAX:
            %   g = MPOLY_GPU.cpulize(p)
            %% INPUTS:
            %  p: a MPOLY_GPU object
            %
            %% OUTPUTS:
            %  g: a MPOLY object
            %
            %% EXAMPLE:
            %       
            %%
            %  See also MPOLY_GPU, MPOLY
            %
            %% COPYRIGHT:
            %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2020/10/11    Initial Coding
            obj= repmat(MPOLY(gpuobj(1).n),size(gpuobj,1),size(gpuobj,2));
            for i=1:numel(obj)
                obj(i)=MPOLY(gpuobj(i).n,gather(gpuobj(i).coef),gather(gpuobj(i).pow));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some local functions used in class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = func2str(obj,format,varname,stringtype)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert polynomial function to string
    if nargin < 3
        varname='x';
    end
    if nargin < 4
        stringtype='matlabstr';
    end
    str = '';
    for i=1:obj.k
        if obj.coef(i)>=0
            sgn = '+';
        else
            sgn = '-';
        end
        abscoef = sprintf(format,full(abs(obj.coef(i))));
        if sum(obj.pow(i,:))==0
            str = sprintf('%s%s%s ',str,sgn,abscoef);
        else
            switch stringtype
                case 'matlabstr'
                    str = sprintf('%s%s%s*%s',str,sgn,abscoef,monostr(obj.pow(i,:),varname));
                case 'latexstr'
                    str = sprintf('%s%s%s %s',str,sgn,abscoef,monolatexstr(obj.pow(i,:),varname));
            end
        end
    end
    if str(1)=='+'
        str(1)=''; % remove first +
    end
end

function str = monostr(pow,varname)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build string for monomial
    str='';
    for i=1:length(pow)
        if pow(i)~=0
            if pow(i)~=1
                str = sprintf('%s%s(%d)^%d*',str,varname,i,full(pow(i)));  % case pow(i)>1
            else
                str = sprintf('%s%s(%d)*',str,varname,i); % case pow(i)=1
            end
        end
    end
    if str(end)=='*'
        str(end)=''; % remove last *
    end
end


function str = monolatexstr(pow,varname)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build latex string for monomial
    str='';
    for i=1:length(pow)
        if pow(i)~=0
            if pow(i)~=1
                str = sprintf('%s%s_{%d}^{%d}',str,varname,i,full(pow(i)));  % case pow(i)>1
            else
                str = sprintf('%s%s_{%d}',str,varname,i); % case pow(i)=1
            end
        end
    end
end

function [newobj1,newobj2] = sizeunified(obj1,obj2)
    %sizeunified
    %% DESCRIPTION:
    %  This function is used for converting obj1 and obj2 into same size polynomial matrix.
    %  Some functions need operands to be same size, e.g., MPOLY_GPU.plus,
    %  MPOLY_GPU.times etc. It can also be used to convert a double scalar or double matrix
    %  operand into polynomial matrix of right size.
    %% SYNTAX:
    %   [newobj1,newobj2] = sizeunified(obj1,obj2)
    %% INPUTS:
    %  obj1: polynomial matrix /double matrix/ double scalar
    %  obj2: polynomial matrix /double matrix/ double scalar
    %% OUTPUTS:
    %  newobj1: polynomial matrix /double matrix/ double scalar
    %  newobj2: polynomial matrix /double matrix/ double scalar
    %% EXAMPLE:
    %
    %%
    %  See also mpoly_plus, mpoly_times
    %
    %% COPYRIGHT:
    %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/08/13    Initial Coding
    
    len1 = numel(obj1);
    len2 = numel(obj2);
    % case 1: polynomial matrix with a scalar
    if isnumeric(obj1) && len1==1
        newobj1=repmpolygpu(obj2(1).n,obj1,size(obj2,1),size(obj2,2));
        newobj2=obj2;
        return;
    elseif isnumeric(obj2) && len2==1
        newobj1=obj1;
        newobj2=repmpolygpu(obj1(1).n,obj2,size(obj1,1),size(obj1,2));
        return;
    end
    % case 2: polynomial matrix and double matrix
    if isnumeric(obj1) && isa(obj2,'MPOLY_GPU') && len1 > 1 && len2 >1
        newobj1=double2mpolygpu(obj2(1).n,obj1);
        newobj2=obj2;
        return;
    elseif isnumeric(obj2) && isa(obj1,'MPOLY_GPU') && len2 > 1 && len1 >1
        newobj1=obj1;
        newobj2=double2mpolygpu(obj1(1).n,obj2);
        return;
    end
    % case 3: polynomial scalar and polynomial matrix
    if isa(obj1,'MPOLY_GPU') && isa(obj2,'MPOLY_GPU') && len1 ==1 && len2 > 1
        newobj1 = repmpolygpu(obj1.n,obj1,size(obj2,1),size(obj2,2));
        newobj2 = obj2;
        return;
    elseif isa(obj1,'MPOLY_GPU') && isa(obj2,'MPOLY_GPU') && len1 >1 && len2 == 1
        newobj1 = obj1;
        newobj2 = repmpolygpu(obj2.n,obj2,size(obj1,1),size(obj1,2));
        return;
    end
    % case 4: polynomial scalar and double matrix
    if isa(obj1,'MPOLY_GPU') && isnumeric(obj2) && len1 ==1 && len2 > 1
        newobj1 = repmpolygpu(obj1.n,obj1,size(obj2,1),size(obj2,2));
        newobj2 =double2mpolygpu(obj1.n,obj2);
        return;
    elseif isnumeric(obj1) && isa(obj2,'MPOLY_GPU') && len1 >1 && len2 == 1
        newobj2 = repmpolygpu(obj2.n,obj2,size(obj1,1),size(obj1,2));
        newobj1 = double2mpolygpu(obj2.n,obj1);
        return;
    end
    % case 5: both polynomial matrices
    if isa(obj1,'MPOLY_GPU') && isa(obj2,'MPOLY_GPU')
        newobj1=obj1;
        newobj2=obj2;
    end
end

function p = repmpolygpu(n,r,mm,nn)
    %repmpolygpu
    %% DESCRIPTION:
    %  repeat (double or mpoly_gpu polynomial) r as a mm x nn block polynomial matrix
    %
    %% SYNTAX:
    %   p = repmpolygpu(n,r,mm,nn)
    %
    %% INPUTS:
    %  n: number of variables
    %  r: scalar (double or mpoly)
    %  mm: number of rows
    %  nn: number of columns
    %
    %% OUTPUTS:
    %  p: polynomial matrix
    %
    %% EXAMPLE:
    %       p = repmpolygpu(3,2,5,6);
    %
    %%
    %  See also mpoly_gpu
    %
    %% COPYRIGHT:
    %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/08/13    Initial Coding
    
    p = MPOLY_GPU.zeros(n,mm,nn);
    for i=1:mm
        for j=1:nn
            if isnumeric(r) && numel(r)==1
                p(i,j).coef=r;
            elseif isa(r,'MPOLY_GPU')
                p(i,j)=r;
            end
        end
    end
end

function p = double2mpolygpu(n,r)
    %DOUBLE2MPOLYGPU - convert a double matrix to polynomial matrix
    %
    %% SYNTAX:
    %   p = double2mpolygpu(n,r)
    %
    %% INPUTS:
    %   n: number of variables
    %   r: double matrix
    %
    %% OUTPUTS:
    %   p: polynomial matrix, see mpolymat for more information
    %
    %% EXAMPLE:
    %   p = double2mpolygpu(3,rand(2,2));
    %
    %%
    % See also mpoly
    %
    %% COPYRIGHT:
    % Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
    % 2019/08/13    Initial Coding
    
    % convert double matrix r to a polynomial matrix with n variables
    [mm,nn]=size(r);
    p = MPOLY_GPU.zeros(n,mm,nn);
    if issparse(r)
        gpur=gpuArray(full(r));
    else
        gpur=gpuArray(r);
    end
    for i=1:numel(p)
        p(i).coef=gpur(i);
    end
end

function x = nextmonopow (n, pow)
    %nextmonopow
    %% DESCRIPTION:
    %  Find next monomial power from the current one
    %% SYNTAX:
    %   x = nextmonopow (n, pow)
    %% INPUTS:
    %  n: number of variables
    %  pow: current power, size of 1 x n
    %% OUTPUTS:
    %  x: target power, size of 1 x n
    %% EXAMPLE:
    %       nextpow = nextmonopow (3,[1 0 0]);
    %%
    %  See also mpoly
    %
    %% COPYRIGHT:
    %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/08/13    Initial Coding
    
    x = pow;
    idx_=find(x,1); % find position of the first nonzero
    if isempty(idx_)
        x(1)=1;
        return;
    else
        idx=idx_(1); % there is in fact only one idex at most, but for coder use, we treat it as a list and choose the first one.
    end
    if idx==1 && idx~=n % first element but not last one
        x(idx)=x(idx)-1;
        x(idx+1)=x(idx+1)+1;
    elseif idx==n % last element
        x(1)=x(n)+1;
        x(n)=0;
    elseif idx~=1 && idx~=n && x(idx)<2 % not first or last one, and the value <=1
        x(idx)=x(idx)-1;
        x(idx+1)=x(idx+1)+1;
    else % not first or last one, and the value >=2
        x(idx+1)=x(idx+1)+1;
        ii=find(x==0,1); % find first zero element
        x(ii)=x(idx)-1;
        x(idx)=0;
    end
end

