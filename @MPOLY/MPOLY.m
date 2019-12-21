classdef MPOLY
    %MPOLY
    %% DESCRIPTION:
    %  Class of multivariate polynomial
    %
    %% SYNTAX:
    %   P = mpoly(n)    return a zero constant polynomial with n variables
    %   P = mpoly(n,coef,pow)   return a polynomial with n variables,
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
    %       P = MPOLY(3)
    %       P = MPOLY(3,1,zeros(1,3))
    %
    %% COPYRIGHT:
    %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/04/10   Initial Coding
    %  2019/08/25   Improvement
    
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
        
        function obj = MPOLY(n,coef,pow)
            switch nargin
                case 0
                    % this case is only for adapting parallel computing e.g., parfor
                    obj.n=0;
                    obj.k=1;
                    obj.coef=sparse(0);
                    obj.pow=sparse(1,0);
                case 1
                    if n<0 || norm(n-round(n),1)~=0
                        error('The number of variables n must be nonnegative integer!');
                    end
                    obj.n=n;
                    obj.k=1;
                    obj.coef=sparse(0);
                    obj.pow=sparse(1,n);
                case 2
                    error('coef and pow must be both provided!');
                case 3
                    if n<0 || norm(n-round(n),1)~=0
                        error('The number of variables n must be nonnegative integer!');
                    end
                    obj.n=n;
                    obj.k=length(coef);
                    obj.coef=sparse(coef(:));
                    obj.pow=sparse(pow);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Print functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function disp(obj)
            %MPOLY.disp
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
            %       p = MPOLY(3,[1;2],[1 0 0;2 1 1]);
            %       p.disp;
            %
            %%
            %  See also MPOLY, disp
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
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
            %MPOLY.sdisp
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
            %       x = MPOLY.mpolyvars(3);
            %       x.sdisp;
            %       x.sdisp('%.2f');
            %%
            %  See also mpoly, disp
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
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
            %MPOLY.eq
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
            %       p = mpoly(3,[1,2],[1 2 3;4 5 6]);
            %       q = p;
            %       p==q
            %%
            %  See also MPOLY, MPOLY.ne
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            b=false;
            if isa(obj1,'MPOLY') && isa(obj2,'MPOLY')
                if size(obj1)==size(obj2)
                    if obj1(1).n == obj2(1).n
                        newobj = obj1-obj2;
                        if MPOLY.iszero(newobj)
                            b=true;
                            return;
                        end
                    end
                end
            end
        end
        
        function b = ne(obj1,obj2)
            %MPOLY.ne
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
            %       p = mpoly(3,[1,2],[1 2 3;4 5 6]);
            %       q = p+1;
            %       p~=q
            %%
            %  See also MPOLY, MPOLY.eq
            b=~eq(obj1,obj2);
        end
        
        function v=degree(obj)
            %MPOLY.degree
            %% DESCRIPTION:
            %  Get degree of polynomial matrix
            %% SYNTAX:
            %   v = obj.degree
            %% INPUTS:
            %  obj: polynomial matrix
            %% OUTPUTS:
            %  v: scalar for degree
            %% EXAMPLE:
            %       p = mpoly(3,[1,2],[1 2 3;4 5 6]);
            %       p.degree
            %%
            %  See also mpoly
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
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
            %       p = MPOLY(3,[1,2],[1 2 3;4 5 6]);
            %       coef = p.coefficients
            %       [coef,monolst] = p.coefficients;
            %       monolst.sdisp;
            %%
            %  See also mpoly, mono
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
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
            %       p = MPOLY(3,[1,2],[1 2 3;4 5 6]);
            %       lst = p.mono(1:2);
            %       lst.sdisp;
            %
            %%
            %  See also mpoly, coefficients
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            if numel(obj)>1
                error('mono is not for matrix');
            end
            if nargin == 1
                idx=1:obj.k;
            end
            nn=length(idx);
            newobj=MPOLY.zeros(obj.n,nn,1);
            for i=1:nn
                if idx(i) > obj.k || idx(i) <=0
                    error('Monomial index must be positive and not exceed the length of monomials in polynomial');
                end
                newobj(i).coef=sparse(1);
                newobj(i).pow=obj.pow(idx(i),:);
            end
        end
        
        function v=eval(obj,x)
            method  = 0; % 0 mex
            switch method
                case 0
                    % mex method
                    %v = mexeval(x,horzcat(obj.coef),horzcat(obj.pow));
                    v = nan(size(obj));
                    %x = x(:);
                    for i=1:numel(v)
                        v(i) = mexeval(x,obj(i).coef,obj(i).pow);
                    end
                    %f=@(coef,pow) mexeval(x,coef,pow);
                    %v = cellfun(f,{obj.coef},{obj.pow});
                    %v = reshape(v,size(obj));
                case 1
                    % first method
                    f=@(coef,pow) coef'*prod(x'.^pow,2);
                    v = cellfun(f,{obj.coef},{obj.pow});
                    v = reshape(v,size(obj));
                case 2
                    % second method
                    i=1:1:numel(obj);
                    f = @(i) obj(i).coef'*prod(x'.^(obj(i).pow),2);
                    v = arrayfun(f,i);
                    v = reshape(v,size(obj));
                case 3
                    % third method
                    v = nan(size(obj));
                    for i=1:numel(v)
                        %v(i)=obj(i).coef'*prod((ones(obj(i).k,1)*x').^(obj(i).pow),2);
                        v(i)=obj(i).coef'*prod(x'.^(obj(i).pow),2);
                    end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Polynomial operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function newobj = plus(obj1,obj2)
            %MPOLY.plus
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
            %       m1 = [MPOLY(3,[1;2],[1 0 0; 0 1 0]); ...
            %               MPOLY(3,[-1;3],[2 1 2; 1 0 1])];
            %       m2 = m1;
            %       m3 = m1 + m2;
            %       m1.sdisp;
            %       m2.sdisp;
            %       m3.sdisp;
            %%
            %  See also MPOLY, sdisp, uplus, minus, uminus, sum
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            % convert all type into polynomial matrices
            [op1,op2] = sizeunified(obj1,obj2);
            % polynomial matrix + polynomial matrix
            if op1(1).n ~= op2(1).n || norm(size(op1)-size(op2),1)~=0
                error('dimension not match!');
            else
                newobj = MPOLY.zeros(op1(1).n,size(op1,1),size(op1,2));
                for i=1:numel(op1)
                    newobj(i).pow=[op1(i).pow; op2(i).pow];
                    newobj(i).coef=[op1(i).coef; op2(i).coef];
                    newobj(i).k=op1(i).k + op2(i).k;
                end
                newobj=newobj.simplify; % comment for speed
            end
        end
        
        function newobj = uplus(obj)
            %MPOLY.uplus
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
            %       m1 = MPOLY(3,[1;2],[1 0 0; 0 1 0]);
            %       m2 = + m1;
            %       m1.sdisp;
            %       m2.sdisp;
            %%
            %  See also MPOLY, sdisp, plus, minus, uminus, sum
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            newobj = obj;
        end
        
        function newobj = sum(obj)
            %MPOLY.sum
            %% DESCRIPTION:
            %  Addition of polynomials.
            %% SYNTAX:
            %   newobj = obj.sum
            %% INPUTS:
            %  obj: polynomial matrix
            %% OUTPUTS:
            %  newobj: polynomial / polynomial matrix
            %% EXAMPLE:
            %%
            %  See also mpoly, simplify, plus
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            [~,nn]=size(obj); % get polynomial matrix size
            
            if isvector(obj) % sum one row or column vector
                newobj=MPOLY(obj(1).n);
                newobj.k=sum([obj.k]);
                newobj.coef=vertcat(obj.coef);
                newobj.pow=vertcat(obj.pow);
            else % sum a matrix by column
                newobj=MPOLY.zeros(obj(1).n,1,nn);
                for j=1:nn
                    newobj(j)=obj(:,j).sum;
                end
            end
            newobj=newobj.simplify;
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
            %       x = MPOLY.mpolyvars(3);
            %       p1 = x(1)+2*x(2);
            %       p2 = x(2);
            %       p3 = p1 - p2;
            %       p1.sdisp;
            %       p2.sdisp;
            %       p3.sdisp;
            %%
            %  See also mpoly, MPOLY.zeros, MPOLY.plus, MPOLY.sdisp,
            %  MPOLY.uminus
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            newobj = plus(obj1,-obj2);
        end
        
        function newobj = uminus(obj)
            %MPOLY.uminus
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
            %       p = MPOLY(3,1,[1 1 1]);
            %       q = - p;
            %       p.sdisp;
            %       q.sdisp;
            %%
            %  See also mpoly, minus, plus, uplus, sdisp
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            newobj = obj;
            for idx=1:numel(obj)
                newobj(idx).coef = -1*newobj(idx).coef;
            end
        end
        
        function newobj = mtimes(obj1,obj2)
            %MPOLY.mtimes
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
            %  See also MPOLY, sdisp, times
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            %  2019/08/23    Improve for speed
            
            len1 = numel(obj1);
            len2 = numel(obj2);
            % case 1: polynomial scalar * polynomial scalar
            if len1==1 && len2 ==1 && isa(obj1,'MPOLY') && isa(obj2,'MPOLY')
                if obj1.n ~= obj2.n
                    error('Dimension not match!');
                else
                    newobj = MPOLY(obj1.n);
                    newobj.k=obj1.k*obj2.k;
                    newobj.coef= reshape(obj2.coef*obj1.coef',newobj.k,1); % kron(obj1.coef,obj2.coef);
                    %newobj.pow=repmat(obj2.pow,obj1.k,1);
                    cc = cell(obj1.k,1);
                    for i=1:obj1.k
                        cc{i} = obj1.pow(i,:)+obj2.pow;
                    end
                    newobj.pow = cell2mat(cc);
                end
                return;
            end
            
            % case 2: double scalar * polynomial matrix (no need simplify)
            if isa(obj1,'double') && len1 == 1
                newobj = obj2;
                for i=1:numel(obj2)
                    newobj(i).coef=obj1*newobj(i).coef;
                end
                return;
            elseif isa(obj2,'double') && len2 == 1
                newobj = obj2*obj1;
                return;
            end
            
            % case 3: polynomial scalar * polynomial matrix
            if isa(obj1,'MPOLY') && len1 == 1 && len2 > 1
                newobj = obj2;
                for i=1:numel(obj2)
                    newobj(i)=obj1*newobj(i);
                end
                return;
            elseif isa(obj2,'MPOLY') && len2 == 1 && len1 > 1
                newobj = obj2*obj1;
                return;
            end
            
            % case 4: double matrix * polynomial matrix
            if isa(obj1,'double') && len1 > 1 % double matrix * polynomial matrix
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
                    newobj = MPOLY.zeros(obj2(1).n,size(obj1,1),size(obj2,2));
                    for i=1:size(obj1,1)
                        for j=1:size(obj2,2)
                            newobj(i,j) = obj1(i,:)*obj2(:,j);
                        end
                    end
                    return;
                end
            elseif isa(obj2,'double') && len2 > 1 % polynomial matrix * double matrix
                newobj = (obj2'*obj1')';
                return;
            end
            
            % case 5: polynomial matrix * polynomial matrix
            if size(obj1,2) ~= size(obj2,1) || obj1(1).n ~= obj2(1).n
                error('Dimension not match!');
            end
            if size(obj1,1)==1 && size(obj2,2)==1 % row vector * col vector
                aux = MPOLY.zeros(obj1(1).n,1,size(obj2,1));
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
                newobj = MPOLY.zeros(obj1(1).n,size(obj1,1),size(obj2,2));
                % vectorize for speed
                for i = 1:size(newobj,1)
                    for j = 1:size(newobj,2)
                        newobj(i,j) = obj1(i,:)*obj2(:,j);
                    end
                end
            end
        end
        
        function newobj = mtimes__(obj1,obj2)
            %MPOLY.mtimes
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
            %  See also MPOLY, sdisp, times
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            %  2019/08/23    Improve for speed
            
            %             len1 = numel(obj1);
            %             len2 = numel(obj2);
            %             % case 1: polynomial scalar * polynomial scalar
            %             if len1==1 && len2 ==1 && isa(obj1,'MPOLY') && isa(obj2,'MPOLY')
            %                 if obj1.n ~= obj2.n
            %                     error('Dimension not match!');
            %                 else
            %                     newobj = MPOLY(obj1.n);
            %                     newobj.k=obj1.k*obj2.k;
            %                     newobj.coef=kron(obj1.coef,obj2.coef); %sparse(newobj.k,1);
            %                     newobj.pow=repmat(obj2.pow,obj1.k,1);
            %                     i=1:obj1.k;
            %                     newobj.pow = newobj.pow + vertcat(repmat(obj1.pow(i,:),obj2.k,1));
            %                     %newobj.pow=sparse(newobj.k,obj1.n);
            %                     % vectorize for speed
            %                     %[i,j]=meshgrid(1:obj2.k,1:obj1.k);
            %                     %cline=i+(j-1).*obj1.k;
            %                     %coef1 = obj1.coef(j);
            %                     %coef2 = obj2.coef(i);
            %                     %newobj.coef(cline) = coef1(:) .* coef2(:);
            %                     %newobj.pow(cline,:) = obj1.pow(j,:) + obj2.pow(i,:);
            %                     %newobj=newobj.simplify; % for speed
            %                 end
            %                 return;
            %             end
            %
            %             % case 2: double scalar * polynomial matrix (no need simplify)
            %             if isa(obj1,'double') && len1 == 1
            %                 newobj = obj2;
            %                 for i=1:numel(obj2)
            %                     newobj(i).coef=obj1*newobj(i).coef;
            %                 end
            %                 return;
            %             elseif isa(obj2,'double') && len2 == 1
            %                 newobj = obj2*obj1;
            %                 return;
            %             end
            %
            %             % case 3: polynomial scalar * polynomial matrix
            %             if isa(obj1,'MPOLY') && len1 == 1 && len2 > 1
            %                 newobj = obj2;
            %                 for i=1:numel(obj2)
            %                     newobj(i)=obj1*newobj(i);
            %                 end
            %                 return;
            %             elseif isa(obj2,'MPOLY') && len2 == 1 && len1 > 1
            %                 newobj = obj2*obj1;
            %                 return;
            %             end
            %
            %             % case 4: polynomial matrix * double/polynomial matrix
            %             if isa(obj1,'double') && len1 > 1 % double matrix * polynomial matrix
            %                 op1 = double2mpoly(obj2(1).n,sparse(obj1));
            %                 op2 = obj2;
            %             elseif isa(obj2,'double') && len2 > 1 % polynomial matrix * double matrix
            %                 op1 = obj1;
            %                 op2 = double2mpoly(obj1(1).n,sparse(obj2));
            %             else
            %                 op1 = obj1;
            %                 op2 = obj2;
            %             end
            %             if size(op1,2) ~= size(op2,1) || op1(1).n ~= op2(1).n
            %                 error('dimension not match!');
            %             elseif size(op1,1)==1 && size(op2,2)==1 % row vector * col vector
            %                 aux = MPOLY.zeros(op1(1).n,1,size(op2,1));
            %                 for ii=1:size(op2,1) % we can not vectorize here since it is elementary to pass to monomial multiplications
            %                     aux(ii) = op1(ii)*op2(ii);
            %                 end
            %                 newobj = aux.sum;
            %             else % matrix * matrix
            %                 newobj = MPOLY.zeros(op1(1).n,size(op1,1),size(op2,2));
            %                 % vectorize for speed
            %                 for i = 1:size(newobj,1)
            %                     for j = 1:size(newobj,2)
            %                         newobj(i,j) = op1(i,:)*op2(:,j);
            %                     end
            %                 end
            %             end
        end
        
        function newobj = times(obj1,obj2)
            %MPOLY.times
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
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
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
            %MPOLY.mpower
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
            %       x = MPOLY.mpolyvars(4);
            %       m = [x(1)^2, x(3)^2; x(2)^2, x(4)^2];
            %       a = 3;
            %       p = m^a;
            %       m.sdisp;
            %       p.sdisp;
            %%
            %  See also mpoly, power
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            if size(obj,1)~=size(obj,2)
                error('obj must be a square matrix!');
            elseif numel(a)~=1
                error('a must be a double scalar!');
            else
                if a==0
                    if MPOLY.iszero(obj)
                        newobj = obj;
                    else
                        % return identity matrix
                        newobj = MPOLY.identity(obj(1).n,size(obj,1));
                    end
                elseif a==1
                    % return obj
                    newobj = obj;
                elseif a>1 && norm(a-round(a),1)==0
                    % return obj^a
                    newobj = obj;
                    for i=2:a
                        newobj = obj*newobj;
                    end
                else
                    error('Exponent can not be fractional or negative!');
                end
            end
        end
        
        function newobj = power(obj,a)
            %MPOLY.power
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
            %       x = MPOLY.mpolyvars(2);
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
            %  See also mpoly, mpower, sdisp, MPOLY.mpolyvars
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            
            newobj=obj;
            if isa(a,'double') && numel(a) == 1
                for i=1:numel(obj)
                    newobj(i)=newobj(i)^a;
                end
            elseif isa(a,'double') && norm(size(obj)-size(a),1)==0
                for i=1:numel(obj)
                    newobj(i)=newobj(i)^a(i);
                end
            else
                error('Dimension not match!');
            end
        end
        
        function newobj = mrdivide(obj,a)
            % /
            if isa(a,'double') && numel(a)==1
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
            if isa(a,'double') && norm(size(obj)-size(a),1)==0
                if any(a==0)
                    error('Can not divide zero');
                else
                    newobj = MPOLY.zeros(obj(1).n,size(obj,1),size(obj,2));
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
            %       p = MPOLY(3,[1;-2],[1 1 1;2 0 1]);
            %       J = p.diff(1);
            %       p.sdisp; J.sdisp;
            %
            %%
            %  See also MPOLY, jacobian
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
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
            %       x = MPOLY.mpolyvars(3); p = x(1)*x(2) + x(3); J = p.jacobian;
            %       J.sdisp;
            %%
            %  See also MPOLY, MPOLY.mpolyvars, MPOLY.zeros, plus, mtimes, sdisp
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved. 2019/08/13
            %  Initial Coding
            
            newobj=MPOLY.zeros(obj(1).n,numel(obj),obj(1).n);
            for i=1:size(newobj,1)
                for j=1:size(newobj,2)
                    newobj(i,j) = obj(i).diff(j);
                end
            end
            newobj = newobj.simplify;
        end
        
        function newobj = simplify(obj,zeroprec)
            %MPOLY.simplify
            %% DESCRIPTION:
            %  simplify a polynomial matrix and eliminate all monomials with zero
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
            %       p = MPOLY(3,[1;2],[1 0 1;1 0 1]);
            %       p.disp;
            %       p.sdisp;
            %       p = p.simplify;
            %       p.disp;
            %       p.sdisp;
            %%
            %  See also MPOLY, disp, sdisp
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
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
                tpow=tpow(repidx==1,:);
                % sum repeated coefs (using sparse matrix)
                tcoef=sparse(cumsum(repidx),sortidx,1)*tcoef;
                
                % Eliminate zero coefs
                idx = find(abs(tcoef)>zeroprec);
                tpow=tpow(idx,:);
                tcoef=tcoef(idx);
                
                % if we get empty list, which means zero polynomial
                if isempty(tcoef)
                    newobj = MPOLY(nvar);
                else
                    % create result polynomial
                    newobj = MPOLY(nvar);
                    newobj.k=length(tcoef);
                    newobj.coef = tcoef;
                    newobj.pow = tpow;
                end
            else % for matrix polynomial
                newobj=MPOLY.zeros(nvar,size(obj,1),size(obj,2));
                for i=1:numel(obj)
                    newobj(i) = obj(i).simplify(zeroprec);
                end
            end
        end
        
        function newobj = diag(obj)
            %MPOLY.diag
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
            %       p = MPOLY.identity(3,2);
            %       p.sdisp;
            %       p = p.diag;
            %       p.sdisp;
            %%
            %  See also MPOLY, disp, sdisp, MPOLY.identity
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/25    Initial Coding
            
            [mm,nn]=size(obj);
            if mm==1 && nn==1
                newobj = obj;
            elseif mm==1 || nn ==1 % row or column vector
                newobj = MPOLY.zeros(obj(1).n,max([mm,nn]));
                for i=1:max([mm,nn])
                    newobj(i,i) = obj(i);
                end
            else % matrix
                newobj = MPOLY.zeros(obj(1).n,min([mm,nn]),1);
                for i=1:min([mm,nn])
                    newobj(i) = obj(i,i);
                end
            end
        end
        
        function newobj = trace(obj)
            %MPOLY.trace
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
            %       p = MPOLY.identity(3,2);
            %       p.trace
            %       p.sdisp;
            %%
            %  See also MPOLY, disp, sdisp, MPOLY.identity, MPOLY.diag
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/25    Initial Coding
            
            [mm,nn]=size(obj);
            if mm~=nn
                error('Must be square matrix!');
            end
            newobj=obj.diag.sum;
        end
        
        function newobj = simplify_bis(obj,prec)
            if nargin<2
                prec = 1e-8;
            end
            newobj = mpoly_matlab.mpolymat(obj.n,size(obj));
            for idx=1:numel(obj)
                newobj(idx).pow=[];
                newobj(idx).coef=[];
                tpow=obj(idx).pow;
                tcoef=obj(idx).coef;
                while ~isempty(tcoef)
                    % ��ȡͬ����
                    lst=find(sum(abs(tpow(1,:)-tpow),2)==0);
                    
                    % �ϲ�ͬ����ϵ��
                    v=sum(tcoef(lst));
                    % �ж�ϵ���Ƿ�Ϊ�� (>prec)
                    if abs(v)>prec
                        newobj(idx).coef = [newobj(idx).coef; sum(tcoef(lst))];
                        newobj(idx).pow = [newobj(idx).pow; tpow(1,:)];
                    end
                    tpow(lst,:)=[];
                    tcoef(lst,:)=[];
                end
                newobj(idx).k=length(newobj(idx).coef);
            end
        end
    end
    
    methods(Static)
        
        function mat = zeros(n,mm,nn)
            %MPOLY.zeros
            %% DESCRIPTION:
            %  Static function for creating a MPOLY polynomial matrix with
            %  all zero polynomials
            %
            %% SYNTAX:
            %   mat = MPOLY.zeros(n)        create a zero polynomial scalar
            %   mat = MPOLY.zeros(n,mm)     create a mm x mm zero
            %     polynomial matrix with n variables
            %   mat = MPOLY.zeros(n,mm,nn)  create a mm x nn zero polynomial
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
            %       mat = MPOLY.zeros(3);
            %       mat.disp;
            %       mat = MPOLY.zeros(3,2);
            %       mat.disp;
            %       mat = MPOLY.zeros(3,2,5);
            %       mat.disp;
            %
            %%
            %  See also MPOLY, disp, sdisp, MPOLY.ones, MPOLY.identity
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            switch nargin
                case 0
                    error('Number of variables n must be given!');
                case 1
                    mat = MPOLY(n);
                    return;
                case 2
                    nn = mm;
                case 3
            end
            a=MPOLY(n); % a is a zero polynomial
            mat=repmat(a,mm,nn);
        end
        
        function obj = identity(n,nn)
            %MPOLY.identity
            %% DESCRIPTION:
            %  Static function for creating an identity square polynomial
            %  matrix
            %% SYNTAX:
            %   obj = MPOLY.identity(n)     create an one constant
            %     polynomial with n variables
            %   obj = MPOLY.identity(n,nn)  create an nn x nn identity
            %     polynomial matrix with n variables
            %% INPUTS:
            %  n: number of variables
            %  nn: number of rows and columns
            %% OUTPUTS:
            %  obj: polynomial matrix
            %% EXAMPLE:
            %       p = MPOLY.identity(3);
            %       p.sdisp;
            %       p = MPOLY.identity(3,2);
            %       p.sdisp;
            %
            %%
            %  See also MPOLY, sdisp, MPOLY.zeros, MPOLY.ones
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            switch nargin
                case 0
                    error('Number of variables n must be given!');
                case 1
                    nn = 1;
            end
            obj = MPOLY.zeros(n,nn,nn);
            for i=1:nn
                obj(i,i).coef = sparse(1);
            end
        end
        
        function obj = ones(n,mm,nn)
            %MPOLY.ones
            %% DESCRIPTION:
            %  Static function for creating an all ones polynomial matrix
            %% SYNTAX:
            %   obj = MPOLY.ones(n)         create a scalar polynomial with
            %     constant one
            %   obj = MPOLY.ones(n,mm)      create an mm x mm
            %     square polynomial matrix with all ones
            %   obj = MPOLY.ones(n,mm,nn)   create an mm x nn polynomial
            %     matrix with all ones
            %% INPUTS:
            %  n: number of variables mm: number of rows nn: number of
            %    columns
            %% OUTPUTS:
            %  obj: polynomial matrix
            %% EXAMPLE:
            %       p = MPOLY.ones(3);
            %       p.sdisp;
            %       p = MPOLY.ones(3,2);
            %       p.sdisp;
            %       p = MPOLY.ones(3,2,3);
            %       p.sdisp;
            %%
            %  See also MPOLY, sdisp, MPOLY.zeros, MPOLY.identity
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            switch nargin
                case 0
                    error('Number of variables n must be given!');
                case 1
                    obj = MPOLY(n);
                    obj.coef = sparse(1);
                    return;
                case 2
                    nn = mm;
            end
            a = MPOLY(n);
            a.coef = sparse(1);
            obj=repmat(a,mm,nn);
        end
        
        function lst = monolist(n,d)
            %MPOLY.monolist
            %% DESCRIPTION:
            %  Static function for creating a list of monomials with n
            %  variables and degree up to d
            %% SYNTAX:
            %   lst = MPOLY.monolist(n,d)
            %% INPUTS:
            %  n: number of variables d: max degree
            %% OUTPUTS:
            %  lst: list of monomials (columnwise polynomial vector)
            %% EXAMPLE:
            %       lst = MPOLY.monolist(3,5);
            %       lst.disp;
            %       lst.sdisp;
            %%
            %  See also MPOLY, nextmonopow, MPOLY.disp, MPOLY.sdisp, nchoosek
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            lst=MPOLY.ones(n,nchoosek(n+d,d),1);
            pow=sparse(1,n);
            for i=1:nchoosek(n+d,d)
                lst(i).pow = pow;
                pow=nextmonopow(n,pow);
            end
        end
        
        function vars = mpolyvars(n)
            %MPOLY.mpolyvars
            %% DESCRIPTION:
            %  Create mpoly variables
            %
            %% SYNTAX:
            %  vars = MPOLY.mpolyvars(n)
            %
            %% INPUTS:
            %  n: number of variables
            %
            %% OUTPUTS:
            %  vars: MPOLY vector
            %
            %% EXAMPLE:
            %       vars = MPOLY.mpolyvars(5);
            %       var.sdisp;
            %
            %%
            %  See also mpoly, MPOLY.zeros
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            if n< 0 && norm(n-round(n),1)~=0
                error('n must be a nonnegative integer!');
            end
            vars = MPOLY.zeros(n,n,1);
            for i=1:n
                vars(i).coef = sparse(1);
                vars(i).pow(i)=1;
            end
        end
        
        function b = iszero(obj)
            %MPOLY.iszero
            %% DESCRIPTION:
            %  Check a zero polynomial or zero matrix
            %% SYNTAX:
            %   b = MPOLY.iszero(obj)
            %% INPUTS:
            %  obj: polynomial matrix/ double matrix
            %% OUTPUTS:
            %  b: true or false
            %% EXAMPLE:
            %       MPOLY.iszero(0)
            %       MPOLY.iszero([0 0])
            %       MPOLY.iszero(1)
            %       MPOLY.iszero(MPOLY(3))
            %       MPOLY.iszero(MPOLY(3,1,[0 0 0]))
            %       MPOLY.iszero(MPOLY.zeros(3,5));
            %%
            %  See also mpoly, MPOLY.zeros
            %
            %% COPYRIGHT:
            %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
            %  2019/08/13    Initial Coding
            
            if isa(obj,'double')
                if norm(obj,1)==0
                    b=true;
                else
                    b=false;
                end
            elseif isa(obj,'MPOLY')
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
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some local functions used in class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = func2str(obj,format)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert polynomial function to string
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
            str = sprintf('%s%s%s*%s',str,sgn,abscoef,monostr(obj.pow(i,:)));
        end
    end
end

function str = monostr(pow)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build string for monomial
    str='';
    for i=1:length(pow)
        if pow(i)~=0
            if pow(i)~=1
                str = sprintf('%sx(%d)^%d*',str,i,full(pow(i)));  % case pow(i)>1
            else
                str = sprintf('%sx(%d)*',str,i); % case pow(i)=1
            end
        end
    end
    if str(end)=='*'
        str(end)=char(0); % replace last * as \0
    end
end

function [newobj1,newobj2] = sizeunified(obj1,obj2)
    %sizeunified
    %% DESCRIPTION:
    %  This function is used for converting obj1 and obj2 into same size polynomial matrix.
    %  Some functions need operands to be same size, e.g., MPOLY.plus,
    %  MPOLY.times etc. It can also be used to convert a double scalar or double matrix
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
    %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/08/13    Initial Coding
    
    len1 = numel(obj1);
    len2 = numel(obj2);
    % case 1: polynomial matrix with a scalar
    if isa(obj1,'double') && len1==1
        newobj1=repmpoly(obj2(1).n,sparse(obj1),size(obj2,1),size(obj2,2));
        newobj2=obj2;
        return;
    elseif isa(obj2,'double') && len2==1
        newobj1=obj1;
        newobj2=repmpoly(obj1(1).n,sparse(obj2),size(obj1,1),size(obj1,2));
        return;
    end
    % case 2: polynomial matrix and double matrix
    if isa(obj1,'double') && isa(obj2,'MPOLY') && len1 > 1 && len2 >1
        newobj1=double2mpoly(obj2(1).n,sparse(obj1));
        newobj2=obj2;
        return;
    elseif isa(obj2,'double') && isa(obj1,'MPOLY') && len2 > 1 && len1 >1
        newobj1=obj1;
        newobj2=double2mpoly(obj1(1).n,sparse(obj2));
        return;
    end
    % case 3: polynomial scalar and polynomial matrix
    if isa(obj1,'MPOLY') && isa(obj2,'MPOLY') && len1 ==1 && len2 > 1
        newobj1 = repmpoly(obj1.n,obj1,size(obj2,1),size(obj2,2));
        newobj2 = obj2;
        return;
    elseif isa(obj1,'MPOLY') && isa(obj2,'MPOLY') && len1 >1 && len2 == 1
        newobj1 = obj1;
        newobj2 = repmpoly(obj2.n,obj2,size(obj1,1),size(obj1,2));
        return;
    end
    % case 4: polynomial scalar and double matrix
    if isa(obj1,'MPOLY') && isa(obj2,'double') && len1 ==1 && len2 > 1
        newobj1 = repmpoly(obj1.n,obj1,size(obj2,1),size(obj2,2));
        newobj2 =double2mpoly(obj1.n,obj2);
        return;
    elseif isa(obj1,'double') && isa(obj2,'MPOLY') && len1 >1 && len2 == 1
        newobj2 = repmpoly(obj2.n,obj2,size(obj1,1),size(obj1,2));
        newobj1 = double2mpoly(obj2.n,obj1);
        return;
    end
    % case 5: both polynomial matrices
    if isa(obj1,'MPOLY') && isa(obj2,'MPOLY')
        newobj1=obj1;
        newobj2=obj2;
    end
end

function p = repmpoly(n,r,mm,nn)
    %repmpoly
    %% DESCRIPTION:
    %  repeat scalar (double or mpoly polynomial) r as a mm x nn polynomial matrix
    %
    %% SYNTAX:
    %   p = repmpoly(n,r,mm,nn)
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
    %       p = repmpoly(3,2,5,6);
    %
    %%
    %  See also mpoly, mpolymat
    %
    %% COPYRIGHT:
    %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2019/08/13    Initial Coding
    
    p = MPOLY.zeros(n,mm,nn);
    for i=1:mm
        for j=1:nn
            if isa(r,'double') && numel(r)==1
                p(i,j).coef=r;
            elseif isa(r,'MPOLY')
                p(i,j)=r;
            end
        end
    end
end

function p = double2mpoly(n,r)
    %DOUBLE2MPOLYMAT - convert a doouble matrix to polynomial matrix
    %
    %% SYNTAX:
    %   p = double2mpoly(n,r)
    %
    %% INPUTS:
    %   n: number of variables
    %   r: double matrix
    %
    %% OUTPUTS:
    %   p: polynomial matrix, see mpolymat for more information
    %
    %% EXAMPLE:
    %   p = double2mpoly(3,rand(2,2));
    %
    %%
    % See also mpoly
    %
    %% COPYRIGHT:
    % Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
    % 2019/08/13    Initial Coding
    
    % convert double matrix r to a polynomial matrix with n variables
    [mm,nn]=size(r);
    p = MPOLY.zeros(n,mm,nn);
    for i=1:numel(p)
        p(i).coef=r(i);
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
    %  Copyright 2019, Yi-Shuai NIU. All Rights Reserved.
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
