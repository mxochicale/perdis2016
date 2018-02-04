function [Ixx] = ami(x,n,m,nb)
% AMI   calculates average mutual information for the input scalar time 
%       series. if no output is requested it plots the AMI curve.
%
% CALL: Ixx=ami(x,n,m,nb);
%
% x  - input scalar time series
% n  - number of points to use in the calculation
% m  - maximum delay parameter to consider
% nb - number of bins to use in calculation (default=16)
% p  - length of continuously sampled data
% Ixx - average mutual information in bits
%
% Copyright and Revised by David Chelidze 2-11-2009
 
if nargin < 3
   error('Requires at least three input arguments!')   
elseif nargin == 3
   nb = 16; % set the number of bins to default
end
 
% check data and make a row vector of scalar data
if min(size(x)) ~= 1
   error('first input should be a vector of scalar data points.')
else
    x=x(:);
end
 
% get the matrix of delayed versions of x
 
nmax = size(x,1) - m; % <- that’s the maximum size possible!!!
 
% get the maximum possible matrix size!!!
if (n == 0) || (n > nmax)
    warning(['Specified data size is zero or exceeds maximum possible'...
        ' value!\n Continuing using maximum possible size!\n'])
    n = nmax;
end
 
y = zeros(n,m+1); % initialize array for embedding full subsets
 
for i = 0:m
    y(:,i+1) = x(1+i:n+i);
end
 
% obtain Px
x1 = y(:,1);
Px = hist(x1, nb); Px = Px/n;
[~, ~, NPx] = find(Px);
 
% initialize Ixx
Ixx = zeros(m+1,1);
 
for ii = 1:m+1 % calculate the average mutual information
   x2 = y(:,ii);
   [Pxy, ~, ~] = joint_prob(x1, x2, nb);
   [~, ~, NPxy] = find(Pxy);
   Ixx(ii) = sum(NPxy.*log2(NPxy)) - 2*sum(NPx.*log2(NPx));
end
 
if nargout == 0
    plot(0:m,Ixx,'o');
    xlabel('Delay (sampling time)');
    ylabel('Average Mutual Information (bits)');
end

%%-------------------------------------------------------------------------
function [Pxy,xx,yy]=joint_prob(x,y,nbx,nby)
%
% Determins the joint Probability P(x,y) based on
% distributed sets of scalar variables x and y
%
% call [Pxy,xx,yy]=joint_prob(x,y,nbx,nby);
% Inputs:
% x & y - vectors (time series)
% nbx   - number of bins for x
% nby   - number of bins for y
% Output: Pxy - (nbx,nby) matrix 
% It plots results if no output is requested
%
% Copyright and Revised by David Chelidze 7-31-2000
% >>>> Needs cond_prob.m <<<<

if (nargin == 0 || nargin == 1)
   error('Requires at least two input arguments')
end
if length(x) ~= length(y)
   error('Both arrays should be the same size')
end
if nargin == 2
   nbx=16; nby=16;
elseif nargin == 3
   nby=nbx;
end
% put in column form
y=y(:);x=x(:);
[Px_y,xx,yy]=cond_prob(x,y,nbx,nby);
Py=histc(y,yy); Py=Py(:)/length(y); Py=Py(1:end-1);
Pxy=Px_y.*(ones(nbx,1)*Py');
end

%%-------------------------------------------------------------------------
function [Px_y, xx, yy] = cond_prob(x, y, nbx, nby)
% COND_PROB     Determins the conditional Probability P(x|y) based on 
%               distributed sets of scalar variables x and y
%
% CALL: [Px_y, xx, yy] = cond_prob(x, y, nbx, nby);
% INPUTS:
% x & y - vectors (scalar time series)
% nbx   - number of bins for x
% nby   - number of bins for y
% Output: Px_y - (nbx,nby) matrix 
% It prolts results if no output is requested
%
% Copyright and Revised by David Chelidze 7-31-2000

if (nargin == 0) || (nargin == 1)
   error('Requires at least two input arguments')
end
if length(x) ~= length(y) % this can be changed for other applications
   error('Both arrays should be the same size')
end
if nargin == 2
   nbx = 16; nby = 16;
elseif nargin == 3
   nby = nbx;
end
% initialize the conditional probability density matrix
Px_y = zeros(nbx, nby);
% put inputs into column form
y = y(:); x = x(:);
% get bin sizes
minx = min(x); maxx = max(x);
binx = (maxx - minx)/nbx;
miny = min(y); maxy = max(y);
biny = (maxy - miny)/nby;
% get bin boundaries (edges)
xx = minx + binx*(0:nbx);
xx(1) = -inf; xx(end) = inf;
yy = miny + biny*(0:nby);
yy(1) = -inf; yy(end) = inf;
for j=1:nby % calculate!
   iy = find((y > yy(j)) & (y <= yy(j+1))); % # points in bin # j
   if ~isempty(iy)
      Hx_y = histc(x(iy), xx);
      Hx_y = Hx_y(:);
      Px_y(:,j) = Hx_y(1:end-1)/length(iy);
   end
end
end

end