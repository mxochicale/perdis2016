function [kd_tree, r] = kd_partition(y, b, c)
% KD_PARTITION  Create a kd-tree and partitioned database for efficiently 
%               finding the nearest neighbors to a point in a 
%               d-dimensional space.
%
% USE: [kd_tree, r] = kd_partition(y, b, c);
%
% INPUT:
%   y: original multivariate data (points arranged columnwise, size(y,1)=d)
%   b: maximum number of distinct points for each bin (default is 100)
%   c: minimum range of point cluster within a final leaf (default is 0)
%
% OUTPUT:
%   kd_tree structure with the following variables:
%       splitdim: dimension used in splitting the node
%       splitval: corresponding cutting point
%       first & last: indices of points in the node
%       left & right: node #s of consequent branches to the current node
%   r: sorted index of points in the original y corresponding to the leafs
%
% to find k-nearest neighbors use kd_search.m
%
% copyrighted (c) and written by David Chelidze, January 28, 2009.
 
% check the inputs
if nargin == 0
    error('Need to input at least the data to partition')
elseif nargin > 3
    error('Too many inputs')
end
 
 
% initializes default variables if needed
if nargin < 2
    b = 100;
end
if nargin < 3
    c = 0;
end
 
[d, n] = size(y); % get the dimension and the number of points in y
 
r = 1:n; % initialize original index of points in y
 
% initializes variables for the number of nodes and the last node
node = 1;
last = 1;
 
% initializes the first node's cut dimension and value in the kd_tree
kd_tree.splitdim = 0;
kd_tree.splitval = 0;
 
% initializes the bounds on the index of all points
kd_tree.first = 1;
kd_tree.last = n;
 
% initializes location of consequent branches in the kd_tree
kd_tree.left = 0;
kd_tree.right = 0;

while node <= last % do until the tree is complete
    
    % specify the index of all the points that are partitioned in this node
    segment = kd_tree.first(node):kd_tree.last(node);
    
    % determines range of data in each dimension and sorts it
    [rng, index] = sort(range(y(:,segment),2));
    
    % now determine if this segment needs splitting (cutting)
    if rng(d) > c && length(segment)>= b % then split
        yt = y(:,segment); 
        rt = r(segment);
        [sorted, sorted_index] = sort(yt(index(d),:));
        % estimate where the cut should go
        lng = size(yt,2);
        cut = (sorted(ceil(lng/2)) + sorted(floor(lng/2+1)))/2;
        L = (sorted <= cut); % points to the left of cut
        if sum(L) == lng % right node is empty
            L = (sorted < cut); % decrease points on the left
            cut = (cut + max(sorted(L)))/2; % adjust the cut
        end
        
        % adjust the order of the data in this node
        y(:,segment) = yt(:,sorted_index); 
        r(segment) = rt(sorted_index);
 
        % assign appropriate split dimension and split value
        kd_tree.splitdim(node) = index(d);
        kd_tree.splitval(node) = cut;
        
        % assign the location of consequent bins and 
        kd_tree.left(node) = last + 1;
        kd_tree.right(node) = last + 2;
        
        % specify which is the last node at this moment
        last = last + 2;
        
        % initialize next nodes cut dimension and value in the kd_tree
        % assuming they are terminal at this point
        kd_tree.splitdim = [kd_tree.splitdim 0 0];
        kd_tree.splitval = [kd_tree.splitval 0 0]; 
 
        % initialize the bounds on the index of the next nodes
        kd_tree.first = [kd_tree.first segment(1) segment(1)+sum(L)];
        kd_tree.last = [kd_tree.last segment(1)+sum(L)-1 segment(lng)];
 
        % initialize location of consequent branches in the kd_tree
        % assuming they are terminal at this point
        kd_tree.left = [kd_tree.left 0 0];
        kd_tree.right = [kd_tree.right 0 0];
        
    end % the splitting process
 
    % increment the node
    node = node + 1;
 
end % the partitioning