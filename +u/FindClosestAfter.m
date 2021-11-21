function [Y, inds, ds] = FindClosestAfter(x, y, uniquey)
% [Y inds] = FindClosest(x, y, uniquey)
% 
% This function finds the value in vector Y that is closest and greater than
% to every value in X. Values that are identical in x and y are discarded.
%
% uniquey, if set to 1, ensure that every x is associated with it's own y.
% if a y is aligned to more than one value of x, only the lowest value of
% x is returned
%
% Returns
%   Y: a vector where Y(i) is taken from Y and minimizes ds = X(i)-Y(i) for ds<0
%   inds: indices in Y that match above
%   ds: a vector == x-Y
%
% arb 11/10/14

if ~exist('uniquey', 'var'), uniquey = 0; end

y = sort(y);

inds = zeros(size(x));

todelete = [];

for i = 1:length(x)
   
    tmp = find(x(i)-y<0, 1, 'first');
        
    if ~isempty(tmp)
        if uniquey % if we cant reuse ys
            if any(tmp==inds) % if y has been used already
                if x(i)<x(tmp==inds) % if this x comes before that x
                    inds(i) = tmp;
                    inds(tmp==inds) = 0;
                else
                    continue
                end
            end
        end
        inds(i) = tmp; 
    end
    
end

todelete = inds==0;

inds(inds==0) = 1;

ds = x-y(inds);

Y = y(inds); % the values in sorted Y

ds(todelete) = 0; % set these to 0 to be deleted

Y(ds>=0) = NaN;
inds(ds>=0) = NaN;
ds(ds>=0) = NaN;

% while elegant, the solution below isnt great for large vectors
% ds = repmat(x(:), 1, length(y)) - repmat(y(:)', length(x), 1);
% 
% ds(ds>=0) = NaN; % set to NaN those values that come after
% 
% [~, inds] = min(abs(ds), [], 2);
% 
% Y = y(inds);
% 
% ds = x-Y;