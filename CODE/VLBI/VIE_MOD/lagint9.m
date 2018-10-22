% lagrange interpolation (degree 9 - 10 points)
% ======================
% takes ten points (five x smaller, five x bigger than xInt) and computes a
% polynomial of nineth degree. If there are less then the needed four 
% values given (e.g. when time series (x,y) is too small), this function 
% produces an error - no checks are done for performance reason.
%
% INPUT
%   x ... given x values [nx1] or [1xn] (must be distinct values especially
%         when you want to interpolate this value, then the return are two
%         values)
%   y ... corresponding y values (same size as x)
%   xInt. vector of x values to be interpolated [1xm] or [mx1]
%
% OUPUT
%   yInt. corresponding (to xInt) interpolated values 
%
% CREATED
%   30 Jan 2012 by Matthias Madzak

function yInt=lagint9(x,y,xInt)

% make one check
if ~isequal(size(x), size(y))
    fprintf('sizes do not fit!\nERROR in function ''lagint9_matthias.m''\n');
    keyboard;
end
% get size of vectors
%ySize=size(y);
xIntSize=size(xInt);


% preallocating output vector (column vector)
yInt=zeros(xIntSize);

for k=1:max(xIntSize)
        
        % a - only get one index and use that
        %===============
        indAll=find(x>xInt(k)); % 1 in original formula is (ind-2)
        ind=indAll(1);          % this index is the first bigger element
        
        yInt(k)=y(ind-5)*( (xInt(k)-x(ind-4))*(xInt(k)-x(ind-3))*(xInt(k)-x(ind-2))*(xInt(k)-x(ind-1))*(xInt(k)-x(ind))  *(xInt(k)-x(ind+1))*(xInt(k)-x(ind+2))*(xInt(k)-x(ind+3))*(xInt(k)-x(ind+4))) / ( (x(ind-5)-x(ind-4))*(x(ind-5)-x(ind-3))*(x(ind-5)-x(ind-2))*(x(ind-5)-x(ind-1))*(x(ind-5)-x(ind)  )*(x(ind-5)-x(ind+1))*(x(ind-5)-x(ind+2))*(x(ind-5)-x(ind+3))*(x(ind-5)-x(ind+4)) ) +...
            y(ind-4)    *( (xInt(k)-x(ind-5))*(xInt(k)-x(ind-3))*(xInt(k)-x(ind-2))*(xInt(k)-x(ind-1))*(xInt(k)-x(ind))  *(xInt(k)-x(ind+1))*(xInt(k)-x(ind+2))*(xInt(k)-x(ind+3))*(xInt(k)-x(ind+4))) / ( (x(ind-4)-x(ind-5))*(x(ind-4)-x(ind-3))*(x(ind-4)-x(ind-2))*(x(ind-4)-x(ind-1))*(x(ind-4)-x(ind)  )*(x(ind-4)-x(ind+1))*(x(ind-4)-x(ind+2))*(x(ind-4)-x(ind+3))*(x(ind-4)-x(ind+4)) ) +...
            y(ind-3)    *( (xInt(k)-x(ind-5))*(xInt(k)-x(ind-4))*(xInt(k)-x(ind-2))*(xInt(k)-x(ind-1))*(xInt(k)-x(ind))  *(xInt(k)-x(ind+1))*(xInt(k)-x(ind+2))*(xInt(k)-x(ind+3))*(xInt(k)-x(ind+4))) / ( (x(ind-3)-x(ind-5))*(x(ind-3)-x(ind-4))*(x(ind-3)-x(ind-2))*(x(ind-3)-x(ind-1))*(x(ind-3)-x(ind)  )*(x(ind-3)-x(ind+1))*(x(ind-3)-x(ind+2))*(x(ind-3)-x(ind+3))*(x(ind-3)-x(ind+4)) ) +...
            y(ind-2)    *( (xInt(k)-x(ind-5))*(xInt(k)-x(ind-4))*(xInt(k)-x(ind-3))*(xInt(k)-x(ind-1))*(xInt(k)-x(ind))  *(xInt(k)-x(ind+1))*(xInt(k)-x(ind+2))*(xInt(k)-x(ind+3))*(xInt(k)-x(ind+4))) / ( (x(ind-2)-x(ind-5))*(x(ind-2)-x(ind-4))*(x(ind-2)-x(ind-3))*(x(ind-2)-x(ind-1))*(x(ind-2)-x(ind)  )*(x(ind-2)-x(ind+1))*(x(ind-2)-x(ind+2))*(x(ind-2)-x(ind+3))*(x(ind-2)-x(ind+4)) ) +...
            y(ind-1)    *( (xInt(k)-x(ind-5))*(xInt(k)-x(ind-4))*(xInt(k)-x(ind-3))*(xInt(k)-x(ind-2))*(xInt(k)-x(ind))  *(xInt(k)-x(ind+1))*(xInt(k)-x(ind+2))*(xInt(k)-x(ind+3))*(xInt(k)-x(ind+4))) / ( (x(ind-1)-x(ind-5))*(x(ind-1)-x(ind-4))*(x(ind-1)-x(ind-3))*(x(ind-1)-x(ind-2))*(x(ind-1)-x(ind)  )*(x(ind-1)-x(ind+1))*(x(ind-1)-x(ind+2))*(x(ind-1)-x(ind+3))*(x(ind-1)-x(ind+4)) ) +...
            y(ind)      *( (xInt(k)-x(ind-5))*(xInt(k)-x(ind-4))*(xInt(k)-x(ind-3))*(xInt(k)-x(ind-2))*(xInt(k)-x(ind-1))*(xInt(k)-x(ind+1))*(xInt(k)-x(ind+2))*(xInt(k)-x(ind+3))*(xInt(k)-x(ind+4))) / ( (x(ind)  -x(ind-5))*(x(ind)  -x(ind-4))*(x(ind)  -x(ind-3))*(x(ind)  -x(ind-2))*(x(ind)  -x(ind-1))*(x(ind)  -x(ind+1))*(x(ind)  -x(ind+2))*(x(ind)  -x(ind+3))*(x(ind)  -x(ind+4)) ) +...
            y(ind+1)    *( (xInt(k)-x(ind-5))*(xInt(k)-x(ind-4))*(xInt(k)-x(ind-3))*(xInt(k)-x(ind-2))*(xInt(k)-x(ind-1))*(xInt(k)-x(ind)  )*(xInt(k)-x(ind+2))*(xInt(k)-x(ind+3))*(xInt(k)-x(ind+4))) / ( (x(ind+1)-x(ind-5))*(x(ind+1)-x(ind-4))*(x(ind+1)-x(ind-3))*(x(ind+1)-x(ind-2))*(x(ind+1)-x(ind-1))*(x(ind+1)-x(ind)  )*(x(ind+1)-x(ind+2))*(x(ind+1)-x(ind+3))*(x(ind+1)-x(ind+4)) ) +...
            y(ind+2)    *( (xInt(k)-x(ind-5))*(xInt(k)-x(ind-4))*(xInt(k)-x(ind-3))*(xInt(k)-x(ind-2))*(xInt(k)-x(ind-1))*(xInt(k)-x(ind)  )*(xInt(k)-x(ind+1))*(xInt(k)-x(ind+3))*(xInt(k)-x(ind+4))) / ( (x(ind+2)-x(ind-5))*(x(ind+2)-x(ind-4))*(x(ind+2)-x(ind-3))*(x(ind+2)-x(ind-2))*(x(ind+2)-x(ind-1))*(x(ind+2)-x(ind)  )*(x(ind+2)-x(ind+1))*(x(ind+2)-x(ind+3))*(x(ind+2)-x(ind+4)) ) +...
            y(ind+3)    *( (xInt(k)-x(ind-5))*(xInt(k)-x(ind-4))*(xInt(k)-x(ind-3))*(xInt(k)-x(ind-2))*(xInt(k)-x(ind-1))*(xInt(k)-x(ind)  )*(xInt(k)-x(ind+1))*(xInt(k)-x(ind+2))*(xInt(k)-x(ind+4))) / ( (x(ind+3)-x(ind-5))*(x(ind+3)-x(ind-4))*(x(ind+3)-x(ind-3))*(x(ind+3)-x(ind-2))*(x(ind+3)-x(ind-1))*(x(ind+3)-x(ind)  )*(x(ind+3)-x(ind+1))*(x(ind+3)-x(ind+2))*(x(ind+3)-x(ind+4)) ) +...
            y(ind+4)    *( (xInt(k)-x(ind-5))*(xInt(k)-x(ind-4))*(xInt(k)-x(ind-3))*(xInt(k)-x(ind-2))*(xInt(k)-x(ind-1))*(xInt(k)-x(ind)  )*(xInt(k)-x(ind+1))*(xInt(k)-x(ind+2))*(xInt(k)-x(ind+3))) / ( (x(ind+4)-x(ind-5))*(x(ind+4)-x(ind-4))*(x(ind+4)-x(ind-3))*(x(ind+4)-x(ind-2))*(x(ind+4)-x(ind-1))*(x(ind+4)-x(ind)  )*(x(ind+4)-x(ind+1))*(x(ind+4)-x(ind+2))*(x(ind+4)-x(ind+3)) );
        

end

