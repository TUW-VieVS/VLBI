function [ yInt ] = lagint9der(x,y, xInt)
    % make one check
    if ~isequal(size(x), size(y))
        fprintf('sizes do not fit!\nERROR in function ''lagint9der.m''\n');
        keyboard;
    end
    % get size of vectors
    %ySize=size(y);
    xIntSize=size(xInt);
    
    % preallocating output vector (column vector)
    %yInt=zeros(xIntSize);
    %l = zeros(x);
    yInt = 0;

    for j=1:length(x)
        q = 0;
        %l = 0;
        for i=1:length(x)
            s=1;
            if(j ~= i)
                s = 1/(x(j)-x(i));                
                p = 1;
                for m=1:length(x)   
                    if(m ~= i && m~= j)
                        p = p * (xInt - x(m))/(x(j)-x(m));
                    end
                end
                q = q + s*p;
            end
        end
        yInt = yInt + y(j)*q;
    end
end