function [ canvas ] = paintContour(contour, dims)
%varargin: BW, size, indexes, arrows,

temp = zeros(dims);
[rl, cl] = size(contour);
for i = 1:rl
    pt = contour(i,:);
    row = floor(pt(1));
    col = floor(pt(2));
%     row = max(1,row);
%     row = min(size(1),row);
%     col = max(1,col);
%     col = min(size(2),col);
     row_low = max(1,row-1);
     col_low = max(1,col-1);
     row_high = min(dims(1), row+1);
     col_high = min(dims(2), col+1);
    temp(row_low:row_high, col_low:col_high) = 1;
end
canvas = temp;
end
