function matVout = interpMat2(matV,matXq,method,extrapval)
%interpMat2 此处显示有关此函数的摘要
%   此处显示详细说明
maxXq = max(round(max(max(matXq(matXq~=Inf)))),size(matV,1));
matVext = zeros(maxXq, size(matV,2));
matVext(1:size(matV,1),:) = matV;
V = vec(matVext);
Xq = vec(matXq + maxXq .* (0:size(matXq,2)-1));
X = 1:size(matV,2)*maxXq;
Vout = interp1(X,V,Xq,method,extrapval);
matVout = reshape(Vout, size(matXq));

    function out = vec(in)
        out = in(:);
    end

end

