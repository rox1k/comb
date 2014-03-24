% Copyright (c) 2013, T. Herr, I. Mirgorodskiy, G. Lihachev, M.L. Gorodetsky
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
% OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function reslist = buildResList(NumberOfModes, CenterFreq, FSR, D2over2pi, D3over2pi)
%returns resonance frequencies in units of Hz according to FSR, D2, D3
global nms_a;
global nms_b;
global linewidth;
list = zeros(NumberOfModes,1);
for k=1:size(list,1)
    list(k) = CenterFreq ...
        + (k-round(size(list,1)/2))*FSR ...
        + (k-round(size(list,1)/2))^2 * D2over2pi ...
        + (k-round(size(list,1)/2))^3 * D3over2pi/3 ...
        + nms_a*linewidth/4/(k-round(size(list,1)/2)-nms_b-0.5);
end

reslist = list;

end