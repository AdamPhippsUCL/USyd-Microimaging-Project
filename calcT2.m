function [T2, S0] = calcT2(vol, TEVec)
% calcT2 Calculate T2 by linear fitting to log of data (multiple TE
% values)
%
% [ADC, S0] = calcT2(vol, TEVec)
%
% vol [ny nx nz nTE]  or  [ny nx nTE]
% TEVec should have length nTE
%
% Solves the following for S0 and ADC using a least squares fit 
% to the log of the data:
%   S=S0 exp(-TE/T2)
% i.e. 
% ln(S) = ln(S0) - TE/T2
% ln(S) = TE.(-1/T2) + ln(S0)
% Y = mX + c
% Y = [X 1] [m c]'
%
% 
% Examples
% dinfo = datparse ;
% [volb, matp] = d2mat(dinfo,{'slice','bv'},'op','fp') ;
% [ADC, S0] = calcADC(volb, matp.bvVec) ;
%  displayADC(ADC*1e6, 2000)
%
% dinfo = dfparse(ffn) ;
% [vb0, mb0] = d2mat(dinfo,{'slice','bv'},'bv',0,'op','fp%) ;
% bvs = unique([dinfo.DiffusionBValue]) ;
% [vbiso, mbiso] = d2mat(dinfo,{'slice','bv','ddty'},'ddty',2,'bv',bvs(2:end),'op','fp%) ;
% volb = cat(4, vb0, vbiso);
% [ADC, S0] = calcADC(volb, bvs) ;
%
%
%
% David Atkinson
%
% See also dwi2ADC

TEUse = TEVec;

% Reshape so that single slice is still 4D (inc bvaues)
ndv = ndims(vol) ;
if ndv==3
    vol = reshape(vol,[size(vol,1) size(vol,2) 1 ndv]) ;
end

vol(~isfinite(vol)) = 0 ;
loc = vol < 0 ;
if ~isempty(loc)
    warning('Negative signal in input to calcT2')
    vol(loc) = 0 ;
end


% % Restrict bvalues if requested.
% if nargin == 3
%    [TEIntersect, ia] = intersect(TEVec,TEUse) ;
%    if length(TEUse)~=length(TEIntersect) 
%        warning(['Not all TE-values in TEUse are in data'])
%    end
% 
%    bvVec = bvVec(ia) ;
%    volb = volb(:,:,:,ia) ;
% end
    
[ny nx nz nTE] = size(vol) ;
if nTE ~= length(TEVec) 
    error(['Last dimension of volb must be same as number of echo times'])
end

A = cat(2, -TEVec(:), ones([length(TEVec) 1])) ;

B = log(vol) ;
B = reshape(B,[ny*nx*nz nTE]) ;

Btest = sum(B,2) ;
[row] = find(~isfinite(Btest)) ;
B(row, :) = 0 ;

B = B.' ;

X = A\B ;

rT2 = X(1,:) ;
S0 = X(2,:) ;

rT2(row) = 0 ;
S0(row) = 0 ;

rT2 = reshape(rT2,[ny nx nz]) ;
T2 = 1./rT2;
S0 = exp(reshape(S0,[ny nx nz])) ;


