function [scalo,f,T,a,wt,wavescaled] = contwtgnmir(x,fmin,fmax,N,wave);

% [scalo,f,T,a,wt,wavescaled] = contwtgnmir(x,fmin,fmax,N,wave);
% Continuous wavelet transform of mirrored 1-D signals
% If x = [a b c e f] is the signal to analyzed, contwtmir.m runs contwt.m
% on the mirrored version XxX = [c b [a  b  c d e f] e d]. The number of
% mirrored samples depends on the analyzed scale and the wavelet length.
% See contwt.m for Inputs/Outputs arguments.
%
% USE AN ORDER "wave" DERIVATIVE OF THE GAUSSIAN

% CHECK INPUT FORMATS

[xr,xc] = size(x) ;
if xr ~= 1 & xc ~= 1
  error('1-D signals only')
elseif xc == 1
  x = x.' ;
end

T = 1 : lenfth(x) ;

% DEFAULT VALUES

nt = size(x,2) ;
if exist('wave') == 0 
  wave = 2 ;
end

if nargin == 1
  XTF = fft(fftshift(x)) ;
  sp = (abs(XTF(1:nt/2))).^2 ;
  f = linspace(0,0.5,nt/2+1) ; f = f(1:nt/2) ;
  plot(f,sp) ; grid ;
  xlabel('frequency');
  title('Analyzed Signal Spectrum') ;
  fmin = input('lower frequency bound = ') ;
  fmax = input('upper frequency bound = ') ;
  N = input('Frequency samples = ') ;
  fmin_s = num2str(fmin) ; fmax_s = num2str(fmax) ; 
  N_s = num2str(N) ;
  disp(['frequency runs from ',fmin_s,' to ',fmax_s,' over ',N_s,' points']) ;
end
if nargin == 4 | nargin == 5
  if fmin >= fmax
    error('fmax must be greater than fmin') ;
  end
end

f = logspace(log10(fmax),log10(fmin),N) ;
a = logspace(log10(1),log10(fmax/fmin),N) ; amax = max(a) ;

for ptr = 1:N
  ha = gaussn(f(ptr),wave) ; nha = (length(ha)-1)/2 ;
  nbmir = min(nt,nha) ;
  x_mir = [x(nbmir:-1:2) x x(nt-1:-1:nt-nbmir+1)] ;
  detail = conv(x_mir,ha) ;
  wt(ptr,1:nt) = detail(nha + nbmir  : nha + nbmir + nt -1 ) ;
end  
wavescaled = wave ;

%%%% pour etre compatible avec le format de donnees de TFTB %%%%%

wt = flipud(wt) ;
a = flipud(a(:)) ;
f = flipud(f(:)) ;
scalo = (wt.*conj(wt)) ;


%%%% pour etre compatible avec le format de donnees de TFTB %%%%%

%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
