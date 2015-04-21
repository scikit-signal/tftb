% time-frequency toolbox for matlab
% written by F. Auger, P. Flandrin, P. Goncalves, O. Lemoine,
% 1996 - 1999, GdR-PRC ISIS, CNRS, France.
%
% signal generation functions:
%  amexpo1s, amexpo2s, amgauss, amrect, amtriang
%  fmconst, fmhyp, fmlin, fmodany, fmpar, fmpower, fmsin, gdpower
%  altes, anaask, anafsk, anabpsk, anafsk, anapulse, anaqpsk, anasing, anastep, atoms,
%  dopnoise, doppler, klauder, mexhat, 
%  noisecg, noisecu
%
% analysis windows generation function
%  window
%  
% signal modification
%  scale, fmt, ifmt, zak, izak
%
% signal characterization
%  instfreq, ifestar2, loctime, locfreq, sgrpdelay
%
% time-frequency representations
%  tfrgabor, tfrstft,
%  tfrbj, tfrbud, tfrcw, tfrgrd, tfrmh, tfrmhs, tfrmmce, tfrpage, tfrpmh, tfrppage, 
%  tfrpwv, tfrri, tfrridb, tfrridbn, tfrridh, tfrridt, tfrsp, tfrspwv, tfrwv, tfrzam,
%  tfrbert, tfrdfla, tfrscalo, tfrspaw, tfrunter, tfrspbk
%
% ambiguity functions
%  ambifunb, ambifuwb
%
% postprocessing functions
%  friedman, holder, htl, margtfr, midpoint, momftfr, momttfr, plotsid; renyi, ridges,
%  tfrsurf, tfrideal
%
% visualization and backup functions
%  plotifl, tfrparam, tfrqview, tfrview, tfrsave
%
% miscellaneous
%  istfr1, istfr2, isaffine, disprog, imextrac, divider, dwindow, integ, integ2d,
%  kaytth, modulo, odd, 
%
% demos
%  movcw4at, movpwdph, movpwjph, movsc2wv, movsp2wv, movwv2at

% F. Auger, oct 1999.
help tftb
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
