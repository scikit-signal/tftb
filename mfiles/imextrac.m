function [Image2,NbDots]=imextrac(Image,trace);
% imextrac(Image) extract and isolate dots in a binary image
%
% example:
% Image=[1 0 0 0 0 0 1 0 0 0 1 0 ;...
%        1 1 0 0 0 1 1 0 0 0 1 1 ;...
%        1 0 0 0 1 0 1 0 0 0 0 1 ;...
%        0 0 0 0 1 1 1 0 0 0 0 0 ;...
%        0 0 0 0 0 0 0 0 0 0 0 0 ];
% image2=imextrac(Image)
%
%	See also tfrsurf.

% 	F. Auger, oct 1999
%	Copyright (c) CNRS - France 1999. 
%
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
if nargin==1, trace=0; end;

Image2=Image;
[Nbrow,Nbcol]=size(Image2);

NbDots=1 ; 
if trace==1, fprintf('extracting dots in the image.\n'); end;

if (Image2(1,1)==1),
 NbDots=NbDots+1; Image2(1,1)=NbDots;
end;

for i=2:Nbcol,
 if (Image2(1,i)==1),
  if (Image2(1,i-1)>1),
   Image2(1,i)=Image2(1,i-1);
  else
   NbDots=NbDots+1; Image2(1,i)=NbDots;
  end;
 end; 
end;
if trace==1, disprog(1,Nbrow,10); end;

for j=1:Nbrow,
 if (Image2(j,1)==1),
  if (Image2(j-1,1)>1),
   Image2(j,1)=Image2(j-1,1);
  else
   NbDots=NbDots+1; Image2(j,1)=NbDots;
  end;
 end;

 for i=2:Nbcol,
  if (Image2(j,i)==1),
   if (Image2(j-1,i)==0)&(Image2(j,i-1)==0),
    NbDots=NbDots+1; Image2(j,i)=NbDots;
   elseif (Image2(j-1,i)==0)&(Image2(j,i-1)>1),
    Image2(j,i)=Image2(j,i-1);
   elseif (Image2(j-1,i)>1)&(Image2(j,i-1)==0),
    Image2(j,i)=Image2(j-1,i);
   elseif (Image2(j-1,i)>1)&(Image2(j,i-1)>1),
    MinDot=min(Image2(j-1,i),Image2(j,i-1));
    MaxDot=max(Image2(j-1,i),Image2(j,i-1));
    Indices=find(Image2==MaxDot); Image2(Indices)=MinDot;
    Image2(j,i)=MinDot;
   else error('should never happen');
   end;
  end; 
 end;
 if trace==1, disprog(j,Nbrow,10); end;
end;

