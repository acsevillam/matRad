function d = Dice3D(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to compute Sorensen-Dice similarity coefficient between binary
% 3D images A and B. Values are in the range [0,1]
% 
% call
%   d = Dice3D(A,B)
%
% input
%   A:      binary image in dimension 3 
%   B:      binary image in dimension 3        
%
% output                                                                                    
%   d:      similirity value                    
%
% References
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % interception size
%     intercepSize = length(find(A+B==2)); % norma de A interceptado B
      intercepSize = nnz(A & B);
%     % image A size
%     imagAsize = length(find(A==1)); 
      imagAsize = nnz(A);
%     % image B size
%     imagBsize = length(find(B==1));
      imagBsize = nnz(B);
    
    d=(2*intercepSize)/(imagAsize + imagBsize); 
%       d = 2 * nnz((A) & (B)) / nnz(A) + nnz(B);
% If=A(:);
% Im=B(:);
% c=length(Im);
% inter=0;
% for i=1:c
%     if If(i) == Im(i)
%         inter=inter+1;
%     end
% end 
% d = inter / c;
% d = 2* nnz((A>0.005) & (B>0.005))/(nnz(A>0.005)+nnz(B>0.005));
end
