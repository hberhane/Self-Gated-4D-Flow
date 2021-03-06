function [ weighted_maxwell_fit ] = maxxCorr2D(xHatb,xHatx )
%MAXXCORR3D Summary of this function goes here
%   Detailed explanation goes here
xHatb = gather(xHatb);
xHatx = gather(xHatx);

b = sum(xHatb,3)/size(xHatb,3);
x = sum(xHatx,3)/size(xHatx,3);

thetaX = angle(x.*conj(b));
mag = abs(x);

% Threshold the DATA to find the support set of proton densities
threshold = 0.03;
threshold =  threshold*max(mag(:));
support = find(mag >= threshold);

% Create a weighting matrix
temporal_preci = 1./var(abs(xHatb),0,3).^0.5;
weights = zeros(size(thetaX));
weights(support) = 1;
weights = weights.*temporal_preci;
weights = weights(:);
W = spdiags(weights,0,size(weights,1),size(weights,1));

% Create a Matrix of Corresponding to a Polynomial of degree N
[y_ind_all, x_ind_all] = meshgrid(1:size(thetaX,2),1:size(thetaX,1));
x_ind_all = x_ind_all(:);
y_ind_all = y_ind_all(:);

degree = 2;
column_ind = 1;
for x_power = 0:degree
    for y_power = 0:degree
        if x_power+y_power<=degree
           A(:,column_ind) = (x_ind_all.^x_power).*(y_ind_all.^y_power);
           column_ind = column_ind +1;
        end
    end
end

%Get the Magnitude of the Raw Maps
y = thetaX(:);

weigthed_coeffs = (A.'*W*A)^-1*(A.'*W)*y;
weighted_maxwell_fit = reshape(A*weigthed_coeffs,size(thetaX));
   
% coeffs = (A.'*A)^-1*(A.')*y;
% maxwell_fit = reshape(A*coeffs,size(thetaX));

end

