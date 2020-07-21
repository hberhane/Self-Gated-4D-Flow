clear;
clc;
close all
% Pseudospiral for 4D flow


%% Phyllotaxis
S   = [96, 96]; %Size of the square grid
TR  = 1500; % Nominal number of samples in each frame
N   = TR*1; % number of samples

d   = 1.5; % variable densisty, d=0.5 uniform, d>0.5 means more samples in the middle
as  = 0.25; % as>0, larger values make the ACS larger along the narrower PE dimension
E   = 1; % Number of encodings 1<E<5
m   = 1*(prod(S)/4)^(1/4); % Nonlinear scaling towards the center to avoid excessive sampling towards the center.
k   = 1/4; % radial shift; k=n implies central n pixels will not be sampled


% gr = (sqrt(5) + 1)/2; % golden ratio for radial increment
% gr = 2-gr;
ga = (sqrt(5) + 1)/2; % golden ratio for radial increment
ga = 2-ga;

% ga = 35^(1/3); % ratio for angular increment; other options: 35, 16, 13
% ga = rem(ga,1);
% gr = 35^(1/3); % ratio for angular increment; other options: 35, 16, 13
gr = sqrt(437);
% gr = rem(gr,1);



figure;
PEInd = zeros(N,2,E);
samp = zeros(S(1), S(2), N/TR, E);
for e=1:E 
    kk=e;
    R  = zeros(N,1);
%     RS = zeros(N,1);
    T  = zeros(N,1); 
    x  = zeros(N,1);%rep(0,N);
    y  = zeros(N,1);%rep(0,N);
    Ri = floor((max(S)-1)/2);

    for n=1:N
        theta=2*pi*ga*(n+(e-1)/E);
    %     theta = theta+fa*e;%(e-1)/E;
        if n==1, r = R(n); 
        else,    r  = rem((R(n-1) + Ri*gr)-0, Ri) + 0;
        end
        R(n) = r;

        x(n) = ((r+m)^(d) - m^(d))                                 *  cos(theta);
        y(n) = ((r+m)^(d*(S(2)/S(1))^as) - m^(d*(S(2)/S(1))^as))   *  sin(theta);

        T(n) = theta;
        blah(e,n) = theta;
    end
    mxX = max(abs(x));
    mxY = max(abs(y));
    x = (S(1)/2-(k+1))*x/mxX + k*cos(T);
    y = (S(2)/2-(k+1))*y/mxY + k*sin(T);
    x=round(x);
    y=round(y);
    PEInd(:, 1, e) = x;
    PEInd(:, 2, e) = y;
    for n=1:N
%       ind = sub2ind([S(1),S(2)], x+S(1)/2+1, y+S(2)/2+1);
        samp(PEInd(n,1,e)+S(1)/2+1, PEInd(n,2,e)+S(2)/2+1, ceil(n/TR), e) = samp(PEInd(n,1,e)+S(1)/2+1, PEInd(n,2,e)+S(2)/2+1, ceil(n/TR), e) + 1;
        imagesc(sqrt(sum(samp(:,:,:,e),3))); axis('off','image');
        pause(1/n);
    end
end

%%
% Time average ACS
sampcum = sum(sum(samp,4),3);
figure; imagesc(logical(sampcum)); axis('image');
%
% Movie of the sampling pattern
figure;
for e=1:E
    for i=1:ceil(N/TR)
        subplot(1,E,e); imagesc(sqrt(samp(:,:,i,e)),[0, max(sqrt(samp(:)))]); axis('image','off');
        pause(0.05);
    end
end

