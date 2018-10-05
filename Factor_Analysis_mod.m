%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Factor Analysis script modified by Jessica Starke 11.03.2017
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

%for further reference see also Hartmann & W�nnemann 2009: Hydrological changes and Holocene climate variations in NW China,
%inferred from lake sediments of Juyanze palaeolake by factor analyses

% in data.mat are all input data saved
% please see Elements.mat to see which variables are used as input
% area [km2] basin area; mean altitude [m] of basin; 10Be blank corrected
% concentration [at/g*yr]; denudation rate [m/Myr]; TRMM derived
% mean precipitation [mm/yr]; Strahler count [Nr.] of river streamorder;
% mean ksn [m09] mean channel steepness of basin; Grain size [mm] used to
% dertermine cosmogenic nuclides; mean slope [�] of basin; mean local relief [m]
% of 5 km radius

% the sript will calculate a cross-corelation matrix (Corr) with significance
% level (p) and the factor loadings (Loadings)
% please see the script for further variable description and the command
% window after calculation

clear
clc

% Data Input
load data.mat % data table
load Elements.mat % heading of data table
X = data(:, 3:end);

%% Select type of transformation

fprintf('Please use the transformation type 2 because 1 and 3 are not suitable for the test dataset! Thanks!\n');
fprintf('1 = Centered log-ratio-tansformation with geometric mean (warning do not use for this data set!)\n');
fprintf('2 = Z-transformation\n');
fprintf('3 = W-transformation to (0,1) after Manson & Imbrie (1964) \n');

M = input('> ');


[n, m] = size(X);
z = zeros(n, m);
switch M
    case 1
        %Centered log-ratio-tansformation with geometric mean
        geom = geomean(X');

        for i = 1:n
            for j = 1:m
                z(i,j) = log(X(i, j)) - log(geom(i));
            end
        end

        clear geom n m
    case 2
        %Z-transformation
        z = zscore(X);
    case 3
        % W-transformation
        for i = 1:n
            for j = 1:m
                z(i, j) = (X(i, j) - min(X(:, j))) / (max(X(:, j)) - min(X(:, j)));
            end
        end

    otherwise
        fprintf('invalid entry \n');
        return;
end


%% Correlation coefficient matrix with 95% confidence interval

[Corr, p, RLO, RUP] = corrcoef(z);
[m, ~] = size(Corr);
l = 1;
for i = 1:m-1
   for j = i+1:m
    k95(l, 1) = RLO(i, j);
    k95(l, 2) = Corr(i, j);
    k95(l, 3) = RUP(i, j);
    p1(l, 1) = p(i, j);
    l = l+1;
   end
end

A = k95(:, 2) > 0.7 | k95(:, 2) < -0.7;
k = k95(A, :);

p2 = p1 > 0.05;
p1 = p1(p2);

%figure
%boxplot(k')

lab_all = nchoosek(Elements, 2);
labelling = cell(length(lab_all), 1);
for i = 1:length(lab_all)
    labelling(i, 1) = strcat(lab_all(i, 1), '-', lab_all(i, 2));
end

labelling_p = labelling(p2); % selects significance which shows no correlation within the confidence interval
labelling = labelling(A); % selects linear relations wihch correlate well with each other
%hold on
for i = 1:length(labelling)
    text(i-.2, k(i, 2), labelling(i))
end
clear k A k95 n m M lab_all i j labelling l p2

%% Criteria for factoranalysis

%1. Number of elements should be > 50 see Winter et al. (2009) oder Wirtz & Nachtigall (2008)
if length (X) < 50
    disp(' ')
    fprintf('The number of elements is below 50.\n')
else
    disp(' ')
    fprintf('The number of elements is above 50.\n')
end

%2. Proof of relationship
if length (X) / length (Elements) < 3
    disp(' ')
    fprintf('The number of Elemnts and variables could be increased.\n')
else
    disp(' ')
    fprintf('The number of elements and variables is perfect.\n')
end

%3. Test of significance
disp(' ')
fprintf('These variables show insignificant correlation within 95 percent confidence interval.\n')
disp('variables    p greater 0.05')
for i = 1:length(p1)
    fprintf('%-12s%8.4f\n', labelling_p{i}, p1(i));
end

%4. test for multicorrelation
detCorr = det(Corr);

if detCorr > 0.00001
    disp(' ')
    fprintf('det(Corr) = %8.6f consequently greater than 0,00001 \n-> the data set shows no multicorrelation\n', detCorr)
else
    disp(' ')
    fprintf('det(Corr) = %8.6f consequently smaller than 0,00001 \n-> the dataset indicates multicorrelation\n', detCorr)
end

%5.structure of correlation
% non diogonal elements differ from 0 with a treshold of 2
Corr_inv = inv(Corr);
Corr_inv_diag = diag(Corr_inv);
Corr_inv_diag_mean = mean(Corr_inv_diag);
Corr_inv_notdiag_mean = sum(sum(Corr_inv - diag(diag(Corr_inv)))) / (length(Corr_inv) * (length(Corr_inv)-1));

% here you can adjust threshold
threshold = 2.0;
fprintf('threshold: %6.2f', threshold)

if abs(Corr_inv_notdiag_mean) < threshold
    disp(' ')
    disp('the mean of the non diagonal elements is smaller than the threshold, the input data matrix is suitable for the factor analysis')
    disp('               meannondiag-meandiag   mean       max       min')
    fprintf('                  %6.2f            %6.2f      %6.2f     %6.2f\n', ...
      Corr_inv_notdiag_mean - Corr_inv_diag_mean, mean(Corr_inv_notdiag_mean), max(max(Corr_inv - diag(diag(Corr_inv)))), min(min(Corr_inv - diag(diag(Corr_inv)))));
else
    disp(' ')
    disp('the mean of the non diagonal elements is greater than the threshold, the input data matrix is not suitable for the factor analysis')
    disp('              meannondiag-meandiag     mean       max       min')
    fprintf('                  %6.2f            %6.2f      %6.2f     %6.2f\n', ...
      Corr_inv_notdiag_mean - Corr_inv_diag_mean, mean(Corr_inv_notdiag_mean), max(max(Corr_inv - diag(diag(Corr_inv)))), min(min(Corr_inv - diag(diag(Corr_inv)))));
end

% 6. Kaiser-Meyer-Olkin test (KMO)

inCorr = inv(Corr);
inCorr_diag = diag(diag((inCorr.^-1)));

AICov = inCorr_diag * inCorr * inCorr_diag;
ICov = Corr + AICov - 2 * inCorr_diag;

Dai = diag(diag(sqrt(AICov)));
ICorr = inv(Dai) * ICov * inv(Dai);
AICor = inv(Dai) * AICov * inv(Dai);

a = sum((AICor - diag(diag(AICor))).^2);
AA = sum(a);
b = sum((Corr - eye(size(Corr))).^2);
BB = sum(b);

MSA = b ./ (b + a); %Measures of sampling adequacy
AICor = AICor - eye(size(AICor)) + diag(MSA);

%7. Examine the anti-image of the correlation matrix. That is the negative of the partial correlations,
%partialling out all other variables.
kmo = BB / (AA + BB);
disp(' ')

fprintf('the Kaiser-Meyer-Olkin (KMO) value is: %3.4f\n', kmo);
if (kmo >= 0.00 && kmo < 0.50);
    disp('the KMO is not acceptable please reconcider your input matrix.')
elseif (kmo >= 0.50 && kmo < 0.60);
    disp('the KMO value is bad')
elseif (kmo >= 0.60 && kmo < 0.70);
    disp('the KMO value could be better')
elseif (kmo >= 0.70 && kmo < 0.80);
    disp('the KMO value is ok')
elseif (kmo >= 0.80 && kmo < 0.90);
    disp('the KMO value is great')
elseif (kmo >= 0.90 && kmo <= 1.00);
    disp('the KMO value is perfect')
end

%8. 25% of non diagonal elements of the anti image correlation is > 0

AICov_d = AICov - diag(diag(AICov));
m = sum(sum(AICov_d > 0.09));

anti_image_covariance = m / length(Corr)^2;

if anti_image_covariance > 0.25
	disp(' ')
    fprintf('%u %% of the non diagonal values of the anti image covarianz are grater than 0.25\n', anti_image_covariance * 100.0);
    fprintf('-> the input data matrix is not suitable for the factor analysis\n');
	disp(' ')
else
    disp(' ')
    fprintf('%u %% of the non diagonal values of the anti image covarianz are between 0 and 0.25\n', anti_image_covariance * 100.0);
    fprintf('\n-> the input data matrix is suitable for the factor analysis\n');
	disp(' ')
end


%% Graphical Interpretion
WM_Z = acos(Corr); %Winkelmatrix in RAD
WM_ZG = WM_Z / pi * 180;%Winkelmatrix in GRAD
r = 1; x = 0; y = 0;% r=length of vectors (1 only if Z-Transformiert), x,y = source
m = length(WM_Z);

%figure if extra interpretation is needed
%for j=1:length(Corr)
    %subplot(ceil(length(Corr)/2), 2,j)
    %for i=1:m
        %u(1,i) = r * cos(WM_Z(j,i)); % convert polar (theta,r) to cartesian in x direction
        %v(1,i) = r * sin(WM_Z(j,i));%in y direction
        %h = quiver(x,y,u(i),v(i),'AutoScale','off');
        %axis([-1 1.2 -0.2 1.2])
        %hold on
    %end
%end

%legend(Elements,'Location','BestOutside')
%hold off

clear r x y Chi2 df m u h v j i WM_Z
%% Calculation of loadings and specific variance (Loadings) and (specVar)

% first solution plot
figure
E = sort(eig(Corr), 'descend');
hold on
plot(E, '-o');

for n = 1:10
    r1 = rand(size(X));
    Rand_Eigen = flipud(sort(eig(corrcoef(r1))));
    plot (Rand_Eigen, 'Color', [0.8 0.8 0.8])
end

% WK: Why question ?
title('Use eigen values above red line','FontSize', 20);
legend('eigen values','random eigen values')
axis([0.9 length(Corr) 0 ceil(max(E))])
line([1, length(Corr)], [1, 1], 'Color', 'r')
for i = 1:length(E)
    text(i+0.2, E(i)+.2, num2str(E(i) * 100 / length(E), '%.2f'));
end

hold off

% calculation of maximum factors
V = size(X, 2);
for i = 1:V
    df = 0.5 * ((V - i)^2 - (V + i));
    if (i > V) || (df < 0)
        break
    end
end
K = i - 1;
clear i df V r1 lambda Rand_Eigen


% Factor analysis matrix rotation
%% Apply varimax rotation see Matlab instruction for further reference

Rot = 'varimax';

m = 1;
for i = 1:K
    [Loadings(:, m+1:m+i), specVar(:, i), T(1:i, m+1:m+i), stats, F(:, m+1:m+i)] = factoran(z, i, 'rotate', Rot);
    eval(['stats_' num2str(i) '=stats;']);
    m = size(Loadings) + 1;
    m = m(2);
end

%% Eigen values (Eigen), declared variance (eVarEig), cummlativ declared variance (KeVarEig), communality (Komm)

[n, m] = size(Loadings);
Eigen = zeros(1, m);
for i = 1:m
    Eigen(1,i) = sum(Loadings(:, i).^2);
end
eVarEig = Eigen * 100 / n;


n = 1;
KeVarEig = zeros(1, K);
for i = 1:K;
    KeVarEig(1,i) = sum(eVarEig(1, n+i:n+i + i-1));
    n = i + n;
end

Komm = 1 - specVar;

for i = 1:length(E)
    if E(i) < 1
        return;
    end
end


%End of factor analysis
%for further help see Starke et al., (hopefully 2017)

disp('The above Warning is due to the selected parameters of Grain Size and Strahler Count. These parameters do not have a unique variance. Consequently this error is expected and minor ^^')
disp('Please continue to work with the variable Corr (Cross-Correlation Matrix) and Loadings (factor loadings indicator of covariance).')
% Warning: Some unique variances are zero: cannot compute significance. is
% not a heavy error due to the parameters selected like Strahler count and grain size.
% Their variance is zero.
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
