function outStruct = mda(inputFile)


%% input required

% accumulation time in hr
accumulationTime_hr = 1;

% false alarm rate per hr
falseAlarmRate_per_hr = 1;

% probability of detection
detectionProbability = 0.95;

% source
source = 'Cs137'; % source
sourceActivity_Bq = 3.7e6; % Bq
sourceStandoffDistance_m = 10; % meters

%% parameters from input

% edge optimization vector
load(inputFile, 'inStruct');
edgeLength_m = inStruct.edgeLength_m;
edgeVelocity_m_per_s = inStruct.edgeVelocity_m_per_s;
edgeIndex = inStruct.edgeIndex;
modeId = inStruct.modeId;

% optimization mode detector info
modeDetConfigFile = 'modeDetConfig.mat';
load(modeDetConfigFile, 'modeDetConfig');
detectorArea_m2 = zeros(size(modeId));
detectorBackgroundRatio = zeros(size(modeId));
detectorId = zeros(size(modeId));
detectorsPerSystem = zeros(size(modeId));
systemsPerMode = zeros(size(modeId));

for i = length(modeId):-1:1
   f = modeId(i) == [modeDetConfig.id];
   detectorName(i) = {modeDetConfig(f).detectorName};
   detector(i) = {sprintf('%s%s', modeDetConfig(f).detectorName, ...
      modeDetConfig(f).detectorDimension)};
   systemsPerMode(i) = modeDetConfig(f).systemsPerMode;
   detectorsPerSystem(i) = modeDetConfig(f).detectorsPerSystem;
   detectorArea_m2(i) = modeDetConfig(f).detectorArea_m2;
   detectorBackgroundRatio(i) = modeDetConfig(f).detectorBackgroundRatio;
   detectorId(i) = modeDetConfig(f).detectorId;
end

inStruct.detector = detector';
inStruct.detectorName = detectorName';
inStruct.detectorsPerSystem = detectorsPerSystem;
inStruct.detectorId = detectorId;
inStruct.modeDetConfigFile = modeDetConfigFile;
inStruct.file = inputFile;
inStruct.systemsPerMode = systemsPerMode;
inStruct.detectorArea_m2 = detectorArea_m2;
inStruct.detectorBackgroundRatio = detectorBackgroundRatio;

% Set source photon energy
expectedSource = {'Cs137', 'Co57', 'Co60'};
validSource = validatestring(source, expectedSource);
f = strcmpi(validSource, expectedSource);
expectedSourcePhotonEnergy_eV = [6.61657e5, 1.2206065e5, 1.332492e6];

sourcePhotonEnergy_eV = expectedSourcePhotonEnergy_eV(f);
inStruct.sourcePhotonEnergy_eV = sourcePhotonEnergy_eV;
inStruct.source = source; % source
inStruct.sourceActivity_Bq = sourceActivity_Bq;
inStruct.sourceStandoffDistance_m = sourceStandoffDistance_m;

% Set background
% systematic
sysBkgDistFile = 'SysBkgDist.mat';
load(sysBkgDistFile, 'BkgDist');
sysBkgDist = BkgDist;

% poisson
for i = length(detectorName):-1:1
   bkgDist(i,1) = {sprintf('BM_%s_%s.mat', detectorName{i}, source)};

   if ~exist(bkgDist{i,1}, 'file')
      bkgDist(i,1) = {sprintf('BM_%s.mat', source)};
   end
   
   bM(i,1) = load(bkgDist{i,1}, 'BM');
end

% set background
inStruct.backgroundSystematicErrorFile = sysBkgDistFile;
inStruct.backgroundSystematicErrorDistribution = sysBkgDist;
inStruct.backgroundFile = bkgDist;
inStruct.backgroundModel = bM;

% Set detector effective area
% Effective detector area in m^2
% DetectorEffectiveArea_m2  = DetectorEfficiency * 
%                          DetectorsPerSystem * 
%                          DetectorArea_m2
% e.g. AEff = .46 * 3 * (4*16)*2.54^2/100^2
% detectorEfficiency is from getDetectorEfficiency
detectorEfficiency = getDetectorEfficiency(detector, ...
      sourcePhotonEnergy_eV);
detectorEffectiveArea_m2 = detectorEfficiency.peak' .* ...
   detectorsPerSystem .* systemsPerMode .* detectorArea_m2;

inStruct.detectorEffectiveArea_m2 = detectorEffectiveArea_m2;

% Set air attenuation for the source photon
% AirAtten.txt file contains various attenuation coefficient for air as a
% function of photon energy. Data was acquired from NIST XCOM website,
% http://www.nist.gov/pml/data/xcom/.
% Total attenuation coefficient (7th column in AirAtten.txt in 
% cm^2/g) in 1/m (factor of 100 is to convert cm to m):
airAttenuationFile = 'AirAtten.txt';
airAtt = load(airAttenuationFile);
airAttenuation = [airAtt(:,1) airAtt(:,7)];

% Density of air in g/cm^3 at 20 degC and 76 cm Hg:
airDensity = 1.2e-3;

% air attenuation in 1/m
sourceAirAttenuationCoef_per_m = ...
   interp1(airAttenuation(:,1),  airAttenuation(:,2), ...
   sourcePhotonEnergy_eV * 1e-6) * airDensity * 100;

inStruct.sourceAirAttenuationCoef_per_m = sourceAirAttenuationCoef_per_m;
inStruct.AirDensity_g_per_cm3 = airDensity;
inStruct.AirAttenuationFile = airAttenuationFile;
inStruct.AirAttenuation = airAttenuation;

% false alarm tolerance
fat = falseAlarmRate_per_hr * accumulationTime_hr;
nModes = length(modeId);
eFAP = fat/nModes;

inStruct.accumulationTime_hr = accumulationTime_hr;
inStruct.falseAlarmRate_per_hr = falseAlarmRate_per_hr;
inStruct.falseAlarmTolerance = fat;
inStruct.edgeFalseAlarmProbability = eFAP;

%% Calculate mda for each edge

UEI = unique(edgeIndex); % unique edge index
Act = sourceActivity_Bq;
mu = sourceAirAttenuationCoef_per_m;

% MDA_Bq = zeros(length(UEI),1);
% wmSNR = zeros(length(UEI),1);
% dwmSNR = zeros(length(UEI),1);
% B = zeros(length(UEI),1);
% X = zeros(length(UEI),1);

for i = 1:length(UEI)
% for i = 1:1
   f1 = UEI(i) == edgeIndex;
   UDID = unique(detectorId(f1)); % unique detector Id
%    bb = zeros(length(UDID),1);
%    xx = zeros(length(UDID),1);
%    SNR = zeros(length(UDID),1);
%    dSNR = zeros(length(UDID),1);
   for j = 1:length(UDID)
%    for j = 3:length(UDID)
      f2 = UDID(j) == detectorId;
      f3 = f1 & f2; % unique systems of detectors Id and unique edge index
      L = edgeLength_m(f3);
      V = edgeVelocity_m_per_s(f3);
      D = sourceStandoffDistance_m * ones(size(modeId(f3)));
      % same systems of detectors therefore using detectorEffectiveArea of
      % one of the modes
      Aeff = detectorEffectiveArea_m2(f3);
      % source counts from a mode
      s = SourceModel(D, V, mu, Act, Aeff(1), L/2);
%      [modeId(f3) detectorId(f3) L V D Aeff s]
      ss = sum(s);
      nPass = length(s); 
      tdt_s = sum(L ./ V); % total dwell time in seconds
      % same systems of detectors therefore using background model of
      % one of the modes
      bMf = bM(f3,1);
      bm0 = bMf(1).BM;
      br = detectorsPerSystem(f3) .* systemsPerMode(f3) .* ...
         detectorBackgroundRatio(f3);
      [pdf, bins] = BackgroundModel(bm0.a, bm0.b, bm0.c);
      res = SampleModel(1, bins, pdf);
      b = br(1) * res .* tdt_s * nPass;
      % background factor
      x = interp2(BkgDist.PFA', BkgDist.nPass', BkgDist.BF', eFAP, nPass);
      x(isnan(x)) = 1;
%      x = sum(ComputeBackgroundFactor(nPass, eFAP, BkgDist));
      [SNR(j,i), dSNR(j,i)] = computeSNR(ss, b, x);
      [nPass tdt_s ss b x SNR(j,i) dSNR(j,i)]
   end
      [SNR dSNR]
      f4 = SNR(:,i) > 0;
      [wmSNR(i,1), dwmSNR(i,1)] = ...
         computeWeightedMean(SNR(f4,i), dSNR(f4,i));
      [wmSNR(i,1), dwmSNR(i,1)]
%    [S(i,1), dS(i,1)] = computeWeightedMean(ss(i,f2), sqrt(ss(i,f2)));
%    [B(i,1), dB(i,1)] = computeWeightedMean(b(i,f2), sqrt(b(i,f2)));
%    [X(i,1), dX(i,1)] = computeWeightedMean(x(i,f2), sqrt(x(i,f2)));
end
% [S B X wmSNR, dwmSNR]
% MDA_Bq = AggregateMDA(Act, wmSNR, B, X, detectionProbability, eFAP);
% outStruct = [UEI MDA_Bq/3.7e4]; % MDA in uCi



