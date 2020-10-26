function [p95_4, p68_2, calprob, medage] = biocal(Fdet, sigdet, calcurve, yeartype, sar, bd, thick, brok, abu, res)
%[p95_4, p68_2, calprob, medage] = biocal(labdet, laberr, calcurve, yeartype, sar, bd, thick, brok, abu, res)
%
% Calibrate a 14C determination using sediment bioturbation priors.
%
% B.C. Lougheed, September 2020
% bryan.lougheed@geo.uu.se
%
% Input parameters
% ========================
% labdet:   laboratory 14C determination in 14C yr BP
%
% laberr:   laboratory 14C measurement uncertainty in 14C yr
%
% calcurve: String specifying calibration curve to use, select from
%           the following (not case sensitive):
%           'IntCal20', 'Marine20', 'SHCal20', 'IntCal13', 'Marine13',
%           'SHCal13, 'IntCal09', 'Marine09', 'IntCal04', 'Marine04',
%           'SHCal04, 'IntCal98', 'Marine98'
%
% yeartype: String specifying how to handle all calendar years:
%           'Cal BP' or 'BCE/CE'
%
% sar:      Sediment accumulation rate (SAR) prior in cm/ka
%
% bd:       Bioturbation depth prior in cm.
%
% thick:    Sediment slice thickness in cm.
%
% brok:     Estimation of fraction of foraminifera population that 
%           that is fragmented, and therefore not picked.
%           Value between 0 and 1.
%
% abu:      Temporal abundance of the measured species. Two possibillities.
%           First possibility:
%           Empty matrix, i.e. [], for constant abundance.
%           Second possiblity:
%           n by 2 matrix, i.e. [age abu] for temporally dynamic abundance.
%           Each row contains an age (in years) in column 1 and relative
%           abundance normalised to between 0 and 1 in column 2.
%           Interpolation will be carried out to assign abundance to intermediate
%           ages. Constant extrapolation of the edge values will be carried out to 
%           assign abudance to ages outside the range provided.
%           Ages (in years) should be entered either as Cal BP or BCE/CE,
%           depending on what is set for yeartype.
%
% res:      Reservoir effect. Three possible types of input.
%           First possibility:
%           Empty matrix, i.e. [], for no reservoir effect.
%           Second possibility:
%           1 by 2 matrix, i.e. [reseff reserr], for a temporally constant
%           reservoir effect (reseff) of reseff ± reserr. Reservoir effect
%           to be entered as 14C yr offset from the desired calibration
%           curve.
%           Third possibility:
%           n by 3 matrix, i.e. [age reseff reserr], for incorporating a temporally
%           dynamic reservoir effect. Each row contains an age (in years) in 
%           column 1 and reservoir effect in column 2 and the reservoir effect
%           uncertainty in column 3. Reservoir effect should be entered in 
%           14C yr BP and ages (in years) should be entered either as Cal BP or 
%           BCE/CE, depending on what is set for yeartype.
%
% Output data
% ========================
% 
% p95_4:     n by 3 matrix containing 95.45% calibrated age probability
%            range interval(s) calculated using highest posterior density.
%            Each row contains a probability range in Cols 1 and 2, and
%            the associated probability for that range in Col 3.
%            Probabilities are normalised to between zero and one.
%
% p68_2:     Same as p95_4, but for the 68.27% calibrated range.
%
% calprob:   n by 2 matrix containing an annualised calibrated age
%            probability density function for implementation in, e.g.,
%            age modelling. Col 1 is a series of annual cal ages,
%            Col 2 contains their associated probability. All probabilities
%            are normalised such that they sum to 1.
%
% medage:    Median age calculated from calprob.

%---- manual function input for testing
% labdet = 2500;
% laberr = 50;
% calcurve = 'IntCal20';
% sar = 5;
% bd = 10;
% brok = 0.1;
% yeartype = 'calbp';
% res = [400 100];
% abu = [];

% Cal curve case and symbols
headerlines = 11;
if strcmpi(calcurve, 'IntCal20') == 1
	calcurve = 'IntCal20';
	cite = '(Reimer et al., 2020)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine20') == 1
	calcurve = 'Marine20';
	cite = '(Heaton et al., 2020)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'SHCal20') == 1
	calcurve = 'SHCal20';
	cite = '(Hogg et al., 2020)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'IntCal13') == 1
	calcurve = 'IntCal13';
	cite = '(Reimer et al., 2013)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine13') == 1
	calcurve = 'Marine13';
	cite = '(Reimer et al., 2013)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'SHCal13') == 1
	calcurve = 'SHCal13';
	cite = '(Hogg et al., 2013)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'IntCal09') == 1
	calcurve = 'IntCal09';
	cite = '(Reimer et al., 2009)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine09') == 1
	calcurve = 'Marine09';
	cite = '(Reimer et al., 2009)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'IntCal04') == 1
	calcurve = 'IntCal04';
	cite = '(Reimer et al., 2004)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine04') == 1
	calcurve = 'Marine04';
	cite = '(Hughen et al., 2004)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'SHCal04') == 1
	calcurve = 'SHCal04';
	cite = '(McCormac et al., 2004)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'IntCal98') == 1
	headerlines = 18;
	calcurve = 'IntCal98';
	cite = '(Stuiver et al., 1998)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine98') == 1
	headerlines = 18;
	calcurve = 'Marine98';
	cite = '(Stuiver et al., 1998)';
	curvetype = 'mar';
else
	error(['Calibration curve "',calcurve,'" unknown. Please specify a valid calibration curve (see help for options)'])
end

% Process some of the the user inputs
sar = sar/1000; % sar from cm/ka to cm/a

% Load cal curve data into worksapce
File = fopen([calcurve,'.14c']);
if File == -1
	error(['Could not find ',calcurve,'.14c in working directory or search paths.'])
end
Contents = textscan(File,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(File);
curvecal = flipud(Contents{1});
curve14c = flipud(Contents{2});
curveerr = flipud(Contents{3});

%process reservoir effect
if isempty(res) % no reservoir effect
	% do nothing
elseif size(res) == [1 2] % constant reservoir effect. Still in 14C yr, so add to cal curve
	curve14c = curve14c + res(1);
	curveerr = sqrt(curveerr.^2 + res(2)^2);
elseif size(res,2) == 3 % temporally dynamic reservoir effect
	try
		if strcmpi(yeartype,'BCE/CE') == 1
			res(:,1) = (1950-res(:,1));
		end
		res = sort(res,1);
		reseff = interp1(res(:,1),res(:,2),curvecal,'linear');
		reseff(curvecal<res(1,1)) = res(1,2); % extrapolate using constant value
		reseff(curvecal>res(end,1)) = res(end,2); % extrapolate using constant value
		curve14c = curve14c + reseff;
		reserr = interp1(res(:,1),res(:,3),curvecal,'linear');
		reserr(curvecal<res(1,1)) = res(1,3); % extrapolate using constant value
		reserr(curvecal>res(end,1)) = res(end,3); % extrapolate using constant value
		curveerr = sqrt(curveerr.^2 + reserr.^2);
	catch
		error('Reservoir effect not entered correctly')
	end
else
	error('Reservoir effect not entered correctly')
end

% linearly interpolate cal curve to 1 cal year resolution
calres = 1; 
curve14c = interp1(curvecal,curve14c,min(curvecal):calres:max(curvecal));
curveerr = interp1(curvecal,curveerr,min(curvecal):calres:max(curvecal));
curvecal = min(curvecal):calres:max(curvecal);
curvef14 = exp(curve14c/-8033);	% Now convert to F14C activity
curveerr = curvef14.*curveerr/8033; % Now convert to F14C activity error

% make prior distribution (ppri) of relative age for the sample
rnga = 0:calres:round((bd*4/sar)); % discrete relative age range of prior (4 bioturbation depths)
ppri = exp(-(rnga.*sar)./bd); % Eq. 1 in manuscript
rk = ( -bd*log(brok) ) / sar; % Eq. 3 in manuscript
ppri = ppri(rnga<rk); % trim prior distribution to whole microfossils only
ppri = ppri-(min(ppri)); % floor the prior distribution to zero
ppri = ppri/sum(ppri); % normalise prior dist such that it sums to 1

% start (ts) and end (te) of sliding windows, based on C14 errors and ppri length
% get estimate (in the future, when computers get faster and have more memory, we can just do
% the entire curve and not need to trim here)
[~,ts] = min(abs( curve14c - (Fdet-(3*sigdet)-numel(ppri)*0.75) ));
ts = ts(1);
[~,te] = min(abs( curve14c - (Fdet+(3*sigdet)+numel(ppri)*0.75) ));
te = te(end);

% initiate prob matrix pmat for all windows t, with space for ppri tail at final t
pmat = zeros(numel(ts:te),numel(ts:te)+numel(ppri)-1); %

% trim cal curve to match prob matrix
ind = find(curvecal >= curvecal(ts) & curvecal <= curvecal(ts)+size(pmat,2)-1);
curvef14 = curvef14(ind);
%curve14c = curve14c(ind); % not used
curveerr = curveerr(ind);
curvecal = curvecal(ind);

% retrim pmat to calcurve (i.e. in case it exceeds end of cal curve)
pmat = pmat(:,1:numel(curvecal));

% populate pmat with ppri. each row (t) of pmat contains ppri placed at a new sliding window starting at t
% where ppri(1) is placed at t
S = size(pmat,2);
N = numel(ppri);
for t = 1:size(pmat,1) % could perhaps be vectorised, but super fast as it is
	if S-t+1 >= N
		pmat(t,t:t+N-1) = ppri;
	else % nearing end of calcurve (and thus pmat)
		pmat(t,t:end) = ppri(1:S-t+1);		
	end
end
pmat = pmat./sum(pmat,2); % normalise all rows of pmat

% apply temporal abundance to all columns of pmat
if isempty(abu)
	abus = 1;
elseif size(abu,2) == 2
	try
		if strcmpi(yeartype,'BCE/CE') == 1
			abu(:,1) = (1950-abu(:,1));
		end
		abus = interp1(abu(:,1),abu(:,2),curvecal);
		abus(curvecal<abu(1,1)) = abu(1,2); % extrapolate using constant value
		abus(curvecal>abu(end,1)) = abu(end,2); % extrapolate using constant value
	catch
		error('Abundance not entered correctly')
	end
else
	error('Abundance not entered correctly')
end

% calculate p14c(T|t) (Eq. 6). Fully vectorised here, each row corresponds to instance of t
p14cTt = sum(   1./(curveerr.*sqrt(2*pi)) .* exp(  -(curvef14-curvef14').^2  ./  (2.*curveerr.^2)  ) , 2 )';

% Equation 7, hdet(t). Once again vectorised, each row of hdet corresponds to an instance of t.
p14cTt = pmat .* p14cTt .* abus; % First part of Eq 7.
p14cTt = p14cTt./sum(p14cTt,2); % normalise all rows. Each row is each instance t
hdet = sum( curvef14 .* p14cTt , 2 ); % the final part of Eq. 7

% probability for each sliding window t based on the closeness of its hdet(t) 
% to labdet as calculated using normal pdf of Fdet and sigdet
Fdet = exp(Fdet/-8033); % first convert Fdet to f14c
sigdet = Fdet*sigdet/8033;
phdet = 1/(sigdet*sqrt(2*pi)) * exp( -(hdet-Fdet).^2 / (2*sigdet^2) );

% Construct final cal age probability distribution
pmat = pmat .* phdet;
pcal = sum(pmat,1);
pcal = pcal/sum(pcal); % normalise
if pcal(1) > 10^-10 || pcal(end) > 10^-10
	warning('Calibration may exceed limits of calibration curve')
end

% output calprob to user
calprob = NaN(numel(curvecal),2);
calprob(:,1) = curvecal;
calprob(:,2) = pcal;

% calculate median cal age
[~, median_ind] = min(abs(cumsum(calprob(:,2))-0.5));
medage = round(mean(calprob(median_ind,1))); % in case more than one

% figure(333)
% plot(calprob(:,1)/1000,calprob(:,2))
% ylabel('p(cal)')
% xlabel('Age (ka)')
% hold on

% find 68.2% and 95.4% cal age credible intervals using highest posterior density (HPD)
% script below is copied from matcal: https://zenodo.org/record/3981591
hpd = calprob(:,1:2);
hpd = sortrows(hpd, 2);
hpd(:,3) = cumsum(hpd(:,2));
% 68.2%
hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)), :);
hpd68_2 = sortrows(hpd68_2,1);
ind1 = find(diff(hpd68_2(:,1)) > 1);
if isempty(ind1)
	p68_2(1,1) = hpd68_2(end,1);
	p68_2(1,2) = hpd68_2(1,1);
	p68_2(1,3) = sum(hpd68_2(1:end,2));
else
	indy1 = NaN(length(ind1)*2,1);
	for i = 1:length(ind1)
		indy1(i*2-1,1) = ind1(i);
		indy1(i*2,1) = ind1(i)+1;
	end
	indy1 = [ 1 ; indy1; length(hpd68_2(:,1)) ];
	p68_2 = NaN(length(2:2:length(indy1)),3);
	for i = 2:2:length(indy1)
		p68_2(i/2,1) = hpd68_2(indy1(i),1);
		p68_2(i/2,2) = hpd68_2(indy1(i-1),1);
		p68_2(i/2,3) = sum(hpd68_2(indy1(i-1):indy1(i),2));
	end
	p68_2 = flipud(p68_2);
end
% 95.4%
hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)), :);
hpd95_4 = sortrows(hpd95_4,1);
ind2 = find(diff(hpd95_4(:,1)) > 1);
if isempty(ind2)
	p95_4(1,1) = hpd95_4(end,1);
	p95_4(1,2) = hpd95_4(1,1);
	p95_4(1,3) = sum(hpd95_4(1:end,2));
else
	indy2 = NaN(length(ind2)*2,1);
	for i = 1:length(ind2)
		indy2(i*2-1,1) = ind2(i);
		indy2(i*2,1) = ind2(i)+1;
	end
	indy2 = [ 1 ; indy2; length(hpd95_4(:,1)) ];
	p95_4 = NaN(length(2:2:length(indy2)),3);
	for i = 2:2:length(indy2)
		p95_4(i/2,1) = hpd95_4(indy2(i),1);
		p95_4(i/2,2) = hpd95_4(indy2(i-1),1);
		p95_4(i/2,3) = sum(hpd95_4(indy2(i-1):indy2(i),2));
	end
	p95_4 = flipud(p95_4);
end

% convert output to BCE/CE if requested
if strcmpi(yeartype,'BCE/CE') == 1
	medage = (medage-1950) * -1;
	calprob(:,1) = (calprob(:,1)-1950) * -1;
	p95_4(:,1:2) = (p95_4(:,1:2)-1950) * -1;
	p68_2(:,1:2) = (p68_2(:,1:2)-1950) * -1;
end




end % end function


