% create a distribution and see if we could successfully calibrate it with prior

File = fopen(['private/Marine20.14c']);
headerlines = 11;
Contents = textscan(File,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(File);
ccal = flipud(Contents{1});
c14c = flipud(Contents{2});
cf14 = exp(c14c/-8033);
cerr = flipud(Contents{3}); 
cerr = cf14.*cerr/8033; % in f14c
cf14 = interp1(ccal,cf14,min(ccal):1:max(ccal));
cerr = interp1(ccal,cerr,min(ccal):1:max(ccal));
ccal = min(ccal):1:max(ccal);




sars = [5  4  6   8  10 12 14 16 18];
bds  = [10 10 10 10  10 10 10 10 10];


[ha, pos] = tightsubplot(3, 3, 0.05, 0.05, 0.05);
delete(ha);

figure(787)
clf

for i = 1:numel(sars)


% sar, bd and foram priors
sarin = (sars(i)); % SAR cm/ka
sar = sarin/1000; % cm/a
bd = bds(i); % BD cm
brok = 0.1; % broken foram fraction
thick = 1;
% make prior distribution of relative age for the sample
% bioturbation pdf with 4 bioturbation depths 
rnga = 0:1:round((bd*4/sar)); % discrete age intervals
rngd = rnga .* sar; % depth domain equivalent of intervals
ppri = exp(-(rngd)./bd); % probabilities for each discrete interval


%error('end')

% get whole foram depth cutoff
brokd = -bd*log(brok); % in depth domain
broka = round(brokd / sar); % in age domain;
% trim pdf to whole forams only
ind = find(rnga<broka);
rnga = rnga(ind);
rngd = rngd(ind);
ppri = ppri(ind);
ppri = ppri - min(ppri);
ppri = ppri/sum(ppri); % normalise such that sums to 1
calrnga = rnga+11000;
rngf14 = cf14(ccal>=min(calrnga) & ccal<=max(calrnga));
meanf14 = sum(rngf14 .* ppri);
mean14c = -8033*log(meanf14);


%---- input variables for function
% lab determination in 14C yr
labdet = mean14c;
laberr = 80;
% chosen cal curve
calcurve = 'Marine20';
% resage and res err
resage = 0;
reserr = 0;
yeartype = 'calbp';

abu = [0    0.5
	   1000 0.5];





hax(i) = axes('position',pos{i});

% ground truth
plot([calrnga(1) calrnga calrnga(end)]/1000,[0 ppri 0]/max(ppri))
hold on

% biocal
tic
% vectorised
[p95_4, p68_2, calprob, medage] = biocal(labdet, laberr, calcurve, yeartype, sarin, bd, thick, brok, abu, []);
disp('vectorised')
toc
plot(calprob(:,1)/1000,calprob(:,2)/max(calprob(:,2)))
hold on



% old-skool
[~, ~, calprob, ~] = matcal(labdet, laberr, calcurve, yeartype, 'plot',0);

plot(calprob(:,1)/1000,calprob(:,2)/max(calprob(:,2)))

xlim([calrnga(1)-1000 calrnga(end)+1000]/1000)

%legend('Ground-truthing distribution','biocal','Traditional 14C calibration')

xlabel('Age (ka)')

hlab(i) = panlabel(hax(i),...
	['SAR: ',num2str(sarin),' cm/ka', newline...
	'BD: ',num2str(bd),' cm', newline...
	'F_{det}: ',num2str(round(labdet)),' \pm',num2str(laberr),' ^1^4C yr'],...
	'top', 'right');

set(hlab(i),'fontsize',6,'fontweight','normal')

set(gca,'Layer','top')


hlet = panlabel(gca, char(96+i), 'bottom', 'left');

end

print2pdf(gcf, '/home/bryan/ble/pubs/customcal/fig1.pdf',15, 10)
