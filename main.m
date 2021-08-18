% main function of lemonade

close all
clear
clc

load input.mat

%% Call muts

data.pvm = ones(nd,1);
data.pvmadj = ones(nd,1);

for s = 1:nst
    for t = 1:nat
        tic
        disp(['Calling ',sampletypelist{s},' samples with ',altypelist{t}, ' type mutations ...'])
        
        alp = Alpha(s,t);
        bet = Beta(s,t);
        noise = ctrlmeanvaf(s,t);
        
        for i = 1:ns
            if strcmp(samplelist.sampleType{i}, sampletypelist{s})
                ix = strcmp(data.sampletype, sampletypelist{s}) ...
                    & strcmp(data.altype, altypelist{t}) ...
                    & strcmp(data.sampleid, samplelist.sampleID{i});
                G1 = data(ix,:);

                data.pvm(ix) = betacdf(noise, alp + G1.altreads, bet + G1.depths - G1.altreads);
            end
        end
        
        % B-H correction by fdr.m, developed by Anderson M. Winkler
        ist = strcmp(data.sampletype, sampletypelist{s}) & strcmp(data.altype, altypelist{t});
        [~,~,data.pvmadj(ist)] = fdr(data.pvm(ist));
        
        toc
    end
end

%% Viz

% Set adj pval cutoff
adjpvalcut = 0.2;

krashots = data(strcmp(data.gene,'KRAS') & data.ishot == 1,:);
nk = size(genehots,1);

vaf0 = zeros(ns,nk);
pv0 = zeros(ns,nk);
prank0 = zeros(ns,nk);
call0 = zeros(ns,nk);
maxvaf = zeros(ns,1);

for i = 1:ns
	vaf0(i,:) = krashots.vaf(strcmp(krashots.sampleid, samplelist.sampleID{i}))';
    pv0(i,:) = krashots.pvmadj(strcmp(krashots.sampleid, samplelist.sampleID{i}))';
    call0(i,pv0(i,:) < adjpvalcut) = 1; % call positive
    if nnz(call0(i,:)) > 0
        vc = vaf0(i,call0(i,:) == 1);
        maxvaf(i) = max(vc);
    end
end

% Fresh
isfresh = strcmp(samplelist.sampleType,'AVM');
ns = nnz(isfresh);
call1 = call0(isfresh,:);
maxvaf1 = maxvaf(isfresh);
vafcut = 10.^(0:(-0.1):-4);
nv = length(vafcut);
nkras1 = zeros(1,nv);

for i = 1:nv
    nkras1(i) = nnz(sum(call1,2) > 0 & maxvaf1 >= vafcut(i))/ns;
end

figure('Position',[0 0 500 500])
hold on
plot(vafcut, nkras1, 'o-', 'linewidth',2, 'color', 'r')
text(vafcut(nv), nkras1(nv) + 0.05, [num2str(round(nkras1(nv)*10000)/100),'%'], ...
    'HorizontalAlignment','right', 'fontsize',20)

% Blood
isblood = strcmp(samplelist.sampleType,'BLO');
ns = nnz(isblood);
call2 = call0(isblood,:);
maxvaf2 = maxvaf(isblood);
nkras2 = zeros(1,nv);

for i = 1:nv
    nkras2(i) = nnz(sum(call2,2) > 0 & maxvaf2 >= vafcut(i))/ns;
end

plot(vafcut, nkras2, 'o-', 'linewidth',2, 'color', 'k')
text(vafcut(end), nkras2(nv) + 0.05, [num2str(round(nkras2(nv)*10000)/100),'%'], ...
    'HorizontalAlignment','right', 'fontsize',20)

legend({'Fresh AVM','Blood'}, 'Location', 'North','Orientation','horizontal','box','off')
set(gca,'tickdir','out','TickLength',[0.015 0.015],'fontsize',18,'box','off',...
    'Xscale','log', 'Xdir', 'reverse','XGrid','on','XMinorTick','on','YGrid','on','linewidth',1.5)
xlabel('VAF Cutoff')
ylabel('Fraction of KRAS Mutant Samples')
ylim([0 1])
axis square
hold off
