close all
clear all

load postpost/patient
load postpost/control

sig_colour=[.9 .9 .9];
pat_colour=[0.8 0 0.2];
hc_colour=[.1 .1 .4];

%% initial blub
frame_taxis = patient_dat.frame_taxis;
datrangelims_nstdvs = 2;
yes_cull_patients = 1;
yes_cull_controls = 1;
yes_2ndcull_controls = 0;

alphalev = 0.05;

tlims4anal_req = [-0.2 1.0];
ilims4anal = dsearchn(frame_taxis',tlims4anal_req')';
idxes4anal = ilims4anal(1):ilims4anal(2);
tlims4anal = frame_taxis(ilims4anal);
tinc4segs = 0.05;
tlims4segs_req = [(tlims4anal(1):tinc4segs:(tlims4anal(2) - tinc4segs))' ((tlims4anal(1) + tinc4segs):tinc4segs:tlims4anal(2))'];
[nsegs,duh] = size(tlims4segs_req);

for iseg = 1:nsegs
  ilims4segs(iseg,:) = dsearchn(frame_taxis',tlims4segs_req(iseg,:)')';
  idxes4segs{iseg}   = ilims4segs(iseg,1):ilims4segs(iseg,2);
  tlims4segs(iseg,:) = frame_taxis(ilims4segs(iseg,:));
end

good_patient_idxes = find(patient_dat.is_good); ngood_patients = length(good_patient_idxes);
good_control_idxes = find(control_dat.is_good); ngood_controls = length(good_control_idxes);

nframes = patient_dat.pert_resp(good_patient_idxes(1)).nframeswin;

ipert_type = [1 2];
nperts = length(ipert_type);

% calculate good patients
clear patient_perttrial
for ipert = 1:nperts
  itrial_last = 0;
  for isubj = 1:ngood_patients
    the_subj = good_patient_idxes(isubj);
    ngood_trials4subj = patient_dat.pert_resp(the_subj).n_good_trials(ipert_type(ipert));
    patient_perttrial{ipert}.dat((itrial_last+1):(itrial_last+ngood_trials4subj),:) = patient_dat.pert_resp(the_subj).cents4comp.pitch_in.dat{ipert_type(ipert)};
    if ipert == 1
      patient_perttrial{ipert}.absdat((itrial_last+1):(itrial_last+ngood_trials4subj),:) = -(patient_dat.pert_resp(the_subj).cents4comp.pitch_in.dat{ipert_type(ipert)});
    else
      patient_perttrial{ipert}.absdat((itrial_last+1):(itrial_last+ngood_trials4subj),:) = +(patient_dat.pert_resp(the_subj).cents4comp.pitch_in.dat{ipert_type(ipert)});
    end
    patient_perttrial{ipert}.subj((itrial_last+1):(itrial_last+ngood_trials4subj),1) = the_subj*ones([ngood_trials4subj 1]);
    itrial_last = itrial_last + ngood_trials4subj;
  end
  patient_perttrial{ipert}.ntrials = itrial_last;
end

%% figure 1 and 2 - the raw data

nhf = 0;
nhf = nhf + 1; hf(nhf) = figure;
nhf_rawdat_patients = nhf;
for ipert = 1:nperts
  patient_perttrial{ipert}.premean = nanmean(patient_perttrial{ipert}.dat,  1);
  for itrial = 1:patient_perttrial{ipert}.ntrials
    patient_perttrial{ipert}.dev_from_premean(itrial,:) = patient_perttrial{ipert}.dat(itrial,:) - patient_perttrial{ipert}.premean;
  end
  patient_perttrial{ipert}.premeanstdv = mean(nanstd( patient_perttrial{ipert}.dat,0,1));
  
  hax(nhf,ipert) = subplot(2,1,ipert);
  plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.dev_from_premean(:,idxes4anal)');
  kids = get(gca,'Children');
  if yes_cull_patients
    datrangelims = datrangelims_nstdvs*patient_perttrial{ipert}.premeanstdv*[-1 1];
%     hl = hline(datrangelims(1),'k');
%     hl = hline(datrangelims(2),'k');
    in_range_dat = (datrangelims(1) < patient_perttrial{ipert}.dev_from_premean(:,idxes4anal)) & (patient_perttrial{ipert}.dev_from_premean(:,idxes4anal) < datrangelims(2));
    out_of_range_dat = ~in_range_dat;
    out_of_range_trials = any(out_of_range_dat')';
    set(kids(out_of_range_trials(end:(-1):1)),'Color',0.7*[1 1 1]);
    axis([tlims4anal -400 400]);
    
    patient_perttrial{ipert}.dat(out_of_range_dat) = NaN;
    patient_perttrial{ipert}.absdat(out_of_range_dat) = NaN;
  end
end

patient_absperttrial.dat = [patient_perttrial{1}.absdat; patient_perttrial{2}.absdat];
patient_absperttrial.ntrials = size(patient_absperttrial.dat,1);

clear control_perttrial
for ipert = 1:nperts
  itrial_last = 0;
  for isubj = 1:ngood_controls
    the_subj = good_control_idxes(isubj);
    ngood_trials4subj = control_dat.pert_resp(the_subj).n_good_trials(ipert_type(ipert));
    control_perttrial{ipert}.dat((itrial_last+1):(itrial_last+ngood_trials4subj),:) = control_dat.pert_resp(the_subj).cents4comp.pitch_in.dat{ipert_type(ipert)};
    if ipert == 1
      control_perttrial{ipert}.absdat((itrial_last+1):(itrial_last+ngood_trials4subj),:) = -(control_dat.pert_resp(the_subj).cents4comp.pitch_in.dat{ipert_type(ipert)});
    else
      control_perttrial{ipert}.absdat((itrial_last+1):(itrial_last+ngood_trials4subj),:) = +(control_dat.pert_resp(the_subj).cents4comp.pitch_in.dat{ipert_type(ipert)});
    end
    control_perttrial{ipert}.subj((itrial_last+1):(itrial_last+ngood_trials4subj),:) = the_subj*ones([ngood_trials4subj 1]);
    itrial_last = itrial_last + ngood_trials4subj;
  end
  control_perttrial{ipert}.ntrials = itrial_last;
end

  
nhf = nhf + 1; hf(nhf) = figure;
nhf_rawdat_controls = nhf;
for ipert = 1:nperts
  control_perttrial{ipert}.premean = nanmean(control_perttrial{ipert}.dat,  1);
  for itrial = 1:control_perttrial{ipert}.ntrials
    control_perttrial{ipert}.dev_from_premean(itrial,:) = control_perttrial{ipert}.dat(itrial,:) - control_perttrial{ipert}.premean;
  end
  control_perttrial{ipert}.premeanstdv = mean(nanstd( control_perttrial{ipert}.dat,0,1));
  
  hax(nhf,ipert) = subplot(2,1,ipert);
  plot(frame_taxis(idxes4anal),control_perttrial{ipert}.dev_from_premean(:,idxes4anal)');
  kids = get(gca,'Children');
  if yes_cull_controls
    datrangelims = datrangelims_nstdvs*control_perttrial{ipert}.premeanstdv*[-1 1];
%     hl = hline(datrangelims(1),'k');
%     hl = hline(datrangelims(2),'k');
    in_range_dat = (datrangelims(1) < control_perttrial{ipert}.dev_from_premean(:,idxes4anal)) & (control_perttrial{ipert}.dev_from_premean(:,idxes4anal) < datrangelims(2));
    out_of_range_dat = ~in_range_dat;
    out_of_range_trials = any(out_of_range_dat')';
    set(kids(out_of_range_trials(end:(-1):1)),'Color',0.7*[1 1 1]);
    axis([tlims4anal -400 400]);
    
    control_perttrial{ipert}.dat(out_of_range_dat) = NaN;
    control_perttrial{ipert}.absdat(out_of_range_dat) = NaN;
  end
end

if yes_cull_controls && yes_2ndcull_controls
  nhf = nhf + 1; hf(nhf) = figure;
  nhf_2ndcull_controls = nhf;
  for ipert = 1:nperts
    hax(nhf,ipert) = subplot(2,1,ipert);
    plot(frame_taxis(idxes4anal),control_perttrial{ipert}.dat(:,idxes4anal)');
    kids = get(gca,'Children');
    hl = hline(datrangelims(1),'k');
    hl = hline(100,'k');
    in_range_trials = (datrangelims(1) < control_perttrial{ipert}.dat(:,idxes4anal(end))) & (control_perttrial{ipert}.dat(:,idxes4anal(end)) < 100);
    out_of_range_trials = ~in_range_trials;
    set(kids(out_of_range_trials(end:(-1):1)),'Color',0.7*[1 1 1]);
    axis([tlims4anal -400 400]);
    
    control_perttrial{ipert}.dat(out_of_range_trials,:) = [];
    control_perttrial{ipert}.absdat(out_of_range_trials,:) = [];
    control_perttrial{ipert}.subj(out_of_range_trials) = [];
    control_perttrial{ipert}.ntrials = length(control_perttrial{ipert}.subj);
  end
end

control_absperttrial.dat = [control_perttrial{1}.absdat; control_perttrial{2}.absdat];
control_absperttrial.ntrials = size(control_absperttrial.dat,1);

for ipert = 1:nperts
  patient_perttrial{ipert}.mean = nanmean(patient_perttrial{ipert}.dat,1);
  patient_perttrial{ipert}.stde = nanstd(patient_perttrial{ipert}.dat,0,1)/sqrt(size(patient_perttrial{ipert}.dat,1));
  control_perttrial{ipert}.mean = nanmean(control_perttrial{ipert}.dat,1);
  control_perttrial{ipert}.stde = nanstd(control_perttrial{ipert}.dat,0,1)/sqrt(size(control_perttrial{ipert}.dat,1));

  patient_perttrial{ipert}.absmean = nanmean(patient_perttrial{ipert}.absdat,1);
  patient_perttrial{ipert}.absstde = nanstd(patient_perttrial{ipert}.absdat,0,1)/sqrt(size(patient_perttrial{ipert}.absdat,1));
  control_perttrial{ipert}.absmean = nanmean(control_perttrial{ipert}.absdat,1);
  control_perttrial{ipert}.absstde = nanstd(control_perttrial{ipert}.absdat,0,1)/sqrt(size(control_perttrial{ipert}.absdat,1));
  all_perttrial{ipert}.mean = nanmean([patient_perttrial{ipert}.mean; control_perttrial{ipert}.mean],1);
end

patient_absperttrial.mean = nanmean(patient_absperttrial.dat,1);
patient_absperttrial.stde = nanstd(patient_absperttrial.dat,0,1)/sqrt(size(patient_absperttrial.dat,1));
control_absperttrial.mean = nanmean(control_absperttrial.dat,1);
control_absperttrial.stde = nanstd(control_absperttrial.dat,0,1)/sqrt(size(control_absperttrial.dat,1));

%% figure 3 - up and down trials separately

nhf = nhf + 1; 
hf(nhf) = figure;
nhf_meandat = nhf;
for ipert = 1:nperts
  hax(nhf,ipert) = subplot(2,1,ipert);
  hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.mean(idxes4anal),'Color', pat_colour, 'LineWidth', 3);
  hold on
  hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.mean(idxes4anal),'Color', hc_colour, 'LineWidth', 3);
  hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.mean(idxes4anal) + patient_perttrial{ipert}.stde(idxes4anal),'Color', pat_colour);
  hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.mean(idxes4anal) - patient_perttrial{ipert}.stde(idxes4anal),'Color', pat_colour);
  hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.mean(idxes4anal) + control_perttrial{ipert}.stde(idxes4anal),'Color', hc_colour);
  hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.mean(idxes4anal) - control_perttrial{ipert}.stde(idxes4anal),'Color', hc_colour);
  goodplot
  a = axis; axis([tlims4anal a(3:4)]);
  axis([-0.2 1 -60 60])
end
legend('patients', 'controls')
%% figure 4 - up and down sep

nhf = nhf + 1; hf(nhf) = figure;
nhf_meandat2 = nhf;
hax(nhf,1) = subplot(1,2,1);
for ipert = 1:nperts
  hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.mean(idxes4anal),'Color', pat_colour, 'LineWidth', 3);
  hold on
  hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.mean(idxes4anal) + patient_perttrial{ipert}.stde(idxes4anal),'Color', pat_colour);
  hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.mean(idxes4anal) - patient_perttrial{ipert}.stde(idxes4anal),'Color', pat_colour);
goodplot
end
a = axis; axis([tlims4anal a(3:4)]);
  axis([-0.2 1 -60 60])
legend('patients')

hax(nhf,2) = subplot(1,2,2);
for ipert = 1:nperts
  hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.mean(idxes4anal),'Color', hc_colour, 'LineWidth', 3);
  hold on
  hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.mean(idxes4anal) + control_perttrial{ipert}.stde(idxes4anal),'Color', hc_colour);
  hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.mean(idxes4anal) - control_perttrial{ipert}.stde(idxes4anal),'Color', hc_colour);
goodplot
end
a = axis; axis([tlims4anal a(3:4)]);
axis([-0.2 1 -60 60])
legend('controls')

print(gcf, '-dpdf', '-r150', 'figures/postpost/up_down_separate_postpost.pdf');


%% figure 5 - up and down flipped and together but not averaged

nhf = nhf + 1; hf(nhf) = figure;
nhf_absmeandat2 = nhf;
hax(nhf,1) = subplot(1,2,1);
for ipert = 1:nperts
  hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.absmean(idxes4anal),'Color', pat_colour, 'LineWidth', 3);
  hold on
  hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.absmean(idxes4anal) + patient_perttrial{ipert}.absstde(idxes4anal),'Color', pat_colour);
  hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.absmean(idxes4anal) - patient_perttrial{ipert}.absstde(idxes4anal),'Color', pat_colour);
goodplot
end
a = axis; axis([tlims4anal a(3:4)]);
axis([-0.2 1 -20 60])
legend('controls')

hax(nhf,2) = subplot(1,2,2);
for ipert = 1:nperts
  hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.absmean(idxes4anal),'Color', hc_colour, 'LineWidth', 3);
  hold on
  hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.absmean(idxes4anal) + control_perttrial{ipert}.absstde(idxes4anal),'Color', hc_colour);
  hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.absmean(idxes4anal) - control_perttrial{ipert}.absstde(idxes4anal),'Color', hc_colour);
goodplot
end
a = axis; axis([tlims4anal a(3:4)]);
axis([-0.2 1 -20 60])
legend('patients')

print(gcf, '-dpdf', '-r150', 'figures/postpost/up_down_together_postpost.pdf');


%% figure 6 - up and down flipped and averaged
%% THIS IS IT
nhf = nhf + 1; hf(nhf) = figure;
nhf_absmeandat3 = nhf;
hax(nhf,1) = axes;
hpl = plot(frame_taxis(idxes4anal),patient_absperttrial.mean(idxes4anal),'Color', pat_colour, 'LineWidth', 3);
hold on
hpl = plot(frame_taxis(idxes4anal),patient_absperttrial.mean(idxes4anal) + patient_absperttrial.stde(idxes4anal),'Color', pat_colour);
hpl = plot(frame_taxis(idxes4anal),patient_absperttrial.mean(idxes4anal) - patient_absperttrial.stde(idxes4anal),'Color', pat_colour);
hpl = plot(frame_taxis(idxes4anal),control_absperttrial.mean(idxes4anal),'Color', hc_colour, 'LineWidth', 3);
hpl = plot(frame_taxis(idxes4anal),control_absperttrial.mean(idxes4anal) + control_absperttrial.stde(idxes4anal),'Color', hc_colour);
hpl = plot(frame_taxis(idxes4anal),control_absperttrial.mean(idxes4anal) - control_absperttrial.stde(idxes4anal),'Color', hc_colour);
a = axis; axis([tlims4anal a(3:4)]);
goodplot
print(gcf, '-dpdf', '-r150', 'figures/postpost/up_down_separate_2.pdf');


for ipert = 1:nperts
  for iseg = 1:nsegs
    for itrial = 1:patient_perttrial{ipert}.ntrials
      patient_perttrial{ipert}.dev(itrial,:) = patient_perttrial{ipert}.dat(itrial,:) - all_perttrial{ipert}.mean;
    end
    for itrial = 1:control_perttrial{ipert}.ntrials
      control_perttrial{ipert}.dev(itrial,:) = control_perttrial{ipert}.dat(itrial,:) - all_perttrial{ipert}.mean;
    end
  end
end

for ipert = 1:nperts
  patient_perttrial{ipert}.devmean = nanmean(patient_perttrial{ipert}.dev,1);
  patient_perttrial{ipert}.devstde = nanstd(patient_perttrial{ipert}.dev,0,1)/sqrt(size(patient_perttrial{ipert}.dev,1));
  control_perttrial{ipert}.devmean = nanmean(control_perttrial{ipert}.dev,1);
  control_perttrial{ipert}.devstde = nanstd(control_perttrial{ipert}.dev,0,1)/sqrt(size(control_perttrial{ipert}.dev,1));
end

yes_plot_devs = 0;
if yes_plot_devs
  nhf = nhf + 1; hf(nhf) = figure;
  nhf_devmeandat = nhf;
  for ipert = 1:nperts
    hax(nhf,ipert) = subplot(2,1,ipert);
    hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.devmean(idxes4anal),'Color', pat_colour, 'LineWidth', 3);
    hold on
    hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.devmean(idxes4anal),'Color', hc_colour, 'LineWidth', 3);
    
    hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.devmean(idxes4anal) + patient_perttrial{ipert}.devstde(idxes4anal),'Color', pat_colour);
    hpl = plot(frame_taxis(idxes4anal),patient_perttrial{ipert}.devmean(idxes4anal) - patient_perttrial{ipert}.devstde(idxes4anal),'Color', pat_colour);
    hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.devmean(idxes4anal) + control_perttrial{ipert}.devstde(idxes4anal),'Color', hc_colour);
    hpl = plot(frame_taxis(idxes4anal),control_perttrial{ipert}.devmean(idxes4anal) - control_perttrial{ipert}.devstde(idxes4anal),'Color', hc_colour);
    a = axis; axis([tlims4anal a(3:4)]);
  end
end

for iseg = 1:nsegs
  for ipert = 1:nperts
    for itrial = 1:patient_perttrial{ipert}.ntrials
      patient_perttrial{ipert}.segmean(itrial,iseg) = nanmean(patient_perttrial{ipert}.dat(itrial,idxes4segs{iseg}),2);
      patient_perttrial{ipert}.segabsmean(itrial,iseg) = nanmean(patient_perttrial{ipert}.absdat(itrial,idxes4segs{iseg}),2);
      patient_perttrial{ipert}.segdevmean(itrial,iseg) = nanmean(patient_perttrial{ipert}.dev(itrial,idxes4segs{iseg}),2);
    end
    for itrial = 1:control_perttrial{ipert}.ntrials
      control_perttrial{ipert}.segmean(itrial,iseg) = nanmean(control_perttrial{ipert}.dat(itrial,idxes4segs{iseg}),2);
      control_perttrial{ipert}.segabsmean(itrial,iseg) = nanmean(control_perttrial{ipert}.absdat(itrial,idxes4segs{iseg}),2);
      control_perttrial{ipert}.segdevmean(itrial,iseg) = nanmean(control_perttrial{ipert}.dev(itrial,idxes4segs{iseg}),2);
    end
  end
  for itrial = 1:patient_absperttrial.ntrials
    patient_absperttrial.segmean(itrial,iseg) = nanmean(patient_absperttrial.dat(itrial,idxes4segs{iseg}),2);
  end
  for itrial = 1:control_absperttrial.ntrials
    control_absperttrial.segmean(itrial,iseg) = nanmean(control_absperttrial.dat(itrial,idxes4segs{iseg}),2);
  end
end

for ipert = 1:nperts
  for iseg = 1:nsegs
    groupseg_anova_idat = 0;
    for itrial = 1:patient_perttrial{ipert}.ntrials
      groupseg_anova_idat = groupseg_anova_idat + 1;
      groupseg_anova_datY(groupseg_anova_idat) = patient_perttrial{ipert}.segmean(itrial,iseg);
      groupseg_anova_devY(groupseg_anova_idat) = patient_perttrial{ipert}.segdevmean(itrial,iseg);
      ng = 0;
      ng = ng + 1; groupseg_anova_G{ng}(groupseg_anova_idat) = 'p';
    end
    for itrial = 1:control_perttrial{ipert}.ntrials
      groupseg_anova_idat = groupseg_anova_idat + 1;
      groupseg_anova_datY(groupseg_anova_idat) = control_perttrial{ipert}.segmean(itrial,iseg);
      groupseg_anova_devY(groupseg_anova_idat) = control_perttrial{ipert}.segdevmean(itrial,iseg);
      ng = 0;
      ng = ng + 1; groupseg_anova_G{ng}(groupseg_anova_idat) = 'c';
    end
    groupseg_varnames = {'group'};
    dims2comp = [1];
    [groupseg_p(ipert,iseg),groupseg_tab{ipert,iseg},groupseg_stats{ipert,iseg},groupseg_terms{ipert,iseg}] = anovan(groupseg_anova_datY,groupseg_anova_G,'varnames',groupseg_varnames,'display','off');
    % [groupseg_p(ipert,iseg),groupseg_tab{ipert,iseg},groupseg_stats{ipert,iseg},groupseg_terms{ipert,iseg}] = anovan(groupseg_anova_devY,groupseg_anova_G,'varnames',groupseg_varnames,'display','off');
  end
end

%% plot significance on final fig 2

for iseg = 1:nsegs
  for ipert = 1:nperts
    if groupseg_p(ipert,iseg) < alphalev
      if yes_plot_devs
        axes(hax(nhf_devmeandat,ipert));
        hp = vpatch(tlims4segs(iseg,1),tlims4segs(iseg,2),sig_colour);
        move2back([],hp);
      end
      axes(hax(nhf_meandat,ipert));
      hp = vpatch(tlims4segs(iseg,1),tlims4segs(iseg,2),sig_colour);
      move2back([],hp);
    end
  end
end


for igroup = 1:2
  clear pertseg_anova_datY pertsegabs_anova_datY pertseg_anova_devY pertseg_anova_G
  for iseg = 1:nsegs
    pertseg_anova_idat = 0;
    for ipert = 1:nperts
      if igroup == 1
        for itrial = 1:patient_perttrial{ipert}.ntrials
          pertseg_anova_idat = pertseg_anova_idat + 1;
          pertseg_anova_datY(pertseg_anova_idat) = patient_perttrial{ipert}.segmean(itrial,iseg);
          pertsegabs_anova_datY(pertseg_anova_idat) = patient_perttrial{ipert}.segabsmean(itrial,iseg);
          pertseg_anova_devY(pertseg_anova_idat) = patient_perttrial{ipert}.segdevmean(itrial,iseg);
          ng = 0;
          ng = ng + 1; pertseg_anova_G{ng}(pertseg_anova_idat) = ipert;
        end
      else
        for itrial = 1:control_perttrial{ipert}.ntrials
          pertseg_anova_idat = pertseg_anova_idat + 1;
          pertseg_anova_datY(pertseg_anova_idat) = control_perttrial{ipert}.segmean(itrial,iseg);
          pertsegabs_anova_datY(pertseg_anova_idat) = control_perttrial{ipert}.segabsmean(itrial,iseg);
          pertseg_anova_devY(pertseg_anova_idat) = control_perttrial{ipert}.segdevmean(itrial,iseg);
          ng = 0;
          ng = ng + 1; pertseg_anova_G{ng}(pertseg_anova_idat) = ipert;
        end
      end
    end
    pertseg_varnames = {'pert'};
    dims2comp = [1];
    [pertseg_p(igroup,iseg),pertseg_tab{igroup,iseg},pertseg_stats{igroup,iseg},pertseg_terms{igroup,iseg}] = anovan(pertseg_anova_datY,pertseg_anova_G,'varnames',pertseg_varnames,'display','off');
    [pertsegabs_p(igroup,iseg),pertsegabs_tab{igroup,iseg},pertsegabs_stats{igroup,iseg},pertsegabs_terms{igroup,iseg}] = anovan(pertsegabs_anova_datY,pertseg_anova_G,'varnames',pertseg_varnames,'display','off');
    % [pertseg_p(igroup,iseg),pertseg_tab{igroup,iseg},pertseg_stats{igroup,iseg},pertseg_terms{igroup,iseg}] = anovan(pertseg_anova_devY,pertseg_anova_G,'varnames',pertseg_varnames,'display','off');
  end
end

for igroup = 1:2
  for iseg = 1:nsegs
    if pertseg_p(igroup,iseg) < alphalev
      axes(hax(nhf_meandat2,igroup));
      hp = vpatch(tlims4segs(iseg,1),tlims4segs(iseg,2),sig_colour);
      move2back([],hp);
    end
    if pertsegabs_p(igroup,iseg) < 0.05 % separate alphalev here
      axes(hax(nhf_absmeandat2,igroup));
      hp = vpatch(tlims4segs(iseg,1),tlims4segs(iseg,2),sig_colour);
      move2back([],hp);
    end
  end
end

%% Final Figures

clear abspertseg_anova_datY abspertseg_anova_devY abspertseg_anova_G abspertseg_p
for iseg = 1:nsegs
  abspertseg_anova_idat = 0;
  for itrial = 1:patient_absperttrial.ntrials
    abspertseg_anova_idat = abspertseg_anova_idat + 1;
    abspertseg_anova_datY(abspertseg_anova_idat) = patient_absperttrial.segmean(itrial,iseg);
    ng = 0;
    ng = ng + 1; abspertseg_anova_G{ng}(abspertseg_anova_idat) = 'p';
  end
  for itrial = 1:control_absperttrial.ntrials
    abspertseg_anova_idat = abspertseg_anova_idat + 1;
    abspertseg_anova_datY(abspertseg_anova_idat) = control_absperttrial.segmean(itrial,iseg);
    ng = 0;
    ng = ng + 1; abspertseg_anova_G{ng}(abspertseg_anova_idat) = 'g';
  end
  abspertseg_varnames = {'pert'};
  dims2comp = [1];
  [abspertseg_p(iseg),abspertseg_tab{iseg},abspertseg_stats{iseg},abspertseg_terms{iseg}] = anovan(abspertseg_anova_datY,abspertseg_anova_G,'varnames',abspertseg_varnames,'display','off');
end

for iseg = 1:nsegs
  if abspertseg_p(iseg) < alphalev
    axes(hax(nhf_absmeandat3,1));
    hp = vpatch(tlims4segs(iseg,1),tlims4segs(iseg,2),sig_colour);
    move2back([],hp);
  end
end
print(gcf, '-dpdf', '-r150', 'figures/postpost/finalfigure.pdf');

save postpost
