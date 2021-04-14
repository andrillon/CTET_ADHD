function trl = CTET_ADHD_trialfun(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr=ft_read_header(cfg.dataset);
trl = [];

% clean events
table=cfg.table;
for i=1:size(table,1)
    % add this to the trl definition
    begsample     = table.Sample(i) - cfg.trialdef.prestim*hdr.Fs;
    endsample     = table.Sample(i) + cfg.trialdef.poststim*hdr.Fs - 1;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl(end+1, :) = [round([begsample endsample offset])];
end
