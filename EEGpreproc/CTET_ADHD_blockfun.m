function trl = CTET_ADHD_blockfun(cfg)

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

for i=1:max(table.BlockN)
    % add this to the trl definition
    subtable=table(table.BlockN==i,:);
    begsample     = subtable.Sample(1) - cfg.trialdef.prestim*hdr.Fs;
    endsample     = subtable.Sample(end) + cfg.trialdef.poststim*hdr.Fs - 1;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl(end+1, :) = [round([begsample endsample offset])];
end
