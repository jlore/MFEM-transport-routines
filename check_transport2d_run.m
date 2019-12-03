clearvars;

root_dir = 'C:\Work\RFSciDAC\Transport2D\dt_e-1';
prefix = 'Transport2D-Parallel';

dir_data = dir(root_dir);

%% Initial settings
fprintf('\nUsing root_dir = %s with prefix = %s\n',root_dir,prefix)


%%  Parse dir to get nt
nt = 0;
for i = 1:length(dir_data)
    if dir_data(i).isdir
        % Check if prefix_######
        if regexp(dir_data(i).name, regexptranslate('wildcard', strcat(prefix,'_******')))
            nt = nt + 1;
        end
    end
end

if nt == 0
    error('Could not find any time steps using prefix: %s',prefix)
else
    fprintf('Found nt = %d timesteps\n',nt)
end

% Parse mfem_root files to get times
time = nan(1,nt);
for i = 0:nt-1
    fname = fullfile(root_dir,strcat(prefix,'_',pad_zeros(i,6),'.mfem_root'));
    [nl,file_lines] = num_lines_file(fname);
    [~,time_loc] = grep('-i -s','"time"',fname);
    npart = time_loc.pcount;
    if npart ~= 1
        error('trouble parsing file %s,fname')
    end
    line = cell2mat(file_lines(time_loc.line));
    coldex = strfind(line,':');
    time(i+1) = sscanf(line(coldex+1:end-1),'%f');    
end




function out_str = pad_zeros(i,n)    
    if i < 0
        error('bad i')
    end
    if n < 1
        error('bad n')
    end
    if i > 10^n - 1
        error('n is too small to fit i')
    end
    
    out_str = num2str(i);
    while length(out_str) < n
        out_str = strcat('0',out_str);
    end   
end