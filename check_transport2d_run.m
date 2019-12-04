clearvars;

get_values = '/home/jjl/MFEM/mfem/miniapps/tools/get-values';
NP = 10;

root_dir = '/home/jjl/MFEM/mfem/miniapps/plasma/dt_e-8';
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

%% Parse mfem_root files to get times
time = nan(1,nt);
for i = 0:nt-1
    fname = fullfile(root_dir,strcat(prefix,'_',pad_zeros(i,6),'.mfem_root'));
    [nl,file_lines] = num_lines_file(fname);
    [~,time_loc] = grep('-i -s','"time"',fname);
    npart = time_loc.pcount;
    if npart ~= 1
        error('trouble parsing file %s,fname')
    end
    line = file_lines{time_loc.line};
    coldex = strfind(line,':');
    time(i+1) = sscanf(line(coldex+1:end-1),'%f');    
end

%% Parse run.log to get convergence info
fname = fullfile(root_dir,'run.log');
[nl,file_lines] = num_lines_file(fname);
[~,newt_loc] = grep('-i -s',"Newton",fname);
for i = 1:newt_loc.pcount
    line = file_lines{newt_loc.line(i)};
%     coldex = strfind(line,':');
%     Newt_iter(i) = sscanf(line(coldex+1:end-1),'%d'); 
    newt_data = textscan(line,'%s %s %d : %s = %f, %s = %f');
    newt_iter(i) = newt_data{3};
    newt_res(i) = newt_data{5};
    if ~isempty(newt_data{7})
        newt_res_norm(i) = newt_data{7};
    else
        newt_res_norm(i) = NaN;
    end
%     aaa=1
end

FS = 12;

figure; hold on; box on; grid on;
subplot(3,1,1); hold on; grid on; box on;
plot(newt_iter,'linew',2)
xlabel('Iteration','fontsize',FS)
title('Newton iterations','fontsize',FS)
subplot(3,1,2); hold on; grid on; box on;
plot(newt_res,'linew',2)
set(gca,'yscale','log')
xlabel('Iteration','fontsize',FS)
title('Newton Residual Norm','fontsize',FS)
subplot(3,1,3); hold on; grid on; box on;
plot(newt_res_norm,'linew',2)
xlabel('Iteration','fontsize',FS)
title('Newton Relative Residual Norm','fontsize',FS)
set(gca,'fontsize',FS)
set(gcf,'color','w')

%%
fname = fullfile(root_dir,'run.log');
[nl,file_lines] = num_lines_file(fname);
[~,iter_loc] = grep('-s',"Iteration",fname);
for i = 1:iter_loc.pcount
    line = file_lines{iter_loc.line(i)};
    iter_data = textscan(line,'%s : %f %s : %f %s %s = %f');
    
    iter_pass(i) = iter_data{2};
    iter_num(i) = iter_data{4};
    iter_res(i) = iter_data{7};

%     aaa=1
end

figure; hold on; box on; grid on;
plot(iter_res)
set(gca,'yscale','log')

%% Call get-values
npoints_want = 100;
Xwant = linspace(0.39,0.65,npoints_want);
Ywant = 0.*Xwant;
infilename = 'points.in';
fid = fopen(fullfile(root_dir,infilename),'w');
fprintf(fid,'%d %d \n',npoints_want,2);
for i = 1:npoints_want
    fprintf(fid,'%f %f\n',Xwant(i),Ywant(i));
end
fclose(fid);
% adsf

% index = 0
% point = [0.5,0.221,0.,0.221,0.5,0];
outfilename = 'myoutput.out';

% command = strcat("mpirun -np ",num2str(NP)," ",get_values," -r ",fullfile(root_dir,prefix)," -c ",num2str(index)," -p ","""",num2str(point),""""," -o ",fullfile(root_dir,outfilename))
for i = 0:nt - 1
    fprintf('Requesting output for index %d of %d\n',i,nt-1)
    index = i;
    command = strcat("mpirun -np ",num2str(NP)," ",get_values," -r ",fullfile(root_dir,prefix)," -c ",num2str(index)," -pf ",fullfile(root_dir,infilename)," -o ",fullfile(root_dir,outfilename));
    [status,output] = system(command);
    
    %% parse output
    % fields: [ B Poloidal, B Toroidal, Electron Temperature, Ion Density, Ion Parallel Velocity, Ion Temperature, Neutral Density, n_e Chi_e Parallel, n_e Chi_e Perpendicular ]
    data = dlmread(fullfile(root_dir,outfilename),'',1,0);
    npoints = size(data,1);
    X(:,i+1) = data(:,2);
    Y(:,i+1) = data(:,3);
    Bpx(:,i+1) = data(:,4);
    Bpy(:,i+1) = data(:,5);
    Bt(:,i+1) = data(:,6);
    Te(:,i+1) = data(:,7);
    ni(:,i+1) = data(:,8);
    vi(:,i+1) = data(:,9);
    Ti(:,i+1) = data(:,10);
    n0(:,i+1) = data(:,11);
    kappa_prl(:,i+1) = data(:,12);
    kappa_perp(:,i+1) = data(:,13);
    
end

LW = 2;
FS = 14;

%% plot
figure; hold on; box on; grid on;
plot(X,Te,'linew',LW)
xlabel('X (m)','fontsize',FS)
ylabel('T_e (eV)','fontsize',FS)
title('Outboard profile at Y = 0')
set(gcf,'color','w')
set(gca,'fontsize',FS)


ind = find(X>0.6,1,'first');

figure; hold on; box on; grid on;
plot(time,Te(ind,:),'o-','linew',LW)
xlabel('time (s)','fontsize',FS)
ylabel('T_e (eV)','fontsize',FS)
title('T_e vs t at X = 0.6, Y = 0')
set(gcf,'color','w')
set(gca,'fontsize',FS)

% figure; hold on; box on; grid on;
% plot(X,n0)
% xlabel('X (m)')
% ylabel('n_0 (m^{-3})')
% title('Outboard profile at Y = 0')
set(gcf,'color','w')


%%
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