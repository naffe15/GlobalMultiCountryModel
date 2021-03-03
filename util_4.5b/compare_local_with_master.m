% utility that compares automatically marginal likelihood and smoothed shocks
% running first locally and afterwards with -Dmypath.

% %To make the local run faster, set in gemc.dyn
%
%         @#define skipinsample = 0
%         @#define do_the_forecast = 0
%         @#define do_the_simulation = 0
%         @#define do_the_estimation=0
%         @#define do_the_slice=0
%         @#define compute_moments = 0
%         @#define realtime = 0
%         @#define plot_the_fit_and_smooth = 0
my_dir = pwd;
duplicate_GM_project([my_dir '_check_dmypath'])

args={'gemc';
    'console';
    'nointeractive';
    '-Dcompare_local_with_master'
    };

if ismac
    args = [args; '-Dmac'];
end

if isempty(dir([my_dir filesep 'gemc_results.mat']))
    dynare(args{:});
    oo0=oo_;
    % store local src
    movefile src src0
else
    oo=load([my_dir filesep 'gemc_results'],'oo_');
    oo0=oo.oo_;
    clear oo;
end
save oo0 oo0
args = [args; '-Dmypath'];
dynare(args{:});
copyfile('gemc_results.mat',[my_dir filesep 'gemc_results_dmypath.mat'])
copyfile('src',[my_dir filesep 'src_dmypath'])
load oo0

struct_compare(1,{oo0.SmoothedShocks,oo_.SmoothedShocks});
struct_compare(0,{oo0.MarginalDensity,oo_.MarginalDensity});

%%
cd(my_dir)
rmdir([my_dir '_check_dmypath'],'s')