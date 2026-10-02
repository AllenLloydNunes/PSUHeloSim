%% run_all_tests.m
% Runs the STEP 1, 2 and 3 verification suites and generates figures.
%
%   matlab -batch "run_all_tests"
%
% Every pass criterion in PLAN.md Sec.9 (steps 1-3) and every safety check in
% PLAN.md Sec.7 is exercised here.
%
% Each test script clears the workspace on entry so it stays runnable
% standalone; results therefore come back via .mat files, not the workspace.

clear; clc;
here = fileparts(mfilename('fullpath'));
cd(here);

fprintf('################################################################\n');
fprintf(' Simplified Landing Model -- verification suite\n');
fprintf(' %s\n', datestr(now));
fprintf('################################################################\n\n');

test_step1_freefall;
fprintf('\n\n');
test_step2_contact;
fprintf('\n\n');
test_step3_friction_slope;

r1 = load('step1_result.mat');
r2 = load('step2_result.mat');
r3 = load('step3_result.mat');

%% Figures
fprintf('\n----------------------------------------------------------------\n');
fprintf(' Generating figures\n');
fprintf('----------------------------------------------------------------\n');
plot_landing_results(r2.S_settle, 'test2A_damped_drop_to_equilibrium', ...
    'Test 2A: damped drop from 0.2 m above contact, settles to static equilibrium')

plot_landing_results(r2.S_econt, 'test2E_undamped_oscillation_in_contact', ...
    'Test 2E: undamped heave oscillation, legs never leave ground (energy conservation, smooth)')

plot_landing_results(r2.S_ebounce, 'test2F_undamped_repeated_bounce', ...
    'Test 2F: undamped repeated bounce through contact on/off (energy conservation, discontinuous)')

plot_landing_results(r3.S3A, 'test3A_hold_on_5deg_slope', ...
    'Test 3A: settle and hold on a 5 deg incline with friction (regularized friction creeps downhill)')

plot_landing_results(r3.S3B{1}, 'test3B_slide_downhill_to_stop', ...
    'Test 3B: slide downhill (tail-first) from 3 m/s to a stop on a 5 deg incline')

plot_landing_results(r3.S3B{2}, 'test3B_slide_uphill_to_stop', ...
    'Test 3B: slide uphill (nose-first) from 3 m/s to a stop on a 5 deg incline')

fprintf('\n################################################################\n');
if r1.step1_pass && r2.step2_pass && r3.step3_pass
    fprintf(' OVERALL: PASS -- steps 1, 2 and 3 complete\n');
else
    fprintf(' OVERALL: FAIL -- step1_pass=%d  step2_pass=%d  step3_pass=%d\n', ...
            r1.step1_pass, r2.step2_pass, r3.step3_pass);
end
fprintf('################################################################\n');
