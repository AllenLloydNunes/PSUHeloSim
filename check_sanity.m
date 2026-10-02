function R = check_sanity(S, label, opts)
%CHECK_SANITY  Applies the PLAN.md Sec.7 safety checks to a completed run.
%
%   R = check_sanity(S, label, opts)
%
%   opts fields (all optional):
%     .checkEnergy   true to enforce an energy-drift bound (only meaningful
%                    when P.c == 0 and P.mu == 0)
%     .energyTol     relative drift tolerance (default 1e-6)
%     .requireAllContact  true to assert every leg stays in contact throughout
%
%   Returns R with .pass plus the measured values, and prints a report.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'checkEnergy'),        opts.checkEnergy = false;       end
if ~isfield(opts,'energyTol'),          opts.energyTol   = 1e-6;        end
if ~isfield(opts,'requireAllContact'),  opts.requireAllContact = false; end

P = S.P;
R = struct('pass', true);

fprintf('  --- sanity checks: %s ---\n', label);

% 1. finite states
R.nFinite = nnz(~isfinite(S.x));
ok = (R.nFinite == 0);
R.pass = R.pass && ok;
report('no NaN/Inf in state history', ok, sprintf('%d non-finite entries', R.nFinite));

% 2. no adhesion: the ground may push but never pull
R.minN = min(S.N(:));
ok = (R.minN >= 0);
R.pass = R.pass && ok;
report('normal force never negative (no adhesion)', ok, sprintf('min(N) = %.3e N', R.minN));

% 3. contact-flag consistency: N > 0 implies the flag is set
R.nInconsistent = nnz( (S.N > 0) & ~S.contact );
ok = (R.nInconsistent == 0);
R.pass = R.pass && ok;
report('contact flag consistent with N > 0', ok, sprintf('%d mismatches', R.nInconsistent));

% 4. penetration is non-negative wherever contact is flagged
dContact = S.delta(S.contact);
if isempty(dContact)
    R.minDeltaInContact = NaN;
    report('penetration positive where in contact', true, 'no contact in this run');
else
    R.minDeltaInContact = min(dContact);
    ok = (R.minDeltaInContact > 0);
    R.pass = R.pass && ok;
    report('penetration positive where in contact', ok, ...
           sprintf('min(delta|contact) = %.3e m', R.minDeltaInContact));
end

% 4b. Coulomb bound: friction magnitude may never exceed mu*N on any leg
if isfield(S, 'Ft')
    FtMag = squeeze(sqrt(sum(S.Ft.^2, 2)));          % (n+1) x nLegs
    if isvector(FtMag), FtMag = reshape(FtMag, size(S.N)); end
    excess = FtMag - P.mu*S.N;
    R.maxCoulombExcess = max(excess(:));
    ok = (R.maxCoulombExcess <= 1e-9*max(1, max(P.mu*S.N(:))));
    R.pass = R.pass && ok;
    report('friction never exceeds mu*N (Coulomb bound)', ok, ...
           sprintf('max(|F_t| - mu*N) = %.3e N', R.maxCoulombExcess));
end

% 5. optional: every leg stays loaded (used by the continuous-contact energy test)
if opts.requireAllContact
    R.nLiftoff = nnz(~S.contact);
    ok = (R.nLiftoff == 0);
    R.pass = R.pass && ok;
    report('all legs remain in contact', ok, sprintf('%d leg-samples off ground', R.nLiftoff));
end

% 6. optional: energy conservation (requires c = 0 and mu = 0)
if opts.checkEnergy
    if P.c ~= 0 || P.mu ~= 0
        R.pass = false;
        report('energy check preconditions (c=0, mu=0)', false, ...
               sprintf('c = %g, mu = %g -- system is dissipative', P.c, P.mu));
    else
        E0        = S.E(1);
        R.Edrift  = max(abs(S.E - E0));
        R.EdriftRel = R.Edrift / max(abs(E0), eps);
        ok = (R.EdriftRel < opts.energyTol);
        R.pass = R.pass && ok;
        report('energy conserved', ok, ...
               sprintf('rel drift = %.3e (tol %.1e), E0 = %.4f J', ...
                       R.EdriftRel, opts.energyTol, E0));
    end
end

fprintf('  --- %s: %s ---\n', label, ternary(R.pass, 'PASS', 'FAIL'));
end

% -------------------------------------------------------------------------
function report(name, ok, detail)
fprintf('    [%s] %-46s %s\n', ternary(ok,'PASS','FAIL'), name, detail);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
