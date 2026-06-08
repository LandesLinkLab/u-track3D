function test_chenouard_v2()
%TEST_CHENOUARD_V2  Smoke test for the rewritten chenouardMetrics.
%  Builds three controlled scenarios and checks that alpha, beta, JSC
%  produce the expected values.
%
%  Scenarios:
%    (1) Perfect: tracked == GT  -> alpha=beta=JSC=1
%    (2) All dummy: trackedTraj empty -> alpha=0, beta=0, JSC=0
%    (3) Half-shifted: tracked is GT with positions shifted by 0.5*gate
%        per frame -> alpha=0.5 (each frame contributes 0.5*gate cost out
%        of gate max)

addpath(fileparts(mfilename('fullpath')));

gate = 100;

% --- Scenario 1: perfect tracking ---
gt = makeGT(5, 50);                    % 5 tracks, 50 frames each
tr = gt;                                % identical
r1 = chenouardMetrics(tr, gt, 'GateDist', gate);
assert(abs(r1.alpha - 1) < 1e-9, 'Scenario 1 alpha != 1: got %.6f', r1.alpha);
assert(abs(r1.beta  - 1) < 1e-9, 'Scenario 1 beta  != 1: got %.6f', r1.beta);
assert(abs(r1.JSC   - 1) < 1e-9, 'Scenario 1 JSC   != 1: got %.6f', r1.JSC);
fprintf('Scenario 1 (perfect):       alpha=%.4f beta=%.4f JSC=%.4f  OK\n', ...
    r1.alpha, r1.beta, r1.JSC);

% --- Scenario 2: empty tracked output ---
r2 = chenouardMetrics([], gt, 'GateDist', gate);
assert(r2.alpha == 0, 'Scenario 2 alpha != 0');
assert(r2.beta  == 0, 'Scenario 2 beta  != 0');
assert(r2.JSC   == 0, 'Scenario 2 JSC   != 0');
fprintf('Scenario 2 (empty tracked): alpha=%.4f beta=%.4f JSC=%.4f  OK\n', ...
    r2.alpha, r2.beta, r2.JSC);

% --- Scenario 3: every detection shifted by 0.5*gate ---
tr3 = gt;
tr3(:, 3) = tr3(:, 3) + 0.5 * gate;
r3 = chenouardMetrics(tr3, gt, 'GateDist', gate);
% Each frame contributes 0.5*gate to d; max possible is gate per frame.
% So d(X,Y)/d(X,empty) = 0.5 -> alpha = 0.5
assert(abs(r3.alpha - 0.5) < 1e-9, 'Scenario 3 alpha != 0.5: got %.6f', r3.alpha);
% No spurious tracks (all GT paired with real est), so beta = alpha.
assert(abs(r3.beta - r3.alpha) < 1e-9, ...
    'Scenario 3 beta should equal alpha (no spurious): got beta=%.6f alpha=%.6f', ...
    r3.beta, r3.alpha);
% All detections are within gate (shift = 0.5*gate < gate) -> all TPs
assert(abs(r3.JSC - 1) < 1e-9, 'Scenario 3 JSC != 1: got %.6f', r3.JSC);
fprintf('Scenario 3 (uniform 0.5g):  alpha=%.4f beta=%.4f JSC=%.4f  OK\n', ...
    r3.alpha, r3.beta, r3.JSC);

% --- Scenario 4: shift > gate -> all detection pairs miss the gate ---
tr4 = gt;
tr4(:, 3) = tr4(:, 3) + 1.5 * gate;
r4 = chenouardMetrics(tr4, gt, 'GateDist', gate);
% Each frame contributes gate (truncated). d(X,Y) = d(X,empty).
assert(abs(r4.alpha) < 1e-9, 'Scenario 4 alpha != 0: got %.6f', r4.alpha);
% No detection pairs within gate -> TP=0
assert(r4.JSC == 0, 'Scenario 4 JSC != 0: got %.6f', r4.JSC);
fprintf('Scenario 4 (uniform 1.5g):  alpha=%.4f beta=%.4f JSC=%.4f  OK\n', ...
    r4.alpha, r4.beta, r4.JSC);

% --- Scenario 5: GT perfectly tracked + ONE spurious track far away ---
% The spurious track has the same frame range as a GT track but is offset
% far in x so it does not match any GT.
spuriousIDnew = max(gt(:,2)) + 1;
% Take a copy of GT track 1's rows and offset position by 10*gate (>> gate)
spurious = gt(gt(:,2) == 1, :);
spurious(:,2) = spuriousIDnew;
spurious(:,3) = spurious(:,3) + 10 * gate;
tr5 = [gt; spurious];
r5 = chenouardMetrics(tr5, gt, 'GateDist', gate);
% Real-match part gives alpha=1; one spurious track of length 50 hurts beta.
% dXempty = sum(|frames(kX)|)*gate = 5*50*gate = 25000
% dYempty = 50*gate = 5000 (the one spurious)
% beta = (25000 - 0) / (25000 + 5000) = 5/6
expectedBeta = 5 / 6;
assert(abs(r5.beta - expectedBeta) < 1e-9, ...
    'Scenario 5 beta != %.6f: got %.6f', expectedBeta, r5.beta);
% JSC: TP = 5*50 = 250, FN = 0, FP = 50 (the spurious 50 detections)
expectedJSC = 250 / (250 + 0 + 50);
assert(abs(r5.JSC - expectedJSC) < 1e-9, ...
    'Scenario 5 JSC != %.6f: got %.6f', expectedJSC, r5.JSC);
fprintf('Scenario 5 (1 spurious):    alpha=%.4f beta=%.4f JSC=%.4f  OK\n', ...
    r5.alpha, r5.beta, r5.JSC);

fprintf('\nAll smoke tests PASSED.\n');
end


function traj = makeGT(nTracks, nFrames)
% Build a [N x 9] trajectory matrix [frame, id, x, y, z, theta, phi, omega, intensity]
rng(0);
traj = [];
for id = 1:nTracks
    for t = 1:nFrames
        x = 100 * id + 0.5 * t;
        y = 100 * id;
        z = 0;
        traj(end+1, :) = [t, id, x, y, z, pi/4, 0, 0.5, 1]; %#ok<AGROW>
    end
end
end
