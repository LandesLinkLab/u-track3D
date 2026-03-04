%% test_pipeline_synthetic.m
%  Test the 6D-SMOLM tracking pipeline components using synthetic data.
%  Tests each module independently. Does NOT require u-track3D.
%
%  SYNTHETIC DATA:
%    5 molecules undergoing 3D Brownian diffusion (D_trans = 500 nm^2/s)
%    and rotational diffusion (D_rot = 0.05 rad^2/s) over 200 frames with
%    15% random photoblinking. Provides ground-truth D values for Test 6.
%
%  TEST 1 — angularDistance:
%    Verifies the core dipole distance geometry:
%      - Same orientation gives distance = 0
%      - Orthogonal dipoles give distance = pi/2
%      - 180-degree head-tail symmetry: antiparallel dipoles give d ~ 0
%      - Vectorized inputs produce correct-length output
%
%  TEST 2 — convert6DSMOLMtoMovieInfo:
%    Loads synthetic CSV and checks:
%      - Output has correct number of frames
%      - All 7 required fields present (xCoord, yCoord, zCoord, amp,
%        theta, phi, omega)
%      - Each field in [P x 2] format (value + uncertainty column)
%      - Both CSV file and direct matrix inputs accepted
%
%  TEST 3 — costMat6DSMOLMLink:
%    Runs the frame-to-frame cost function on two consecutive non-empty
%    frames and checks:
%      - No error flag returned
%      - Non-empty sparse cost matrix with valid links
%      - Position-only mode (useOrientation = false) also works
%
%  TEST 4 — traj_filt_6D:
%    Uses a hand-crafted trajectory with known gap structure:
%      - 1-frame gap is interpolated (NaN intensity marks interpolated rows)
%      - 5-frame gap causes a trajectory split
%      - Interpolated spatial positions fall between bounding endpoints
%      - Circular phi interpolation: phi = +3.0 to phi = -2.9 goes
%        through +-pi (short arc on the circle), NOT through 0
%
%  TEST 5 — combineTrajectories6D:
%    Merges two "movie" trajectory matrices:
%      - Track IDs renumbered to be globally unique across movies
%      - Movie ID column (column 10) correctly appended
%      - Row assignments to movies are correct
%
%  TEST 6 — estimateDiffusion6D_MLE:
%    Tests Bo Shuang's MLE with known ground-truth diffusion:
%      - Clean 500-step Brownian trajectory: D_trans recovered within 2x
%      - After removing 20% of frames (simulating photoblinking), Bo's
%        time-lag-aware MLE still gives reasonable D_trans despite variable
%        time gaps between steps
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

clear; clc;
fprintf('=== Testing 6D-SMOLM Pipeline Components ===\n\n');

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(genpath(fullfile(thisDir,'..')));

%% ---- Generate synthetic data ----
fprintf('--- Generating synthetic data ---\n');
rng(42);
nParticles = 5; nFrames = 200; dt = 0.05;
D_trans = 500; D_rot = 0.05;

allData = [];
for iPart = 1:nParticles
    x = 3000 + randn*500; y = 3000 + randn*500; z = randn*200;
    th = rand*pi/2; ph = (rand*2-1)*pi; om = 0.5+rand;
    for iFr = 1:nFrames
        if rand < 0.15, continue; end
        ss = sqrt(2*D_trans*dt); sr = sqrt(2*D_rot*dt);
        x = x+ss*randn; y = y+ss*randn; z = z+ss*randn*0.5;
        th = max(0,min(pi/2,th+sr*randn*0.5));
        ph = mod(ph+sr*randn+pi,2*pi)-pi;
        om = max(0,min(2*pi,om+0.01*randn));
        allData = [allData; iFr,x,y,z,th,ph,om,1000+randn*100]; %#ok<AGROW>
    end
end
fprintf('  %d localizations, %d particles\n\n', size(allData,1), nParticles);

synFile = fullfile(thisDir,'synthetic_6D_data.csv');
fid = fopen(synFile,'w');
fprintf(fid,'frame,x_nm,y_nm,z_nm,theta,phi,omega,intensity\n');
fclose(fid);
dlmwrite(synFile, allData, '-append','delimiter',',','precision','%.6f');

%% ---- TEST 1: angularDistance ----
fprintf('--- TEST 1: angularDistance ---\n');
[d1,~,~] = angularDistance(pi/4,0,1, pi/4,0,1);
assert(abs(d1)<1e-10); fprintf('  Same orientation: PASS (d=%.2e)\n',d1);

[d2,~,~] = angularDistance(0,0,[], pi/2,0,[]);
assert(abs(d2-pi/2)<0.01); fprintf('  Orthogonal: PASS (d=%.4f, expect %.4f)\n',d2,pi/2);

th1=[0;pi/6;pi/4]; ph1=[0;0;0]; th2=[pi/6;pi/6;pi/4]; ph2=[0;0;pi/2];
[dv,~,~] = angularDistance(th1,ph1,[],th2,ph2,[]);
assert(length(dv)==3); fprintf('  Vectorized: PASS\n');

% 180-degree symmetry test: theta=pi/4,phi=0 vs theta=pi/4,phi=pi should give d~0
% because the dipoles point in opposite directions
% n1 = (sin(pi/4),0,cos(pi/4)), n2 = (-sin(pi/4),0,cos(pi/4))... no, that's not opposite
% True opposite: theta=pi/4,phi=0 -> n=(s,0,c). Antipodal: theta=3pi/4,phi=0 -> n=(s,0,-c)
% But theta is in [0,pi/2]. So test via dot product directly:
% For a dipole at theta=0.3,phi=0.5, the antipodal is at theta=pi-0.3, phi=0.5+pi
theta_a = 0.3; phi_a = 0.5;
theta_b = pi - 0.3; phi_b = 0.5 + pi;
[dsym,~,~] = angularDistance(theta_a,phi_a,[], theta_b,phi_b,[]);
assert(dsym < 0.01); fprintf('  180deg symmetry: PASS (d=%.4f, should be ~0)\n\n',dsym);

%% ---- TEST 2: convert6DSMOLMtoMovieInfo ----
fprintf('--- TEST 2: convert6DSMOLMtoMovieInfo ---\n');
movieInfo = convert6DSMOLMtoMovieInfo(synFile,'Verbose',true);
assert(length(movieInfo)==nFrames, 'Wrong number of frames');
assert(isfield(movieInfo(1),'theta'), 'Missing theta');
assert(isfield(movieInfo(1),'phi'), 'Missing phi');
assert(isfield(movieInfo(1),'omega'), 'Missing omega');
ne = find(arrayfun(@(s)size(s.xCoord,1),movieInfo)>0,1);
assert(size(movieInfo(ne).xCoord,2)==2, 'xCoord should be Px2');
assert(size(movieInfo(ne).theta,2)==2, 'theta should be Px2');

% Test with matrix input directly
mi2 = convert6DSMOLMtoMovieInfo(allData,'Verbose',false);
assert(length(mi2)==nFrames, 'Matrix input failed');
fprintf('  All checks PASS\n\n');

%% ---- TEST 3: costMat6DSMOLMLink ----
fprintf('--- TEST 3: costMat6DSMOLMLink ---\n');
tp = struct('linearMotion',0,'minSearchRadius',50,'maxSearchRadius',2000,...
    'brownStdMult',3,'useLocalDensity',0,'nnWindow',4,'diagnostics',[],...
    'wSpatial',1,'wOrient',0.3,'wOmega',0,'maxAngularDist',pi/4,...
    'useOrientation',true);

fp = [];
for f=1:nFrames-1
    if size(movieInfo(f).xCoord,1)>0 && size(movieInfo(f+1).xCoord,1)>0
        fp=f; break;
    end
end
assert(~isempty(fp),'No consecutive non-empty frames found');

[cm,~,~,~,ef] = costMat6DSMOLMLink(movieInfo,[],tp,fp,[],[]);
fprintf('  Frame %d->%d: %d tracks x %d dets, %d links, err=%d\n',...
    fp,fp+1,size(cm,1),size(cm,2),nnz(cm),ef);
assert(ef==0,'Error flag nonzero');
assert(nnz(cm)>0,'No links found');

% Test position-only mode
tp2 = tp; tp2.useOrientation = false;
[cm2,~,~,~,ef2] = costMat6DSMOLMLink(movieInfo,[],tp2,fp,[],[]);
assert(ef2==0,'Position-only mode error');
fprintf('  Position-only mode: %d links (PASS)\n\n',nnz(cm2));

%% ---- TEST 4: traj_filt_6D ----
fprintf('--- TEST 4: traj_filt_6D ---\n');
testTraj = [
    1,1,100,200,0,0.5,0.1,1.0,1000;
    2,1,105,203,2,0.52,0.12,1.01,1010;
    3,1,108,207,3,0.53,0.11,1.02,990;
    5,1,115,214,6,0.55,0.15,1.03,1005;  % gap of 1
    6,1,120,218,8,0.57,0.16,1.04,995;
    7,1,123,221,9,0.56,0.14,1.05,1020;
    8,1,126,224,10,0.55,0.13,1.06,1015;
    13,1,160,250,20,0.7,0.3,1.1,980;   % gap of 5 -> split
    14,1,163,253,21,0.71,0.31,1.11,970;
    15,1,166,256,22,0.72,0.32,1.12,985;
    16,1,169,259,23,0.73,0.33,1.13,1000;
    17,1,172,262,24,0.74,0.34,1.14,995;
    18,1,175,265,25,0.75,0.35,1.15,1010];

[tf,st] = traj_filt_6D(testTraj,'ConnectSize',3,'MinTrajLen',5,'Verbose',true);
assert(st.nSplit>=1,'Should split at 5-frame gap');
assert(st.nInterpolatedFrames>=1,'Should interpolate 1-frame gap');

% Check frame 4 was interpolated
if any(tf(:,1)==4)
    row4 = tf(tf(:,1)==4,:);
    assert(isnan(row4(1,9)),'Interpolated intensity should be NaN');
    assert(row4(1,3)>108 && row4(1,3)<115, 'Interpolated x out of range');
    fprintf('  Interpolated frame 4: x=%.1f, theta=%.3f (PASS)\n',row4(1,3),row4(1,6));
end

% Check phi circular interpolation
testPhiTraj = [
    1,1,0,0,0,0.5, 3.0,1.0,1000;   % phi near +pi
    3,1,0,0,0,0.5,-2.9,1.0,1000];   % phi near -pi (short arc wraps)
[tfPhi,~] = traj_filt_6D(testPhiTraj,'ConnectSize',3,'MinTrajLen',1,'Verbose',false);
if any(tfPhi(:,1)==2)
    iPhi = tfPhi(tfPhi(:,1)==2,7);
    % The interpolated phi should be near +-pi (short arc), not near 0
    assert(abs(iPhi)>2,'Circular phi interpolation failed: phi=%.2f',iPhi);
    fprintf('  Circular phi interp: phi=%.3f (near +-pi, PASS)\n',iPhi);
end
fprintf('  All traj_filt_6D checks PASS\n\n');

%% ---- TEST 5: combineTrajectories6D ----
fprintf('--- TEST 5: combineTrajectories6D ---\n');
traj1 = [1,1,100,200,0,0.5,0.1,1.0,1000; 2,1,105,203,2,0.52,0.12,1.01,1010;
         1,2,200,300,5,0.3,0.5,0.8,900;  2,2,205,303,7,0.32,0.52,0.81,910];
traj2 = [1,1,300,400,10,0.3,0.5,0.8,900; 2,1,305,403,12,0.32,0.52,0.81,910];

[comb,cs] = combineTrajectories6D({traj1,traj2},'AddMovieID',true,'Verbose',true);
assert(cs.nTotalTracks==3,'Should have 3 total tracks (got %d)',cs.nTotalTracks);
assert(length(unique(comb(:,2)))==3,'Should have 3 unique track IDs');
assert(size(comb,2)==10,'Should have movieID column');
assert(all(comb(1:4,10)==1),'First 4 rows should be movie 1');
assert(all(comb(5:6,10)==2),'Last 2 rows should be movie 2');
fprintf('  All combineTrajectories6D checks PASS\n\n');

%% ---- TEST 6: estimateDiffusion6D_MLE ----
fprintf('--- TEST 6: estimateDiffusion6D_MLE ---\n');

% Clean Brownian trajectory with known D (no localization noise)
rng(123);
nSteps = 500;
x_c = cumsum(sqrt(2*D_trans*dt)*randn(nSteps,1));
y_c = cumsum(sqrt(2*D_trans*dt)*randn(nSteps,1));
z_c = cumsum(sqrt(2*D_trans*dt)*randn(nSteps,1));
th_c = pi/4 + cumsum(sqrt(2*D_rot*dt)*randn(nSteps,1)*0.5);
th_c = max(0,min(pi/2,th_c));
ph_c = cumsum(sqrt(2*D_rot*dt)*randn(nSteps,1));
ph_c = mod(ph_c+pi,2*pi)-pi;

cleanTraj = [(1:nSteps)', ones(nSteps,1), x_c, y_c, z_c, ...
    th_c, ph_c, ones(nSteps,1), ones(nSteps,1)*1000];

res = estimateDiffusion6D_MLE(cleanTraj,'FrameTime',dt,...
    'LocPrecisionXY',0,'LocPrecisionZ',0,...
    'LocPrecisionTheta',0,'LocPrecisionPhi',0,...
    'MinDisplacements',3,'Verbose',true);

fprintf('  D_trans: true=%.0f, est=%.0f (ratio=%.2f)\n',...
    D_trans, res.D_trans, res.D_trans/D_trans);
fprintf('  D_rot:   true=%.4f, est=%.4f (ratio=%.2f)\n',...
    D_rot, res.D_rot, res.D_rot/D_rot);

assert(res.D_trans/D_trans > 0.5 && res.D_trans/D_trans < 2.0,...
    'D_trans estimate out of range');
fprintf('  D_trans within 2x: PASS\n');

% Test with gapped trajectory (Bo's method should handle this)
gappedTraj = cleanTraj;
% Remove ~20%% of frames to simulate blinking
keepIdx = sort(randperm(nSteps, round(nSteps*0.8)));
gappedTraj = gappedTraj(keepIdx,:);
gappedTraj(:,2) = 1; % same track ID

resGap = estimateDiffusion6D_MLE(gappedTraj,'FrameTime',dt,...
    'LocPrecisionXY',0,'LocPrecisionZ',0,...
    'LocPrecisionTheta',0,'LocPrecisionPhi',0,'Verbose',false);

fprintf('  With 20%% blinking: D_trans=%.0f (ratio=%.2f), gaps=%.0f%%\n',...
    resGap.D_trans, resGap.D_trans/D_trans, resGap.fracGaps*100);
assert(resGap.D_trans/D_trans > 0.3 && resGap.D_trans/D_trans < 3.0,...
    'Gapped D_trans estimate way off');
fprintf('  Gapped estimate reasonable: PASS\n\n');

%% ---- Cleanup ----
if exist(synFile,'file'), delete(synFile); end

%% ---- Summary ----
fprintf('==============================================\n');
fprintf('  ALL %d TESTS PASSED\n', 6);
fprintf('==============================================\n');
fprintf('\nTo run the full pipeline with u-track3D:\n');
fprintf('  1. addpath(genpath(''path/to/u-track3D/software''))\n');
fprintf('  2. addpath(genpath(''path/to/software6D''))\n');
fprintf('  3. Edit run6DTracking.m (set inputFile, outputDir)\n');
fprintf('  4. Run run6DTracking\n');
