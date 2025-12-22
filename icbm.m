clc; clear; close all;

%% ================= CONSTANTS =================
Re  = 6371;                 % Earth radius (km)
mu  = 398600;               % Earth GM (km^3/s^2)

%% ================= USER PARAMETERS =================
Nmc        = 40;
dt         = 1.0;
Tmax       = 4200;
hit_radius = 5;
h_reentry  = 100;
Rmax_I     = 4000;
Npn        = 4;

launch_ECI = [Re; 0];

%% ================= STORAGE =================
success = 0;
miss_all = NaN(Nmc,1);
ICBM_paths = cell(Nmc,1);
INT_paths  = cell(Nmc,1);

%% ================= MONTE CARLO =================
for mc = 1:Nmc

    %% ---- ICBM INITIAL POSITION (SCATTERED) ----
    r0 = [Re+7000+500*randn; 6000+500*randn];

    %% ---- LAMBERT TARGETING (STRICT) ----
    Timpact = 3600 + 300*randn;     % seconds
    [v0, ~] = lambert2D(r0, launch_ECI, Timpact, mu);

    rT = r0;
    vT = v0;

    %% ---- INTERCEPTOR ----
    rI = launch_ECI;
    vI_mag = 9.8;
    vI = [0;0];
    launched = false;

    rT_hist = [];
    rI_hist = [];

    %% ---- TIME LOOP ----
    for t = 0:dt:Tmax

        %% === TRUE 2‑BODY ICBM ===
        aT = -mu * rT / norm(rT)^3;
        vT = vT + aT*dt;
        rT = rT + vT*dt;
        rT_hist(:,end+1) = rT;

        alt = norm(rT) - Re;
        if alt < 0
            break
        end

        %% === RELATIVE GEOMETRY ===
        r_rel = rT - rI;
        Rrel = norm(r_rel);
        v_rel = vT - vI;

        %% === LAUNCH GATE ===
        if ~launched && alt > h_reentry && Rrel < Rmax_I
            launched = true;
            vI = vI_mag * r_rel/Rrel;
        end

        %% === PN GUIDANCE ===
        if launched
            lambda_dot = (r_rel(1)*v_rel(2) - r_rel(2)*v_rel(1))/(Rrel^2);
            Vc = -dot(r_rel,v_rel)/Rrel;
            a_lat = Npn * Vc * lambda_dot;

            v_hat = vI/norm(vI);
            n_hat = [-v_hat(2); v_hat(1)];
            vI = vI + a_lat*n_hat*dt;
            vI = vI_mag * vI/norm(vI);
            rI = rI + vI*dt;
        end
        rI_hist(:,end+1) = rI;

        %% === INTERCEPT ===
        if launched && Rrel < hit_radius
            success = success + 1;
            break
        end
    end

    %% === MISS DISTANCE ===
    L = min(size(rT_hist,2),size(rI_hist,2));
    miss_all(mc) = min(vecnorm(rT_hist(:,1:L) - rI_hist(:,1:L)));

    %% === LOCAL TANGENT FRAME ===
    ICBM_paths{mc} = rT_hist - launch_ECI;
    INT_paths{mc}  = rI_hist - launch_ECI;
end

%% ================= STATS =================
fprintf('\nENGAGEMENT STATISTICS ===\n');
fprintf('POI                : %.2f\n',success/Nmc);
fprintf('Mean miss distance : %.2f km\n',mean(miss_all,'omitnan'));
fprintf('Min miss distance  : %.2f km\n',min(miss_all,[],'omitnan'));

%% ================= PLOT =================
figure; hold on; grid on; axis equal;
for i = 1:Nmc
    plot(ICBM_paths{i}(1,:),ICBM_paths{i}(2,:),'b','Color',[0.2 0.4 1 0.3]);
    plot(INT_paths{i}(1,:),INT_paths{i}(2,:),'r','Color',[1 0.3 0.3 0.3]);
end

theta = linspace(0,2*pi,400);
plot(Rmax_I*cos(theta),Rmax_I*sin(theta),'k--','LineWidth',1.5);
yline(h_reentry,'k:','Re‑entry boundary');
scatter(0,0,80,'k','filled');

xlabel('Downrange (km)');
ylabel('Altitude (km)');
title('Strict Keplerian ICBM Intercept (Lambert‑Targeted)');
legend('ICBM','Interceptor','Interceptor range','Re‑entry','Launch');

function [v1, v2] = lambert2D(r1, r2, tof, mu)

    r1n = norm(r1); r2n = norm(r2);
    cosd = dot(r1,r2)/(r1n*r2n);
    dtheta = acos(max(min(cosd,1),-1));

    A = sin(dtheta)*sqrt(r1n*r2n/(1-cosd));

    z = 0;
    for k = 1:1000
        [S,C] = stumpff(z);
        y = r1n + r2n + A*(z*S-1)/sqrt(C);
        F = (y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*tof;
        dF = (y/C)^1.5*(0.5/C*(C-1.5*S/C)+1.5*S^2/C) + A*(0.5/y)*sqrt(y);
        z = z - F/dF;
    end

    f = 1 - y/r1n;
    g = A*sqrt(y/mu);
    gdot = 1 - y/r2n;

    v1 = (r2 - f*r1)/g;
    v2 = (gdot*r2 - r1)/g;
end

function [S,C] = stumpff(z)
    if z > 0
        S = (sqrt(z)-sin(sqrt(z)))/sqrt(z^3);
        C = (1-cos(sqrt(z)))/z;
    elseif z < 0
        S = (sinh(sqrt(-z))-sqrt(-z))/sqrt((-z)^3);
        C = (cosh(sqrt(-z))-1)/(-z);
    else
        S = 1/6;
        C = 1/2;
    end
end