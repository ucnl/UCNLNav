

function [ result ] = Nav_build_ellipsoid(mjsa_m, ifltn)

if (mjsa_m <= 0)
    error('Major semi-axis should be greater than zero');
end

if (ifltn <= 0)
    error('Inverse flattening should be greater than zero');
end

fltn = 1 / ifltn;
mnsa_m = mjsa_m * (1 - fltn);
ecnt_sq = (mjsa_m^2 - mnsa_m^2) / mjsa_m^2;
ecnt = sqrt(ecnt_sq);

result.mjsa_m = mjsa_m;
result.ifltn = ifltn;
result.fltn = fltn;
result.mnsa_m = mnsa_m;
result.ecnt = ecnt;
result.ecnt_sq = ecnt_sq;

end
