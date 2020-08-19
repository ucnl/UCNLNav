

function [ result ] = Nav_build_standard_ellipsoid(el_name)

if (strcmp(el_name, 'WGS72') ~= 0)
   result = Nav_build_ellipsoid(6378135.0, 298.26);
else
    if (strcmp(el_name, 'WGS84') ~= 0)
        result = Nav_build_ellipsoid(6378137.0, 298.257223563);
    else
        if (strcmp(el_name, 'GRS80') ~= 0)
            result = Nav_build_ellipsoid(6378137.0, 298.257222100882711);
        else
            if (strcmp(el_name, 'PZ90') ~= 0)
                result = Nav_build_ellipsoid(6378136.0, 298.257839303);
            else
                if (strcmp(el_name, 'IERS') ~= 0)
                    result = Nav_build_ellipsoid(6378136.49, 298.25645);
                else
                    if (strcmp(el_name, 'KRSY') ~= 0)
                        result = Nav_build_ellipsoid(6378245.0, 298.3);
                    else
                        if (strcmp(el_name, 'sphere') ~= 0)
                            result = Nav_build_ellipsoid(6378137.0, 1E16);
                        else
                            error('Unknown ellipsoid name "%s"', el_name);
                        end
                    end
                end
            end
        end
    end
end

end