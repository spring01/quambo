function cartesian = ZMatrixToCartesian(zMatrix)
xyz = zeros(size(zMatrix,1) + 2, 3);
xyz(1:2, [3 1]) = [1 1; 0 1];
for i = 2:size(zMatrix,1)
    ref = zMatrix(i, [2 4 6]);
    bond = zMatrix(i, 3);
    if(i == 2)
        ref(2) = i - 2;
        angle = 90;
    else
        angle = zMatrix(i, 5);
    end
    if(i <= 3)
        ref(3) = i - 3;
        dihedral = 0;
    else
        dihedral = zMatrix(i, 7);
    end
    xyz(i+2, :) = xyz(ref(1)+2, :);
    xyzRef21 = xyz(ref(2)+2, :) - xyz(ref(1)+2, :);
    xyzRef21 = xyzRef21 ./ norm(xyzRef21);
    xyzRef32 = xyz(ref(3)+2, :) - xyz(ref(2)+2, :);
    xyzRef32 = xyzRef32 - xyzRef21 * xyzRef32' .* xyzRef21;
    if(norm(xyzRef32) > 0)
        xyzRef32 = xyzRef32 ./ norm(xyzRef32);
    end
    xyz(i+2, :) = xyz(i+2, :) + bond .* ( ...
        cosd(angle) .* xyzRef21 + ...
        sind(angle) .* (cosd(dihedral).*xyzRef32 ...
        - sind(dihedral) .* cross(xyzRef21, xyzRef32)));
end
cartesian = [zMatrix(:, 1) xyz(3:end, :)];
end