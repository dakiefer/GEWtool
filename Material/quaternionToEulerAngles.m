function [theta1, theta2, theta3] = quaternionToEulerAngles(q, i, j, k)
% quaternionToEulerAngles - Convert Cartesian coordinates on a (hpyer-)sphere to Euler angles.
%
% Arguments: 
% - q          (Nx4 real) N points in Cartesian coordinates on the unit-hypersphere in R4.
%              It represents the final position after rotation and can be interpreted as 
%              quaternions.
% - i, j, k    (one of {1, 2, 3}, where i ≠ j, j ≠ k) three axes of rotation
% 
% Output: 
% - theta1, theta2, theta3:   (scalar real) Euler angles for rotation sequence
%                             indicated by i,j and k.
%
% Algorithm from 
% E. Bernardes and S. Viollet, “Quaternion to Euler angles conversion: A direct,
% general and computationally efficient method,” PLOS ONE, vol. 17, no. 11, p.
% e0276302, Oct. 2022, doi: 10.1371/journal.pone.0276302.
% 
% 2024 - Gatien Clement, Institut Langevin, Paris, France

    if i == k
        notproper = false;
        k = 6 - i - j; % because i + j + k = 1 + 2 + 3 = 6
    else
        notproper = true;
    end
    
    epsilon = (i - j) * (j - k) * (k - i) / 2; % equivalent to the Levi-Civita symbol
    
    if notproper
        a = q(:,1) - q(:,j + 1);
        b = q(:,i + 1) + q(:,k + 1) * epsilon;
        c = q(:,j + 1) + q(:,1);
        d = q(:,k + 1) * epsilon - q(:,i + 1);
    else
        a = q(:,1);
        b = q(:,i + 1);
        c = q(:,j + 1);
        d = q(:,k + 1) * epsilon;
    end
    
    theta2 = acos(2 * (a.^2 + b.^2) ./ (a.^2 + b.^2 + c.^2 + d.^2) - 1);
    theta_plus = atan2(b, a);
    theta_minus = atan2(d, c);%-pi;
    
  

    theta1 = 0.* (theta2 == 0)  +   0.*(theta2 == pi / 2)   +   (theta_plus - theta_minus).*not((theta2 == pi / 2)+(theta2 == 0));
    theta3 = (2 * theta_plus - theta1).*(theta2 == 0)   +   (2 * theta_minus + theta1).*(theta2 == pi / 2)  +   (theta_plus + theta_minus).*not((theta2 == pi / 2)+(theta2 == 0));
    
    if notproper
        theta3 = epsilon * theta3;
        theta2 = theta2 - pi / 2;
    end
    
    [theta1, theta3] = wrap(theta1, theta3); % "wrap" assures theta1, theta3 ∈ (-pi, pi] or theta1, theta3 ∈ [0, 2*pi)
end

function [theta1_wrapped, theta3_wrapped] = wrap(theta1, theta3)
    % Wraps angles to be within (-pi, pi] or [0, 2*pi)
    theta1_wrapped = mod(theta1 + pi, 2 * pi) - pi;
    theta3_wrapped = mod(theta3 + pi, 2 * pi) - pi;
end
