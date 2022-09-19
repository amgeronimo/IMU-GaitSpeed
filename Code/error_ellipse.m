function [handle] = error_ellipse(data1,data2,plot_data)
data = [data1 data2];

% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);

% Get the 95% confidence interval error ellipse
 chisquare_val = 2.4477;


theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);


x = data(:,1); y = data(:,2);
handle.outlier = (((x-X0).*cos(phi)+(y-Y0).*sin(phi)).^2/a^2) + (((x-X0).*sin(phi)-(y-Y0).*cos(phi)).^2/b^2)>1;


C = spring;
if plot_data > 0
    handle.e = ellipse(a,b,phi,X0,Y0,C);
    handle.e.LineWidth = .5;
    if plot_data > 1
        hold on;
        % Plot the original data - all individual points, input plot_data=2
        handle.d = plot(data(:,1), data(:,2), '.','Color',handle.e.Color);

    end
end
