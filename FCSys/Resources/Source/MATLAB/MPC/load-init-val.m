% Requires tnindex.m from the Dymola installation.

dos(['alist.exe -b ', cd, '\dsfinal.txt ', cd, '\dsfinal.mat']); % Convert final value data file from ASCII to MATLAB binary file.
load 'dsfinal.mat' % Load final value data file.

% Get the value of each model state.
x0 = zeros(n, 1);
for i=1:n % Step through all states.
    name = deblank(dslin.xuyName(i,:));
    ind = tnindex(initName, name);
    disp([name, '(0) = ', num2str(initVal(ind,2))])
    x0(i) = initVal(ind,2);
end
disp(' ')

% Get the value of each input.
u0 = zeros(m, 1);
for i=1:m % Step through all states.
    name = deblank(dslin.xuyName(n+i,:));
    ind = tnindex(initName, name);
    disp([name, '(0) = ', num2str(initVal(ind,2))])
    x0(i) = initVal(ind,2);
end
disp(' ')

% Get the value of each output.
y0 = zeros(p, 1);
for i=1:p % Step through all states.
    name = deblank(dslin.xuyName(n+m+i,:));
    ind = tnindex(initName, name);
    disp([name, '(0) = ', num2str(initVal(ind,2))])
    x0(i) = initVal(ind,2);
end
