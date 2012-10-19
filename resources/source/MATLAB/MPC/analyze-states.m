[modalfc,T] = canon(fc, 'modal');


disp(' ')
disp('Eigenvalues of model in modal form:')
e = eig(modalfc);
for i=1:n
    display(num2str(e(i)))
end

surf(T)
title('State Transformation Matrix from Original States into States of Modal Form')
xlabel('Index of old states')
ylabel('Index of new states')
zlabel('Correlation')