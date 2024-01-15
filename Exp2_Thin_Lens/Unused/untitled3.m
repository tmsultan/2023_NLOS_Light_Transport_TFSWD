z = [0.1, 1,5,10]
x = -5:0.01:5

[X, Z] = meshgrid(x,z); 
rsq = X.^2 + Z.^2;


I = 1./rsq; 

% figure; imagesc(x, z, I)
% figure; plot(I(100,:))

figure; plot(x, (I'))
xlabel('\rho (m)'); ylabel('1/r^2')
legend([ repmat('z=', [4,1]), num2str(z'), repmat('m', [4,1])])

figure; plot(x, log(I'))
xlabel('\rho (m)'); ylabel('log(1/r^2)')
legend([ repmat('z=', [4,1]), num2str(z'), repmat('m', [4,1])])