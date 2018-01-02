for i=1:10
    A = dlmread(['D:\Science\Nonlinear refractive index meaurement by temporal phase reconstruction\FROG\FROG Extended\output_svd\'  num2str(i)  '.txt']);
    plot(A(:,2))
    hold on
end