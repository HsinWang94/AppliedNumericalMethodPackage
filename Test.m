clc;clear
obj = GA_Optimzer( );
obj.LowBound = -1 * ones(1,30);
obj.HighBound = 1 * ones(1,30);
obj.obj_function = @ackley;
obj.N_Itr = 500;
tic
[x,xset]= obj.Normal_Optim();
toc

yset = zeros(obj.N_Itr,1);
for j = 1 : obj.N_Itr
    yset(j) = ackley(xset(j,:));
end
 plot(yset)