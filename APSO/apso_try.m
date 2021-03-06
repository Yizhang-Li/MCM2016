% 改进的快速粒子群优化算法 (APSO):
function apso

Lb=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %下边界
Ub=[100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100]; %上边界
% 默认参数
para=[30 250 0.95]; %[粒子数，迭代次数，gama参数]

% APSO 优化求解函数
[gbest,fmin]=pso_mincon(@cost,@constraint,Lb,Ub,para);

% 输出结果
Bestsolution=gbest % 全局最优个体
fmin

%% 目标函数
function f=cost(x)
f=log((x(1)+x(2)+x(3)+x(4)+x(5)+51)*(x(6)+x(7)+x(8)+x(9)+x(10)+8)*(x(11)+x(12)+x(13)+x(14)+x(15)+53)*(x(16)+x(17)+x(18)+x(19)+x(20)+66));

% 非线性约束
function [g,geq]=constraint(x)
% 不等式限制条件
g=[];
% 如果没有等式约束，则置geq=[]；
geq(1)=100-x(1)-x(6)-x(11)-x(16);
geq(2)=100-x(2)-x(7)-x(12)-x(17);
geq(3)=100-x(3)-x(8)-x(13)-x(18);
geq(4)=100-x(4)-x(9)-x(14)-x(19);
geq(5)=100-x(5)-x(10)-x(15)-x(20);
%%  APSO Solver
function [gbest,fbest]=pso_mincon(fhandle,fnonlin,Lb,Ub,para)
if nargin<=10,
    para=[20 150 0.95];
end
n=para(1);% 粒子种群大小
time=para(2); %时间步长，迭代次数
gamma=para(3); %gama参数
 scale=abs(Ub-Lb); %取值区间
% 验证约束条件是否合乎条件
if abs(length(Lb)-length(Ub))>0,
    disp('Constraints must have equal size');
    return
end

  alpha=0.2; % alpha=[0,1]粒子随机衰减因子
  beta=0.5;  % 收敛速度(0->1)=(slow->fast);

% 初始化粒子群
best=init_pso(n,Lb,Ub);

fbest=1.0e+100;
% 迭代开始
for t=1:time,     
   
%寻找全局最优个体
  for i=1:n,   
    fval=Fun(fhandle,fnonlin,best(i,:)); 
    % 更新最有个体
    if fval<=fbest, 
        gbest=best(i,:);
        fbest=fval;
    end
        
  end

% 随机性衰减因子
 alpha=newPara(alpha,gamma);

% 更新粒子位置 
  best=pso_move(best,gbest,alpha,beta,Lb,Ub);  
 
% 结果显示
	str=strcat('Best estimates: gbest=',num2str(gbest));
	str=strcat(str,'  iteration='); str=strcat(str,num2str(t));
	disp(str);

    fitness1(t)=fbest;
    plot(fitness1,'r','Linewidth',2)
    grid on
    hold on
    title('适应度')
end

% 初始化粒子函数
function [guess]=init_pso(n,Lb,Ub)
ndim=length(Lb);
for i=1:n,
    guess(i,1:ndim)=Lb+rand(1,ndim).*(Ub-Lb); 
end

%更新所有的粒子 toward (xo,yo)
function ns=pso_move(best,gbest,alpha,beta,Lb,Ub)
% 增加粒子在上下边界区间内的随机性
n=size(best,1); ndim=size(best,2);
scale=(Ub-Lb);
for i=1:n,
    ns(i,:)=best(i,:)+beta*(gbest-best(i,:))+alpha.*randn(1,ndim).*scale;
end
ns=findrange(ns,Lb,Ub);

% 边界函数
function ns=findrange(ns,Lb,Ub)
n=length(ns);
for i=1:n,
  % 下边界约束
  ns_tmp=ns(i,:);
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % 上边界约束 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  
  %更新粒子
  ns(i,:)=ns_tmp; 
end

% 随机性衰减因子
function alpha=newPara(alpha,gamma);
alpha=alpha*gamma;

% 带约束的d维目标函数的求解
function z=Fun(fhandle,fnonlin,u)
% 目标
z=fhandle(u);

z=z+getconstraints(fnonlin,u); % 非线性约束

function Z=getconstraints(fnonlin,u)
% 罚常数 >> 1
PEN=10^15;
lam=PEN; lameq=PEN;

Z=0;
% 非线性约束
[g,geq]=fnonlin(u);

%通过不等式约束建立罚函数
for k=1:length(g),
    Z=Z+ lam*g(k)^2*getH(g(k));
end
% 等式条件约束
for k=1:length(geq),
   Z=Z+lameq*geq(k)^2*geteqH(geq(k));
end

% Test if inequalities 
function H=getH(g)
if g<=0, 
    H=0; 
else
    H=1; 
end

% Test if equalities hold
function H=geteqH(g)
if g==0,
    H=0;
else
    H=1; 
end