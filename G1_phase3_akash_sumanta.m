clc;
clear all;

global q;
global epsilon1;
% global x2;
global delta;
global d;
global R;
global fun_eval;
% global x0;
% global s0;

%----------------------------------------------
%reading data from xls files
xp0=input('enter the initial guess column vector in format [1;2;3] : ');
q=readmatrix('input_file_phase3.xlsx','Sheet','Sheet1','Range','C3:C3'); %input('enter the function number: ');
epsilonp2=readmatrix('input_file_phase3.xlsx','Sheet','Sheet1','Range','C11:C11');
c=readmatrix('input_file_phase3.xlsx','Sheet','Sheet1','Range','C12:C12'); % r update parameter
epsilon1=readmatrix('input_file_phase3.xlsx','Sheet','Sheet1','Range','C4:C4'); %input('enter the 1st termination parameter: ');
epsilon2=readmatrix('input_file_phase3.xlsx','Sheet','Sheet1','Range','C5:C5'); %input('enter the 2nd termination parameter: ');
epsilon3=readmatrix('input_file_phase3.xlsx','Sheet','Sheet1','Range','C6:C6'); %input('enter the 3rd termination parameter: ');
% x2 =readmatrix('input_file_phase3.xlsx','Sheet','Sheet1','Range','C8:C8'); % Get the initial guess from the user for bounding phase
delta = readmatrix('input_file_phase3.xlsx','Sheet','Sheet1','Range','C9:C9');  % Get the increment value from the user for bounding phase
M = readmatrix('input_file_phase3.xlsx','Sheet','Sheet1','Range','C7:C7'); %input('enter maximun no. of iteration: ');
%----------------------------------------------------------------------------

%initializing required variables

d=size(xp0,1);
tp=1; % setting initial counter for penalty sequence.
R=0.1; % initalizing penalty parameter
fun_eval=0;
%-----------------------------------------------------------------------------

% penalty constrained optimizaton algorithm main code-------------------------------
while(1)
   fprintf('$$$$$$$$$$$$$$$$  Penalty sequencce no. %g  $$$$$$$$$$$$$$$$$$$$\n',tp);
   fprintf('\n here we have R=%g, epsilonp2=%g, c=%g, tp=%g \n xp0=\n,',R,epsilonp2,c,tp);disp(xp0');
   penalty=penfun(xp0,R,q); 
   fprintf('Sequence no. %g penalty is %g\n',tp,penalty);
   x0=xp0;
   xp_old=x0;

% conjugate gradient method code-------------------------------------
    disp('-------------- staring conjugate gradient menthod for new x0 ------');

    % gradient_x0=gradient_eval(x0,q);
    % disp('gradient at x0 is: '); disp(gradient_x0);
   
    s0=-1*gradient_eval(x0,q);   % initial direction of search vector
    s0=s0/norm(s0);
    alpha=unidirectional(x0,s0,q);
    fprintf('\n alpha value is: '); disp(alpha);
    
    x_new=x0+alpha*s0;
    gradient_x_new = gradient_eval(x_new,q);
    % disp('gradient_x_new is: '); disp(gradient_x_new);
    k=1; % setting initial counter as 1 for conjugate gradient.
    while(1)
        fprintf('\n-----------------------iternation %g of cgm-------------------------\n',k);
        s_old=s0/norm(s0);
        gradient_x0=gradient_eval(x0,q);
        s0 = -1*gradient_eval(x_new,q) + (norm(gradient_x_new)/norm(gradient_x0))^2*s0;
        s0=s0/norm(s0);
        
        alpha =unidirectional(x_new,s0,q);
        % fprintf('\n new alpha value is: '); disp(alpha);
        x0=x_new; 
        
        x_new=x0+alpha*s0;
        % checking linearly dependance here
        theta = acos(dot(s_old/norm(s_old), s0/norm(s0)));
        % fprintf('\n theta angle between search directions is: %g degree\n',theta*180/pi);
        
        gradient_x_new = gradient_eval(x_new,q);
        % fprintf('\n gradient_x_new is: '); disp(gradient_x_new);

        if ((norm(x_new-x0)/(norm(x0))) <= epsilon2)
            % disp('breaking in norm 1 of cgm');
            break;
        elseif(norm(gradient_x_new) <= epsilon3)
            % disp('breaking in 2 nd gradeitn of cgm');
            break;
        elseif(k >= M)
            % disp('breaking in exceeding of k value 3');
            break;
        else
            k=k+1;
            
        end
        fprintf('penalty value at the end of this cgm iteration is %g\n',penfun(x_new,R,q));
        % disp('gradient vlaue s0 is'); disp(s0');
        disp('x_new value in this iteration is '); disp(x_new');
    end
  
fprintf('\n The optimum point using cgm is : \n');disp(x_new');
fprintf('optimum function value using cgm is:');disp(penfun(x_new,R,q));

 % ------- conjugate gradient mehtod end ------------

   xp0=x_new;
   R_old=R;
   if tp>0
       if abs(penfun(xp0,R,q)-penfun(xp_old,R_old,q))<=epsilonp2
            % disp('we got constrained problem optimum point as'); disp(xp0');
            fprintf(' in sequence %g we got constrained problem optimum point as ',tp); disp(xp0');
            fprintf('the constrained problem optimum vlaue is %g ',funvalue(xp0,q));
            fprintf('\n Total no. of functions evaluations are %g ',fun_eval);
            break;
       end 
   end
   R=c*R;
   tp=tp+1;
end

% Penalty functions definitions ---------------------------------------- 
function f = funvalue(x,q)
    global fun_eval;
    fun_eval=fun_eval+1;
    % global q;
    if (q==1)  % Problem 1 
        f = (x(1)-10)^3 + (x(2)-20)^3;
    end
    if(q==2) %problem 2
        f= (sin(2*pi*x(1)))^3 * sin(2*pi*x(2)) / (x(1)^3 * (x(1) + x(2)));
    end
    if(q==3) %problem 3
        f= x(1)+x(2)+x(3);
    end
     if (q==4) % himmelblau function
        f= (x(1)^2+x(2)-11)^2 + (x(1)+x(2)^2-7)^2;
     end
end

function P = penfun(x,R,q)
    % global q
    global fun_eval;
    fun_eval=fun_eval+1;
    if (q==1)  % Problem 1 
        f = (x(1)-10)^3 + (x(2)-20)^3;
        P = f + R*(ohm_fun(x,q));
        
    end
    if(q==2) %problem 2
        f= -(sin(2*pi*x(1)))^3 * sin(2*pi*x(2)) / (x(1)^3 * (x(1) + x(2)));
        P = f + R*(ohm_fun(x,q));
    end
    if(q==3) %problem 3
       
        f= x(1)+x(2)+x(3);
        
        P = f + R*(ohm_fun(x,q));

    end
    if (q==4) % himmelblau function
        f= (x(1)^2+x(2)-11)^2 + (x(1)+x(2)^2-7)^2;
        
        P= f + R*(ohm_fun(x,q));

    end
end
function ohm = ohm_fun(x,q)
    % global fun_eval;
    if q==1 % Problem 1
    
        g(1) = min((x(1)-5)^2 + (x(2)-5)^2 - 100,0);
        g(2) = min(82.81 - ((x(1)-6)^2 + (x(2)-5)^2),0);
        g(3) = min((x(1)-13),0);
        g(4) = min((20-x(1)),0);
        g(5) = min(x(2),0);
        g(6) = min(4-x(2),0);
        ohm = sum(g.^2);
        % fun_eval=fun_eval+7;
        
    end
    if q==2  % Problem 2
    
        g(1) = min(-(x(1)^2 - x(2) + 1),0);
        g(2) = min(-(1 - x(1) + (x(2)-4)^2),0);
        g(3) = min(10-x(1),0);
        g(4) = min(10-x(2),0);
        g(5) = min(x(1),0);
        g(6) = min(x(2),0);
        ohm = sum(g.^2);
        % fun_eval=fun_eval+7;
    end
    
    if q==3 % Problem 3
        g(1) = min(1-0.0025*(x(4) + x(6)),0);
        g(2) = min(1-0.0025*(-x(4) + x(5) + x(7)),0);
        g(3) = min(1-0.01*(-x(6) + x(8)),0);
        g(4) = min(-(100*x(1) + (-x(6)*x(1)) + (833.33252*x(4)) + (-83333.333)),0);
        g(5) = min(-(x(2)*x(4) + (-x(2)*x(7)) + (-1250*x(4)) + (1250*x(5))),0);
        g(6) = min(-(x(3)*x(5) + (-x(3)*x(8)) + (-2500*x(5)) + 1250000),0);
        g(7) = min(10000-x(1),0);
        g(8) = min(10000-x(2),0);
        g(9) = min(10000-x(3),0);
        g(10) = min(x(1)-100,0);
        g(11) = min(x(2)-1000,0);
        g(12) = min(x(3)-1000,0);
        g(13) = min(1000-x(4),0);
        g(14) = min(1000-x(5),0);
        g(15) = min(1000-x(6),0);
        g(16) = min(1000-x(7),0);
        g(17) = min(1000-x(8),0);
        g(18) = min(x(4)-10,0);
        g(19) = min(x(5)-10,0);
        g(20) = min(x(6)-10,0);
        g(21) = min(x(7)-10,0);
        g(22) = min(x(8)-10,0);
        ohm = sum(g.^2);
        % fun_eval=fun_eval+23;
    end
end

% uppder and lower bound value picking
function x_lb = lb(q)
    if q==1
        x_lb= [13;0];
    end
    if q==2
        x_lb= [0;0];
    end
    if q==3
        x_lb=[100;1000;1000;10;10;10;10;10];
    end
    if q==4
        x_lb=[0;0];
    end
end

function x_ub = ub(q)
    if q==1
        x_ub= [20;4];
    end
    if q==2
        x_ub= [10;10];
    end
    if q==3
        x_ub=[10000;10000;10000;1000;1000;1000;1000;1000];
    end
    if q==4
        x_ub=[5;5];
    end
end
 %---x2 value for unidirectional and getting alpha-----------------------
        function x2 = uni_x2(x,s)
           
            global q;
            global d;
            
            x_lb=lb(q);
            x_ub=ub(q);
            alpha=zeros(2*d,1);
            for j=1:d
                alpha(j,1) = (x_ub(j,1)-x(j,1))/s(j,1);
                alpha(d+j,1) = (x_lb(j,1)-x(j,1))/s(j,1);
            end
            alpha=sort(alpha);
            alpha1 = alpha(d);
            alpha2 = alpha(d+1);
         
            x2 = alpha1 + (alpha2-alpha1)*rand(1);
        end
       % --------------------------------------------------------
% unidirectional search algorithm-----bounding phase----interval halving--
function alpha = unidirectional(x,s,q)
while(1)
    global fun_eval;
    x2 = uni_x2(x,s);
    delta=0.5;
    k1=0;
    fx1 = uni_f(x2-abs(delta),x,s,q);
    fx2 = uni_f(x2,x,s,q);
    fx3 = uni_f(x2+abs(delta),x,s,q);
    fun_eval=fun_eval+3;
    % fprintf('x2=%g,function values at fx1(x2-abs(delta))=%g, fx1(x2) =%g, fx1(x2+abs(delta))=%g: ',x2,fx1,fx2,fx3);
    if (fx1>=fx2 && fx2>=fx3)
        delta=1*delta;
        break;
    elseif (fx1<=fx2 && fx2<=fx3)
        delta=-1*delta;
        break;
    else
        disp('setting initial parameters agian');
        

    end
end
while(1)
    % fprintf('\n========iteration %g===================\n',k1+1);
    a10= x2-(2^(k1-1))*delta;
    x_old = x2;
    x2=x2+(2^k1)*delta;
    % fprintf('x2=%g and fx1(x2)=%g',x2,uni_f(x2,x,s,q));
    if uni_f(x2,x,s,q)<uni_f(x_old,x,s,q) %fx1(x2)<fx1(x_old)
        k1=k1+1;
        fun_eval=fun_eval+2;
    else
        b10=x2;
        a1=min(a10,b10);
        b1=max(a10,b10);
        fprintf('\nminima bracketed between %g to %g in iteration %g: ',a1,b1,k1);
        break;
    end
end

   
    global epsilon1;
    epsilon=epsilon1; %input('enter the solution accuracy for interval halving : ');
    p=0; %set interation number as zero
    while(1)  % 1 means our while will always run and the breaking of while loop depends on further conditions
        % fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
        p=p+1; % increase interation number by 1 after running each times while loop
        % fprintf('\n starting iteration number : %g',p);
        xm=(a1+b1)/2; %xm represents mean value
        L=b1-a1;
        fxm=uni_f(xm,x,s,q);
        x1=a1+L/4;
        x2=b1-L/4;
        fx11= uni_f(x1,x,s,q);
        fx22= uni_f(x2,x,s,q); 
        fun_eval=fun_eval+3;
        % fprintf('\n we have the values as a1:%g b1:%g xm: %g L: %g x1: %g x2:%g ',a1,b1,xm,L,x1,x2);
        % fprintf('\n we got the vlaues of function as fx1: %g  fx2:%g  fx:%g  ',fx1,fx2,fxm);
        if fx11<fxm 
            % fprintf('left side value is minimun so eliminating right of xm by choosing upper bound as xm')
            b1=xm;
            xm=x1;
        elseif fx22<fxm
            % fprintf('right side vlaue is minimum so eliminating left of xm by choosing lower bound as xm')
            a1=xm;   %choosing lower bound as xm
            xm=x2;  % storing x2 in xm  
        else
            a1=x1;
            b1=x2;     
        end
         if abs(L)<epsilon %breaking while loop when we meet the solution as per our requried approximation
             break;
         end
    end
    alpha=xm;
    fprintf('\n uisng interval halving optimum(alpha) as: %g in %g iterations \n',alpha,p);
end


% for using in unidirectional search making function f_alpha as function
% of alpha.---------------------------------------------
function f_alpha = uni_f(y,x,s,q)    %y is alpha value // we are making function in term of alpha
   
    global d; % d was number of elements in initial guess vector
    global R;
    global fun_eval;
    fun_eval=fun_eval+1;
    for i = 1:d
        x(i) = x(i) + y*s(i);
    end
    f_alpha = penfun(x,R,q);
end

% ----------------------gradient evaluation------------------
function gradient = gradient_eval(x,q) 
global d;
global R;
global fun_eval;
fun_eval=fun_eval+2;
gradient = zeros(d,1);
h = 0.001;
for i = 1:d
    y = x;
 y(i) = y(i)+h;
 a = penfun(y,R,q);
 y(i) = y(i)-2*h;
 b = penfun(y,R,q);
 gradient(i) = (a - b)/(2*h);
end
end
