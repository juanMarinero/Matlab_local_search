function [xopt,fopt,iter,cont,vec_x_opt,vec_f_opt]=checkMethods(Fopt,x0,method,par,N_1cont)
tic
if 0>1
    
    clc,clear all,close all
    Fopt=@(x) 10*2+(x(:,1).^2-10*cos(2*pi*x(:,1)))+(x(:,2).^2-10*cos(2*pi*x(:,2)));  %equivale a: 10*2+sum(x(:).^2-10*cos(2*pi*x(:))); pero en el mesh de esta manera da problemas
    
    % method fminsearch (method 1)
    % built-in method
    method=1;
    x0=[0.4665,2.5];
    [xopt,fopt]=checkMethods(Fopt,x0,method);
    
    % method Climbing_Hill (method 2)
    method=2;
    x0=[0.4665,2.5];
    N=1e2;Delta=0.01;
    par=struct('N', N, 'Delta', Delta);% par es opcional
    [xopt,fopt,iter,cont,vec_x_opt,vec_f_opt]=checkMethods(Fopt,x0,method,par);
   
    % method Multi_Start (method 3)
    method=3;
    x0=[-5,5];
    N=1e2;Delta=0.01;
    par=struct('N', N, 'Delta', Delta);% par es opcional e inidica lo mismo q en Climbing_Hill
    N_1cont=1e2; % nº de ptos aleatorios de 1º vez
    [xopt,fopt,iter,cont,vec_x_opt,vec_f_opt]=checkMethods(Fopt,x0,method,par,N_1cont);
    
    % method Montecarlo
    method=4;
    x0=[-5,5];
    N=1e3;
    par=struct('N', N);
    [xopt,fopt]=checkMethods(Fopt,x0,method,par);

end

% VARIABLES DE ENTRADA
% method: 'Climbing_Hill','f_minsearch' or 'Multi_Start'

% Fopt: handle fcn

% Climbing_Hill y fminsearch
% x0: es un pto de D-1 dimensiones
% si Fopt es 2D, x0 es vector(1x1) q indica abcisas(eje1) de cada pto, y Fopt(xopt) son ordenadas(eje2)
% si Fopt es 3D, x0 es matriz(1x2) q indica eje1 y eje2 de cada pto, y Fopt(xopt) eje3
% ...
% ej. si Fopt es 4D, x0 puede ser x0=[0,0,0,0]'

% Multi_Start
% x0 es para formar los aleatorios de 1º iteración, NO inidica 1º xopt,
% sino los límites para la aleatoriedad de ptos de 1º iteración
% x0(1) es lim inferior de x e y,
%x0(1) es lim inferior de x e y,

% par: estrucutra opcional con campos de los valores de N y Delta
% N: cantidad de puntos para evaluar alrededor del xopt de cada iter
% Delta: participa en el espaciado de ptos aleatorios, q estén mas o menos
% alejados del xopt anterior


% VARIABLES DE SALIDA
% xopt es un pto de D-1 dimensiones
% fopt es un pto de 1x1
% cont: nº de iteraciones (veces q se ejecuta for) necesarias para obtener fopt
% iter: nº de iteraciones totales (sum x evaluados totales={cont}*{nº ptos evaluados en cada cont}) necesarias para obtener fopt
% vec_x_opt: vector con todos los xopt obtenidos
% vec_f_opt: vector con todos los fopt obtenidos

if method==1%'f_minsearch'
    disp(['Method f_minsearch (method=',num2str(method),').']);
elseif method==2%'Climbing_Hill'
    disp(['Method Climbing_Hill (method=',num2str(method),').']);
elseif method==3%'Multi_Start'
    disp(['Method Multi_Start (method=',num2str(method),').']);
elseif method==4%'Multi_Start'
    disp(['Method Montercarlo (method=',num2str(method),').']);    
end



% 1º iteración (todo method tiene al menos una iteración):

    xopt=x0;
    fopt=Fopt(xopt);
    vec_x_opt(:,1)=xopt;
    vec_f_opt(1)=fopt;



if method==1
    
    [xopt,fopt,output]=fminsearch(Fopt,x0);% fminsearch starts at point x0
    
elseif method==4
    %Parametros por defecto
    if isfield(par,'N'),N=par.N; end
    x=rand(N,2);
    
    x=-5+x*10;
    f=Fopt(x);
    
    [fmin,ind]=min(f);
    xmin=x(ind,:);
    xopt=xmin;
    fopt=fmin;
    vec_x_opt(:,1)=xmin;
    vec_f_opt(1)=fmin;
    
elseif method==2 | method==3
    
    % PARÁMETROS
    %Comprobar entradas y marcar errores
    if nargin <2,
        error ('[xopt,fopt]=hill(Fopt,x0,opt)');
    elseif nargin==2,
        par=[];
    end
    
    %Parametros por defecto
    N=1e2;
    Delta=0.01;
    
    % ó parametros del usuario que sustituyen a los de defecto
    % tomamos valor de variables de fields de la estructura par homónima a variable
    % isfiled(strucuture,'fieldname') may be 0 or 1
    if isfield(par,'N'),N=par.N; end
    if isfield(par,'Delta'),Delta=par.Delta; end
    
    
    %% Inicializacion de variables
    iter=0;
    
    dim= length (x0);
    % como expliqué x0 es un vector (1 x dim), donde dim es el número de
    % dimensiones de las que depende Fopt, ej. si Fopt es una curva (2D), el
    % eje x se lo da x0, luego, x0 sería un vector(1 x dim) con dim=1, así dim=D-1
    % si Fopt es una plano (3D), x0 debe dar eje x e y, x0(1 x dim) con dim=2
    
    
    
    
    % bucle evaluador
    cont=1;
    
    f= zeros (N,1);
    while 1,
        
        % BREVE REPASO DE STRUCTURE
        % 1: mira strucutre de BasicoMatlab
        % 2: v(a,b).filed1=k; v(a,b).filed2='k'; v(c,d).field1='k1'; y así hasta
        % que quiera, v sería una matriz de estrucutras <row x columns struct>,
        % o sea, q es formato structure, en que cada celda contiene una
        % estructura <1x1 struct>
        % 3: "v.fieldx" muestra los valores del campo filedx de cada celda
        % de matriz x, si una posicion de matriz, ej (m,n), no tiene ese field
        % saldría simplemente v(m,n).fieldx >> []
        
        cont=cont+1;
        
        %Nuevos puntos para evaluar 
        
        if method==2 | (method==3 && cont>2)
            
            x=xopt(ones(N,1),:).*(1+Delta*(rand(N,dim)-0.5));
            % 1: xopt(ones(N,1),:)es para tener todas filas iguales, ej.
            % a=1:5,n=3,repmat(a,[n,1]),ones(n,1)*a,a(ones(n,1),:) % son 3 métodos equivalentes
            % 2:
        elseif method==3 % &&cont==2
            if cont==2
                
                M=rand(N_1cont,dim);
                y1=x0(1);
                y2=x0(2);
                
                v1=M(:,1);
                x1=min(v1);
                x2=max(v1);
                
                x(:,1)=y1+(y2-y1)/(x2-x1)*v1;% escalado
                
                v2=M(:,2);
                x1=min(v2);
                x2=max(v2);
                
                x(:,2)=y1+(y2-y1)/(x2-x1)*v2;% escalado
            end
        end
        
        
        
        %Valor de la funcion en puntos (vectorial)
        for k=1:N,f(k)=Fopt(x(k,:)); end
        
        iter=iter+N; %Estadistica del algoritmo
        
        [fmin,I]= min (f); %Mejor obtenido
        
        %Cambio a mejor o salida
        if fmin<fopt,
            fopt= fmin ;xopt=x(I,:);
            
            vec_x_opt(:,cont)=xopt;
            vec_f_opt(cont)=fopt;
        else
            break ;
        end
    end
    if method==3
        vec_x_opt(:,1)=vec_x_opt(:,2);
        vec_f_opt(1)=vec_f_opt(2);
    end
    % plot de cada xopt obtenido, y cada fopt obtenido, y comparaciones de su evolución
    
    n=length(vec_f_opt)-1;
    m=n*.4; %pa los vectores var solo comparo m veces (separadas igualmente, linspace)
    z=round(linspace(n,2,m));
    if length(z)>1
        h=z(1)-z(2);
        var_x=zeros(1,length(z));var_y=var_x;var_xy=var_x;var_f_opt=var_x;
        cont=0;
        for k=round(linspace(n,2,m))
            cont=cont+1;
            if k-h>0
                var_x(cont)=(vec_x_opt(1,k)-vec_x_opt(1,k-h));
                var_y(cont)=(vec_x_opt(2,k)-vec_x_opt(2,k-h));
                var_xy(cont)=sqrt((var_x(cont))^2+(var_x(cont))^2);
                var_f_opt(cont)=(vec_f_opt(k)-vec_f_opt(k-h));
            end
        end
        % invertir posiciones del vector
        var_x=fliplr(var_x);
        var_y=fliplr(var_y);
        var_xy=fliplr(var_xy);
        var_f_opt=fliplr(var_f_opt);
        
        figure(3)
        % ejex
        subplot(4,2,1),plot(vec_x_opt(1,:),'.'),title('x','fontweight','bold','fontsize',10);grid on
        ylim([min(vec_x_opt(1,:)),max(vec_x_opt(1,:))])
        % varaiación eje x
        subplot(4,2,2),plot(var_x,'.'),title('var_x','fontweight','bold','fontsize',10);grid on
        if min(var_x)<max(var_x),ylim([min(var_x),max(var_x)]);end
        % eje y
        subplot(4,2,3),plot(vec_x_opt(2,:),'.'),title('y','fontweight','bold','fontsize',10);grid on
        ylim([min(vec_x_opt(2,:)),max(vec_x_opt(2,:))])
        % variación del eje y
        subplot(4,2,4),plot(var_y,'.'),title('var_y','fontweight','bold','fontsize',10);grid on
        if min(var_y)<max(var_y), ylim([min(var_y),max(var_y)]);end
        % eje x e eje y
        subplot(4,2,5)
        hold on
        plot(vec_x_opt(1,:),vec_x_opt(2,:),'.'),title('x,y','fontweight','bold','fontsize',10);grid on
        plot(vec_x_opt(1,end),vec_x_opt(2,end),'r.','MarkerSize',10)
        text(vec_x_opt(1,end),vec_x_opt(2,end),'\leftarrow[xopt,yopt]','fontweight','bold','Rotation',45)
        hold off
        xlim([min(vec_x_opt(1,:)),max(vec_x_opt(1,:))])
        ylim([min(vec_x_opt(2,:)),max(vec_x_opt(2,:))])
        % variación de eje x e y
        subplot(4,2,6),plot(var_xy,'.'),title('var_{xy} (abs value)','fontweight','bold','fontsize',10);grid on
        if 0.95*min(var_xy)<1.05*max(var_xy),ylim([0.95*min(var_xy),1.05*max(var_xy)]);end
        % fopt
        subplot(4,2,7),plot(vec_f_opt,'.'),title('f_{opt}','fontweight','bold','fontsize',10);grid on
        ylim([min(vec_f_opt),max(vec_f_opt)])
        % variación de fopt
        subplot(4,2,8),plot(var_f_opt,'.'),title('var_{f_{opt}}','fontweight','bold','fontsize',10);grid on
        if min(var_f_opt)<max(var_f_opt),ylim([min(var_f_opt),max(var_f_opt)]);end
    end
    
end

%%
figure(1)
hold on
plot3_mesh_or_surf(Fopt)

% plot3 fopt final
figure(1)
hold on
n=50;a=linspace(vec_f_opt(1,1),vec_f_opt(1,1)+100,n);
b=ones(1,n)*vec_x_opt(1,1);c=ones(1,n)*vec_x_opt(2,1);
plot3(b,c,a,...
    '*b','MarkerSize',2,'LineWidth',2);
plot3(xopt(1),xopt(2),fopt*1.5,...
    'o','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot3([xopt(1),xopt(1)],[xopt(2),xopt(2)],[fopt,fopt+70],...
    'color','k','LineWidth',2);
%     mArrow3([xopt(1),xopt(2),3*h+fopt],[xopt(1),xopt(2),fopt],...
%         'color','k','stemWidth',var_stemWidth,'facealpha',var_facealpha);

titulo={sprintf('Final resault in black:\n\nx_{\ropt}=[ %d , %d ]\n\nf_{\ropt}= %d      aprox %1.3d ',xopt(1),xopt(2),fopt,round(fopt*1e5)/1e5)};
title(titulo,'fontweight','bold','fontsize',20);
view([85,65])
hold off

figure(2)
hold on
if method==1 | method==2
    plot(vec_x_opt(1,1),vec_x_opt(2,1),...
        'o','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k')
    
end
plot(xopt(1),xopt(2),...
    'o','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor','k')

titulo={sprintf('Final resault in black:\n\nx_{\ropt}=[ %d , %d ]\n\nf_{\ropt}= %d      aprox %1.3d ',xopt(1),xopt(2),fopt,round(fopt*1e5)/1e5)};
title(titulo,'fontweight','bold','fontsize',20);
hold off

if method==2 | method==3
    flechas(vec_x_opt,vec_f_opt,cont)
end

toc
return


function plot3_mesh_or_surf(Fopt)
% Fopt=@(x) 10*2+(x(:,1).^2-10*cos(2*pi*x(:,1)))+(x(:,2).^2-10*cos(2*pi*x(:,2)));  %equivale a: 10*2+sum(x(:).^2-10*cos(2*pi*x(:))); pero en el mesh de esta manera da problemas

figure(1)
hold on
% newStr = regexprep(str,expression,replace)
Fopt=func2str(Fopt);
Fopt= regexprep(Fopt,':',':,:');
Fopt=str2func(Fopt);

lim=5.12;
x(:,1)=linspace(-lim,lim,1e2);
x(:,2)=x(:,1);

v(:,:,1)=repmat(x(:,1),[1,length(x(:,1))])';
v(:,:,2)=repmat(x(:,2),[1,length(x(:,2))]);
a=x(:,1);
b=x(:,2);
x=v;
% plot(x(:,:,1),x(:,:,2),'.')

f=Fopt(x);
mesh(a,b,f),hold on,plot3(x(:,:,1),x(:,:,2),f,'.','MarkerSize',4.55),grid on,hold off
xlim([-lim,lim]),ylim([-lim,lim])
hold off

%%
figure(2)
contour(a,b,f)
colorbar


% está pa q yo lo vea mediante "Run Section", y el resto del codigo
% funciona perfectamente sin esta parte
if 0>1
    
    % algoritmo Climbing-Hill
    clc,clear all,close all
    Fopt=@(x) 10*2+...
        (x(:,:,1).^2-10*cos(2*pi*x(:,:,1)))+...
        (x(:,:,2).^2-10*cos(2*pi*x(:,:,2)));
    x0=[1,1];
    lim=5.12;
    x(:,1)=linspace(-lim,lim,1e2);
    x(:,2)=x(:,1);
    
    v(:,:,1)=repmat(x(:,1),[1,length(x(:,1))])';
    v(:,:,2)=repmat(x(:,2),[1,length(x(:,2))]);
    a=x(:,1);
    b=x(:,2);
    x=v;
    %     plot(x(:,:,1),x(:,:,2),'.')
    
    f=Fopt(x);
    plot3(x(:,:,1),x(:,:,2),f,'.','MarkerSize',4.55),grid on
    xlim([-lim,lim]),ylim([-lim,lim])
    
    %%
    % surf y mesh    % aqui griddata no es necesario al ser f ya matriz
    
    figure(2)
    surf(a,b,f),grid on
    xlim([-lim,lim]),ylim([-lim,lim])
    figure(3)
    mesh(a,b,f),hold on,plot3(x(:,:,1),x(:,:,2),f,'.','MarkerSize',4.55),grid on,hold off
    xlim([-lim,lim]),ylim([-lim,lim])
    
end

return

function flechas(vec_x_opt,vec_f_opt,cont)
%% plot3 & mArrow3 de cada fopt que kiera
% necesario el fichero mArrow3.m
% http://www.mathworks.com/matlabcentral/fileexchange/25372-marrow3-m-easy-to-use-3d-arrow

%path pal mArrow3.m
addpath C:\Users\trini\Desktop\juan\Matlab_Mathematica_LaTeX\MATLAB\MATLAB_trabajos\functions,

R=.3;
h=25;var_stemWidth=0.02;% grosor de flechas
var_facealpha=.5;% grado de transparencia flechas
n=20;%(length(vec_x_opt)-1)
n=min(n,length(vec_x_opt)-1);
col=hsv(n);
k_col=0;
figure(1)
for k=round(linspace(1,max(.99*cont,n),n));% otra opcion podría ser concentrarse en los últimos puntos: k=(cont-n):(cont-1);
    
    if k>length(vec_x_opt);break;end
    k_col=k_col+1;
    xi=vec_x_opt(1,k);
    yi=vec_x_opt(2,k);
    zi=vec_f_opt(k);
    
    
    grad=repmat(linspace(0,2*pi,n),[1,ceil(cont/n)]);
    xf=xi+R*cos(grad(k));
    yf=yi+R*sin(grad(k));
    zf=zi+h;
    
    figure(2)
    hold on
    plot(xi,yi,...
        'o','MarkerSize',10,'MarkerFaceColor',col(k_col,:),'MarkerEdgeColor',col(k_col,:))
    hold off
    
end


return

