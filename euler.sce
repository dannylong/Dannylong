clear; 
function y = f_analitica(x)
    //y = 70/9 * exp(-0.3*x) - 43/9 * exp(-1.2*x);
    y = -.5*x.^4 + 4*x.^3 -10*x.^2 + 8.5*x + 1;
endfunction

function dy = f(x, y)
    //dy = -1.2*y + 7*exp(-0.3*x)
    dy = -2*x.^3 + 12*x.^2 - 20*x + 8.5;
endfunction

function [x, y] = euler_basic(f, a, b, h, x0, y0)
    k = 1;
    x(k) = x0;
    y(k) = y0;
    for i = (a+h):h:(b)
        x(k+1) = x(k) + h;
        y(k+1) = y(k) + f(x(k), y(k))*h;
        k = k + 1;
    end
endfunction

function [x, y] = euler_modified(f, a, b, h, x0, y0)
    k = 1;
    x(k) = x0;
    y(k) = y0;
    for i = (a+h):h:(b)
        
        x(k+1) = x(k) + h;
        
        slope_euler     = f(x(k), y(k));
        y_KP1_euler     = y(k) + slope_euler * h;
        slope_KP1_euler = f(x(k+1), y_KP1_euler);
        slope           = (slope_euler + slope_KP1_euler)/2;
         
        y(k+1) = y(k) + slope*h;
        k = k + 1;
    end
endfunction

function [x, y] = euler_midpoint(f, a, b, h, x0, y0)
    k = 1;
    x(k) = x0;
    y(k) = y0;
    for i = (a+h):h:(b)
        xH      = x(k) + h/2;
        yH      = y(k) + f(x(k), y(k)) * h/2;
        slope   = f(xH,yH);
        x(k+1)  = x(k) + h;
        y(k+1)  = y(k) + slope * h;
        k       = k + 1;
    end
endfunction

function [x, y] = euler_RK3(f, a, b, h, x0, y0)
    k = 1;
    x(k) = x0;
    y(k) = y0;
    for i = (a+h):h:(b)
        
        k1 = f(x(k)      , y(k)              );
        k2 = f(x(k) + h/2, y(k) + k1 * h/2   );
        k3 = f(x(k) + h  , y(k) - (k1-2*k2)*h);
        
        x(k+1)  = x(k) + h;
        y(k+1)  = y(k) + (1/6) * (k1+4*k2+k3) * h;
        k       = k + 1;
    end
endfunction

function [x, y] = euler_RK4(f, a, b, h, x0, y0)
    k = 1;
    x(k) = x0;
    y(k) = y0;
    for i = (a+h):h:(b)
        
        k1 = f(x(k)      , y(k)         );
        k2 = f(x(k) + h/2, y(k) + k1*h/2);
        k3 = f(x(k) + h/2, y(k) + k2*h/2);
        k4 = f(x(k) + h  , y(k) + k3*h  );
        
        x(k+1)  = x(k) + h;
        y(k+1)  = y(k) + (1/6) * (k1+2*k2+2*k3+k4) * h;
        k       = k + 1;
    end
endfunction

// Opções para a variável "method": euler_basic, euler_modified, euler_midpoint, euler_RK3 e euler_RK4
// Exemplos de uso:
//plot_EDO(euler_midpoint, f, 0, 4, .5, 0, 1, f_analitica, .1)    
//plot_EDO(euler_RK4, f, 0, 4, .5, 0, 1, f_analitica, .1)
function plot_EDO(method, f, a, b, h, x0, y0, f_analytics, h_analytics)
    [x, y] = method(f, a, b, h, x0, y0);
    t = a:h_analytics:b;
    plot(t, f_analytics(t))
    plot(x, y, 'r*')
    xgrid;
endfunction

function plot_EDO_all_methods(f, a, b, h, x0, y0, f_analytics, h_analytics)
    
    methods             = list(euler_basic, euler_modified, euler_midpoint, euler_RK3, euler_RK4);
    styleOfDots         = list('bo', 'ro', 'g*', 'bx', 'r*');
    legendDescriptions  = ['Solução exata'; 'Euler'; 'Euler Modificado'; 'Ponto Médio'; 'Runge-Kutta3'; 'Runge-Kutta4'];
    
    if (length(methods) == length(styleOfDots)) then
        t = a:h_analytics:b;
        plot(t, f_analytics(t))
        for (i=1:length(methods))
            [x, y] = methods(i)(f, a, b, h, x0, y0);
            plot(x, y, styleOfDots(i))
        end
        hl=legend(legendDescriptions);
        xgrid;
    end
    
endfunction

function dy = f_especifica(x, y)
    // Definir a função derivada específica para o problema
    // dy = ...;
    dy = x+y; // A expressão ( <-- )deve ser alterada!!!
endfunction

function y = f_analitica_especifica(x)
    // Definir a função analítica específica para o problema
    // y = ...;
    y = x;// A expressão deve ser alterada!!!
endfunction

// USO da função plot_EDO(f, a, b, h, x0, y0, f_analytics, h_analytics)
// Explicação dos parâmetros:
// f=f(x, y)=dy/dx; a=x_min; b=x_max; h=step; x0=initionDot; y0 = y(x0)
// f            : função derivada (idealmente, deve-se ser criada uma função própria para cada problema específico, modificando o modelo da função f_especifica)
// a            : limite inferior do intervalo do DOMÍNIO da função resposta (y(x))
// b            : limite superior do intervalo do DOMÍNIO da função resposta (y(x))
// h            : passo para cálculo dos pontos da função resposta (y(x)) no DOMÍNIO do intervalo [a, b]
// x0           : condição inicial para o valor de x (coordenada x do ponto de partida para plotagem do gráfico de y(x))
// y0           : condição inicial para o valor de y(x0) (coordenada y do ponto de partida para plotagem do gráfico de y(x))
// f_analytics  : função resposta analítica (idealmente, deve-se ser criada uma função própria para cada problema específico, modificando o modelo da função f_especifica) 


//SISTEMAS DE EQUAÇÕES DIFERENCIAIS
//Computar as EDOs
function dx1 = f1(t, x1, x2)
    //dy = (3/10) * x2 - (8/5) * x1; 
    dx1 = -.5*x1;
endfunction

function dx2 = f2(t, x1, x2)
    //dz = (8/5) * x1 - (4/5) * x2; 
    dx2 = 4 - .3*x2 - .1*x1;
endfunction

//Exemplo de uso da função resolveSystemODE:
// ICs '==' [x1(startPoint) x2(startPoint) ... xn(startPoint)]
function output = resolveSystemODERK4(listODEs, startPoint, endingPoint, step, ICs)
    systemOrder = length(listODEs);
    nRows       = length(startPoint:step:endingPoint);
    nColumns    = systemOrder + 1;
    output      = zeros(nRows, nColumns);
    
    // Assign the first row of the OUTPUT matriz
    output(1, 1)    = startPoint;
    for (i = 2:nColumns)
        output(1, i) = ICs(i-1);
    end
    // #Assign the first row of the OUTPUT matriz
    
    k1 = []; // k1s das n EDOs em ordem de ocorrência
    k2 = []; // k2s das n EDOs em ordem de ocorrência
    k3 = []; // k3s das n EDOs em ordem de ocorrência
    k4 = []; // k4s das n EDOs em ordem de ocorrência
    
    allKs = [k1; k2; k3; k4];
    k = 1;
    
    aux_K1 = [', '          , '0'               ,' * 1'        ];
    aux_K2 = [' + step/2, ' , 'allKs(1, iK)'    ,' * step/2'   ];
    aux_K3 = [' + step/2, ' , 'allKs(2, iK)'    ,' * step/2'   ];
    aux_K4 = [' + step, '   , 'allKs(3, iK)'    ,' * step'     ];
    aux_Kn = [aux_K1; aux_K2; aux_K3; aux_K4];
    
    for i = (startPoint + step):step:endingPoint // Número de iterações
        
        // Calculando todos os Ks
        for (b = 1:4)
            for (l = 1:systemOrder)
                text_aux = 'listODEs(l)( output(k, 1)' + aux_Kn(b, 1);
                iK = 1;
                for (m = 2:nColumns)
                    iK = m-1;
                    text_aux = text_aux + 'output(k, ' + string(m) + ') +' + string( evstr(aux_Kn(b, 2)) ) + aux_Kn(b, 3);
                    if (m <> nColumns)
                        text_aux = text_aux + ', ';
                    end
                end
                text_aux = text_aux + ')';
                allKs(b, l) = evstr(text_aux);
                //disp(allKs)
            end
        end
        
        // Calculando os Ts e XNs
        output(k+1, 1) = output(k, 1) + step;
        for (p = 1:systemOrder)
            slope6              = allKs(1, p) + 2*( allKs(2, p) + allKs(3, p) ) + allKs(4, p);
            output(k+1, p+1)    = output(k, p+1) + (1/6) * slope6 * step;
        end
        k = k + 1;
        
    end // # Número de iterações
    
endfunction

//Exemplo de uso da função resolveSystemODE:
//out = resolveSystemODEEulerB(list(f1, f2), 0, 2, .5, [4 6])
function output = resolveSystemODEEulerB(listODEs, startPoint,endingPoint, step, ICs)
    systemOrder = length(listODEs);
    nRows       = length(startPoint:step:endingPoint);
    nColumns    = systemOrder + 1;
    output      = zeros(nRows, nColumns);
    
    // Assign the first row of the OUTPUT matriz
    output(1, 1)    = startPoint;
    for (i = 2:nColumns)
        output(1, i) = ICs(i-1);
    end
    
    k = 1;
    
    for i = (startPoint + step):step:endingPoint
        output(k+1, 1) = output(k, 1) + step;

        for j = 2:nColumns
            text_aux = 'listODEs(j-1)(';
            for m = 1:nColumns
                text_aux = text_aux + 'output(k, ' + string(m) + ')';
                if (m <> nColumns)
                    text_aux = text_aux + ',';
                end
            end
            text_aux = text_aux + ')';
            
            output(k+1, j) = output(k, j) + evstr(text_aux) * step;
        end
        k = k + 1;
    end
endfunction

//Para o exemplo do último slide, usar no console o seguinte
//out1 = resolveSystemODEEulerB(list(f1, f2), 0, 2, .5, [4 6])
//out2 = resolveSystemODERK4(list(f1, f2), 0, 2, .5, [4 6])
