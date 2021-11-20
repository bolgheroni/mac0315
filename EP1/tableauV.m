
0;
% format long;
function [ind v ] = simplex_tab_intern(indB, tableau, iteracao)
    disp(["Iterando: ", int2str( iteracao)]);

    #display("Tableau:");
    display(tableau) 

    n = size(tableau)(2)-1;
    m = size(tableau)(1)-1;
    
    disp("Variaveis basicas:");
    for i = 1:m
        disp ([indB(i), tableau(1+i, 1)]);
    endfor
    #display("Instancia de x_b:");
    x_b = tableau(2:end, 1)

    disp(["Valor da funcao objetivo: ", num2str(-tableau(1, 1))]);

    c_t = tableau(1, 2:end);
    disp(["Vetor de custos reduzidos: ", mat2str(c_t)]);

    hasNegativeCost = false;
    for j = 1:n
        c_red = c_t(j);
        if (c_red < 0 && !any( indB == j) )
            hasNegativeCost = true;
            break
        endif
    endfor
    if(!hasNegativeCost)
        custo = -tableau(1, 1);
        ind = 1;
        v = zeros(1, n);
        for i = 1:m
            v( indB(i) ) = tableau(i+1, 1);
        endfor

        disp(strcat("Solucao otima encontrada com custo= ", num2str(custo), ":"));
        for i = 1:n
            disp([i, v(i)]);
        endfor
        return
    endif
    display("Ha custo reduzido negativo");
    disp(["Entra na base: ", int2str(j)]);

    u = tableau(2:end, j+1);
    disp(["Direcao a ser percorrida para proxima instancia de x: u=", mat2str(u)]);
    
    t = -1;
    teta = 0;
    hasPositive = false;
    
    for i = 1:m
        if (u(i) > 0)
            val = x_b(i)/u(i);
            if (!hasPositive)
                t = i;
                teta = val;
                hasPositive = true;
            endif
            if (val < teta)
                t = i;
                teta = val;
            endif
        endif
    endfor
    if (hasPositive)
        display("Valor de u positivo encontrado.");
        disp(["Sai da base: ", int2str(indB(t))]);
        disp(["Teta*: ", num2str(teta)]);
        

        indB(t) = j;
        % calculando o proximo tableau
        for i = 1:m+1
            if (i != t+1)
                tableau(i, :) += tableau(t+1, :)*(-tableau(i, j + 1)/tableau(t +1, j+ 1));
            endif
        endfor
        tableau(t+1, :) = tableau(t+1, :)/(u(t));

        display("Proximo tableau encontrado:");
        tableau
        [ind v] = simplex_tab_intern(indB, tableau, iteracao + 1);
        return
    else 
        disp("Problema de solucao ilimitada");
        ind =-1 ;
        d = zeros(1, n);
        d_B = u;
        for i = 1:m
            d(indB(i)) = d_B(i);
        endfor
        d(j) = 1;
        v=d;
        disp(["direcao em que decresce infinitamente: ", mat2str(d)]);
        return
    endif
endfunction
function [ind v ] = simplex_tab(indB, tableau)
    iteracao = 1;
    [ind v ] = simplex_tab_intern(indB, tableau, iteracao)
    return
endfunction 

% b = [20 20 20];
% c = [-10  -12 -12 0 0 0];
% A = [ [1; 2; 3] [2; 1; 2] [2; 2; 4] [1; 0; 0] [0; 1; 0] [0; 0; 1] ];
% m = 3;
% n = 6;
% x = [0, 0, 0, 20, 20, 20];
% indB = [4, 5, 6];

% indB = [4, 5, 6];
% tableau = [ [0; 20; 20; 20] [-10; 1; 2; 2 ] [-12; 2; 1; 2] [-12; 2; 2; 1] [0; 1; 0; 0] [0; 0; 1; 0] [0; 0; 0; 1]];

% indB = [5, 6, 7];
% tableau = [ [0; 15; 120; 100] [-4; 1; 7; 3] [-5; 1; 5; 5] [-9; 1; 3; 10] [-11; 1; 2; 15] [0; 1; 0; 0] [0; 0; 1; 0] [0; 0; 0; 1] ]

indB = [1 2];
tableau = [ [-31; 4; 3] [0; 1; 0] [0; 0; 1] [4; -2; 1] [-2; -1; 1]];

% inf
% indB = [1 2];
% tableau = [ [-31; 4; 3] [0; 1; 0] [0; 0; 1] [4; -2; 1] [-2; -1; 1] [-1; 0; 0]];


simplex_tab(indB, tableau)