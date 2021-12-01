
0;
warning('off');
function [ind v ] = simplex_res_intern(A,b,c,m,n,x,indB, Binv, iteracao)
    disp(["Iterando: ", int2str( iteracao)]);
    B = [];
    disp("Variaveis basicas:");
    for b_i = indB
        disp ([b_i, x(b_i)])
        B = [B A(:, b_i)];
    endfor
    disp(["Valor da funcao objetivo: ", mat2str(c*x')]);
    
    display(B);


    if (mod(iteracao, 10) == 0)
        display("Binv recalculado");
        Binv = inv(B);
    endif
    display(Binv) 

    c_B = c(indB);
    disp(["Vetor de custos da base: ", mat2str(c_B)]);
    p = c_B*Binv; # deveria ser c_B' mas os vetores aqui ja sao deitados por padrao
    disp(["p: ", mat2str(p)]);

    % calcular custos reduzidos para as variaveis nao basicas usando
    % c_j = c(j) - p * A(:, j) 
    disp(["Custos reduzidos das nao basicas (ate primeiro negativo): "]);
    hasNegativeCost = false;
    for j = 1:n
        if !any( indB == j)
            c_red = c(j) - p * A(:, j);
            disp([j, c_red]);
            if (c_red < 0)
                hasNegativeCost = true;
                break
            endif
        endif
    endfor
    if(!hasNegativeCost)
        disp(strcat("Solucao otima encontrada com custo= ", mat2str(c*x'), ":"));
        for k = 1:n
            disp([k, x(k)]);
        endfor
        ind = 1;
        v = x;
        return
    endif
    display("Ha custo reduzido negativo");
    disp(["Entra na base: ", int2str(j)]);

    u = Binv * A(:, j);
    disp(["Direcao a ser percorrida para proxima instancia de x: u=", mat2str(u)]);
    t = -1;
    teta = 0;
    hasPositive = false;
    for i = 1:m
        if (u(i) > 0)
            val = x(indB(i))/u(i);
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
        
        x(j) = teta;
        for i = 1:m
            if (i != t)
                x(indB(i)) -= teta*u(i);
            else
                x(indB(i)) = 0;
            endif
        endfor
        indB(t) = j;
        disp("\n");
        % calculando a proxima base inversa   
        nextBinv = Binv;
        for i = 1:m
            if (i != t)
                nextBinv(i, :) += nextBinv(t, :)*(-u(i)/u(t));
            endif
        endfor
        nextBinv(t, :) = nextBinv(t, :)/(u(t));

        display("Proxima Base de indices (invertida) encontrada:");
        nextBinv
        [ind v] = simplex_res_intern(A,b,c,m,n,x,indB, nextBinv, iteracao + 1);
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
function [ind v ] = simplex_res(A,b,c,m,n,x,indB, Binv)
    iteracao = 1;
    [ind v ] = simplex_res_intern(A,b,c,m,n,x,indB, Binv, iteracao)
    return
endfunction 

% EXEMPLO 1
% b = [20 20 20];
% c = [-10  -12 -12 0 0 0];
% A = [ [1; 2; 3] [2; 1; 2] [2; 2; 4] [1; 0; 0] [0; 1; 0] [0; 0; 1] ];
% m = 3;
% n = 6;
% x = [0, 0, 0, 20, 20, 20];
% indB = [4, 5, 6];

% EXEMPLO 2
% b = [20 20 20];
% c = [-10  -12 -12 0 0 0];
% A = [ [1; 2; 2] [2; 1; 2] [2; 2; 1] [1; 0; 0] [0; 1; 0] [0; 0; 1] ];
% m = 3;
% n = 6;
% x = [0, 0, 0, 20, 20, 20];
% indB = [4, 5, 6];
% B = [ [1; 0; 0] [0; 1; 0] [0; 0; 1]];
% Binv = inv(B);

% EXEMPLO 3
% b = [15 120 100];
% c = [-4  -5 -9 -11 0 0 0];
% A = [ [1; 7; 3] [1; 5; 5] [1; 3; 10] [1; 2; 15] [1; 0; 0] [0; 1; 0] [0; 0; 1] ];
% m = 3;
% n = 7;
% x = [0, 0, 0, 0, 15, 120, 100];
% indB = [5, 6, 7];
% B = [ [1; 0; 0] [0; 1; 0] [0; 0; 1]];
% Binv = inv(B);

% EXEMPLO 4
% A = [ [1; 0] [2; 1] [0; 1] [1; 1] ];
% b = [10 3];
% c = [4  5  1 -1];
% m = 2;
% n = 4;
% x = [4 3 0 0];
% indB = [1 2];
% B = [ [1; 0] [2; 1] ];
% Binv = inv(B);

% EXEMPLO 4
% solucao ilimitada
% b = [10 3];
% c = [4  5  1 -1 -1];
% A = [ [1; 0] [2; 1] [0; 1] [1; 1] [0; 0] ];
% m = 2;
% n = 5;
% x = [4 3 0 0 0];
% indB = [1 2];
% B = [ [1; 0] [2; 1] ];
% Binv = inv(B);

% Para executar um exemplo, descomente a secao do exemplo a sua escolha e execute o arquivo
simplex_res(A,b,c,m,n,x,indB, Binv)