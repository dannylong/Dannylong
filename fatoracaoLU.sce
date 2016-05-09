clear;

A = [3 2 4; 1 1 2; 4 3 2];

// Fatoração LU sem pivotamento parcial
function [L, U] = fatoracaoLU(A)
    
    [linhas colunas] = size(A);
    L = eye(linhas, colunas);
    for j = 1:colunas
        pivo = A(j,j);
        for i = (j+1):linhas
            //zerar todos os termos abaixo do pivo
            //zerar o elemento Ab(i,j) 
            lambda = A(i,j)/pivo;
            A(i,:) = A(i, :) - lambda * A(j,:);
            L(i, j) = lambda;
        end
    end
    // In the end(for a LS of order 3, e. g.), L = [1 0 0; lambda21 1 0; lambda31 lambda32 1];
    
    // ----------------------------------------------------------------------------------------------------------------------------    
    
    // U recebe todos os termos da matriz A modificada (matriz triangular superior), mas deve ser corrigida posteriormente, pois...
    // erros de arredondamente podem deixar termos abaixo dos pivôs muito pequenos, mas ainda assim não nulos
    U = A;
    // Garantir que a parte à esquerda(ou abaixo, sei lá xD) da diagonal principal tenha todos os termos nulos
    for j = 1:colunas
        for i = (j+1):linhas
            U(i, j) = 0;
        end
    end 
    
endfunction
