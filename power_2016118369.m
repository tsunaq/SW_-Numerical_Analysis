function [eig_val, eig_vec] = power_2016118369(in_A,in_es,in_maxit)
% in_A: 입력받을 행렬
% in_es: 입력받을 허용 오차
% in_maxit: 입력받을 최대 반복수

    n = input('Eigen value: the 1. largest 2. smallest ? ');
% 변수 n에 1이나 2를 입력받아 저장
    
    [mat_a, mat_b] = size(in_A); % mat_a, mat_b에 행렬의 크기를 저장
    temp = mat_a & mat_b; % 정방행렬인지 검사하기 위한 변수 temp

    if n == 2 % 입력받은 n이 2라면 알고리즘에 의해 in_A에 역행렬을 저장
        in_A = inv(in_A);        
    end

    if ~(issymmetric(in_A) && temp && mat_a > 2) % 주어진 조건이 아니라면 에러 반환
        error('Wrong Operation');        
    end
    
    X = ones(mat_a,1); % 초기값 설정
    lamda = max(abs(X)); % lamda의 값을 X의 절대값 중 가장 큰 값으로 설정
    lamda_old = lamda;
    
    for i = 1:in_maxit % 1부터 maxit까지 반복
        lamda = max(abs(in_A*X));
        X = in_A*X./lamda;
        ea = 100 * abs((lamda-lamda_old)/lamda); % 오차 계산

        if ea <= in_es % 오차가 허용 오차보다 작으면 반복문 종료
           break;
        end
        
        lamda_old = lamda; % 허용오차 계산을 위한 lamda_old 변수
    end
    
    eig_val = lamda;
    eig_vec = X;
    if n == 2 % 입력받은 n이 2라면 알고리즘에 의해 eig_val의 값을 역수로 저장
            eig_val = 1/eig_val;
    end
end