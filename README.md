# SW_Numerical_Analysis

## 과제 #1) Taylor Series를 이용하여, x가 구간 [0, 1]에 대해  파형과 함수값들을 나타내는 function 파일 작성

![func](https://user-images.githubusercontent.com/58457978/70209778-7c766300-1774-11ea-8772-77c4835b4b74.jpg)

> exp_2015XXXXXX(0,1,21,0.5,6)

![KakaoTalk_20191205_154749328](https://user-images.githubusercontent.com/58457978/70210701-ae88c480-1776-11ea-9adf-8686818a55a2.png)

→ 0에서 1까지의  파형과 함수값들이 나타남

<pre><code>
function y = exp_2016118369(min,max,data,a,deg)
% min: 분석할 값의 최소값
% max: 분석할 값의 최대값
% data: 데이터 개수
% a: a의 값
% deg: 차수

result = zeros(deg+1,data); % (차수+1) * (data) 크기의 영벡터 생성
x = linspace(min,max,data); % 최소값부터 최대값까지 data만큼 쪼갬
result(1,:) = exp(-a);

for i = 1:deg % 테일러 시리즈를 계산하는 알고리즘
    if mod(i,2)
    result(i+1, :) = result(i, :) - ((x-a).^i .* exp(-a)) ./ factorial(i);
    else
    result(i+1, :) = result(i, :) + ((x-a).^i .* exp(-a)) ./ factorial(i);
    end
end

y = result(deg+1,:); % 결과 벡터를 변수 y에 저장
plot(x,y,'g:o'),grid; % 결괏값 그래프로 그리기

end
</code></pre>

----------

## 과제 #2) power method로 대칭 행렬의 가장 큰 eigenvalue와 가장 작은 eigenvalue를 구하고, 그에 따른 eigenvector 구하기

> A = [ 40 –20 0; -20 40 –20; 0 –20 40 ];  
> es = 0.05; maxit = 100;  
> [eig_val, eig_vec ] = power_2015XXXXXX(A, es, maxit)  
Eigen value: the 1. largest 2. smallest ? 1  

→ 가장 큰 eigenvalue를 선택할 것인지(선택1), 가장 작은 eigenvalue를 선택할 것인지(선택2)에 따라 eigenvalue와 그에 따른 eigenvector가 나타남

<pre><code>
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
</code></pre>
