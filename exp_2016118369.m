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

