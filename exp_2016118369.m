function y = exp_2016118369(min,max,data,a,deg)
% min: �м��� ���� �ּҰ�
% max: �м��� ���� �ִ밪
% data: ������ ����
% a: a�� ��
% deg: ����

result = zeros(deg+1,data); % (����+1) * (data) ũ���� ������ ����
x = linspace(min,max,data); % �ּҰ����� �ִ밪���� data��ŭ �ɰ�
result(1,:) = exp(-a);

for i = 1:deg % ���Ϸ� �ø�� ����ϴ� �˰���
    if mod(i,2)
    result(i+1, :) = result(i, :) - ((x-a).^i .* exp(-a)) ./ factorial(i);
    else
    result(i+1, :) = result(i, :) + ((x-a).^i .* exp(-a)) ./ factorial(i);
    end
end

y = result(deg+1,:); % ��� ���͸� ���� y�� ����
plot(x,y,'g:o'),grid; % �ᱣ�� �׷����� �׸���

end

