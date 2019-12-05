function [eig_val, eig_vec] = power_2016118369(in_A,in_es,in_maxit)
% in_A: �Է¹��� ���
% in_es: �Է¹��� ��� ����
% in_maxit: �Է¹��� �ִ� �ݺ���

    n = input('Eigen value: the 1. largest 2. smallest ? ');
% ���� n�� 1�̳� 2�� �Է¹޾� ����
    
    [mat_a, mat_b] = size(in_A); % mat_a, mat_b�� ����� ũ�⸦ ����
    temp = mat_a & mat_b; % ����������� �˻��ϱ� ���� ���� temp

    if n == 2 % �Է¹��� n�� 2��� �˰��� ���� in_A�� ������� ����
        in_A = inv(in_A);        
    end

    if ~(issymmetric(in_A) && temp && mat_a > 2) % �־��� ������ �ƴ϶�� ���� ��ȯ
        error('Wrong Operation');        
    end
    
    X = ones(mat_a,1); % �ʱⰪ ����
    lamda = max(abs(X)); % lamda�� ���� X�� ���밪 �� ���� ū ������ ����
    lamda_old = lamda;
    
    for i = 1:in_maxit % 1���� maxit���� �ݺ�
        lamda = max(abs(in_A*X));
        X = in_A*X./lamda;
        ea = 100 * abs((lamda-lamda_old)/lamda); % ���� ���

        if ea <= in_es % ������ ��� �������� ������ �ݺ��� ����
           break;
        end
        
        lamda_old = lamda; % ������ ����� ���� lamda_old ����
    end
    
    eig_val = lamda;
    eig_vec = X;
    if n == 2 % �Է¹��� n�� 2��� �˰��� ���� eig_val�� ���� ������ ����
            eig_val = 1/eig_val;
    end
end