function C = detect_error(c, cutoff_k, range)
%********************************************%
% �жϺ���
% input:c:finite-difference coefficient
%       cutoff_k:��ֹ����
%       range:������Χ
% output:C:0,����Ҫ��1,������
%********************************************%
k=0:0.001:1;
N=length(c);
st=0;
C=0;

for n=1:N
    st=st+c(n)*sin(n*k*pi);
end
Et=(2.*st)-(k*pi); 

for i=1:cutoff_k*1000+1
    if(abs(Et(i))>=range)
        C=C+1;
        
    end
end

end