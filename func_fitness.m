function E=func_fitness(c, cutoff_k)

k=0:0.001:1;
N=length(c);
st=0;
cutoff_K=cutoff_k*1000;

for n=1:N
    st=st+c(n)*sin(n*k*pi);
end
Et=(2.*st)-(k*pi);

E=sum(abs(Et(1:cutoff_K)))/cutoff_K;
% for i=1:cutoff_k*1000+1
%     if(abs(Et(i))>=range)
%         C=C+1;
%     end
% end

end