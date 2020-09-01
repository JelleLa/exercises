function [ips]=myfind_multiple(S,str)
n=length(S);
m=length(str);
for j=1:m
    ip = [];
    for i=1:n
        if (strcmp(S(i),str(j)))
            ip=[ip i];
        end
    end
    if isempty(ip)
        fprintf('%s NOT FOUND\n',str(j));
    end
    ips(j).ip = ip;
end