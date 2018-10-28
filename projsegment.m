function proj=projsegment(a,b)
f=-dot(a,b-a);
s=norm(b-a)^2;

if f<=0
    proj=a;
else
    if (f>0) && (f<s)
        proj=a +(f/s)*(b-a);
    else
        proj=b;
    end
end
end
