function center_point = find_center(u, x)
L = x(end) - x(1); u_2 = dot(u,u);
xx = x; center_point = 0;
for ii=1:10
    c = dot(u, xx.*u)/u_2;
    if abs(c)<1e-3
        break;
    else
        center_point = center_point + c;
        xx = xx - c;
        xx(xx<-L/2) = xx(xx<-L/2) + L;
        xx(xx>L/2) = xx(xx>L/2) - L;
    end
end