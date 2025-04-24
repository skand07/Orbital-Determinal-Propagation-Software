function [r0,v0] = laplace_cal(lat,lst,alt,ra,dec,JD)

mu = 398600.4354; %km^3/s^2
R_earth = 6378.1366; %km
wE = 7.2921159*10^(-5); %rad/s


% calculations:
L = cell(1,3);
for i=1:1:3
    L{i}=[cosd(dec(i))*cosd(ra(i)); cosd(dec(i))*sind(ra(i)); sin(dec(i))];
end

T = JD*86400; %time in seconds

% differentiating to get first and second derivatives of L 
L_td1 = ((T(2)-T(3)) / ((T(1)-T(2))*(T(1)-T(3))) * L{1} + ...
           (2*T(2)-T(1)-T(3)) / ((T(2)-T(1))*(T(2)-T(3))) * L{2} + ...
           (T(2)-T(1)) / ((T(3)-T(1))*(T(3)-T(2))) * L{3});

L_td2 = (2 / ((T(1)-T(2))*(T(1)-T(3))) * L{1} + ...
             2 / ((T(2)-T(1))*(T(2)-T(3))) * L{2} + ...
             2 / ((T(3)-T(1))*(T(3)-T(2))) * L{3});

% finding rsite vectors
f_E = 1/298.257223563; 
R_site = R_earth + alt;
sl = sind(lat);
e2 = f_E * (2 - f_E);


p = (R_earth/sqrt(1-e2*sl*sl) + alt) * cosd(lat);
%p = (R_site) * cosd(lat);

q = (R_earth*(1-e2)/sqrt(1-e2*sl*sl) + alt) * sl;
%q = (R_site) * sl;

ct = cosd(lst(2)); 
st = sind(lst(2));
rsite = [p*ct; p*st; q];
rs_d1 = [cross([0 0 wE],rsite)]';
rs_d2 = [cross([0 0 wE],rs_d1)]';


% initial determinant:
D = [L{2} L_td1 L_td2];
D = 2*det(D);

%determinants D1, D2:
D1 = [L{2} L_td1 rs_d2];
D1 = det(D1);
D2 = [L{2} L_td1 rsite];
D2 = det(D2);


% 8th degree poly for r2:
C = dot(L{2},rsite);
R=norm(rsite);
Z = roots([1 0 (4*C*D1/D-4*D1^2/D^2-R^2) 0 0 mu*(4*C*D2/D-8*D1*D2/D^2) 0 0 -4*mu^2*D2^2/D^2]);
for (i = 1:numel(Z))
    if (isreal(Z(i)) && Z(i) > 0)
        r2 = Z(i);
        break;
    end
end

% calculate range rho2
rho2 = -2*D1/D-2*mu/r2^3*D2/D;

%establish determinants D3, D4:
D3 = [L{2} rs_d2 L_td2];
D3 = det(D3);

D4 = [L{2} rsite L_td2];
D4 = det(D4);


%calculate range rate rhodot2:
rhodot2 = D3/D-mu/r2^3*D4/D;

% solve for r0 in ECI 
r0 = rho2*L{2}+rsite;

%solve for v0 in ECI
v0 = rhodot2*L{2}+rho2*L_td1+rs_d1;

end
