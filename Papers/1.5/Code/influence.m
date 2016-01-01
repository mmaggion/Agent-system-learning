function [a]=influence(distance,choice)

switch choice
    
    case 1
        
        a = heaviside(distance - 5).*-100.*((10).^2./(distance).^3 - (10).^1.5./(distance).^2.5);
        a = a + heaviside(5 - distance).*-100.*((10).^2./(5).^3 - (10).^1.5./(5).^2.5);
         
    case 2
        
        %initial conditions:
        %x00(i,1)=4*cos(i+sqrt(2)) + rand;
        %x00(i,2)=4*cos(i+2*sqrt(2)) + rand;
        a = (heaviside(0.3 - distance).*log(distance));
        a = distance.*(a + sin(distance.*3).*(heaviside(distance - 0.3).*heaviside(10 - distance).*log(distance) + heaviside(distance - 10).*log(10)));
        
    case 3
        
        %initial conditions:
        %x00(i,1)=2*cos(i+sqrt(2)) + rand;
        %x00(i,2)=2*cos(i+2*sqrt(2)) + rand;
        a = (heaviside(0.3 - distance).*log(distance));
        a = distance.*(a + sin(distance.*3).*heaviside(distance - 0.3).*heaviside(10 - distance).*log(distance) + sin(30).*heaviside(distance - 10).*log(10));
        
    case 4
        
        %initial conditions:
        %x00(i,1)=0.1*cos(i+sqrt(2)) + rand;
        %x00(i,2)=0.1*cos(i+2*sqrt(2)) + rand;
        a = heaviside(0.3 - distance).*log(distance);
        a = a + heaviside(distance - 0.3).*heaviside(5 - distance).*cos((distance - 0.3).*7)./cos(7).*log(0.3).^2;
        a = a - heaviside(distance - 5).*cos(4.7*7)./cos(7).*log(0.3).*(distance - 5).^2;
        
    case 5
        
        a = heaviside(1 - distance).*(-2./(distance+0.3).^2 + 4);
        a = a + heaviside(distance - 1).*(-2./(1.3).^2 + 4).*sin(distance.*pi./2).*cos(distance - 1).*sin(distance.^(1.5) + pi/2 - 1);
        
    case 6
        
        a = -2./(distance+0.3).^2 + 4;
        a = a - sin(distance./5).*100./(3 + (distance - 5).^2) + sin(distance.*20);
        b = -2./(15+0.3).^2 + 4 - sin(15./5).*100./(3 + (15 - 5).^2) + sin(15.*20);
        a = heaviside(15 - distance).*a + heaviside(distance - 15).*(b - (distance - 15).^2);
        a = 0.1*a;
        
    case 7
        
        a = -10*exp(-(distance./10).^2).*cos(4*distance);
        
    case 8
        
        r_m = 0.7;
        eps = 100;
        gamma = 5;
        x_f = 0.6;
        y_f = -2.*r_m.*eps.*(1./(x_f).^2 - r_m./(x_f).^3);
        y_F = -2.*r_m.*eps.*(3*r_m./(x_f).^4 - 2./(x_f).^3);
        B = -y_F/(gamma*x_f^(gamma-1)*y_f);
        A = y_f*exp(B*x_f^gamma);
        
        a = -2.*r_m.*eps.*(1./(distance).^2 - r_m./(distance).^3);        
        a = A*exp(-B*(distance).^gamma).*heaviside(x_f - distance) + a.*heaviside(distance - x_f);
        a = -0.1*a;
        
end



end